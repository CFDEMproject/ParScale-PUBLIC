/*------------------------------------------------------------------------------------*\

                                      /$$$$$$                      /$$          
                                     /$$__  $$                    | $$          
        /$$$$$$   /$$$$$$   /$$$$$$ | $$  \__/  /$$$$$$$  /$$$$$$ | $$  /$$$$$$ 
       /$$__  $$ |____  $$ /$$__  $$|  $$$$$$  /$$_____/ |____  $$| $$ /$$__  $$
      | $$  \ $$  /$$$$$$$| $$  \__/ \____  $$| $$        /$$$$$$$| $$| $$$$$$$$
      | $$  | $$ /$$__  $$| $$       /$$  \ $$| $$       /$$__  $$| $$| $$_____/
      | $$$$$$$/|  $$$$$$$| $$      |  $$$$$$/|  $$$$$$$|  $$$$$$$| $$|  $$$$$$$
      | $$____/  \_______/|__/       \______/  \_______/ \_______/|__/ \_______/
      | $$                                                                      
      | $$                                                                      
      |__/        A Compilation of Particle Scale Models

   Copyright (C): 2014 DCS Computing GmbH (www.dcs-computing.com), Linz, Austria
                  2014 Graz University of Technology (ippt.tugraz.at), Graz, Austria
---------------------------------------------------------------------------------------
License
    ParScale is licensed under the GNU LESSER GENERAL PUBLIC LICENSE (LGPL).

    Everyone is permitted to copy and distribute verbatim copies of this license
    document, but changing it is not allowed.

    This version of the GNU Lesser General Public License incorporates the terms
    and conditions of version 3 of the GNU General Public License, supplemented
    by the additional permissions listed below.

    You should have received a copy of the GNU Lesser General Public License
    along with ParScale. If not, see <http://www.gnu.org/licenses/lgpl.html>.

	This code is designed to simulate transport processes (e.g., for heat and
	mass) within porous and no-porous particles, eventually undergoing
	chemical reactions.

	Parts of the code were developed in the frame of the NanoSim project funded
	by the European Commission through FP7 Grant agreement no. 604656.
\*-----------------------------------------------------------------------------------*/

#include "comm.h"
#include "input.h"
#include "output.h"
#include "control.h"
#include "particle_data.h"
#include "simulation_state.h"
#include "coupling_model.h"
#include "coupling.h"
#include "stdlib.h"
#include "vector_pascal.h"
#include "mpi_pascal.h"

using namespace PASCAL_NS;

/* ----------------------------------------------------------------------
   Constructor / Destructor
------------------------------------------------------------------------- */

Comm::Comm(ParScale *ptr,MPI_Comm &communicator) : ParScaleBaseAccessible(ptr),
    me_(-1),
    nprocs_(-1),
    world_(communicator),
    buf_send_(0),
    buf_recv_(0),
    maxsend_(0),
    maxrecv_(0)
{
    vectorZeroize3D(boxlo_);
    vectorZeroize3D(boxhi_);
    vectorZeroize3D(subboxlo_);
    vectorZeroize3D(subboxhi_);
    vectorZeroize3D(myloc_);
    vectorZeroize2D(procneigh_[0]);
    vectorZeroize2D(procneigh_[1]);
    vectorZeroize2D(procneigh_[2]);

   // do some MPI stuff

    MPI_Comm_rank(world_,&me_);
    MPI_Comm_size(world_,&nprocs_);

    maxsend_ = BUFMIN();
    create(buf_send_,maxsend_+BUFEXTRA());
    maxrecv_ = BUFMIN();
    create(buf_recv_,maxrecv_);
}

Comm::~Comm()
{
//    MPI_Comm_free(&world_);
}

/* ----------------------------------------------------------------------
   Pull box properties from LIGGGHTS and exchange atoms
------------------------------------------------------------------------- */

void Comm::pull()
{
    // no action required for serial
    if(1 == nprocs_)
        return;
    // this is executed before Coupling::pull() ParticleData::pull(), so this is before
    // all other pull operations
    if(coupling().external_code_in_control())
    {
        CouplingModel &cm = coupling().couplingModel();
        // pull simulation box properties and 3d proc grid infos
        error().throw_error_one(FLERR,"call coupling fct that pulls box extents from LIGGGHTS\n");
        cm.pull_box(boxlo_,boxhi_,subboxlo_,subboxhi_);
        cm.pull_proc_info(procgrid_,myloc_,procneigh_);
        // exchange particles to neighboring procs
        exchange();
    }
}

/* ----------------------------------------------------------------------
   settings
------------------------------------------------------------------------- */

void Comm::init()
{
    //TODO check if CommModel allocated
}

/* ----------------------------------------------------------------------
   settings
------------------------------------------------------------------------- */

void Comm::parse_command(int narg,char const* const* arg)
{
    output().write_screen_one("Comm is parsing something:");
    output().write_screen_one(arg[0]);
}

/* ----------------------------------------------------------------------
   abort simulation, called by one proc (crash) or
   finalize simulation, called by all proc (clean termination)
------------------------------------------------------------------------- */

void Comm::abort_one() const
{
    //MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Abort(world_,1);
}

void Comm::finalize_all() const
{
    //TODO: need barrier here?!

    MPI_Finalize();
    exit(1);
}

/* ----------------------------------------------------------------------
   broadcast data so all procs have it
------------------------------------------------------------------------- */

void Comm::bcast()
{
    // no action required for serial
    if(1 == nprocs_)
        return;

    int nsend = 0;
    MPI_Request request;
    MPI_Status status;

    // scale translate rotate not needed here
    // OPERATION_COMM_BCAST means that only those properties
    // which have been read from file are bcasted
    OperationProperties op(OPERATION_COMM_BCAST,false,false,false);

    for (int dim = 0; dim < 3; dim++)
    {
        // no action if no neigh procs in this dim
        if (procgrid_[dim] == 1)
            continue;

        // push data to buffer
        // fill buffer with elements leaving my box, using < and >=
        // when elements is deleted, fill it in with last atom

        if(0 == me_)
            nsend = pushExchangeBcast(dim,op);

        // bcast number of elements so can grow buf if necessary
        MPI_Bcast(&nsend,1,MPI_INT,0,world_);
        if (nsend > maxsend_)
            grow_send(nsend,1);

        // bcast data
        MPI_Bcast(buf_send_,nsend,MPI_DOUBLE,0,world_);

        // pop data if in my subbox
        popExchangeBcast(nsend,dim, buf_send_,op);
    }

    // re-calculate nbody_all, error if lost particle
    particleData().data().recalc_nbody_all(true);
}

/* ----------------------------------------------------------------------
   exchange elements with nearby processors
------------------------------------------------------------------------- */

void Comm::exchange()
{
    // no action required for serial
    if(1 == nprocs_)
        return;

    int nrecv, nsend = 0;
    int nrecv1,nrecv2;
    double *buf;
    MPI_Request request;
    MPI_Status status;

    // scale translate rotate not needed here
    OperationProperties op(OPERATION_COMM_EXCHANGE,false,false,false);

    for (int dim = 0; dim < 3; dim++)
    {
        // push data to buffer
        // fill buffer with elements leaving my box, using < and >=
        // when elements is deleted, fill it in with last atom

        nsend = pushExchangeBcast(dim,op);

        // send/recv in both directions
        // if 1 proc in dimension, no send/recv, set recv buf to send buf
        // if 2 procs in dimension, single send/recv
        // if more than 2 procs in dimension, send/recv to both neighbors

        if (procgrid_[dim] == 1)
        {
            nrecv = nsend;
            buf = buf_send_;
        }
        else
        {
            MPI_Sendrecv(&nsend,1,MPI_INT,procneigh_[dim][0],0,&nrecv1,1,MPI_INT,procneigh_[dim][1],0,world_,&status);
            nrecv = nrecv1;

            if (procgrid_[dim] > 2)
            {
                MPI_Sendrecv(&nsend,1,MPI_INT,procneigh_[dim][1],0,&nrecv2,1,MPI_INT,procneigh_[dim][0],0,world_,&status);
                nrecv += nrecv2;
            }

            if (nrecv > maxrecv_) grow_recv(nrecv);

            MPI_Irecv(buf_recv_,nrecv1,MPI_DOUBLE,procneigh_[dim][1],0,world_,&request);
            MPI_Send(buf_send_,nsend,MPI_DOUBLE,procneigh_[dim][0],0,world_);
            MPI_Wait(&request,&status);

            if (procgrid_[dim] > 2)
            {
                MPI_Irecv(&buf_recv_[nrecv1],nrecv2,MPI_DOUBLE,procneigh_[dim][0],0,world_,&request);
                MPI_Send(buf_send_,nsend,MPI_DOUBLE,procneigh_[dim][1],0,world_);
                MPI_Wait(&request,&status);
            }

            buf = buf_recv_;
        }

        // check incoming elements to see if they are in my box
        // if so, add on this proc

        popExchangeBcast(nrecv,dim, buf,op);

    }

    // re-calculate nbody_all, no error if lost particle
    // (could be by exit of simulation domain in liggghts)
    particleData().data().recalc_nbody_all(false);
}

/* ----------------------------------------------------------------------
   push data to buffer for atom exchange or for container broadcast
------------------------------------------------------------------------- */

int Comm::pushExchangeBcast(int dim,OperationProperties &op)
{
      // do NOT make a local copy of buf_send_ as this fct calls re-allocation!!

      double checklo,checkhi;
      CustomValueTracker &data = particleData().data();
      ContainerVector<double,3>& pos = *data.getElementProperty<ContainerVector<double,3> >("pos");

      checklo = subboxlo_[dim];
      if(subboxhi_[dim] == boxhi_[dim])
        checkhi = boxhi_[dim] + SMALL_COMM();
      else
        checkhi = subboxhi_[dim];

      int nsend = 0, nsend_this = 0;
      int i = 0;
      while(i < particleData().nbody())
      {
          if(!(pos(i)[dim] >= checklo && pos(i)[dim] < checkhi))
          {
              nsend_this = data.pushElemToBuffer(i,&(buf_send_[nsend+1]),op);
              buf_send_[nsend] = static_cast<double>(nsend_this+1);
              nsend += (nsend_this+1);

              if (nsend > maxsend_)
              {
                  grow_send(nsend,1);
              }
              data.deleteElement(i); //NP deleteElement() decreases nLocal

              char bufout[100];
              sprintf(bufout,"nBody %d, nbody_all %d",particleData().nbody(),particleData().nbody_all());
              error().throw_error_one(FLERR,"check if nLocal correct:",bufout);
              // maybe have to decrease bnlocal

          }
          else i++;
      }
      return nsend;
}

/* ----------------------------------------------------------------------
   pop data from buffer for atom exchange
------------------------------------------------------------------------- */

int Comm::popExchangeBcast(int nrecv,int dim,double *buf,OperationProperties &op)
{
    double center_elem[3];
    double checklo,checkhi;
    int m = 0, nrecv_this;

    CustomValueTracker &data = particleData().data();

    checklo = subboxlo_[dim];
    if(subboxhi_[dim] == boxhi_[dim])
        checkhi = boxhi_[dim] + SMALL_COMM();
    else
        checkhi = subboxhi_[dim];

    while (m < nrecv)
    {
        // number of values is first in buffer
        nrecv_this = static_cast<int>(buf[m]);

        // center is next in buffer, test it
        vectorCopy3D(&(buf[m+1]),center_elem);

        if(center_elem[dim] >= checklo && center_elem[dim] < checkhi)
        {
            // inc = true increases nlocal
            data.popElemFromBuffer(&(buf[m+1]),op,true);
        }

        m += nrecv_this;
    }

    return m;
}

  /* ----------------------------------------------------------------------
   realloc the size of the send buffer as needed with BUFFACTOR & BUFEXTRA
   if flag = 1, realloc
   if flag = 0, don't need to realloc with copy, just free/malloc
  ------------------------------------------------------------------------- */

  void Comm::grow_send(int n, int flag)
  {
      maxsend_ = static_cast<int> (BUFFACTOR() * n);

      if (flag)
        grow(buf_send_,maxsend_+BUFEXTRA());
      else {
        destroy(buf_send_);
        create(buf_send_,maxsend_+BUFEXTRA());
      }
  }

  /* ----------------------------------------------------------------------
   free/malloc the size of the recv buffer as needed with BUFFACTOR
  ------------------------------------------------------------------------- */

  void Comm::grow_recv(int n)
  {
      maxrecv_ = static_cast<int> (BUFFACTOR() * n);
      destroy(buf_recv_);
      create(buf_recv_,maxrecv_);
  }

