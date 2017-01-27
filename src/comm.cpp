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
    maxrecv_(0),
    data_(NULL)
{
    vectorZeroize3D(boxlo_);
    vectorZeroize3D(boxhi_);
    vectorZeroize3D(subboxlo_);
    vectorZeroize3D(subboxhi_);
    vectorZeroize3D(myloc_);
    vectorZeroize2D(procneigh_[0]);
    vectorZeroize2D(procneigh_[1]);
    vectorZeroize2D(procneigh_[2]);
    neighAgoCaller_ = -1;
    timeStepFromRun_ = -1;

    // do some MPI stuff
    MPI_Comm_rank(MPI_COMM_WORLD,&me_);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs_);

    maxsend_ = BUFMIN();
    create(buf_send_,maxsend_+BUFEXTRA());
    maxrecv_ = BUFMIN();
    create(buf_recv_,maxrecv_);

    create<int>(exchangeCounts_,       nprocs_);
    create<int>(exchangeDisplacements_,nprocs_);
    verbose_ = false; //true; //Developer to set here if needed for debugging

}

Comm::~Comm()
{
    destroy(buf_send_);
    destroy(buf_recv_);
    destroy<int>(exchangeCounts_);
    destroy<int>(exchangeDisplacements_);
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
        cm.pull_box(boxlo_,boxhi_,subboxlo_,subboxhi_);
        cm.pull_proc_info(procgrid_,myloc_,procneigh_);
        double timeStep;
        cm.pull_timeStepping_info(timeStep,neighAgoCaller_,timeStepFromRun_);

        exchangeEventsLocalId_          = cm.exchangeEventsLocalId();
        exchangeEventsReceivingProcess_ = cm.exchangeEventsReceivingProcess();

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

    if(verbose_)
    {
        output().write_screen_one("Comm is parsing something:");
        output().write_screen_one(arg[0]);
    }
}

/* ----------------------------------------------------------------------
   abort simulation, called by one proc (crash) or
   finalize simulation, called by all proc (clean termination)
------------------------------------------------------------------------- */

void Comm::abort_one() const
{
    //MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Abort(MPI_COMM_WORLD,1);
}

void Comm::finalize_all() const
{
    //TODO: need barrier here?!

    MPI_Finalize();
    exit(1);
}

/* ----------------------------------------------------------------------
   exchange elements with nearby processors
------------------------------------------------------------------------- */
void Comm::exchange() //this is similar the class MultiNodeMeshParallel in LIGGGHTS
{
    // no action required for serial
    if(1 == nprocs_)
        return;

    // exchange particles to neighboring procs only after neighboring, or during very first step
    if(neighAgoCaller_!=0 && timeStepFromRun_ != 1)
        return;

    int mpi_err;

    // scale translate rotate not needed here
    OperationProperties op(OPERATION_COMM_EXCHANGE,false,false,false);
    sizeExchangeEvents_   = exchangeEventsLocalId_->size();
    data_ = &(particleData().data());

    //Loop all processes to gather data to currProcess
    for (int currProcess=0; currProcess<nprocs_; currProcess++)
    {
        int nrecv=0; int nsend=0;

        //fill puffers with data for sending to currProcess
        if( currProcess!=me_ )
        {
            nsend = pushExchangeBcast(currProcess,op);
            if(verbose_ && nsend > 0)
              printf("[%d/%d]: **exchange: will send %d doubles to process %d \n",
                 me_, nprocs_,
                 nsend, currProcess
                );
        }

        wait(); //be sure we are in sync
        mpi_err = MPI_Gather(&nsend,           1,MPI_INT,
                              exchangeCounts_, 1,MPI_INT,
                              currProcess,MPI_COMM_WORLD
                            );

        //Calculate displacements and the size of the recv array
        if( currProcess==me_ )
        {
            exchangeDisplacements_[0] = 0;
            for(int iPro=1; iPro<nprocs_; iPro++)
                exchangeDisplacements_[iPro] = exchangeCounts_[iPro-1]
                                             + exchangeDisplacements_[iPro-1];

            nrecv=0;
            for(int iPro=0; iPro<nprocs_; iPro++)
                nrecv += exchangeCounts_[iPro];

            if (nrecv > maxrecv_) grow_recv(nrecv);

            if(verbose_ && nrecv > 0)
            {
              printf("[%d/%d]: **exchange: MPI_Gather yields the following exchangeCounts_/displacements:",
                      me_, nprocs_
                    );
              for(int k=0; k<nprocs_; k++)
                  printf(" %d/%d, ",
                          exchangeCounts_[k], exchangeDisplacements_[k]
                        );
              printf("\n");
             }
        }

        wait(); //be sure we are in sync
        mpi_err = MPI_Gatherv(buf_send_,   nsend,                                  MPI_DOUBLE,
                              buf_recv_,   exchangeCounts_, exchangeDisplacements_,MPI_DOUBLE,
                              currProcess, MPI_COMM_WORLD);

        if(verbose_ && nrecv > 0 )
          printf("[%d/%d]: **exchange: MPI_Gatherv to %d complete (error: %d)! Poping buffer with size %d... \n",
                 me_, nprocs_,
                 currProcess, 
                 mpi_err,
                 nrecv
                );

        if( currProcess==me_ && nrecv > 0 )
            popExchangeBcast(nrecv, buf_recv_,op);

    } //loop over all processes

    //Delete all elements that have been sent
    for(int event=(sizeExchangeEvents_-1); event>=0; event--)
        data_->deleteElement( (*exchangeEventsLocalId_)[event] );

    // re-calculate nbody_all, no error if lost particle
    // (could be by exit of simulation domain in liggghts)
    data_->recalc_nbody_all(false);

    exchangeEventsLocalId_->clear();
    exchangeEventsReceivingProcess_->clear();

}

/* ----------------------------------------------------------------------
   push data to buffer for atom exchange or for container broadcast
------------------------------------------------------------------------- */
int Comm::pushExchangeBcast(int currProcess,OperationProperties &op)
{
      int nsend=0; int nsend_this = 0;
      for(int event=0; event<sizeExchangeEvents_; event++)
      {
        if( currProcess==(*exchangeEventsReceivingProcess_)[event] )
        {
            int iToPush = (*exchangeEventsLocalId_)[event];
            if( iToPush >= particleData().nbody() )
            {
                char errorMsg[500];
                sprintf(errorMsg,"pushExchangeBcast: attempting to push body %d, but only have %d bodies.",
                        iToPush,particleData().nbody());
                error().throw_error_all(FLERR,errorMsg);
            }

            nsend_this = data_->pushElemToBuffer(iToPush,&(buf_send_[nsend+1]),op);
            buf_send_[nsend] = static_cast<double>(nsend_this+1);
            nsend += (nsend_this+1);
        }
      }

      if (nsend > maxsend_)
      {
            grow_send(nsend,1);
      }

      if( verbose_ && nsend>0 )
      printf("[%d/%d]: **pushExchangeBcast: size of buffer (nsend): %d\n",
              me_, nprocs_,
              nsend);

      return nsend;
}

/* ----------------------------------------------------------------------
   pop data from buffer for atom exchange
------------------------------------------------------------------------- */

int Comm::popExchangeBcast(int nrecv, double *buf, OperationProperties &op)
{
    int m = 0, nrecv_this;

    while (m < nrecv)
    {
        // number of values is first in buffer
        nrecv_this = static_cast<int>(buf[m]);

        if(verbose_)
           printf("[%d/%d]: **popExchangeBcast: receiving %d elements, start reading buffer at %d current nbody/nbody_all: %d/%d \n",
                   me_, nprocs_,
                   nrecv_this, m+1,
                   particleData().nbody(),particleData().nbody_all()
                  );

        data_->popElemFromBuffer(&(buf[m+1]),op,true);
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


/* ----------------------------------------------------------------------
   broadcast data so all procs have it -
------------------------------------------------------------------------- */
void Comm::bcast() //this is similar the class MultiNodeMeshParallel in LIGGGHTS
{
// DONT NEED THIS FUNCTION RIGHT NOW, SINCE
/*
    // no action required for serial
    if(1 == nprocs_)
        return;

    //can only execute if particle positions are present
    if(!particleData().haveParticlePositions())
        return;

    if(verbose_)
        output().write_screen_all("**Comm::bcast");

    int nsend = 0;
    MPI_Request request;
    MPI_Status status;

    // scale translate rotate not needed here
    // OPERATION_COMM_BCAST means that only those properties
    // which have been read from file are bcasted
    OperationProperties op(OPERATION_COMM_BCAST,false,false,false);

    for (int dim = 0; dim < 3; dim++)
    {
        nsend = pushExchangeBcast(dim,op);

        // no action if no neigh procs in this dim
        if (procgrid_[dim] == 1)
            continue;

        // push data to buffer
        // fill buffer with elements leaving my box, using < and >=
        // when elements is deleted, fill it in with last atom


        // bcast number of elements so can grow buf if necessary
        MPI_Bcast(&nsend,1,MPI_INT,0,world_); //DEFECT-TODO
        if (nsend > maxsend_)
            grow_send(nsend,1);

        if( nsend == 0)
            continue;

        if(verbose_)
            output().write_screen_all("**MPI_Bcasting now...");



//        // bcast data
//        MPI_Bcast(buf_send_,nsend,MPI_DOUBLE,0,world_);

//        // pop data if in my subbox
//        popExchangeBcast(nsend,dim, buf_send_,op);
    }

    // re-calculate nbody_all, error if lost particle
//    particleData().data().recalc_nbody_all(true); //DEFECT-TODO

*/
}
