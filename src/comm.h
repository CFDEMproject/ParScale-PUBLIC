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

#ifndef PASC_COMM_H
#define PASC_COMM_H

#include "stdio.h"
#include "pascal_base_accessible.h"
#include "pascal_base_interface.h"
#include "mpi.h"

namespace PASCAL_NS
{

class Comm : public ParScaleBaseAccessible, public ParScaleBaseInterface
{
    public:

      Comm(ParScale *ptr,MPI_Comm &communicator);
      ~Comm();

      // inherited from ParScaleBaseInterface

      void pull();
      void init();
      void bcast();
      void parse_command(int narg,char const* const* arg);

      // inline access

      inline int me()           const { return me_;}
      inline int nprocs()       const { return nprocs_;}
      inline bool is_proc_0()   const { return (0 == me_) ? true: false; }
      inline bool is_parallel() const { return (nprocs_ > 1) ? true: false; }

      void abort_one() const;
      void finalize_all() const;

      inline MPI_Comm world() const { return world_; }

      //TODO route this to CommModel
      void exchange();

    private:

      inline double SMALL_COMM()
      { return 1.0e-8; }

      inline double BUFFACTOR()
      { return 1.5; }

      inline int BUFEXTRA()
      { return 2000; }

      inline int BUFMIN()
      { return 2000; }

      int popExchangeBcast(int nrecv,int dim,double *buf,class OperationProperties &op);
      int pushExchangeBcast(int dim,class OperationProperties &op);

      void grow_send(int n, int flag);
      void grow_recv(int n);

      // TODO
      // class CommModel *commModel_;
      // can be CommModelCartesian or CommModelMany2Many

      int me_, nprocs_;        // number of MPI procs and my rank
      MPI_Comm &world_;        // MPI communicator

      // domain extent for global box and my box
      double boxlo_[3],boxhi_[3];
      double subboxlo_[3],subboxhi_[3];

      // 3D grid of processors

      int procgrid_[3];                  // procs assigned in each dim of 3d grid
      int myloc_[3];                     // which proc I am in each dim
      int procneigh_[3][2];              // my 6 neighboring procs, 0/1 = left/right

      // communication buffers

      double *buf_send_;                 // send buffer for all comm
      double *buf_recv_;                 // recv buffer for all comm
      int maxsend_,maxrecv_;              // current size of send/recv buffer

};

} //end PASCAL_NS

#endif
