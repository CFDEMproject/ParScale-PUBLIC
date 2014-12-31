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

     Copyright (C): 2012 - 2014 DCS Computing GmbH (www.dcs-computing.com), Linz, Austria
                    2012 - 2014 Department of Particulate Flow Modelling, JKU Linz
                              (www.jku.at/pfm), Linz, Austria

   This file was originally part of LIGGGHTS (www.cfdem.com), and is now re-distributed
   under LGPL as part of ParScale with the permission of the copyright holders
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

#ifndef LMP_MPI_PASCAL_H
#define LMP_MPI_PASCAL_H

#include "mpi.h"
#include "stdio.h"

/* ---------------------------------------------------------------------- */
// a poor man's inline MPI wrappers for LIGGGHTS/ParScale
/* ---------------------------------------------------------------------- */

namespace PASCAL_NS
{

/* ----------------------------------------------------------------------
   Helper function to be able to templetize wrappers
------------------------------------------------------------------------- */

template<typename T>
inline MPI_Datatype mpi_type()
{
  printf("**************ILLEGAL CALL TO mpi_type()*************");
  return 0;
}

template<>
inline MPI_Datatype mpi_type<double>()
{
  return MPI_DOUBLE;
}

template<>
inline MPI_Datatype mpi_type<int>()
{
  return MPI_INT;
}

/* ---------------------------------------------------------------------- */

template<typename T>
inline void MPI_Sum_Vector(T* vector, int len, MPI_Comm comm)
{
  MPI_Allreduce(MPI_IN_PLACE, vector, len, mpi_type<T>(), MPI_SUM, comm);
}

/* ---------------------------------------------------------------------- */

template<typename T>
inline void MPI_Sum_Scalar(T& scalar, MPI_Comm comm)
{
  MPI_Allreduce(MPI_IN_PLACE, &scalar, 1, mpi_type<T>(), MPI_SUM, comm);
}

/* ---------------------------------------------------------------------- */

template<typename T>
inline void MPI_Sum_Scalar(T& scalar, T& scalar_all, MPI_Comm comm)
{
  MPI_Allreduce(&scalar, &scalar_all, 1, mpi_type<T>(), MPI_SUM, comm);
}

/* ---------------------------------------------------------------------- */

template<typename T>
inline void MPI_Min_Scalar(T& scalar, MPI_Comm comm)
{
  MPI_Allreduce(MPI_IN_PLACE, &scalar, 1, mpi_type<T>(), MPI_MIN, comm);
}

/* ---------------------------------------------------------------------- */

template<typename T>
inline void MPI_Min_Scalar(T scalar, T& scalar_all, MPI_Comm comm)
{
  MPI_Allreduce(&scalar, &scalar_all, 1, mpi_type<T>(), MPI_MIN, comm);
}

/* ---------------------------------------------------------------------- */

template<typename T>
inline void MPI_Max_Scalar(T& scalar, MPI_Comm comm)
{
  MPI_Allreduce(MPI_IN_PLACE, &scalar, 1, mpi_type<T>(), MPI_MAX, comm);
}

/* ---------------------------------------------------------------------- */

template<typename T>
inline void MPI_Max_Scalar(T scalar, T& scalar_all, MPI_Comm comm)
{
  MPI_Allreduce(&scalar, &scalar_all, 1, mpi_type<T>(), MPI_MAX, comm);
}

/* ---------------------------------------------------------------------- */

template<typename T>
inline void MPI_Max_Vector(T *vector, int len, MPI_Comm comm)
{
  MPI_Allreduce(MPI_IN_PLACE, vector, len, mpi_type<T>(), MPI_MAX, comm);
}


/* ---------------------------------------------------------------------- */

template<typename T>
inline void MPI_Min_Vector(T* vector, int len, MPI_Comm comm)
{
  MPI_Allreduce(MPI_IN_PLACE, vector, len, mpi_type<T>(), MPI_MIN, comm);
}

/* ---------------------------------------------------------------------- */

inline void MPI_Allgather_Sum_Scalar(int scalar,int &scalar_acc,MPI_Comm comm)
{
    int rank,size, *allg;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    allg = new int[size];

    MPI_Allgather(&scalar,1,MPI_INT,allg,1,MPI_INT,comm);

    scalar_acc = 0;
    for (int iproc = 1; iproc < rank; iproc++)
       scalar_acc = scalar_acc + allg[iproc-1];

    delete []allg;
}

/* ----------------------------------------------------------------------
   Gather vector data from all processors at proc 0
   returns allocated and populated array vector0 to caller
------------------------------------------------------------------------- */

inline int MPI_Gather0_Vector(double *vector, int size ,double *&vector_0,MPI_Comm comm)
{
    int me,nprocs, *recvcnts, *displs;
    int size_0;

    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &me);

    recvcnts = new int[nprocs];
    displs = new int[nprocs];

    //NP gather recvcount for each processor
    //NP MPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm)
    MPI_Allgather(&size,1,MPI_INT,recvcnts,1,MPI_INT,comm);

    size_0 = 0;
    displs[0] = 0;
    for (int iproc = 1; iproc < nprocs; iproc++)
    {
        size_0 += recvcnts[iproc-1];
        displs[iproc] = displs[iproc-1] + recvcnts[iproc-1];
    }
    size_0 += recvcnts[nprocs-1];


    //NP allocate
    if(me == 0)
        vector_0 = new double[size_0];
    else
        vector_0 = 0;

    //NP use MPI_Gatherv to gather vector data at proc 0
    //NP MPI_Gatherv(sendbuf, sendcnt, sendtype, recvbuf, recvcnts, displs, recvtype, root, comm)

    MPI_Gatherv(vector,size,MPI_DOUBLE,vector_0, recvcnts, displs, MPI_DOUBLE,0, comm);

    delete []recvcnts;
    delete []displs;

    return size_0;
}

/* ----------------------------------------------------------------------
   Allgather vector data from all processors
   returns allocated and populated array vector_all to caller
------------------------------------------------------------------------- */

template<typename T>
inline int MPI_Allgather_Vector(T *vector, int size ,T *&vector_all,MPI_Comm comm)
{
    int me,nprocs, *recvcnts, *displs;
    int size_all;

    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &me);

    recvcnts = new int[nprocs];
    displs = new int[nprocs];

    //NP gather recvcount for each processor
    MPI_Allgather(&size,1,MPI_INT,recvcnts,1,MPI_INT,comm);

    size_all = 0;
    displs[0] = 0;
    for (int iproc = 1; iproc < nprocs; iproc++)
    {
        size_all += recvcnts[iproc-1];
        displs[iproc] = displs[iproc-1] + recvcnts[iproc-1];
    }
    size_all += recvcnts[nprocs-1];

    //NP allocate
    vector_all = new T[size_all];

    //NP use MPI_Allgatherv to gather vector data at each proc
    MPI_Allgatherv(vector,size,mpi_type<T>(),vector_all, recvcnts, displs, mpi_type<T>(), comm);

    delete []recvcnts;
    delete []displs;

    return size_all;
}

}; // end namespace PASCAL_NS


#endif
