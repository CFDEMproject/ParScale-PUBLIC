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

#ifndef LMP_CONTAINER_GENERAL
#define LMP_CONTAINER_GENERAL

#include "error.h"
#include "container_base.h"
#include "memory_ns.h"
#include "input_properties.h"
#include "coupling_base.h"
#include "math_extra_pascal.h"
#include <string.h>
#include <fstream>

#define GROW 100

using namespace PASCAL_MEMORY_NS;
using namespace PASCAL_NS;

namespace PASCAL_NS
{

  template<typename T, int NUM_VEC, int LEN_VEC>
  class ContainerGeneral : public ContainerBase
  {
      public:

          bool isDoubleData();
          bool isIntData();
          bool isBoolData();

          void add(T** elem);
          void addZero();

          void copy(int from,int to);
          void del(int n);
          void delForward(int n,OperationProperties &op);
          void delRestart(int n,OperationProperties &op);
          void clearReverse(OperationProperties &op);
          void read(OperationProperties &op,InputBase const* input_base);
          void pull(OperationProperties &op,CouplingBase const* coupling_base);
          void write(OperationProperties &op,InputBase const* input_base);
          void push(OperationProperties &op,CouplingBase const* coupling_base);

          void get(int n, T** elem);

          void setAll(T def);
          void setAll(int to, T def);
          void set(int i, T** elem);
          void set(int i, int j, T* elem);

          bool setFromContainer(ContainerBase *cont);

          T max_scalar();
          T min_scalar();

          T**& operator()(int n);
          T** const& operator()(int n) const;
          T*** begin();
          virtual void* begin_slow_dirty();

          inline void scale(double factor);
          inline void move(double *dx);
          inline void moveElement(int i,double *dx);
          inline void rotate(double *dQ);

          // all push and pop functions return number of bytes taken from / added to buf
          // all push and pop functions expect buf to point to first element with usable data

          // push / pop all elements
          //NP used for restart
          inline int bufSize(OperationProperties &op);
          inline int pushToBuffer(double *buf, OperationProperties &op);
          inline int popFromBuffer(double *buf, OperationProperties &op);

          // push / pop a list elements
          //NP used for borders and forward comm
          //NP better to call this than calling for each element since cannot
          //NP be inlined

          //NP ASSUMPTION: all elements have the same properties / same communication need
          //NP so is ok to just have first and n (otherwise would mess up ghost data)
          inline int elemListBufSize(int n, OperationProperties &op);
          inline int pushElemListToBuffer(int n, int *list, double *buf, OperationProperties &op);
          inline int popElemListFromBuffer(int first, int n, double *buf, OperationProperties &op);
          inline int pushElemListToBufferReverse(int first, int n, double *buf, OperationProperties &op);
          inline int popElemListFromBufferReverse(int n, int *list, double *buf, OperationProperties &op);

          // push / pop one single element
          //NP used for exchange - not performance critical to call this for each element
          //NP since number of elements to exchange is typically low
          inline int elemBufSize(OperationProperties &op);
          inline int pushElemToBuffer(int i, double *buf, OperationProperties &op);
          inline int popElemFromBuffer(double *buf, OperationProperties &op);

          void addUninitialized(int n);

          inline int size()
          { return numElem_; }

          inline int nVec()
          { return NUM_VEC; }

          inline int lenVec()
          { return LEN_VEC; }

          inline int lenVecUsed()
          { return lenVecUsed_; }

          void setLenVecUsed(int x)
          {
            if(x>LEN_VEC)
            {
                lenVecUsed_ = LEN_VEC;
                isFull_     = true;
            }
            lenVecUsed_ = x;
            return;
           }

          inline int capacity()
          { return maxElem_; }

          inline void empty()
          { numElem_ = 0; }

          inline bool isFull()
          { return isFull_;}

          inline bool isFilled()
          { return isFilled_;}

      protected:

          ContainerGeneral(ParticleDataContainerProperties &cp);
          ContainerGeneral(ContainerGeneral<T,NUM_VEC,LEN_VEC> const &orig);
          virtual ~ContainerGeneral();

          // shall return the size of an entry in bytes
          int getElemSize();

          int numElem_, maxElem_;

          int lenVecUsed_; //used number of vector elements

          bool isFilled_;

          bool isDumped_;

          bool isFull_;

          bool verbose_; //TODO: read from JSON

          T*** arr_;
  };

  // *************************************
  #include "container_general_I.h"
  // *************************************

} /* LAMPPS_NS */
#endif /* CONTAINER_GENERAL_ */
