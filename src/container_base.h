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

#ifndef LMP_CONTAINER_BASE_H
#define LMP_CONTAINER_BASE_H

#include "string.h"
#include "input_base.h"
#include "coupling_base.h"
#include "particle_data_container_properties.h"

namespace PASCAL_NS
{


  class ContainerBase
  {
      public:

          ContainerBase(class ParticleDataContainerProperties &cp);

          //NP need to make this virtual to have the correct destructor of derived class
          //NP called via delete in AssociativePointerArray class
          virtual ~ContainerBase();

          void setProperties(class ParticleDataContainerProperties &cp);
          bool propertiesSetCorrectly();

          inline const ParticleDataContainerProperties& prop() const
          { return containerProperties_; }

          inline void id(const char *_id)
          { _id = containerProperties_.id(); }

          inline bool matches_id(const char *_id)
          { return containerProperties_.matches_id(_id); }

          virtual bool isDoubleData() = 0;
          virtual bool isIntData() = 0;
          virtual bool isBoolData() = 0;

          virtual void addZero() = 0;
          virtual void addUninitialized(int n) = 0;
          virtual int size() = 0;
          virtual int nVec() = 0;
          virtual int lenVec() = 0;
          virtual void* begin_slow_dirty() = 0;

          virtual void copy(int from,int to) = 0;
          virtual void del(int n) = 0;
          virtual void delForward(int n,OperationProperties &op) = 0;
          virtual void delRestart(int n,OperationProperties &op) = 0;
          virtual void clearReverse(OperationProperties &op) = 0;

          virtual void read(OperationProperties &op,InputBase const* input_base) = 0;
          virtual void pull(OperationProperties &op,CouplingBase const* coupling_base) = 0;
          virtual void write(OperationProperties &op,InputBase const* input_base) = 0;
          virtual void push(OperationProperties &op,CouplingBase const* coupling_base) = 0;

          virtual bool setFromContainer(ContainerBase *cont) = 0;

          virtual void scale(double factor) = 0;
          virtual void move(double *dx) = 0;
          virtual void moveElement(int i,double *dx) = 0;
          virtual void rotate(double *dQ) = 0;

          virtual inline int lenVecUsed() { return -1; };
          virtual void setLenVecUsed(int x) = 0;

          // buffer functions for parallelization

          virtual int bufSize(OperationProperties &op) = 0;
          virtual int popFromBuffer(double *buf,OperationProperties &op) = 0;
          virtual int pushToBuffer(double *buf,OperationProperties &op) = 0;

          virtual int elemListBufSize(int n,OperationProperties &op) = 0;
          virtual int pushElemListToBuffer(int n, int *list, double *buf, OperationProperties &op) = 0;
          virtual int popElemListFromBuffer(int first, int n, double *buf, OperationProperties &op) = 0;
          virtual int pushElemListToBufferReverse(int first, int n, double *buf, OperationProperties &op) = 0;
          virtual int popElemListFromBufferReverse(int n, int *list, double *buf, OperationProperties &op) = 0;

          virtual int elemBufSize(OperationProperties &op) = 0;
          virtual int pushElemToBuffer(int n, double *buf,OperationProperties &op) = 0;
          virtual int popElemFromBuffer(double *buf,OperationProperties &op) = 0;

     protected:

          ContainerBase(ContainerBase const &orig);

          ParticleDataContainerProperties containerProperties_;

     private:

         ContainerBase();
  };

} /* PASCAL_NS */
#endif /* CONTAINERBASE_H_ */
