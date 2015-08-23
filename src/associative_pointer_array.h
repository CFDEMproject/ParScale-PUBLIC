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

#ifndef LMP_ASSOCIATIVE_POINTER_ARRAY_H
#define LMP_ASSOCIATIVE_POINTER_ARRAY_H

#include <string.h>
#include "memory.h"
#include "stdio.h"
#include "input_base.h"
#include "coupling_base.h"
#include "error_base.h"
#include "psctype.h"
#include "particle_data_container_properties.h"

namespace PASCAL_NS
{
  #define ID_LEN 100

template<typename T>
class AssociativePointerArray
{
      public:
        AssociativePointerArray();
        ~AssociativePointerArray();

        template <typename U>
        U* add(ParticleDataContainerProperties &cp);

        void remove(const char *_id);

        template <typename U>
        U* getPointerById(const char *_id);

        T* getBasePointerById(const char *_id);

        template <typename U>
        U* getPointerByIndex(int i);

        T* getBasePointerByIndex(int i);

        void grow(int to);

        int size();

        inline void copyElement(int from, int to);
        inline void addZeroElement();
        inline void deleteElement(int n);
        inline void deleteForwardElement(int n,OperationProperties &op);
        inline void deleteRestartElement(int n,OperationProperties &op);

        inline void clearReverse(OperationProperties &op);

        inline void read(OperationProperties &op,InputBase const* input_base);
        inline void pull(OperationProperties &op,CouplingBase const* coupling_base);
        inline void write(OperationProperties &op,InputBase const* input_base);
        inline void push(OperationProperties &op,CouplingBase const* coupling_base);

        inline void storeOrig(class AssociativePointerArray &orig);
        inline void storeOrig(const char *_id,class AssociativePointerArray &orig);
        inline bool reset(class AssociativePointerArray &orig);
        inline bool reset(const char *_id,class AssociativePointerArray &orig);

        void rotate(double *dQ);
        void move(double *delta);
        void moveElement(int i,double *delta);
        void scale(double factor);

        void sortPropsByExtMap(int *_id, int _nlocal, int &_len_id, int *_map,int _len_map,class ErrorBase const*err, bool verbose, int me);

        inline int bufSize(OperationProperties &op);
        inline int pushToBuffer(double *buf, OperationProperties &op);
        inline int popFromBuffer(double *buf, OperationProperties &op);

        inline int elemListBufSize(int n,OperationProperties &op);
        inline int pushElemListToBuffer(int n, int *list, double *buf, OperationProperties &op);
        inline int popElemListFromBuffer(int first, int n, double *buf, OperationProperties &op);
        inline int pushElemListToBufferReverse(int first, int n, double *buf, OperationProperties &op);
        inline int popElemListFromBufferReverse(int n, int *list, double *buf, OperationProperties &op);

        inline int elemBufSize(OperationProperties &op);
        inline int pushElemToBuffer(int n, double *buf, OperationProperties &op);
        inline int popElemFromBuffer(double *buf, OperationProperties &op);

        int idToIndex(const char *_id);
        void indexToId(int index, char *_id);

      private:

        T **content_;
        int numElem_, maxElem_;

        void growArrays();
};

  // *************************************
  #include "associative_pointer_array_I.h"
  // *************************************

} /* PASCAL_NS */
#endif /* ASSOCIATIVEPOINTERARRAY_H_ */
