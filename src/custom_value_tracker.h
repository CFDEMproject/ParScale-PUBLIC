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

#ifndef PSC_CUSTOM_VALUE_TRACKER_H
#define PSC_CUSTOM_VALUE_TRACKER_H

#include "associative_pointer_array.h"
#include "pascal_base.h"
#include "container.h"

namespace PASCAL_NS
{
  class CustomValueTracker : protected ParScaleBase
  {
      public:
        CustomValueTracker(ParScale *ptr, ParticleData *pdata);
        ~CustomValueTracker();

        // per-element properties

        template<typename T>
        T* addElementProperty(const char *_id, const char* _comm, const char* _ref, const char *_restart,
                              const char* _coupling ,const char* _do_read, const char* _do_output,
                              const char* _scope = 0,int _scalePower = 1);

        template<typename T>
        T* getElementProperty(const char *_id);

        inline ContainerBase* getElementPropertyBase(const char *_id);

        inline int getElementPropertyIndex(const char *_id);

        template<typename T, typename U>
        void setElementProperty(const char *_id, U def);

        void removeElementProperty(const char *_id);

        // global (e.g. mesh) properties

        template<typename T>
        T* addGlobalProperty(const char *_id, const char* _comm, const char* _ref, const char *_restart,
                             const char* _coupling ,const char* _do_read, const char* _do_output,
                             const char* _scope = 0,int _scalePower = 1);

        template<typename T>
        T* getGlobalProperty(const char *_id);

        template<typename T, typename U>
        void setGlobalProperty(const char *_id, U def);

        void removeGlobalProperty(const char *_id);

        // initial allocation
        void set_n_body(int _nlocal,int _nglobal);
        void allocate();

        void read(OperationProperties &op);
        void pull(OperationProperties &op);
        void write(OperationProperties &op);
        void push(OperationProperties &op);

        // particle # management

        inline int nbody()     const {return nbody_;}
        inline int nbody_all() const {return nbody_all_;}
        void recalc_nbody_all(bool errflag);

        // operation with
        // per-element properties

        inline void copyElement(int from, int to);
        inline void addZeroElement();
        inline void deleteElement(int i);
        inline void deleteForwardElement(int i,OperationProperties &op);
        inline void deleteRestartElement(int i,OperationProperties &op);
        void clearReverse(OperationProperties &op);

        void storeOrig();
        void resetToOrig();

        inline void storeGlobalPropOrig(const char *_id);
        inline void resetGlobalPropToOrig(const char *_id);

        inline void moveElement(int i, double *delta);
        void move(double *vecTotal, double *vecIncremental);
        void move(double *vecIncremental);
        void rotate(double *totalQ, double *dQ);
        void rotate(double *dQ);
        void scale(double factor);

        void sortPropsByExtMap(int *_id, int _nlocal, int &_len_id,int *_map,int _len_map, bool verbose, int me);

        // buffer operations

        inline int elemListBufSize(int n,OperationProperties &op);
        inline int pushElemListToBuffer(int n, int *list, double *buf, OperationProperties &op);
        inline int popElemListFromBuffer(int first, int n, double *buf, OperationProperties &op);
        inline int pushElemListToBufferReverse(int first, int n, double *buf, OperationProperties &op);
        inline int popElemListFromBufferReverse(int n, int *list, double *buf, OperationProperties &op);

        inline int elemBufSize(OperationProperties &op);
        inline int pushElemToBuffer(int n, double *buf, OperationProperties &op);
        inline int popElemFromBuffer(double *buf, OperationProperties &op,bool inc = false);

        inline int globalPropsBufSize(OperationProperties &op);
        inline int pushGlobalPropsToBuffer(double *buf, OperationProperties &op);
        inline int popGlobalPropsFromBuffer(double *buf, OperationProperties &op);

        // mem managenement

        int getCapacity();
        inline void grow_nbody(int _nbody,int _nbody_all);
        inline void grow(int to);

      private:

        // # of local (owned) particles, # particles on all procs
        int nbody_, nbody_all_;

        class ParticleData &owner_;

        int capacityElement_; //NP only up-to-date at start, push/pop may change
        class AssociativePointerArray<ContainerBase> elementProperties_;
        class AssociativePointerArray<ContainerBase> globalProperties_;
        class AssociativePointerArray<ContainerBase> globalProperties_orig_;
  };

  // *************************************
  #include "custom_value_tracker_I.h"
  // *************************************

} /* PASCAL_NS */
#endif
