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

#ifndef PSC_CUSTOM_VALUE_TRACKER_I_H
#define PSC_CUSTOM_VALUE_TRACKER_I_H

  /* ----------------------------------------------------------------------
   add property
  ------------------------------------------------------------------------- */

  template<typename T>
  T* CustomValueTracker::addElementProperty(const char *_id, const char* _comm, const char* _ref,
                                           const char *_restart,const char* _coupling, const char* _do_read,
                                           const char* _do_output, const char* _scope, int _scalePower)
  {
     ParticleDataContainerProperties cp(_id,_comm,_ref,_restart,_coupling,_do_read,_do_output,
                                         "element_property",_scope,_scalePower);

     // error if property exists already
     if(elementProperties_.getPointerById<T>(_id))
     {
         char *errmsg = new char[strlen(_id)+200];
         sprintf(errmsg,"Illegal command, features are incompatible - element property '%s' exists already",_id);
         printf("%s\n",errmsg);
         delete []errmsg;
     }

     // add property
     elementProperties_.add<T>(cp);

     // check if properties were set correctly
     // error here since ContainerBase not derived from Pointers
     if(!elementProperties_.getPointerById<T>(_id)->propertiesSetCorrectly())
     {
         char *errmsg = new char[strlen(_id)+200];
         sprintf(errmsg,"Illegal element property, comm or frame property not set correctly for property '%s'",_id);
         printf("%s\n",errmsg);
         delete []errmsg;
     }

     // return pointer
     return elementProperties_.getPointerById<T>(_id);
  }

  template<typename T>
  T* CustomValueTracker::addGlobalProperty(const char *_id, const char* _comm, const char* _ref,
                                           const char *_restart,const char* _coupling, const char* _do_read,
                                           const char* _do_output, const char* _scope,int _scalePower)
  {
     ParticleDataContainerProperties cp(_id,_comm,_ref,_restart,_coupling,_do_read,_do_output,
                                         "global_property",_scope,_scalePower);

     // error if property exists already
     if(globalProperties_.getPointerById<T>(_id))
     {
         char *errmsg = new char[strlen(_id)+200];
         sprintf(errmsg,"Illegal command, features are incompatible - global property '%s' already exists",_id);
         printf(errmsg);
         delete []errmsg;
     }

     // add property
     globalProperties_.add<T>(cp);
     globalProperties_orig_.add<T>(cp);

     // check if properties were set correctly
     // error here since ContainerBase not derived from Pointers
     if(!globalProperties_.getPointerById<T>(_id)->propertiesSetCorrectly())
     {
         char *errmsg = new char[strlen(_id)+200];
         sprintf(errmsg,"Illegal global property, comm or frame property not set correctly for property '%s'",_id);
         printf(errmsg);
         delete []errmsg;
     }

     // allocate memory
     //globalProperties_.getPointerById<T>(_id)->addUninitialized(capacityElement_);

     // return pointer
     return globalProperties_.getPointerById<T>(_id);
  }

  /* ----------------------------------------------------------------------
   mem management
  ------------------------------------------------------------------------- */

  void CustomValueTracker::grow_nbody(int _nbody,int _nbody_all)
  {
     grow(_nbody);
     set_n_body(_nbody,_nbody_all);
  }

  void CustomValueTracker::grow(int to)
  {
      elementProperties_.grow(to);
      capacityElement_ = to;
  }

  /* ----------------------------------------------------------------------
   get reference
  ------------------------------------------------------------------------- */

  template<typename T>
  T* CustomValueTracker::getElementProperty(const char *_id)
  {
     return elementProperties_.getPointerById<T>(_id);
  }

  inline ContainerBase* CustomValueTracker::getElementPropertyBase(const char *_id)
  {
     return elementProperties_.getBasePointerById(_id);
  }

  inline int CustomValueTracker::getElementPropertyIndex(const char *_id)
  {
     return elementProperties_.idToIndex(_id);
  }

  template<typename T>
  T* CustomValueTracker::getGlobalProperty(const char *_id)
  {
     return globalProperties_.getPointerById<T>(_id);
  }

  /* ----------------------------------------------------------------------
   set property
  ------------------------------------------------------------------------- */

  template<typename T, typename U>
  void CustomValueTracker::setElementProperty(const char *_id, U def)
  {
     elementProperties_.getPointerById<T>(_id)->set(def);
  }

  template<typename T, typename U>
  void CustomValueTracker::setGlobalProperty(const char *_id, U def)
  {
     //NP global properties are container classes with just one element contained
     if(globalProperties_.getPointerById<T>(_id)->size() == 0)
        globalProperties_.getPointerById<T>(_id)->addUninitialized(1);
     globalProperties_.getPointerById<T>(_id)->set(0,def);

     //NP global properties are container classes with just one element contained
     if(globalProperties_orig_.getPointerById<T>(_id)->size() == 0)
        globalProperties_orig_.getPointerById<T>(_id)->addUninitialized(1);
     globalProperties_orig_.getPointerById<T>(_id)->set(0,def);
  }

  /* ----------------------------------------------------------------------
   store global property orig - only needs to be done manually for
   special cases, eg moving mesh ref points
  ------------------------------------------------------------------------- */

  inline void CustomValueTracker::storeGlobalPropOrig(const char *_id)
  {
      globalProperties_.storeOrig(_id,globalProperties_orig_);
  }

  /* ----------------------------------------------------------------------
   reset global property to orig - only needs to be done manually for
   special cases, eg moving mesh ref points
  ------------------------------------------------------------------------- */

  inline void CustomValueTracker::resetGlobalPropToOrig(const char *_id)
  {
      globalProperties_.reset(_id,globalProperties_orig_);
  }

  /* ----------------------------------------------------------------------
   copy data from element from to element to
  ------------------------------------------------------------------------- */

  void CustomValueTracker::copyElement(int from, int to)
  {
      elementProperties_.copyElement(from,to);
  }

  /* ----------------------------------------------------------------------
   add an element and initialize its properties with 0
  ------------------------------------------------------------------------- */

  void CustomValueTracker::addZeroElement()
  {
      elementProperties_.addZeroElement();
  }

  /* ----------------------------------------------------------------------
   delete element i
  ------------------------------------------------------------------------- */

  void CustomValueTracker::deleteElement(int i)
  {
      elementProperties_.deleteElement(i);
      nbody_--;
  }

  /* ----------------------------------------------------------------------
   delete forward comm properties of element i
  ------------------------------------------------------------------------- */

  void CustomValueTracker::deleteForwardElement(int i,OperationProperties &op)
  {
      elementProperties_.deleteForwardElement(i,op);
  }

  /* ----------------------------------------------------------------------
   delete restart properties of element i
  ------------------------------------------------------------------------- */

  void CustomValueTracker::deleteRestartElement(int i,OperationProperties &op)
  {
      elementProperties_.deleteRestartElement(i,op);
  }

  /* ----------------------------------------------------------------------
   move element i
  ------------------------------------------------------------------------- */

  void CustomValueTracker::moveElement(int i, double *delta)
  {
      elementProperties_.moveElement(i,delta);
  }

  /* ----------------------------------------------------------------------
   push / pop for list of elements
  ------------------------------------------------------------------------- */

  int CustomValueTracker::elemListBufSize(int n,OperationProperties &op)
  {
    return elementProperties_.elemListBufSize(n,op);
  }

  int CustomValueTracker::pushElemListToBuffer(int n, int *list, double *buf, OperationProperties &op)
  {
    return elementProperties_.pushElemListToBuffer(n,list,buf,op);
  }

  int CustomValueTracker::popElemListFromBuffer(int first, int n, double *buf, OperationProperties &op)
  {
    return elementProperties_.popElemListFromBuffer(first,n,buf,op);
  }

  int CustomValueTracker::pushElemListToBufferReverse(int first, int n, double *buf, OperationProperties &op)
  {
    return elementProperties_.pushElemListToBufferReverse(first,n,buf,op);
  }

  int CustomValueTracker::popElemListFromBufferReverse(int n, int *list, double *buf, OperationProperties &op)
  {
    return elementProperties_.popElemListFromBufferReverse(n,list,buf,op);
  }

  /* ----------------------------------------------------------------------
   push / pop for element i
  ------------------------------------------------------------------------- */

  int CustomValueTracker::elemBufSize(OperationProperties &op)
  {
    /*NL*///char _id[300];
    /*NL*///for(int i=0;i<elementProperties_.numElem_;i++)
    /*NL*///{
    /*NL*///  elementProperties_.getBasePointerByIndex(i)->id(_id);
    /*NL*///  fprintf(this->screen,"prop %s size %d\n",_id,elementProperties_.getBasePointerByIndex(i)->elemBufSize(op));
    /*NL*///}
    return elementProperties_.elemBufSize(op);
  }

  int CustomValueTracker::pushElemToBuffer(int i, double *buf, OperationProperties &op)
  {
    return elementProperties_.pushElemToBuffer(i,buf,op);
  }

  int CustomValueTracker::popElemFromBuffer(double *buf, OperationProperties &op,bool inc)
  {
    int size = elementProperties_.popElemFromBuffer(buf,op);
    if(inc)
      nbody_++;
    return size;
  }

  /* ----------------------------------------------------------------------
   push / pop for global properties
  ------------------------------------------------------------------------- */

  int CustomValueTracker::globalPropsBufSize(OperationProperties &op)
  {
    int n = 0;
    n += globalProperties_.bufSize(op);
    n += globalProperties_orig_.bufSize(op);
    return n;
  }

  int CustomValueTracker::pushGlobalPropsToBuffer(double *buf, OperationProperties &op)
  {
    int n = 0;
    n += globalProperties_.pushToBuffer(&(buf[n]),op);
    n += globalProperties_orig_.pushToBuffer(&(buf[n]),op);
    return n;
  }

  int CustomValueTracker::popGlobalPropsFromBuffer(double *buf, OperationProperties &op)
  {
    int n = 0;
    n += globalProperties_.popFromBuffer(&(buf[n]),op);
    n += globalProperties_orig_.popFromBuffer(&(buf[n]),op);
    return n;
  }

#endif
