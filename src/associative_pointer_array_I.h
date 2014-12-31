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


#ifndef LMP_ASSOCIATIVE_POINTER_ARRAY_I_H
#define LMP_ASSOCIATIVE_POINTER_ARRAY_I_H

  /* ----------------------------------------------------------------------
   constructors, destructor
  ------------------------------------------------------------------------- */

  template<typename T>
  AssociativePointerArray<T>::AssociativePointerArray()
   : content_(0), numElem_(0), maxElem_(1)
  {
    content_ = new T*[1];
    content_[0] = 0;
  }

  template<typename T>
  AssociativePointerArray<T>::~AssociativePointerArray()
  {
    for(int i = 0; i < numElem_; i++)
      delete content_[i];

    delete[] content_;
  }

  /* ----------------------------------------------------------------------
   add for per-element and per-mesh properties
  ------------------------------------------------------------------------- */

  template<typename T> template<typename U>
  U* AssociativePointerArray<T>::add(ParticleDataContainerProperties &cp)
  {
    if(numElem_ == maxElem_)
      growArrays();

    content_[numElem_] = static_cast<T*>(new U(cp));
    numElem_++;
    /*NP
    printf("numElem_ %d\n",numElem_);
    for(int i=0;i<numElem_;i++)
      printf("  %d %s\n",i,content_[i]->id_);
    */
    return static_cast<U*>(content_[numElem_-1]);
  }

  /* ----------------------------------------------------------------------
   delete properties
  ------------------------------------------------------------------------- */

  template<typename T>
  void AssociativePointerArray<T>::remove(const char *_id)
  {
    int index = idToIndex(_id);
    if(index == -1) return;

    numElem_--;

    delete content_[index];

    if(numElem_ > 0)
        content_[index] = content_[numElem_];
  }

  /* ----------------------------------------------------------------------
   get pointer to property
  ------------------------------------------------------------------------- */

  template<typename T> template<typename U>
  U* AssociativePointerArray<T>::getPointerById(const char *_id)
  {
    int ind = idToIndex(_id);
    if(ind != -1)
      return getPointerByIndex<U>(ind);
    else
      return 0;
  }

  template<typename T>
  T* AssociativePointerArray<T>::getBasePointerById(const char *_id)
  {
    int ind = idToIndex(_id);
    if(ind != -1)
      return getBasePointerByIndex(ind);
    else
      return 0;
  }

  template<typename T> template<typename U>
  U* AssociativePointerArray<T>::getPointerByIndex(int i)
  {
    if(i >= size() || i < 0) return 0;
    else return dynamic_cast<U*>(content_[i]);
  }

  template<typename T>
  T* AssociativePointerArray<T>::getBasePointerByIndex(int i)
  {
    if(i >= size() || i < 0) return 0;
    else return content_[i];
  }

  /* ----------------------------------------------------------------------
   memory management
  ------------------------------------------------------------------------- */

  template<typename T>
  void AssociativePointerArray<T>::growArrays()
  {

    // for(int i=0;i<numElem_+1;i++)
    //  printf("%d %s %d\n",i,id_[i], strcmp(id_[i],"v"));

    T ** tmp = new T*[maxElem_];

    for(int i = 0; i < maxElem_; i++)
        tmp[i] = content_[i];

    delete[] content_;

    maxElem_++;
    content_ = new T*[maxElem_];

    for(int i = 0; i < numElem_; i++)
        content_[i] = tmp[i];

    delete[] tmp;

    //for(int i=0;i<numElem_+1;i++)
    //  printf("%d %s %d\n",i,id_[i], strcmp(id_[i],"v"));
  }

  template<typename T>
  void AssociativePointerArray<T>::grow(int to)
   {
      int by;
      for(int i = 0; i < maxElem_; i++)
      {
          by = to - getBasePointerByIndex(i)->size();
          if(by > 0)
            getBasePointerByIndex(i)->addUninitialized(by);
      }
  }

  template<typename T>
  int AssociativePointerArray<T>::size()
  {
    return numElem_;
  }

  /* ----------------------------------------------------------------------
   copy data from element from to element to
  ------------------------------------------------------------------------- */

  template<typename T>
  void AssociativePointerArray<T>::copyElement(int from, int to)
  {
      for(int i=0;i<numElem_;i++)
        content_[i]->copy(from,to);
  }

  /* ----------------------------------------------------------------------
   add an element and initialize its properties with 0
  ------------------------------------------------------------------------- */

  template<typename T>
  void AssociativePointerArray<T>::addZeroElement()
  {
      for(int i=0;i<numElem_;i++)
        content_[i]->addZero();
  }

  /* ----------------------------------------------------------------------
   delete element n
  ------------------------------------------------------------------------- */

  template<typename T>
  void AssociativePointerArray<T>::deleteElement(int n)
  {
      for(int i=0;i<numElem_;i++)
        content_[i]->del(n);
  }

  /* ----------------------------------------------------------------------
   delete forward properties of element i
  ------------------------------------------------------------------------- */

  template<typename T>
  void AssociativePointerArray<T>::deleteForwardElement(int n,OperationProperties &op)
  {
      for(int i=0;i<numElem_;i++)
        content_[i]->delForward(n,op);
  }

  /* ----------------------------------------------------------------------
   delete restart properties of element i
  ------------------------------------------------------------------------- */

  template<typename T>
  void AssociativePointerArray<T>::deleteRestartElement(int n,OperationProperties &op)
  {
      for(int i=0;i<numElem_;i++)
        content_[i]->delRestart(n,op);
  }

  /* ----------------------------------------------------------------------
   clear reverse properties, i.e. reset all of them to 0
  ------------------------------------------------------------------------- */

  template<typename T>
  void AssociativePointerArray<T>::clearReverse(OperationProperties &op)
  {
      for(int i=0;i<numElem_;i++)
        content_[i]->clearReverse(op);
  }

  /* ----------------------------------------------------------------------
   read input from JSON file
  ------------------------------------------------------------------------- */

  template<typename T>
  void AssociativePointerArray<T>::read(OperationProperties &op,
                                        InputBase const* input_base)
  {
      for(int i=0;i<numElem_;i++)
        content_[i]->read(op,input_base);
  }

  /* ----------------------------------------------------------------------
   pull from LIGGGHTS
  ------------------------------------------------------------------------- */

  template<typename T>
  void AssociativePointerArray<T>::pull(OperationProperties &op,
                                        CouplingBase const* coupling_base)
  {
      for(int i=0;i<numElem_;i++)
        content_[i]->pull(op,coupling_base);
  }

  /* ----------------------------------------------------------------------
   write input to JSON file
  ------------------------------------------------------------------------- */

  template<typename T>
  void AssociativePointerArray<T>::write(OperationProperties &op,
                                        InputBase const* input_base)
  {
      for(int i=0;i<numElem_;i++)
        content_[i]->write(op,input_base);
  }

  /* ----------------------------------------------------------------------
   push to LIGGGHTS
  ------------------------------------------------------------------------- */

  template<typename T>
  void AssociativePointerArray<T>::push(OperationProperties &op,
                                        CouplingBase const* coupling_base)
  {
      for(int i=0;i<numElem_;i++)
        content_[i]->push(op,coupling_base);
  }

  /* ----------------------------------------------------------------------
   id 2 index
  ------------------------------------------------------------------------- */

  template<typename T>
  int AssociativePointerArray<T>::idToIndex(const char *_id)
  {
    for(int i=0;i<numElem_;i++)
      if(content_[i]->matches_id(_id))
        return i;
    return -1;
  }

  template<typename T>
  void AssociativePointerArray<T>::indexToId(int index, char *_id)
  {
      content_[index]->id(_id);
  }

  /* ----------------------------------------------------------------------
   store original value for reset
  ------------------------------------------------------------------------- */

  template<typename T>
  void AssociativePointerArray<T>::storeOrig(AssociativePointerArray &orig)
  {
      for(int i = 0; i < numElem_; i++)
          orig.content_[i]->setFromContainer(content_[i]);
  }

  template<typename T>
  void AssociativePointerArray<T>::storeOrig(const char *_id, AssociativePointerArray &orig)
  {
      /*NL*/ //printf("storeOrig called for %s \n",_id);
      for(int i = 0; i < numElem_; i++)
          if(content_[i]->matches_id(_id))
            orig.content_[i]->setFromContainer(content_[i]);
  }

  /* ----------------------------------------------------------------------
   reset to original value
  ------------------------------------------------------------------------- */

  template<typename T>
  bool AssociativePointerArray<T>::reset(AssociativePointerArray &orig)
  {
      /*NL*/ //printf("numElem_ %d\n",numElem_);

      for(int i = 0; i < numElem_; i++)
          content_[i]->setFromContainer(orig.content_[i]);

      return true;
  }

  template<typename T>
  bool AssociativePointerArray<T>::reset(const char *_id, AssociativePointerArray &orig)
  {
      /*NL*/ //printf("numElem_ %d\n",numElem_);

      for(int i = 0; i < numElem_; i++)
          if(content_[i]->matches_id(_id))
            content_[i]->setFromContainer(orig.content_[i]);

      return true;
  }

  /* ----------------------------------------------------------------------
   move, rotate scale all properties
  ------------------------------------------------------------------------- */

  template<typename T>
  void AssociativePointerArray<T>::rotate(double *dQ)
  {
      for(int i = 0; i < numElem_; i++)
        content_[i]->rotate(dQ);
  }

  template<typename T>
  void AssociativePointerArray<T>::scale(double factor)
  {
      for(int i = 0; i < numElem_;i++)
        content_[i]->scale(factor);
  }

  template<typename T>
  void AssociativePointerArray<T>::move(double *delta)
  {
      for(int i = 0; i < numElem_;i++)
        content_[i]->move(delta);
  }

  template<typename T>
  void AssociativePointerArray<T>::moveElement(int n,double *delta)
  {
      for(int i = 0; i < numElem_;i++)
        content_[i]->moveElement(n,delta);
  }

  /* ----------------------------------------------------------------------
   sort list
  ------------------------------------------------------------------------- */

  template<typename T>
  void AssociativePointerArray<T>::sortPropsByExtMap(int *_id,int _len_id,  //Pascal's LOCAL ->GLOBAL ids to be set
                                                     int *_map,int _len_map,//external GLOBAL->LOCAL  map with data for setting id
                                                     ErrorBase const* err)
  {

    printf("sortPropsByExtMap is sorting your data... \n");

    bool verbose_=false; 

    // check consistency of number of elements in containers
    // new nbody might have been pulled before already
    int nbody_old = content_[0]->size();
    for(int icontainer=1;icontainer<numElem_;icontainer++)
        if(nbody_old != content_[icontainer]->size() || _len_id != nbody_old)
            err->throw_error_all(FLERR,"internal error; assertion failed because containers have different size");

    // ensure there is enough space in the containers for swap, so add 1
    grow(nbody_old+1);

    int  ibody = 0,ibody_new;
    int  globalCurrID;              //current GLOBAL ID of the atom
    int  globalCurrMem;             //position in the global MEM
    
    
    //Check if all bodies have a valid global ID
    for(ibody=0; ibody<nbody_old; ibody++)
    {
        if(_id[ibody]==-1) //global ID was not set before or is invalid: Assume consecutive ordering
        {
            _id[ibody] = ibody+1; //ids start with 1!
            if(verbose_)
            {
                printf("ibody: %d, resetting global ID to %d\n",
                        ibody, _id[ibody]);
            }
        }            
    }
    

    // go through all bodies/elements available in ParScale
    ibody=0;
    while(ibody < nbody_old)    //loop all GLOBAL values and get values of the map
    {
        //A - Get the current global ID 
        globalCurrID = _id[ibody];
        globalCurrMem = globalCurrID-1;

        //B - Get the data location from the global map
        int dataLocation = _map[globalCurrMem]; //this is the LOCAL data location
        
        if(verbose_)
        {
            printf("sortPropsByExtMap::nbody_old: %d, ibody: %d, globalCurrID: %d, dataLocation: %d\n",
                    nbody_old, ibody, globalCurrID, dataLocation);
        }

        //Nothing to do if data location corresponds to index
        if(dataLocation == ibody)
        {
            ibody++;
            continue;
        }
        // body not in tag list, delete from ParScale
        else if(dataLocation < 0) //local index not found --> particle not on this machine
        {
            if(verbose_)
                printf("sortPropsByExtMap::data location not found. Will delete the particle...\n");

            // delete data of ibody
            for(int icontainer=0;icontainer<numElem_;icontainer++)
                content_[icontainer]->del(ibody);

            // copy GLOBAL ID from last to current (same way as done by del() call)
            _id[ibody] = _id[nbody_old-1];

            nbody_old--;
        }
        // swap
        else
        {
            if(verbose_)
            {
              printf("sortPropsByExtMap::will swap data for global id: %d. \n",
                     globalCurrID);
              printf("sortPropsByExtMap::spare new location at %d by putting data to end (= %d) \n",
                     dataLocation, nbody_old);
              printf("sortPropsByExtMap::copy data from %d to %d \n",
                     ibody, dataLocation);
              printf("sortPropsByExtMap::copy spare from %d to %d \n",
                     nbody_old, ibody);
            }

            // copy data at inew to the end of containers
            copyElement(dataLocation,nbody_old);

            // copy data at i to inew
            copyElement(ibody,dataLocation);

            // copy from spare place to i
            copyElement(nbody_old,ibody);

            // update the IDs
            int id_tmp = _id[dataLocation];  
            _id[dataLocation] = globalCurrID; //set GLOBAL ID at new data location
            _id[ibody] = id_tmp;
            
            if(verbose_)
              printf("new GLOBAL ID at ibody: %d, and at dataLocation: %d \n",
                     _id[ibody], _id[dataLocation]);
        }
    }
    
    nbody_old = content_[0]->size();
    deleteElement(nbody_old-1);
    if(verbose_)
    {
         printf("Deleting last element in array with size: %d, size is now: %d \n",
                 nbody_old, content_[0]->size());
    }
    
    for(ibody=0; ibody<(content_[0]->size()-1); ibody++)
    {
       if(verbose_)
       {
                printf("final check: ibody: %d, global ID: %d, dataLocation: %d\n",
                        ibody, _id[ibody], _map[_id[ibody]]);
       }
    }
  }

  /* ----------------------------------------------------------------------
   buf size, push, pop for all elements
  ------------------------------------------------------------------------- */

  template<typename T>
  int AssociativePointerArray<T>::bufSize(OperationProperties &op)
  {
    int buf_size = 0;
    for(int i=0;i<numElem_;i++)
      buf_size += getBasePointerByIndex(i)->bufSize(op);
    return buf_size;
  }

  template<typename T>
  int AssociativePointerArray<T>::pushToBuffer(double *buf, OperationProperties &op)
  {
    int nsend = 0;
    for(int i=0;i<numElem_;i++)
      nsend += getBasePointerByIndex(i)->pushToBuffer(&(buf[nsend]),op);
    return nsend;
  }

  template<typename T>
  int AssociativePointerArray<T>::popFromBuffer(double *buf, OperationProperties &op)
  {
    int nrecv = 0;
    for(int i=0;i<numElem_;i++)
      nrecv += getBasePointerByIndex(i)->popFromBuffer(&(buf[nrecv]),op);
    return nrecv;
  }

  /* ----------------------------------------------------------------------
   buf size, push, pop for list of elements
  ------------------------------------------------------------------------- */

  template<typename T>
  int AssociativePointerArray<T>::elemListBufSize(int n,OperationProperties &op)
  {
    int buf_size = 0;
    for(int i=0;i<numElem_;i++)
      buf_size += getBasePointerByIndex(i)->elemListBufSize(n,op);
    return buf_size;
  }

  template<typename T>
  int AssociativePointerArray<T>::pushElemListToBuffer(int n, int *list, double *buf, OperationProperties &op)
  {
    int nsend = 0;
    for(int i=0;i<numElem_;i++)
      nsend += getBasePointerByIndex(i)->pushElemListToBuffer(n,list,&buf[nsend],op);
    return nsend;
  }

  template<typename T>
  int AssociativePointerArray<T>::popElemListFromBuffer(int first, int n, double *buf, OperationProperties &op)
  {
    int nrecv = 0;
    for(int i=0;i<numElem_;i++)
      nrecv += getBasePointerByIndex(i)->popElemListFromBuffer(first,n,&buf[nrecv],op);
    return nrecv;
  }

  template<typename T>
  int AssociativePointerArray<T>::pushElemListToBufferReverse(int first, int n, double *buf, OperationProperties &op)
  {
    int nrecv = 0;
    for(int i=0;i<numElem_;i++)
      nrecv += getBasePointerByIndex(i)->pushElemListToBufferReverse(first,n,&buf[nrecv],op);
    return nrecv;
  }

  template<typename T>
  int AssociativePointerArray<T>::popElemListFromBufferReverse(int n, int *list, double *buf, OperationProperties &op)
  {
    int nsend = 0;
    for(int i=0;i<numElem_;i++)
      nsend += getBasePointerByIndex(i)->popElemListFromBufferReverse(n,list,&buf[nsend],op);
    return nsend;
  }

  /* ----------------------------------------------------------------------
   buf size, push, pop for single element
  ------------------------------------------------------------------------- */

  template<typename T>
  int AssociativePointerArray<T>::elemBufSize(OperationProperties &op)
  {
    int buf_size = 0;
    for(int i=0;i<numElem_;i++)
      buf_size += getBasePointerByIndex(i)->elemBufSize(op);
    return buf_size;
  }

  template<typename T>
  int AssociativePointerArray<T>::pushElemToBuffer(int n, double *buf, OperationProperties &op)
  {
    int nsend = 0;
    for(int i=0;i<numElem_;i++)
      nsend += getBasePointerByIndex(i)->pushElemToBuffer(n,&buf[nsend],op);
    return nsend;
  }

  template<typename T>
  int AssociativePointerArray<T>::popElemFromBuffer(double *buf, OperationProperties &op)
  {
    int nrecv = 0;
    for(int i=0;i<numElem_;i++)
      nrecv += getBasePointerByIndex(i)->popElemFromBuffer(&buf[nrecv],op);
    return nrecv;
  }

#endif

