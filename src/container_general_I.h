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

#ifndef LMP_CONTAINER_GENERAL_I_H
#define LMP_CONTAINER_GENERAL_I_H

  /* ----------------------------------------------------------------------
   constructors
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  ContainerGeneral<T,NUM_VEC,LEN_VEC>::ContainerGeneral(ParticleDataContainerProperties &cp)
  : ContainerBase(cp),
    numElem_(0),
    maxElem_(GROW),
    lenVecUsed_(LEN_VEC),
    isFilled_(false),
    isDumped_(false),
    isFull_(false),
    verbose_(true)
  {
          create<T>(arr_,GROW,NUM_VEC,LEN_VEC);
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  ContainerGeneral<T,NUM_VEC,LEN_VEC>::ContainerGeneral(ContainerGeneral<T,NUM_VEC,LEN_VEC> const &orig)
  : ContainerBase(orig),
    numElem_(orig.numElem_),
    maxElem_(orig.numElem_)
  {
          create<T>(arr_,maxElem_,NUM_VEC,LEN_VEC);
          for(int i=0;i<maxElem_;i++)
                  for(int ii=0;ii<NUM_VEC;ii++)
                          for(int jj=0;jj<LEN_VEC;jj++)
                                  arr_[i][ii][jj] = orig.arr_[i][ii][jj];
  }

  /* ----------------------------------------------------------------------
   destructor
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  ContainerGeneral<T,NUM_VEC,LEN_VEC>::~ContainerGeneral()
  {
          destroy<T>(arr_);
  }

  /* ----------------------------------------------------------------------
   check if data is of type double
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  bool ContainerGeneral<T,NUM_VEC,LEN_VEC>::isDoubleData()
  {
      // partial templatization does not work
      // std::is_same<T,double>::value is from C++11
      // this is work-around

      if(sizeof(T) == sizeof(double))
        return true;
      else
        return false;
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  bool ContainerGeneral<T,NUM_VEC,LEN_VEC>::isIntData()
  {
      // partial templatization does not work
      // std::is_same<T,double>::value is from C++11
      // this is work-around

      if(sizeof(T) == sizeof(int))
        return true;
      else
        return false;
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  bool ContainerGeneral<T,NUM_VEC,LEN_VEC>::isBoolData()
  {
      // partial templatization does not work
      // std::is_same<T,double>::value is from C++11
      // this is work-around

      if(sizeof(T) == sizeof(bool))
        return true;
      else
        return false;
  }

  /* ----------------------------------------------------------------------
   add element(s)
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  void ContainerGeneral<T,NUM_VEC,LEN_VEC>::add(T** elem)
  {
          if(numElem_ == maxElem_)
          {
                  grow<T>(arr_,maxElem_+GROW,NUM_VEC,LEN_VEC);
                  maxElem_ += GROW;
          }
          for(int i=0;i<NUM_VEC;i++)
                  for(int j=0;j<LEN_VEC;j++)
                          arr_[numElem_][i][j] = elem[i][j];
          numElem_++;
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  void ContainerGeneral<T,NUM_VEC,LEN_VEC>::addZero()
  {
          if(numElem_ == maxElem_)
          {
                  grow<T>(arr_,maxElem_+GROW,NUM_VEC,LEN_VEC);
                  maxElem_ += GROW;
          }
          for(int i=0;i<NUM_VEC;i++)
                  for(int j=0;j<LEN_VEC;j++)
                          arr_[numElem_][i][j] = static_cast<T>(0);
          numElem_++;
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  void ContainerGeneral<T,NUM_VEC,LEN_VEC>::addUninitialized(int n)
  {
        numElem_ += n;
        if(numElem_ >= maxElem_)
        {
            grow(arr_,numElem_+GROW,NUM_VEC,LEN_VEC);
            maxElem_ = numElem_ + GROW;
        }
  }

  /* ----------------------------------------------------------------------
   delete an element
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  void ContainerGeneral<T,NUM_VEC,LEN_VEC>::del(int n)
  {
          numElem_--;
          if(numElem_ == n) return;
          for(int i=0;i<NUM_VEC;i++)
                  for(int j=0;j<LEN_VEC;j++)
                          arr_[n][i][j] = arr_[numElem_][i][j];
  }


  /* ----------------------------------------------------------------------
   copy element data
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  void ContainerGeneral<T,NUM_VEC,LEN_VEC>::copy(int from,int to)
  {
          for(int i=0;i<NUM_VEC;i++)
                  for(int j=0;j<LEN_VEC;j++)
                          arr_[to][i][j] = arr_[from][i][j];
  }

  /* ----------------------------------------------------------------------
   delete an element
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  void ContainerGeneral<T,NUM_VEC,LEN_VEC>::delForward(int n,OperationProperties &op)
  {
          op.set_operation(OPERATION_COMM_FORWARD);
          // do only delete property if it is a forward comm property
          if(!this->containerProperties_.decidePackUnpackOperation(op))
            return;

          numElem_--;
          if(numElem_ == n) return;
          for(int i=0;i<NUM_VEC;i++)
                  for(int j=0;j<LEN_VEC;j++)
                          arr_[n][i][j] = arr_[numElem_][i][j];
  }

  /* ----------------------------------------------------------------------
   clear reverse properties, i.e. reset all of them to 0
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  void ContainerGeneral<T,NUM_VEC,LEN_VEC>::clearReverse(OperationProperties &op)
  {
      op.set_operation(OPERATION_COMM_REVERSE);
      // do only reset property if it is a reverse comm property
      if(!this->containerProperties_.decidePackUnpackOperation(op))
        return;

      int len = size();
      for(int i = 0; i < len; i++)
            for(int j = 0; j < NUM_VEC; j++)
                for(int k = 0; k < LEN_VEC; k++)
                    arr_[i][j][k] = 0.;
  }

  /* ----------------------------------------------------------------------
    read properties from JSON file
    ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  void ContainerGeneral<T,NUM_VEC,LEN_VEC>::read(OperationProperties &op,
                                                 InputBase const* input_base)
  {
      op.set_operation(OPERATION_READ);

      // do only read property if (a) it is a must read property, (b) container is empty
      // this statement also ensures that properties that are pulled via pull() upon init()
      // do not get overwritten by read()
      // also ensures that container are over-written by reads of successive "run" commands
    
      if(verbose_)
        printf("...reading container with id '%s': isFilled %s, decidePackUnpack %s\n",
                this->containerProperties_.id(),
                isFilled_?"true":"false",
                this->containerProperties_.decidePackUnpackOperation(op)?"true":"false"
              );

      if(isFilled_ || !this->containerProperties_.decidePackUnpackOperation(op))
            return;

      //do not read containers that require pulling and are not read
      if(containerProperties_.needsPull() && !containerProperties_.needsRead() )
      {
          if(verbose_)
            printf( "Container with id '%s': skipping read because container content will be pulled and does not require read.\n",
                    containerProperties_.id()
                  );
          isFilled_ = true;
          return;
      }

      input_base->fill_container_from_json(static_cast<ContainerBase&>(*this));

      // flag that container is filled
      isFilled_ = true;
     
      // flag that container filled from file, i.e. has to be bcasted before using
      containerProperties_.needBCast();
  }

  /* ----------------------------------------------------------------------
    pull from LIGGGHTS
    ------------------------------------------------------------------------- */
  template<typename T, int NUM_VEC, int LEN_VEC>
  void ContainerGeneral<T,NUM_VEC,LEN_VEC>::pull(OperationProperties &op,
                                                 CouplingBase const* coupling_base)
  {
      op.set_operation(OPERATION_PULL);

      //only perform opertation for a pull property
      if(!this->containerProperties_.decidePackUnpackOperation(op) && coupling_base->forceFill() )
            return;

      // set to true only if actually filled, i.e. if coupling really pulled it
      // coupling none will return false
      isFilled_ = coupling_base->fill_container_from_coupling(static_cast<ContainerBase&>(*this));

      // flag that need no bcast
      containerProperties_.needNoBCast();

//      if( !isFilled_ ) //TODO: print WARNING just in case external_code_in_control() s
//        printf("\nWARNING: ContainerGeneral::pull could not fill from coupling.");
      
  }

  /* ----------------------------------------------------------------------
    write containers
    ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  void ContainerGeneral<T,NUM_VEC,LEN_VEC>::write(OperationProperties &op,
                                                 InputBase const* input_base)
  {
      op.set_operation(OPERATION_OUTPUT);

      input_base->write_containersJSON(static_cast<ContainerBase&>(*this));

      input_base->write_containersHDF5(static_cast<ContainerBase&>(*this));
  }

  /* ----------------------------------------------------------------------
    push to LIGGGHTS
    ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  void ContainerGeneral<T,NUM_VEC,LEN_VEC>::push(OperationProperties &op,
                                                 CouplingBase const* coupling_base)
  {
      op.set_operation(OPERATION_PUSH);

      //Check when not to dump!
      if(!containerProperties_.needsPush())
        return;
//      if(!this->containerProperties_.decidePackUnpackOperation(op))
//            return;
      isDumped_ = coupling_base->dump_container_to_coupling(static_cast<ContainerBase&>(*this));
        
      return;

  }

  /* ----------------------------------------------------------------------
   delete an element
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  void ContainerGeneral<T,NUM_VEC,LEN_VEC>::delRestart(int n,OperationProperties &op)
  {
          op.set_operation(OPERATION_RESTART);
          // do only delete property if it is a restart property
          if(!this->containerProperties_.decidePackUnpackOperation(op))
            return;

          /*NL*/ //printf("del restart for %s, numElem_ %d, n %d\n",this->id_,numElem_,n);

          numElem_--;
          if(numElem_ == n) return;
          for(int i=0;i<NUM_VEC;i++)
                  for(int j=0;j<LEN_VEC;j++)
                          arr_[n][i][j] = arr_[numElem_][i][j];
  }

  /* ----------------------------------------------------------------------
   get an element
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  void ContainerGeneral<T,NUM_VEC,LEN_VEC>::get(int n, T** elem)
  {
          for(int i=0;i<NUM_VEC;i++)
                  for(int j=0;j<LEN_VEC;j++)
                          elem[i][j] = arr_[n][i][j];
  }

  /* ----------------------------------------------------------------------
   operator()
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  T**& ContainerGeneral<T,NUM_VEC,LEN_VEC>::operator() (int n)
  {
          return arr_[n];
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  T** const& ContainerGeneral<T,NUM_VEC,LEN_VEC>::operator() (int n) const
  {
          return arr_[n];
  }

  /* ----------------------------------------------------------------------
   set all data by copy from other container
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  bool ContainerGeneral<T,NUM_VEC,LEN_VEC>::setFromContainer(ContainerBase *cont)
  {
      ContainerGeneral<T,NUM_VEC,LEN_VEC> *gcont = static_cast<ContainerGeneral<T,NUM_VEC,LEN_VEC>* >(cont);

      /*NL*/// printf("container %s sizes %d %d nvec %d %d lenvec %d %d\n",
      /*NL*///        id_,size(),gcont->size(),nVec(),gcont->nVec(),lenVec(),gcont->lenVec());

      /*NL*/// printf("TRYING set container %s, sizes %d %d \n",id_,size(), gcont->size());

      //NP only copy if identical
      if(size() != gcont->size() || nVec() != gcont->nVec() || lenVec() != gcont->lenVec())
        return false;

      /*NL*/// printf("SETTING container %s\n",id_);

      int len = size();
      for(int n = 0; n < len; n++)
          for(int i=0;i<NUM_VEC;i++)
                  for(int j=0;j<LEN_VEC;j++)
                  {
                          arr_[n][i][j] = gcont->arr_[n][i][j];
                          /*NL*/ //printf("   gcont->arr_[n][i][j] %f \n",static_cast<double>(gcont->arr_[n][i][j]));
                  }

      return true;
  }

  /* ---------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  void ContainerGeneral<T,NUM_VEC,LEN_VEC>::set(int n, T** elem)
  {
          for(int i = 0; i < NUM_VEC; i++)
                          for(int j = 0; j < LEN_VEC; j++)
                                  arr_[n][i][j] = elem[i][j];
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  void ContainerGeneral<T,NUM_VEC,LEN_VEC>::set(int n, int m, T* elem)
  {
      for(int j = 0; j < LEN_VEC; j++)
          arr_[n][m][j] = elem[j];
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  void ContainerGeneral<T,NUM_VEC,LEN_VEC>::setAll(T def)
  {
      int len = size();
      for(int n = 0; n < len; n++)
          for(int i = 0; i < NUM_VEC; i++)
                          for(int j = 0; j < LEN_VEC; j++)
                                  arr_[n][i][j] = def;
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  void ContainerGeneral<T,NUM_VEC,LEN_VEC>::setAll(int to,T def)
  {
      int len = MathExtraPascal::min(to,size());
      for(int n = 0; n < len; n++)
          for(int i = 0; i < NUM_VEC; i++)
                          for(int j = 0; j < LEN_VEC; j++)
                                  arr_[n][i][j] = def;
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  T*** ContainerGeneral<T,NUM_VEC,LEN_VEC>::begin()
  {
          return arr_;
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  void* ContainerGeneral<T,NUM_VEC,LEN_VEC>::begin_slow_dirty()
  {
          return (void*) arr_;
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  int ContainerGeneral<T,NUM_VEC,LEN_VEC>::getElemSize()
  {
          return NUM_VEC*LEN_VEC*sizeof(T);
  }

  /* ----------------------------------------------------------------------
   min,max
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  T ContainerGeneral<T,NUM_VEC,LEN_VEC>::max_scalar()
  {
      T max = arr_[0][0][0];

      int len = size();
      for(int i = 0; i < len; i++)
            for(int j = 0; j < NUM_VEC; j++)
                for(int k = 0; k < LEN_VEC; k++)
                    if(arr_[i][j][k] > max)
                        max = arr_[i][j][k];

      return max;
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  T ContainerGeneral<T,NUM_VEC,LEN_VEC>::min_scalar()
  {
      T min = arr_[0][0][0];

      int len = size();
      for(int i = 0; i < len; i++)
            for(int j = 0; j < NUM_VEC; j++)
                for(int k = 0; k < LEN_VEC; k++)
                    if(arr_[i][j][k] < min)
                        min = arr_[i][j][k];

      return min;
  }

  /* ----------------------------------------------------------------------
   translate, rotate, scale
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  void ContainerGeneral<T,NUM_VEC,LEN_VEC>::scale(double factor)
  {
      if(this->containerProperties_.isScaleInvariant()) return;

      double factorApplied = 1.;
      for(int i = 0; i < containerProperties_.scalePower(); i++)
        factorApplied *= factor;

      int len = size();
      for(int i = 0; i < len; i++)
            for(int j = 0; j < NUM_VEC;j++)
                for(int k = 0; k < LEN_VEC; k++)
                    arr_[i][j][k] *= factorApplied;
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  void ContainerGeneral<T,NUM_VEC,LEN_VEC>::move(double *delta)
  {
      if(this->containerProperties_.isTranslationInvariant()) return;

      int len = size();

      for(int i = 0; i < len; i++)
            for(int j = 0; j < NUM_VEC; j++)
                for(int k = 0; k < LEN_VEC; k++)
                    arr_[i][j][k] += delta[k];
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  void ContainerGeneral<T,NUM_VEC,LEN_VEC>::moveElement(int i,double *delta)
  {
      if(this->containerProperties_.isTranslationInvariant()) return;

            for(int j = 0; j < NUM_VEC; j++)
                for(int k = 0; k < LEN_VEC; k++)
                    arr_[i][j][k] += delta[k];
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  void ContainerGeneral<T,NUM_VEC,LEN_VEC>::rotate(double *dQ)
  {
      if(this->containerProperties_.isRotationInvariant()) return;

      // ATTENTION: only correct for 3D vectors
      int len = size();
      for(int i = 0; i < len; i++)
            for(int j = 0; j < NUM_VEC; j++)
              MathExtraPascal::vec_quat_rotate(arr_[i][j],dQ);
  }

  /* ----------------------------------------------------------------------
   buffer size for all elements, push / pop for all elements
   used for global properties
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  int ContainerGeneral<T,NUM_VEC,LEN_VEC>::bufSize(OperationProperties &op)
  {
      if(!this->containerProperties_.decidePackUnpackOperation(op))
            return 0;

      if(!this->containerProperties_.decideCommOperation(op))
            return 0;

      return (1 + size()*NUM_VEC*LEN_VEC);
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  int ContainerGeneral<T,NUM_VEC,LEN_VEC>::pushToBuffer(double *buf,OperationProperties &op)
  {
          //TODO throw error if sizeof(T) > sizeof(double)

          int m = 0;

          if(!this->containerProperties_.decidePackUnpackOperation(op))
            return 0;

          int len = size();

          buf[m++] = static_cast<double>(len);

          for(int i = 0; i < len; i++)
            for(int j = 0; j < NUM_VEC; j++)
                for(int k = 0; k < LEN_VEC; k++)
                    buf[m++] = static_cast<double>(arr_[i][j][k]);

          return (1 + len*NUM_VEC*LEN_VEC);
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  int ContainerGeneral<T,NUM_VEC,LEN_VEC>::popFromBuffer(double *buf,OperationProperties &op)
  {
          int nNew, m = 0;

          if(!this->containerProperties_.decidePackUnpackOperation(op))
            return 0;

          //NP always uses buffer data

          if(this->containerProperties_.decideCreateNewElements(op))
          {
              T** tmp;
              create<T>(tmp,NUM_VEC,LEN_VEC);

              nNew = static_cast<int>(buf[m++]);

              for(int i = 0; i < nNew; i++)
              {
                for(int j = 0; j < NUM_VEC; j++)
                    for(int k = 0; k < LEN_VEC; k++)
                        tmp[j][k] = static_cast<T>(buf[m++]);
                add(tmp);
              }

              destroy<T>(tmp);

              return (1 + nNew*NUM_VEC*LEN_VEC);
          }
          else return 0;
  }

  /* ----------------------------------------------------------------------
   buffer size for a list of elements, push / pop a list of elements
   used for borders, fw and rev comm for element properties
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  int ContainerGeneral<T,NUM_VEC,LEN_VEC>::elemListBufSize(int n,OperationProperties &op)
  {
      if(!this->containerProperties_.decidePackUnpackOperation(op))
            return 0;

      if(!this->containerProperties_.decideCommOperation(op))
            return 0;

      return (n*NUM_VEC*LEN_VEC);
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  int ContainerGeneral<T,NUM_VEC,LEN_VEC>::pushElemListToBuffer(int n, int *list,double *buf,OperationProperties &op)
  {
        int i,m = 0;

        if(!this->containerProperties_.decidePackUnpackOperation(op))
            return 0;

        if(!this->containerProperties_.decideCommOperation(op))
            return 0;

        for(int ii = 0; ii < n; ii++)
        {
            i = list[ii];
            for(int j = 0; j < NUM_VEC; j++)
                for(int k = 0; k < LEN_VEC; k++)
                    buf[m++] = static_cast<double>(arr_[i][j][k]);
        }

        return (n*NUM_VEC*LEN_VEC);
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  int ContainerGeneral<T,NUM_VEC,LEN_VEC>::popElemListFromBuffer(int first, int n, double *buf,OperationProperties &op)
  {
        int m = 0;

        if(!this->containerProperties_.decidePackUnpackOperation(op))
            return 0;

        //NP check if uses communicates data
        //NP both cases possible for borders
        //NP fw comm always uses buffer data
        bool pullBuf = this->containerProperties_.decideCommOperation(op);

        //NP borders and exchange create new elements
        //NP fw comm overwrites existing data
        bool createElem = this->containerProperties_.decideCreateNewElements(op);

        T** tmp;
        create<T>(tmp,NUM_VEC,LEN_VEC);

        for(int i = first; i < first+n; i++)
        {
            for(int j = 0; j < NUM_VEC; j++)
                for(int k = 0; k < LEN_VEC; k++)
                    (createElem ? tmp[j][k] : arr_[i][j][k]) = (pullBuf ? static_cast<T>(buf[m++]) : static_cast<T>(0));

            if(createElem) add(tmp);
        }

        destroy<T>(tmp);

        return m;
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  int ContainerGeneral<T,NUM_VEC,LEN_VEC>::pushElemListToBufferReverse(int first, int n, double *buf,OperationProperties &op)
  {
        int m = 0;

        if(!this->containerProperties_.decidePackUnpackOperation(op))
            return 0;

        //NP always uses buffer data

        for(int i = first; i < first+n; i++)
        {
            for(int j = 0; j < NUM_VEC; j++)
                for(int k = 0; k < LEN_VEC; k++)
                    buf[m++] = static_cast<double>(arr_[i][j][k]);
        }

        return (n*NUM_VEC*LEN_VEC);
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  int ContainerGeneral<T,NUM_VEC,LEN_VEC>::popElemListFromBufferReverse(int n, int *list,double *buf,OperationProperties &op)
  {
        int i,m = 0;

        if(!this->containerProperties_.decidePackUnpackOperation(op))
            return 0;

        //NP always uses buffer data

        //NP never creates new elements, always unpacks at existing ones
        for(int ii = 0; ii < n; ii++)
        {
            i = list[ii];
            for(int j = 0; j < NUM_VEC; j++)
                for(int k = 0; k < LEN_VEC; k++)
                    arr_[i][j][k] += static_cast<T>(buf[m++]);
        }

        return (n*NUM_VEC*LEN_VEC);
  }

  /* ----------------------------------------------------------------------
   buffer size for a single element, push / pop a single element
   used for exchange of single elements
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  int ContainerGeneral<T,NUM_VEC,LEN_VEC>::elemBufSize(OperationProperties &op)
  {
      /*NL*/ //if(OPERATION_RESTART == operation) printf("Container ID %s called\n",id_);

      if(!this->containerProperties_.decidePackUnpackOperation(op))
            return 0;

      if(!this->containerProperties_.decideCommOperation(op))
            return 0;
      /*NL*/ //if(OPERATION_RESTART == operation) printf("   (size is %d)\n",NUM_VEC*LEN_VEC);
      return (NUM_VEC*LEN_VEC);
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  int ContainerGeneral<T,NUM_VEC,LEN_VEC>::pushElemToBuffer(int i, double *buf,OperationProperties &op)
  {
        int m = 0;

        if(!this->containerProperties_.decidePackUnpackOperation(op))
            return 0;

        if(!this->containerProperties_.decideCommOperation(op))
            return 0;

        for(int j = 0; j < NUM_VEC; j++)
            for(int k = 0; k < LEN_VEC; k++)
                buf[m++] = static_cast<double>(arr_[i][j][k]);

        return m;
  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  int ContainerGeneral<T,NUM_VEC,LEN_VEC>::popElemFromBuffer(double *buf,OperationProperties &op)
  {
        int m = 0;

        if(!this->containerProperties_.decidePackUnpackOperation(op))
            return 0;

        //NP pop for a single element always creates a new element

        //NP for comm_none and comm_reverse types, do not use buffer data
        bool pullBuf = this->containerProperties_.decideCommOperation(op);

        T** tmp;
        create<T>(tmp,NUM_VEC,LEN_VEC);

        for(int j = 0; j < NUM_VEC; j++)
            for(int k = 0; k < LEN_VEC; k++)
                tmp[j][k] = pullBuf ? static_cast<T>(buf[m++]) : static_cast<T>(0);

        add(tmp);
        destroy<T>(tmp);

        return m;
  }

#endif
