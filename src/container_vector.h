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

#ifndef LMP_CONTAINER_VECTOR
#define LMP_CONTAINER_VECTOR

#include "container_general.h"
#include "memory.h"

namespace PASCAL_NS
{
  template<typename T, int LEN_VEC>
  class ContainerVector : public ContainerGeneral <T, 1, LEN_VEC>
  {
    public:
          ContainerVector(ParticleDataContainerProperties &cp);
          ContainerVector(ContainerVector<T,LEN_VEC> const &orig);
          virtual ~ContainerVector();

          void add(T* elem);
          void get(int n, T* elem);
          void set(int n, T* elem);

          T max_elem(int n);
          T min_elem(int n);

          //void setAll(T def);
          T*& operator() (int n);
          T* const& operator() (int n) const;

          T** begin();
          void* begin_slow_dirty();
  };

  /* ----------------------------------------------------------------------
   constructors
  ------------------------------------------------------------------------- */

  template<typename T, int LEN_VEC>
  ContainerVector<T,LEN_VEC>::ContainerVector(ParticleDataContainerProperties &cp)
  : ContainerGeneral<T,1,LEN_VEC>(cp)
  {

  }

  template<typename T, int LEN_VEC>
  ContainerVector<T,LEN_VEC>::ContainerVector(ContainerVector<T,LEN_VEC> const &orig)
  : ContainerGeneral<T,1,LEN_VEC>(orig)
  {

  }

  /* ----------------------------------------------------------------------
   destructor
  ------------------------------------------------------------------------- */

  template<typename T, int LEN_VEC>
  ContainerVector<T,LEN_VEC>::~ContainerVector()
  {

  }

  /* ----------------------------------------------------------------------
   add element
  ------------------------------------------------------------------------- */

  template<typename T, int LEN_VEC>
  void ContainerVector<T,LEN_VEC>::add(T* elem)
  {
          if(ContainerGeneral<T,1,LEN_VEC>::numElem_ == ContainerGeneral<T,1,LEN_VEC>::maxElem_)
          {
                  grow(ContainerGeneral<T,1,LEN_VEC>::arr_,ContainerGeneral<T,1,LEN_VEC>::maxElem_+GROW,1,LEN_VEC);
                  ContainerGeneral<T,1,LEN_VEC>::maxElem_ += GROW;
          }
          for(int i=0;i<LEN_VEC;i++)
                  ContainerGeneral<T,1,LEN_VEC>::arr_[ContainerGeneral<T,1,LEN_VEC>::numElem_][0][i] = elem[i];

          ContainerGeneral<T,1,LEN_VEC>::numElem_++;
  }

  /* ----------------------------------------------------------------------
   access
  ------------------------------------------------------------------------- */

  template<typename T, int LEN_VEC>
  T*& ContainerVector<T,LEN_VEC>::operator() (int n)
  {
          return ContainerGeneral<T,1,LEN_VEC>::arr_[n][0];
  }

  template<typename T, int LEN_VEC>
  T* const& ContainerVector<T,LEN_VEC>::operator() (int n) const
  {
          return ContainerGeneral<T,1,LEN_VEC>::arr_[n][0];
  }

  template<typename T, int LEN_VEC>
  void ContainerVector<T,LEN_VEC>::get(int n, T* elem)
  {
          for(int i = 0; i < LEN_VEC; i++)
                  elem[i] = ContainerGeneral<T,1,LEN_VEC>::arr_[n][0][i];
  }

  template<typename T, int LEN_VEC>
  void ContainerVector<T,LEN_VEC>::set(int n, T* elem)
  {
          for(int i = 0; i < LEN_VEC; i++)
                  ContainerGeneral<T,1,LEN_VEC>::arr_[n][0][i] = elem[i];
  }

  template<typename T, int LEN_VEC>
  T ContainerVector<T,LEN_VEC>::max_elem(int n)
  {
          T max = ContainerGeneral<T,1,LEN_VEC>::arr_[n][0][0];

          for(int i = 1; i < LEN_VEC; i++)
                  if(ContainerGeneral<T,1,LEN_VEC>::arr_[n][0][i] > max)
                    max = ContainerGeneral<T,1,LEN_VEC>::arr_[n][0][i];

          return max;
  }

  template<typename T, int LEN_VEC>
  T ContainerVector<T,LEN_VEC>::min_elem(int n)
  {
          T min = ContainerGeneral<T,1,LEN_VEC>::arr_[n][0][0];

          for(int i = 1; i < LEN_VEC; i++)
                  if(ContainerGeneral<T,1,LEN_VEC>::arr_[n][0][i] < min)
                    min = ContainerGeneral<T,1,LEN_VEC>::arr_[n][0][i];

          return min;
  }
/*
  template<typename T, int LEN_VEC>
  void ContainerVector<T,LEN_VEC>::setAll(T def)
  {
      int len = this->size();
      for(int n = 0; n < len; n++)
          for(int i = 0; i < LEN_VEC; i++)
                  ContainerGeneral<T,1,LEN_VEC>::arr_[n][0][i] = def;
  }
*/
  template<typename T, int LEN_VEC>
  T** ContainerVector<T,LEN_VEC>::begin()
  {
          return &(ContainerGeneral<T,1,LEN_VEC>::arr_[0][0]);
  }

  template<typename T, int LEN_VEC>
  void* ContainerVector<T,LEN_VEC>::begin_slow_dirty()
  {
          return (void*) &(ContainerGeneral<T,1,LEN_VEC>::arr_[0][0]);
  }

} /* PASCAL_NS */
#endif /* CONTAINER_VECTOR */
