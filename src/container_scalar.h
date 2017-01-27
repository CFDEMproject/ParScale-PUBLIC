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

#ifndef LMP_CONTAINER_SCALAR
#define LMP_CONTAINER_SCALAR

#include "container_general.h"
#include "memory.h"
#include <limits>

namespace PASCAL_NS
{
  template<typename T>
  class ContainerScalar : public ContainerGeneral <T, 1, 1>
  {
    public:
          ContainerScalar(ParticleDataContainerProperties &cp);
          ContainerScalar(ContainerScalar<T> const &orig);
          virtual ~ContainerScalar();

          void add(T elem);
          T get(int n);
          void set(int n, T elem);
          void setAll(T def);
          T* begin();
          void* begin_slow_dirty();
          T& operator() (int n);
          T const& operator() (int n) const;
          T max();
          T max(int);
  };

  /* ----------------------------------------------------------------------
   constructors
  ------------------------------------------------------------------------- */

  template<typename T>
  ContainerScalar<T>::ContainerScalar(ParticleDataContainerProperties &cp)
  : ContainerGeneral<T,1,1>(cp)
  {

  }

  template<typename T>
  ContainerScalar<T>::ContainerScalar(ContainerScalar<T> const &orig)
  : ContainerGeneral<T,1,1>(orig)
  {

  }

  /* ----------------------------------------------------------------------
   destructor
  ------------------------------------------------------------------------- */

  template<typename T>
  ContainerScalar<T>::~ContainerScalar()
  {

  }

  /* ----------------------------------------------------------------------
   add element
  ------------------------------------------------------------------------- */

  template<typename T>
  void ContainerScalar<T>::add(T elem)
  {
          if(this->numElem_ == this->maxElem_)
          {
                  grow<T>(this->arr_,this->maxElem_+GROW,1,1);
                  this->maxElem_ += GROW;
          }
          this->arr_[this->numElem_][0][0] = elem;
          this->numElem_++;
  }

  /* ----------------------------------------------------------------------
   access
  ------------------------------------------------------------------------- */

  template<typename T>
  T& ContainerScalar<T>::operator() (int n)
  {
          return this->arr_[n][0][0];
  }

  template<typename T>
  T const& ContainerScalar<T>::operator() (int n) const
  {
          return this->arr_[n][0][0];
  }

  template<typename T>
  T ContainerScalar<T>::get(int n)
  {
          return this->arr_[n][0][0];
  }

  template<typename T>
  void ContainerScalar<T>::set(int n, T elem)
  {
          this->arr_[n][0][0] = elem;
  }

  template<typename T>
  void ContainerScalar<T>::setAll(T def)
  {
      int len = this->size();
      for(int n = 0; n < len; n++)
          this->arr_[n][0][0] = def;
  }

  template<typename T>
  T* ContainerScalar<T>::begin()
  {
          return &(this->arr_[0][0][0]);
  }

  template<typename T>
  void* ContainerScalar<T>::begin_slow_dirty()
  {
          return (void*) &(this->arr_[0][0][0]);
  }

  template<typename T>
  T ContainerScalar<T>::max()
  {
      int len = this->size();

      if(len == 0)
        return  (std::numeric_limits<T>::min)();

      T maxim = this->arr_[0][0][0];

      for(int n = 0; n < len; n++)
          if(this->arr_[n][0][0] > maxim)
            maxim = this->arr_[n][0][0];

      return maxim;
  }

  template<typename T>
  T ContainerScalar<T>::max(int to)
  {
      T maxim;

      int nn = MathExtraPascal::min(to,this->size());

      if(nn == 0)
        return  (std::numeric_limits<T>::min)();

      maxim = this->arr_[0][0][0];

      for(int n = 1; n < nn; n++)
          if(this->arr_[n][0][0] > maxim)
            maxim = this->arr_[n][0][0];

      return maxim;
  }
} /* PASCAL_NS */
#endif /* CONTAINER_SCALAR */
