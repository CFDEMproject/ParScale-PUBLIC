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

#ifndef LMP_CONTAINER_MULTI_VECTOR
#define LMP_CONTAINER_MULTI_VECTOR

#include "container_general.h"

namespace PASCAL_NS{
  template<typename T, int NUM_VEC, int LEN_VEC>
  class ContainerMultiVector : public ContainerGeneral <T, NUM_VEC, LEN_VEC>
  {
      public:
          ContainerMultiVector(ParticleDataContainerProperties &cp);
          ContainerMultiVector(ContainerMultiVector<T,NUM_VEC,LEN_VEC> const &orig);
          virtual ~ContainerMultiVector();
  };

  /* ----------------------------------------------------------------------
   constructors
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  ContainerMultiVector<T,NUM_VEC,LEN_VEC>::ContainerMultiVector(ParticleDataContainerProperties &cp)
  : ContainerGeneral<T,NUM_VEC,LEN_VEC>(cp)
  {

  }

  template<typename T, int NUM_VEC, int LEN_VEC>
  ContainerMultiVector<T,NUM_VEC,LEN_VEC>::ContainerMultiVector(ContainerMultiVector<T,NUM_VEC,LEN_VEC> const &orig)
  : ContainerGeneral<T,NUM_VEC,LEN_VEC>(orig)
  {

  }

  /* ----------------------------------------------------------------------
   destructor
  ------------------------------------------------------------------------- */

  template<typename T, int NUM_VEC, int LEN_VEC>
  ContainerMultiVector<T,NUM_VEC,LEN_VEC>::~ContainerMultiVector()
  {

  }

} /* PASCAL_NS */
#endif /* CONTAINER_MULTI_VECTOR */
