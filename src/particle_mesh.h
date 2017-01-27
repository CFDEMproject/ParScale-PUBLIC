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

   Copyright (C): 2014 DCS Computing GmbH (www.dcs-computing.com), Linz, Austria
                  2014 Graz University of Technology (ippt.tugraz.at), Graz, Austria
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
/*-----------------------------------------------------------------------------------
Description
    This class is the base class for meshes used for spatial discretization of each
    particle.
-----------------------------------------------------------------------------------*/

#ifndef PASC_PARTICLE_MESH_H
#define PASC_PARTICLE_MESH_H

#include "stdio.h"
#include "pascal_base.h"
#include "pascal_base_interface.h"

namespace PASCAL_NS
{

class ParticleMesh : public ParScaleBase, public ParScaleBaseInterface
{
    public:

      ParticleMesh(ParScale *ptr);

      int  nGridPoints() const {return nGridPoints_; };
      void setGridPoints(int x) {nGridPoints_ = x; };

      void parse_command(int narg,char const* const* arg);

    private:

      //TODO data structure must reflect discretization scheme
      // eg 27 points, shell model,...
      int       nGridPoints_;

      //TODO: properties to implement
      //double**  geometryBounds       //overall bounds of the particle, push-pull
      //bool      geometryIsConstant   //true in case particle can change its size, push-pulll

};

} //end PASCAL_NS

#endif
