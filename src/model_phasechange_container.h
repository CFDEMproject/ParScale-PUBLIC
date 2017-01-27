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
/*
    Contributing author:    Stefan Radl, TU Graz
    Container for phase change models

\*-----------------------------------------------------------------------------------*/


#ifndef PASC_MODEL_PHASECHANGE_CONTAINER_H
#define PASC_MODEL_PHASECHANGE_CONTAINER_H

#include "stdio.h"
#include "pascal_base.h"
#include "pascal_base_interface.h"
#include "model_phasechange.h"
#include <map>
#include <string>
#include <boost/ptr_container/ptr_vector.hpp>

namespace PASCAL_NS
{


class ModelPhaseChangeContainer : public ParScaleBase, public ParScaleBaseInterface
{
    public:

      ModelPhaseChangeContainer(ParScale *ptr);
      virtual ~ModelPhaseChangeContainer();

      void parse_command(int narg,char const* const* arg);

      void postParticleDataBuild();

      void begin_of_step();
      void pre_middle_of_step();
      void post_middle_of_step();
      void end_of_step();

      int nrPhaseChangeEqns()    const      {return modelPhaseChangeEqns_.size();};

      //Access model eqns
      const ModelPhaseChange* modelPhaseChangeEqn(int i)    const {return &modelPhaseChangeEqns_[i];};

    private:

      typedef ModelPhaseChange *(*ModelPhaseChangeCreator)(ParScale*, char *name);
      std::map<std::string,ModelPhaseChangeCreator> *model_map_;

      template <typename T> static ModelPhaseChange *model_creator(ParScale *ptr, char *name);

      boost::ptr_vector<ModelPhaseChange> modelPhaseChangeEqns_;


};

} //end PASCAL_NS

#endif
