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


#ifndef PASC_MODEL_CHEMISTRY_CONTAINER_H
#define PASC_MODEL_CHEMISTRY_CONTAINER_H

#include "stdio.h"
#include "pascal_base.h"
#include "pascal_base_interface.h"
#include "model_chemistry.h"
#include <map>
#include <string>

namespace PASCAL_NS
{


class ModelChemistryContainer : public ParScaleBase, public ParScaleBaseInterface
{
    public:

      ModelChemistryContainer(ParScale *ptr);
      virtual ~ModelChemistryContainer();

      void parse_command(int narg,char const* const* arg);

      void begin_of_step();
      void pre_middle_of_step();
      void post_middle_of_step();
      void end_of_step();
      
      int nrChemistryEqns()    const      {return modelChemistryEqns_.size();};
      
      //Access model eqns
      ModelChemistry* modelChemistryEqn(int i)    const {return modelChemistryEqns_[i];};
      
    private:

      typedef ModelChemistry *(*ModelChemistryCreator)(ParScale*, char *name);
      std::map<std::string,ModelChemistryCreator> *model_map_;

      template <typename T> static ModelChemistry *model_creator(ParScale *ptr, char *name);

      vector<ModelChemistry*> modelChemistryEqns_;


};

} //end PASCAL_NS

#endif
