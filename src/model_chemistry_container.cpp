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



#include "model_chemistry_container.h"
#include "input.h"
#include "style_model_chemistry.h"

using namespace PASCAL_NS;

/* ----------------------------------------------------------------------
   ModelEqnContainer Constructor
------------------------------------------------------------------------- */

ModelChemistryContainer::ModelChemistryContainer(ParScale *ptr) : ParScaleBase(ptr),
    model_map_(new std::map<std::string,ModelChemistryCreator>())
{
  // fill map with all models listed in style_model_eqn.h

#define MODEL_CHEMISTRY_CLASS
#define ModelChemistryStyle(key,Class) \
  (*model_map_)[#key] = &model_creator<Class>;
#include "style_model_chemistry.h"
#undef ModelChemistryStyle
#undef MODEL_CHEMISTRY_CLASS

}

/* ---------------------------------------------------------------------- */

ModelChemistryContainer::~ModelChemistryContainer()
{
    delete model_map_;

    // TODO destroy all members of modelEqns_
}

/* ----------------------------------------------------------------------
   one instance per model in style_model_eqn.h
------------------------------------------------------------------------- */

template <typename T>
ModelChemistry *ModelChemistryContainer::model_creator(ParScale *ptr, char *name)
{
  return new T(ptr, name);
}

/* ----------------------------------------------------------------------
   settings
------------------------------------------------------------------------- */

void ModelChemistryContainer::parse_command(int narg, char const* const* arg)
{
    int n = strlen(arg[1]) + 1;   
    char *modelName = new char[n];
    strcpy(modelName,arg[1]);

    if (model_map_->find(arg[0]) != model_map_->end())
    {
        ModelChemistryCreator model_creator = (*model_map_)[arg[0]];
        modelChemistryEqns_.push_back(model_creator(pascal_ptr(), modelName));

        printf("...this chemistry model of type %s is registered with ID %lu\n", 
                arg[0],
                modelChemistryEqns_.size()-1);

        modelChemistryEqns_[modelChemistryEqns_.size()-1]->init(narg, arg, REACTION, modelChemistryEqns_.size()-1);
    }
    else
        printf("FAIL: ModelChemistryContainer PARSING: model name not found\n");

}


// ----------------------------------------------------------------------

void ModelChemistryContainer::begin_of_step()
{
	
    for(uint iEqn=0; iEqn<modelChemistryEqns_.size(); iEqn++)
            modelChemistryEqns_[iEqn]->begin_of_step();

}

// ----------------------------------------------------------------------
void ModelChemistryContainer::pre_middle_of_step()
{
}
// ----------------------------------------------------------------------
void ModelChemistryContainer::post_middle_of_step()
{
}// ----------------------------------------------------------------------
void ModelChemistryContainer::end_of_step()
{
}


