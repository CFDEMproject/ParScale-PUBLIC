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



#include "model_eqn_container.h"
#include "input.h"
#include "style_model_eqn.h"

using namespace PASCAL_NS;

/* ----------------------------------------------------------------------
   ModelEqnContainer Constructor
------------------------------------------------------------------------- */

ModelEqnContainer::ModelEqnContainer(ParScale *ptr) : ParScaleBase(ptr),
    model_map_(new std::map<std::string,ModelEqnCreator>())
{
  // fill map with all models listed in style_model_eqn.h

#define MODEL_EQN_CLASS
#define ModelEqnStyle(key,Class) \
  (*model_map_)[#key] = &model_creator<Class>;
#include "style_model_eqn.h"
#undef ModelEqnStyle
#undef MODEL_EQN_CLASS

}

/* ---------------------------------------------------------------------- */

ModelEqnContainer::~ModelEqnContainer()
{
    delete model_map_;

    // TODO destroy all members of modelEqns_
}

/* ----------------------------------------------------------------------
   one instance per model in style_model_eqn.h
------------------------------------------------------------------------- */

template <typename T>
ModelEqn *ModelEqnContainer::model_creator(ParScale *ptr, char *name)
{
  return new T(ptr, name);
}

/* ----------------------------------------------------------------------
   settings
------------------------------------------------------------------------- */

void ModelEqnContainer::parse_command(int narg, char const* const* arg)
{
    int n = strlen(arg[1]) + 1;   
    char *modelName = new char[n];
    strcpy(modelName,arg[1]);

    if (model_map_->find(arg[0]) != model_map_->end())
    {
        ModelEqnCreator model_creator = (*model_map_)[arg[0]];
        
        if(strstr(modelName, "heat") != NULL)
        {
            modelHeatEqns_.push_back(model_creator(pascal_ptr(), modelName));
            printf("...this heat transport Eqn of type %s is registered with ID %d\n", 
                     arg[0],
                     modelHeatEqns_.size()-1);
            modelHeatEqns_[modelHeatEqns_.size()-1]->init(narg, arg, HEAT, modelHeatEqns_.size()-1);
        }
        else if(strstr(modelName, "species") != NULL)
        {
            modelSpeciesEqns_.push_back(model_creator(pascal_ptr(), modelName));
            printf("...this species transport Eqn of type %s is registered with ID %d\n",
                    arg[0],
                    modelSpeciesEqns_.size()-1);
            modelSpeciesEqns_[modelSpeciesEqns_.size()-1]->init(narg, arg, SPECIES, modelSpeciesEqns_.size()-1);
        }
        else
        {
            modelOtherEqns_.push_back(model_creator(pascal_ptr(), modelName));
            printf("...this other transport Eqn of type %s is registered with ID %d\n", 
                     arg[0],
                     modelOtherEqns_.size()-1);
            modelOtherEqns_[modelOtherEqns_.size()-1]->init(narg, arg, SPECIES, modelOtherEqns_.size()-1);
        }
        
    }
    else
        printf("FAIL: ModelEqnContainer PARSING: model name not found\n");

}


// ----------------------------------------------------------------------
void ModelEqnContainer::setupParticle()
{
    for(uint iEqn=0; iEqn<modelHeatEqns_.size(); iEqn++)
            modelHeatEqns_[iEqn]->setupParticle();

    for(uint iEqn=0; iEqn<modelSpeciesEqns_.size(); iEqn++)
            modelSpeciesEqns_[iEqn]->setupParticle();

    for(uint iEqn=0; iEqn<modelOtherEqns_.size(); iEqn++)
            modelOtherEqns_[iEqn]->setupParticle();

}

// ----------------------------------------------------------------------
void ModelEqnContainer::computeParticleProps()
{
  
		for(uint iEqn=0; iEqn<modelHeatEqns_.size(); iEqn++)
		        modelHeatEqns_[iEqn]->computeParticleProps();

		for(uint iEqn=0; iEqn<modelSpeciesEqns_.size(); iEqn++)
		        modelSpeciesEqns_[iEqn]->computeParticleProps();

		for(uint iEqn=0; iEqn<modelOtherEqns_.size(); iEqn++)
		        modelOtherEqns_[iEqn]->computeParticleProps();

}
// ----------------------------------------------------------------------

void ModelEqnContainer::begin_of_step()
{
	
    for(uint iEqn=0; iEqn<modelHeatEqns_.size(); iEqn++)
            modelHeatEqns_[iEqn]->begin_of_step();

    for(uint iEqn=0; iEqn<modelSpeciesEqns_.size(); iEqn++)
            modelSpeciesEqns_[iEqn]->begin_of_step();

    for(uint iEqn=0; iEqn<modelOtherEqns_.size(); iEqn++)
            modelOtherEqns_[iEqn]->begin_of_step();
   
}

// ----------------------------------------------------------------------
void ModelEqnContainer::pre_middle_of_step()
{
}
// ----------------------------------------------------------------------
void ModelEqnContainer::post_middle_of_step()
{
}// ----------------------------------------------------------------------
void ModelEqnContainer::end_of_step()
{
}


