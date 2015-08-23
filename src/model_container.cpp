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



#include "model_container.h"
#include "input.h"
#include "style_model.h"

using namespace PASCAL_NS;

/* ----------------------------------------------------------------------
   ModelContainer Constructor
------------------------------------------------------------------------- */

ModelContainer::ModelContainer(ParScale *ptr) : ParScaleBase(ptr),
    model_map_(new std::map<std::string,ModelCreator>())
{
  // fill map with all models listed in style_model.h

#define MODEL_CLASS
#define ModelStyle(key,Class) \
  (*model_map_)[#key] = &model_creator<Class>;
#include "style_model.h"
#undef ModelStyle
#undef MODEL_CLASS

}

/* ---------------------------------------------------------------------- */

ModelContainer::~ModelContainer()
{
    delete model_map_;
}

/* ----------------------------------------------------------------------
   one instance per model in style_model.h
------------------------------------------------------------------------- */

template <typename T>
ModelBase *ModelContainer::model_creator(ParScale *ptr, char *name)
{
  return new T(ptr,name);
}

/* ----------------------------------------------------------------------
   settings
------------------------------------------------------------------------- */

void ModelContainer::parse_command(int narg,char const* const* arg)
{
    int n = strlen(arg[1]) + 1;   
    char *modelName = new char[n];
    strcpy(modelName,arg[1]);

    if (model_map_->find(arg[0]) != model_map_->end())
    {
        ModelCreator model_creator = (*model_map_)[arg[0]];
        models_.push_back(model_creator(pascal_ptr(), modelName ) );

        printf("...this model is registered with ID %d\n", models_.size()-1);
        //fix[ifix] = fix_creator(lmp,narg,arg);

        //map_[map_string]->parse_command(narg,arg);
        //intialize this model equation

        models_[models_.size()-1].init(narg, arg);
    }
    else
        printf("ModelContainer PARSING: model name not found\n");

    delete [] modelName;
}


// ----------------------------------------------------------------------
void ModelContainer::begin_of_step()
{

    for(int iModel=0; iModel<models_.size(); iModel++)
            models_[iModel].begin_of_step();
   
}

// ----------------------------------------------------------------------
void ModelContainer::pre_middle_of_step()
{

    for(int iModel=0; iModel<models_.size(); iModel++)
            models_[iModel].pre_middle_of_step();

}
// ----------------------------------------------------------------------
void ModelContainer::post_middle_of_step()
{

    for(int iModel=0; iModel<models_.size(); iModel++)
            models_[iModel].post_middle_of_step();

}// ----------------------------------------------------------------------
void ModelContainer::end_of_step()
{

    for(int iModel=0; iModel<models_.size(); iModel++)
            models_[iModel].end_of_step();

}


