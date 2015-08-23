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


#ifndef PASC_MODEL_EQN_CONTAINER_H
#define PASC_MODEL_EQN_CONTAINER_H

#include "stdio.h"
#include "pascal_base.h"
#include "pascal_base_interface.h"
#include "model_eqn.h"
#include <map>
#include <string>
#include <boost/ptr_container/ptr_vector.hpp>

namespace PASCAL_NS
{

class ModelEqnContainer : public ParScaleBase, public ParScaleBaseInterface
{
    public:

      ModelEqnContainer(ParScale *ptr);
      virtual ~ModelEqnContainer();

      void parse_command(int narg,char const* const* arg);
      void setupParticle();
      void computeParticleProps();

      void begin_of_step();
      void pre_middle_of_step();
      void post_middle_of_step();
      void end_of_step();
      
      int nrHeatEqns()    const      {return modelHeatEqns_.size();};
      int nrSpeciesEqns() const      {return modelSpeciesEqns_.size();};
      int nrOtherEqns()   const      {return modelOtherEqns_.size();};
      int nrEqns()        const      {return modelHeatEqns_.size() + modelSpeciesEqns_.size() + modelOtherEqns_.size();};

      int lookupEqn(const char *name) const ;

      //Access model eqns
      const ModelEqn* modelHeatEqn(int i)    const {return &modelHeatEqns_[i];};
      const ModelEqn* modelSpeciesEqn(int i) const {return &modelSpeciesEqns_[i];};
      const ModelEqn* modelOtherEqn(int i)   const {return &modelOtherEqns_[i];};

      const ModelEqn* modelEqn(int i)   const //return global model Eqn
      { 
        if(i<nrHeatEqns())
            return &modelHeatEqns_[i];

        else if(i < (nrHeatEqns()+nrSpeciesEqns()) )
            return &modelSpeciesEqns_[i-nrHeatEqns()];

        else if(i < nrEqns() )
            return &modelOtherEqns_[i-nrHeatEqns()-nrSpeciesEqns()];

        else
            return NULL;
      };
 
    private:

      typedef ModelEqn *(*ModelEqnCreator)(ParScale*, char *name);
      std::map<std::string,ModelEqnCreator> *model_map_;

      template <typename T> static ModelEqn *model_creator(ParScale *ptr, char *name);
    
      boost::ptr_vector<ModelEqn> modelHeatEqns_;
      boost::ptr_vector<ModelEqn> modelSpeciesEqns_;
      boost::ptr_vector<ModelEqn> modelOtherEqns_;

};

} //end PASCAL_NS

#endif
