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
	This class is the base class for all chemical models (e.g., for reaction rates,
	heat release due to reactions, sintering, etc.).
-----------------------------------------------------------------------------------*/

#ifndef PASC_MODEL_CHEMISTRY_H
#define PASC_MODEL_CHEMISTRY_H

#include "model_base.h"
#include "particle_mesh.h"
#include "integrator.h"
#include <sundials/sundials_direct.h>
#include "particle_data.h"
#include "memory_ns.h"
#include "qjson_includes.h"

using namespace PASCAL_MEMORY_NS;

namespace PASCAL_NS
{

class ModelChemistry : public ModelBase
{
    friend class ModelChemistryMultiReact;
    friend class ModelChemistrySingleReaction;

    public:

      ModelChemistry(ParScale *ptr,  char *name);
      virtual ~ModelChemistry();
      
      virtual void    init(int narg, char const* const* arg, int eqnType, int modelChemistryID);

      virtual void    begin_of_step();

      virtual inline  double  kSurface_arrhenius(int particleID) {return kSurface_arrhenius_;}; //return the surface reaction rate
      
      

    private:
	
      mutable QJsonObject   properties_;
      QJsonObject           rates_;

      vector<QString>  speciesNames_;
      int       heatEqnID_;     
      int       speciesEqnID_;     
      int       speciesGridPoints_;     
      bool      species1DModel_;        //Indicate if we have a 1D model or other (e.g., shrinking core model)

      //TODO: Need to be cleaned
      bool      eqnReactionActive_;     //switch to activate or deactive reactions

      double**  speciesSource_;         //source for each position and species
      
      double*   heatSourceJac_;         //diagonal elements of the Jacobian from chemical source terms, heat
      double**  speciesSourceJac_;      //diagonal elements of the Jacobian from chemical source terms, species
      
	  
	  double  * tempIntraDataSpecies_;
      double  * tempIntraDataHeat_; 
 
      double    kSurface_arrhenius_;    //surface reation rate for arrhenius eqn 
 
        
      bool verbose_;
    
      bool elementary_reaction_;      

};

} //end PASCAL_NS


#endif
