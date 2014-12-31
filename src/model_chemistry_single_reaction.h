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

/* Model chemistry single reaction 
User must explicitly specify all reactaion parameters
Always a
    - irreversible
    - non-elementary 
reaction!

\*-----------------------------------------------------------------------------------*/

#ifdef MODEL_CHEMISTRY_CLASS

ModelChemistryStyle(SingleReaction, ModelChemistrySingleReaction)

#else
 
#ifndef PASC_CHEM_SIMPLE_H
#define PASC_CHEM_SIMPLE_H

#include "model_chemistry.h"
#include "model_eqn_container.h"
#include "chemistry_reaction_single.h"
#include "chemistry_grainmodel.h"
#include "input.h"


namespace PASCAL_NS
{

class ModelChemistrySingleReaction : public ModelChemistry
{
    public:

      ModelChemistrySingleReaction(ParScale *ptr, char *name);

     void begin_of_step();
     void init(int narg, char const* const* arg, int eqnType, int modelChemistryID);
     void read_chemisty_JSON();

    private:

    ParScale *ptr_;

    QJsonObject single_reaction_file_;
    
    ChemistryReactionSingle* reaction_;
    ChemistryGrainModel*     grainModel_;

    //needed for calculating of reaction rate, to be read in by chemistry_reader TODO: link to parser/chemkin reader
    bool                  isIsoThermal_;              //specify if a homogeneous temperature is assumed
    double                temperature_;               //specify a homogeneous temperature FOR ALL PARTICLES
    vector<std::string>   species_names_;             //name of all species in the reaction
    vector<double>        species_stoich_;            //stoichometry of all species in the reaction
    vector<double>        species_order_;             //reaction order of each species

    double arrhenius_A_;                        //pre-exponential factor
    double arrhenius_beta_;                     //exponential feactor 
    double arrhenius_E_A_;                      //activation energy 


    //needed for calculating of chemical source term
    double         deltaH_r;                    //heat of reaction
    vector<int>    original_IDs_reactant_;      //orginal ids of the reactants in the particle mem
    vector<double> k_f_i_;                      //reaction rate constants, needed for shrinking core model
    vector<double> q_i_dot_;                    //rates of the reaction at each grid point
    vector<double> molar_mass_;                 //molar mass of species

    QJsonObject     global_properties_;
};

} //end PASCAL_NS

#endif

#endif

