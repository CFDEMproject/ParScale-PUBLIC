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

#ifndef PASC_CHEM_REAC_SINGLE_H
#define PASC_CHEM_REAC_SINGLE_H


#include "stdio.h"
#include "error.h"
#include "input_base.h"
#include "pascal_base_accessible.h"
#include "container.h"
#include "pascal_base_interface.h"
#include "particle_mesh.h"
#include "model_eqn.h"
#include "model_eqn_container.h"
#include "particle_data.h"
#include "input.h"

#define UNIVERSAL_GAS_CONSTANT 8314.4621                //SI UNITS: J/kMol/K !!
#define INVERSE_UNIVERSAL_GAS_CONSTANT 0.0001202723625  //SI UNITS: kMol*K/J !!
#define MINIMUM_CONCENTRATION  1e-16                    //below which 0-order reaction will be linearly reduced

namespace PASCAL_NS
{

class ChemistryReactionSingle : public ParScaleBaseAccessible, ParScaleBaseInterface 
{
    public:

      ChemistryReactionSingle(ParScale *ptr);
      ~ChemistryReactionSingle();

    void calculate_rate_constants(int particleID_);
    void calculate_delta_S_i_o(int grid_point_);
    void calculate_delta_H_i_o(int grid_point_);
    void calculate_K_P_i(int grid_point_);
    void calculate_K_c_i(int grid_point_);
    void calculate_k_r_i(int grid_point_);
    void calculate_k_f_i(int grid_point_);
    void calculate_q_i(int grid_point_,int particleID_);
    void correct_q_i(double factor, int grid_point_,int particleID_);
    double integrateVolume_q_i(int particleID, double &particleVolume);
    
    inline double arrheniusRate(double& temperature)
    {
        double rate = MathExtraPascal::fastPow(temperature, arrhenius_beta_i_);
        rate *=  arrhenius_A_i_ 
               * exp(
                          -arrhenius_E_i_
                       / (UNIVERSAL_GAS_CONSTANT*temperature)
                    );
        return rate;
    }

    //This is the prefactor to the arrheniusRate that will give the 
    //derivative of the arrheniusRate with respect to the temperature
    //this is needed for the Jacobi calculation
    inline double arrheniusRateDerivativePrefactor(double& temperature)
    {
        double inverseTemperature = 1.0 / temperature;
        return  arrhenius_beta_i_ * inverseTemperature
              + arrhenius_E_i_    * inverseTemperature * inverseTemperature * INVERSE_UNIVERSAL_GAS_CONSTANT;
    }
    inline double computeArrheniusRateDerivativePrefactor(int _grid_point)
    {
        return arrheniusRateDerivativePrefactor(actual_temp_[_grid_point]);
    }

    inline double reactionProduct(vector<double>& concentration, vector<double>& exponent)
    {
        double product = 1.0;
        for (uint i_reactant = 0; i_reactant < concentration.size(); i_reactant++)
        {
            //compute product, however,
            //if reactant species is fully depleted, reaction stops
            if( exponent[i_reactant]==0 && concentration[i_reactant]<cMinimum_  ) 
                product *= concentration[i_reactant]/cMinimum_; //linear approximate for very small concentrations
            else
                product *= MathExtraPascal::fastPow(concentration[i_reactant],
                                                    exponent[i_reactant]);
                                               
        if (verbose_)
          printf("***reactant %d, conc: %g, exponent: %g, product: %g  \n",
                  i_reactant, concentration[i_reactant], exponent[i_reactant], product
                );
        }
        
        return product;
    }

    //Since chemistry_reaction_single has all properties connected of the local ID, the global corresponding species ID
    //for every reaction and reactant/product has to be connected with the global ID of that species.    
    void getSpeciesEqnIDs();  

    //function for empty all calculated vectors
    void empty_vectors();
    
    double returnArrhenius_A()      {return arrhenius_A_i_;}; 
    double returnArrhenius_beta()   {return arrhenius_beta_i_;};
    double returnArrhenius_E()      {return arrhenius_E_i_;};    


    void returnKfi(vector<double> & k_f_i) 
    {
        k_f_i=k_f_i_;
    }
    
    void returnQDoti(vector<double> & q_Dot_i) 
    {
        q_Dot_i=q_dot_i_;
    }

    void resetOriginalIDs()
    {
      if(original_IDs_reactant_.size()!=species_names_reactant_.size())
      {
        original_IDs_reactant_.resize(species_names_reactant_.size(),0);
        
        for(uint i =0; i<original_IDs_reactant_.size(); i++)
        {
            original_IDs_reactant_[i]= i + modelEqnContainer().nrHeatEqns();
            if(original_IDs_reactant_[i]>=modelEqnContainer().nrEqns())
                error().throw_error_one(FLERR,"ERROR: failed to reset original IDs for reactants    . Most likely, you have not enough species equations in your input file. \n");
            
            if (species_names_reactant_[i].compare(modelEqnContainer().modelEqn(original_IDs_reactant_[i])->nameSpecies()) !=0 )
            {
                printf("detected problem with reactant id %d (name: %s) when comparing to modelEqn dataId=%d, name: %s (nrHeatEqns: %d).\n", 
                        i, species_names_reactant_[i].c_str(),
                        modelEqnContainer().modelEqn(original_IDs_reactant_[i])->particleDataID(),
                        modelEqnContainer().modelEqn(original_IDs_reactant_[i])->nameSpecies(),
                        modelEqnContainer().nrHeatEqns()
                      );
                error().throw_error_one(FLERR,"ERROR: Mismatch of names in chemistry and in modelEqns detected. Re-order you modelEqns to match the names! \n");
            }
        }
      }
      
      if(original_IDs_product_.size()==species_names_product_.size()) 
            return;
            
      original_IDs_product_.resize(species_names_product_.size(),0);
      for(uint i =0; i<original_IDs_product_.size(); i++)
      {
            original_IDs_product_[i]= i + modelEqnContainer().nrHeatEqns();
            if(original_IDs_product_[i]>=modelEqnContainer().nrEqns())
                error().throw_error_one(FLERR,"ERROR: failed to reset original IDs for products. Most likely, you have not enough species equations in your input file. \n");
            if (species_names_product_[i].compare(modelEqnContainer().modelEqn(original_IDs_product_[i])->nameSpecies()) !=0 )
            {
                printf("detected problem with product id %d (name: %s) when comparing to modelEqn dataId=%d, name: %s (nrHeatEqns: %d). \n", 
                        i, species_names_product_[i].c_str(),
                        modelEqnContainer().modelEqn(original_IDs_product_[i])->particleDataID(),
                        modelEqnContainer().modelEqn(original_IDs_product_[i])->nameSpecies(),
                        modelEqnContainer().nrHeatEqns()
                      );
                error().throw_error_one(FLERR,"ERROR: Mismatch of names in chemistry and in modelEqns detected. Re-order you modelEqns to match the names! \n");
            }
            
      }
    }

    void returnOriginalIDsreactant(vector<int> & original_IDs_reactant) 
    {
        original_IDs_reactant=original_IDs_reactant_;
    }

    void returnOriginalIDsproduct(vector<int> & original_IDs_product) 
    {
        original_IDs_product=original_IDs_product_;
        
    }

    void returnSpeciesstoichreactant(vector<double> & species_stoich_reactant) 
    {
        species_stoich_reactant=species_stoich_reactant_;
    }

    void returnSpeciesstoichproduct(vector<double> & species_stoich_product) 
    {
        species_stoich_product=species_stoich_product_;
    }

    double returnStoichReactants(int n_) 
    {
        return species_stoich_reactant_.at(n_);
    }
    
    double returnStoichProducts(int m_) 
    {
        return species_stoich_product_.at(m_);
    }


    //Setting Arrhenius Constants for reaction
    void setArrheniusConstants(double x,double y, double z) 
    {
                arrhenius_A_i_ = x; 
                arrhenius_beta_i_ = y;
                arrhenius_E_i_ = z;
    }; 


    //Setting names of reactants and products for reaction
    void SetSpeciesNamesReactant(std::string name_reactant_) 
    {
                species_names_reactant_.push_back(name_reactant_);
    };

    void SetSpeciesNamesProduct(std::string name_product_) 
         {
                species_names_product_.push_back(name_product_);
         };
 

    //Setting stoichometric factors of reactants and products for reaction
    void SetSpeciesStoichReactant(double stoich_reactant_) 
    {
                species_stoich_reactant_.push_back(stoich_reactant_);
    };

    void SetSpeciesStoichProduct(double stoich_product_) 
    {
                species_stoich_product_.push_back(stoich_product_);
    };

    //Setting bools for dependencies //TODO: implement this
    void SetReversible(bool rev)                   {reversible_ = rev;};
    void SetTempDepend(bool temp_depend)           {temperature_dependend_ = temp_depend;};
    void SetPressureDepend(bool pressure_depend)   {pressure_dependend_ = pressure_depend;};
    void SetThirdBody(bool third_body)             {third_bodies_ = third_body;};             //TODO: implement this
    void SetElementaryReaction(bool elementayReac) {elementary_reaction_=elementayReac;};
   
    //setting the actual temperature to particular time step
    void SetActualTemp (double* actual_temp)         {actual_temp_ = actual_temp;};

    void setCMinimum(double _x)    
    {
         cMinimum_ = _x;
    }

    //Setting the Enthalpy, vector includes all species Enthalpies in order of species ID
    void SetEnthalpySpecies(vector<double> Enthalpy) 
        {
            delta_H_k_0_full_=Enthalpy;
        };

     //Setting the Enthalpy, vector includes all species Entropies in order of species ID
    void SetEntropySpecies(vector<double> Entropy)
        {
            delta_S_k_0_full_=Entropy;
        }; 
    
    void SetReactionOrders(double ReactionOrder)
        {
            reaction_order_full_.push_back(ReactionOrder);
        };

    void SetReactionOrderReactant(vector<double> ReactionOrder)
        {
            reaction_order_reactant_=ReactionOrder;
        };
    
    void printStatus()
    {
        if(  (species_names_reactant_.size() != reaction_order_reactant_.size()) 
           ||(species_names_reactant_.size() != species_stoich_reactant_.size())
          )
            error().throw_error_one(FLERR,"ERROR: names, order, or stoichiometry not set for all speicies. \n");

        printf("\n *** ChemistryReactionSingle - Status *** \n");
        printf("reversible: %d, temperature_dependend: %d, elementary_reaction_: %d, involved species: %lu \n", 
               reversible_, temperature_dependend_,
               elementary_reaction_,
               species_names_reactant_.size()
              );
        printf("Arrhenius-Props: A_i: %g,  beta: %g,  E_i: %g \n", 
               arrhenius_A_i_, arrhenius_beta_i_, arrhenius_E_i_
              );
        for(uint i=0; i<species_names_reactant_.size();i++)
        {
            printf("*   reactant: %s, \n", species_names_reactant_[i].c_str());
            printf("    order   : %g, \n", reaction_order_reactant_[i]);
            printf("    stoichm.: %g, \n", species_stoich_reactant_[i]);
            printf("  \n");
        }
        printf("\n *** ChemistryReactionSingle - Status *** \n \n");

    }

    private:

    //********Containers with size of grid points***********//
    int n_grip_points_;
    double* actual_temp_;

    vector<double> k_f_i_;                   //Forward rate constant of the ith reaction [dependend on order of reaction]
    vector<double> k_r_i_;                   //Backward rate constant of the ith reaction [dependend on order of reaction]
    vector<double> K_c_i_;                    //equilibrium constant
    vector<double> K_p_i_;  

    vector<double> q_dot_i_;                  //Rate of progress of the ith reaction 
    double         q_dot_i_Integral_;         //Volume integral of Rate of ith reaction 

    vector<double> delta_S_i_0_;              //molar entropy used for the calculation of K_p_i                 - to be calculated
    vector<double> delta_H_i_0_;              //molar enthalpy used for the calculation of K_p_i                - to be calculated


    //********Containers with size of number of species***********//

    vector<double> delta_S_k_0_full_;     //molar entropy of the kth species - includes all species
    vector<double> delta_H_k_0_full_;     //molar enthalpy of the kth species - includes all species
    vector<double> reaction_order_full_;  //Reaction order of kth species

   //********Containers with size of 1***********//

    bool reversible_;               //determines whether a reaction is reversible or not                - to be read in from chemkin
    bool temperature_dependend_;    //determines whether a reaction is temperature dependend or not     - to be read in from chemkin
    bool pressure_dependend_;       //determines whether a reaction is pressure dependend or not        - to be read in from chemkin
    bool third_bodies_;             //determines whether a third body reaction is influencing           - not supported yet- to be read in from chemkin

    double cMinimum_;
    double arrhenius_A_i_;          //pre-exponential factor                                            - to be read in from chemkin 
    double arrhenius_beta_i_;       //exponential factor for temperature dependency                    - to be read in from chemkin
    double arrhenius_E_i_;          //activation energy                                                 - to be read in from chemkin
    //Arrhenius k_f_i = A_i * T^(beta_i) * exp ((-E_i)/(R_c*T))

    double P_atm_;           //[Pa] WARNING: SI Units!


    //********Containers with size of number reactants/products***********//

    vector<std::string> species_names_reactant_;        //name of all reactant species in the reaction          - to be read in from chemkin
                                                        //needed for looking up enthalpy and entropy
    vector<std::string> species_names_product_;         //name of all product species in the reaction           - to be read in from chemkin
                                                        //needed for looking up enthalpy and entropy
       
    vector<double> species_stoich_reactant_;    //stoichometry of all reactant species in the reaction  - to be read in from chemkin
    vector<double> species_stoich_product_;     //stoichometry of all product species in the reaction   - to be read in from chemkin

    vector<double> delta_S_k_0_reactant_;     //molar entropy of the kth species - includes reactant species
    vector<double> delta_S_k_0_product_;      //molar entropy of the kth species - includes product species
 
    vector<double> delta_H_k_0_product_;      //molar enthalpy of the kth species - includes reactant species
    vector<double> delta_H_k_0_reactant_;     //molar enthalpy of the kth species - includes product species

    vector<double> reaction_order_reactant_;  //reaction order of the kth reactant - manually set
    vector<double> reaction_order_product_;   //reaction order of the kth product - manually set
 
    vector<int>    original_IDs_reactant_;    //vector with original IDs for reactants - corresponding speciesID from model_eqn
    vector<int>    original_IDs_product_;     //vector with original IDs for products- corresponding speciesID from model_eqn

    //*****other variables****//

    bool verbose_; 
    mutable bool    elementary_reaction_;

    QJsonObject     global_properties_;

};

} //end PASCAL_NS

#endif
