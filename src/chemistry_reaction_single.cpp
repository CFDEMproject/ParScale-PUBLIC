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

#include "chemistry_reaction_single.h"
#include "input.h"

using namespace PASCAL_NS;

/* ----------------------------------------------------------------------
   Constructor / Destructor
------------------------------------------------------------------------- */

ChemistryReactionSingle::ChemistryReactionSingle(ParScale *ptr) 
 : 
   ParScaleBaseAccessible(ptr),
   cMinimum_(MINIMUM_CONCENTRATION),
   arrhenius_A_i_(-1.),
   arrhenius_beta_i_(-1.),
   arrhenius_E_i_(-1.),
   reversible_(0),
   pressure_dependend_(0),
   temperature_dependend_(1),
   P_atm_(1.013e5),
   K_c_i_(particleMesh().nGridPoints()),
   k_r_i_(particleMesh().nGridPoints()),
   k_f_i_(particleMesh().nGridPoints()),
   K_p_i_(particleMesh().nGridPoints()),
   q_dot_i_(particleMesh().nGridPoints()),
   delta_S_i_0_(particleMesh().nGridPoints()),
   delta_H_i_0_(particleMesh().nGridPoints())   
{
    //read verbose file
    input().openJsonFile("settings", "verbose", "verbose", global_properties_ );
    if(global_properties_["chemistry_single_reaction"].isNull())
        error().throw_error_one(FLERR,"ERROR: could not read verbose settings for this class. \n");
    verbose_   = global_properties_["chemistry_single_reaction"].toBool();    

}

// ----------------------------------------------------------------------
void ChemistryReactionSingle::calculate_rate_constants(int particleID_)
{

    for (int grid_point_ = 0; grid_point_< particleMesh().nGridPoints(); grid_point_++)
    {
        //getSpeciesEqnIDs();
        calculate_delta_S_i_o(grid_point_);
        calculate_delta_H_i_o(grid_point_);
        calculate_K_c_i(grid_point_);
        calculate_k_f_i(grid_point_);
        calculate_k_r_i(grid_point_);
        calculate_q_i(grid_point_,particleID_);  
    }
}

// ----------------------------------------------------------------------
void ChemistryReactionSingle::getSpeciesEqnIDs()
{
    int species_counter_;
    std::string str1="species"; 
    species_counter_ = species_names_reactant_.size() + species_names_product_.size();
//    printf("species_counter_1 = %i \n",species_counter_);
   
    for(int iEqn=0; iEqn < modelEqnContainer().nrSpeciesEqns(); iEqn++)
    {
          for(uint i_=0; i_ < species_names_reactant_.size(); i_++)
          {     
//               printf("...comparing strings '%s' and '%s' \n",
//                       species_names_reactant_[i_].insert(0,str1).c_str(), 
//                       modelEqnContainer().modelSpeciesEqn(iEqn)->name()
//                     );        
               //species_names_reactant_[i_].erase(0,7);

               if(strcmp(
                         species_names_reactant_[i_].insert(0,str1).c_str(), 
                         modelEqnContainer().modelSpeciesEqn(iEqn)->name()
                        )==0
                 )
               {
//                   printf(" ...found correspondence. \n",modelEqnContainer().modelSpeciesEqn(iEqn)->name(),modelEqnContainer().modelEqn(iEqn + modelEqnContainer().nrHeatEqns())->name());
                   delta_H_k_0_reactant_.push_back(delta_H_k_0_full_[iEqn]);
                   delta_S_k_0_reactant_.push_back(delta_S_k_0_full_[iEqn]);

                   if (elementary_reaction_ != true)
                   {
                        reaction_order_reactant_.push_back(reaction_order_full_[iEqn]);
                   }

                   original_IDs_reactant_.push_back(iEqn + modelEqnContainer().nrHeatEqns());
//                   printf("Push back for %d ractants done. \n", original_IDs_reactant_.size()); 
                    species_counter_--;  
//                    printf("species_counter_2 = %i \n",species_counter_);            
               }   
               species_names_reactant_[i_].erase(0,7);       
          }

          for(uint i2_=0; i2_ < species_names_product_.size(); i2_++)
          {     
                //printf("You want to compare %s == %s \n",species_names_product_[i2_].insert(0,str1).c_str(), modelEqnContainer().modelSpeciesEqn(iEqn)->name()); 
                //species_names_product_[i2_].erase(0,7);

               if(strcmp(species_names_product_[i2_].insert(0,str1).c_str(), modelEqnContainer().modelSpeciesEqn(iEqn)->name())==0)
               {
                    //printf("Is %s and %s the same?\n",modelEqnContainer().modelSpeciesEqn(iEqn)->name(),modelEqnContainer().modelEqn(iEqn + modelEqnContainer().nrHeatEqns())->name());
                    delta_H_k_0_product_.push_back(delta_H_k_0_full_[iEqn]);
                    delta_S_k_0_product_.push_back(delta_S_k_0_full_[iEqn]); 

                    if (elementary_reaction_ != true)
                    {
                        reaction_order_product_.push_back(reaction_order_full_[iEqn]);
                    }

                    original_IDs_product_.push_back(iEqn + modelEqnContainer().nrHeatEqns());
                    //printf("Push back for products done \n");  
                    species_counter_--;
//                    printf("species_counter_3 = %i \n",species_counter_);
               }  
               species_names_product_[i2_].erase(0,7);        
          }
    }
    
    if (species_counter_ != 0)
    {
        printf("It seems like you have a problem in this reaction. Please make sure every species has a corresponding species eqn.\n");
        error().throw_error_one(FLERR," Check consitency of species eqn and species mentioned in chemkin.inp under SPECIES.\n");
    }
}

// ----------------------------------------------------------------------
void ChemistryReactionSingle::calculate_delta_S_i_o(int grid_point_)
{
    
    double Si_sum = 0.0;

    //Loop reactants
    for (uint i_reactant = 0; i_reactant < species_stoich_reactant_.size(); i_reactant++)
    {
        double K_th_reactant_ = species_stoich_reactant_[i_reactant]*delta_S_k_0_reactant_[i_reactant];
        Si_sum += K_th_reactant_;
    } 
    //Loop products
    for (uint i_product_ = 0; i_product_ < species_stoich_product_.size(); i_product_++)
    {
        double K_th_product_ = species_stoich_product_[i_product_]*delta_S_k_0_product_[i_product_];
        Si_sum += K_th_product_;
    }  
    
    delta_S_i_0_[grid_point_] = Si_sum;
    
    if (verbose_)
        printf("Sum Si = %g \n",Si_sum);
    

}

// ----------------------------------------------------------------------
void ChemistryReactionSingle::calculate_delta_H_i_o(int grid_point_)
{
    double Hi_sum = 0.0;

    //Loop reactants
    //printf("size of reactant vector in single_reaction 1 = %i \n",species_stoich_reactant_.size());
    for (uint i_reactant = 0; i_reactant < species_stoich_reactant_.size(); i_reactant++)
    {
        double H_th_reactant_ = species_stoich_reactant_[i_reactant]*delta_H_k_0_reactant_[i_reactant];
        Hi_sum += H_th_reactant_;
    } 
    //Loop products
    for (uint i_product_ = 0; i_product_ < species_stoich_product_.size(); i_product_++)
    {
        double H_th_product_ = species_stoich_product_[i_product_]*delta_H_k_0_product_[i_product_];
        Hi_sum += H_th_product_;
    }  

    delta_H_i_0_[grid_point_] = Hi_sum;

    if (verbose_)
        printf("Sum Hi = %g \n",Hi_sum);

}

// ----------------------------------------------------------------------
void ChemistryReactionSingle::calculate_K_c_i(int grid_point_)
{
    calculate_K_P_i(grid_point_);
    
    //Summation over all stoichometric factors
    double sum_stoich_total_ = 0.;
    //Loop reactants
    for (uint i_reactant = 0; i_reactant < species_stoich_reactant_.size(); i_reactant++)
    {
        sum_stoich_total_ += species_stoich_reactant_[i_reactant];
    } 
    //Loop products
    for (uint i_product_ = 0; i_product_ < species_stoich_product_.size(); i_product_++)
    {
        sum_stoich_total_ += species_stoich_product_[i_product_];
    }    


    double K_c_i_element = K_p_i_[grid_point_]*((P_atm_/(UNIVERSAL_GAS_CONSTANT*actual_temp_[grid_point_]))*exp(sum_stoich_total_));
    K_c_i_[grid_point_] = K_c_i_element;

    if (verbose_)
        printf("K_c_i_element = %g \n",K_c_i_element);
    
}

// ----------------------------------------------------------------------
void ChemistryReactionSingle::calculate_K_P_i(int grid_point_)
{
    double K_p_i_element = exp((delta_S_i_0_[grid_point_]/actual_temp_[grid_point_])-(delta_H_i_0_[grid_point_]/((UNIVERSAL_GAS_CONSTANT*actual_temp_[grid_point_]))));
    K_p_i_[grid_point_] = K_p_i_element;

    if (verbose_)
        printf("K_p_i_element = %g \n",K_p_i_element);
}

// ----------------------------------------------------------------------
void ChemistryReactionSingle::calculate_k_f_i(int _grid_point)
{
    //Arrhenius
    if(verbose_)
        printf("Grid point = %i, Acctual temp = %g, arrhenius_A_i = %g, arrhenius_beta_i_ = %g, -arrhenius_E_i_ = %g \n",
               _grid_point,actual_temp_[_grid_point],arrhenius_A_i_,arrhenius_beta_i_,arrhenius_E_i_ );

    k_f_i_[_grid_point] =  arrheniusRate(actual_temp_[_grid_point]);

    if (verbose_)
        printf("k_f_i_element = %g \n", k_f_i_[_grid_point]);

}


// ----------------------------------------------------------------------
void ChemistryReactionSingle::calculate_k_r_i(int grid_point_)
{
    double k_r_i_element = k_f_i_[grid_point_]/K_c_i_[grid_point_];
    k_r_i_[grid_point_]= k_r_i_element;

    if (verbose_)
        printf("k_r_i_element = %g \n",k_r_i_element);
}

// ----------------------------------------------------------------------
void ChemistryReactionSingle::correct_q_i(double factor, int grid_point_,int particleID_)
{
    q_dot_i_[grid_point_] *= factor;
}

// ----------------------------------------------------------------------
void ChemistryReactionSingle::calculate_q_i(int grid_point_,int particleID_)
{
    double species_value_     = 0.0;
    double pi_value_reactant_ = 1.0;
    double pi_value_product_  = 1.0;
    double exponent           = 0.0;
    vector<double>* exponentsReactant;
    vector<double>* exponentsProduct;

    if (verbose_)
        printf("*** ...computing reaction rates *** \n");

    if (elementary_reaction_)  
    {
        exponentsReactant = &species_stoich_reactant_;
        exponentsProduct  = &species_stoich_product_;  
    }
    else 
    {
        exponentsReactant = &reaction_order_reactant_;
        exponentsProduct  = &reaction_order_product_; 
    }

    for (uint i_reactant = 0; i_reactant < exponentsReactant->size(); i_reactant++)
    {
            particleData().setParticleIDPointer(original_IDs_reactant_[i_reactant],particleID_);
            particleData().returnDataPoint(species_value_,grid_point_);

            exponent            = (*exponentsReactant)[i_reactant];
            //compute product, however,
            //if reactant species is fully depleted and exponent is 0,
            // reaction must be forced to stop
            if( exponent==0 && species_value_< MINIMUM_CONCENTRATION ) 
                pi_value_reactant_ = 0.0;
            else
                pi_value_reactant_ *= MathExtraPascal::fastPow(species_value_,exponent);

            if (verbose_)
              printf("species_value_ = %g, pi_value_reactant_ =  %g, species_stoich_reactant_[%u] = %g \n",
                    species_value_, pi_value_reactant_,
                    i_reactant,(*exponentsReactant)[i_reactant]
                    );
     }

     for (uint i_product_ = 0; i_product_ < exponentsProduct->size(); i_product_++)
     {
            particleData().setParticleIDPointer(original_IDs_product_[i_product_],particleID_);
            particleData().returnDataPoint(species_value_,grid_point_);
            
            exponent           = (*exponentsProduct)[i_product_];
            pi_value_product_ *= MathExtraPascal::fastPow(species_value_,exponent);
           
            if (verbose_)
                printf("species_value_ = %g, pi_value_product_ =  %g, species_stoich_product_[%i] = %g \n",
                        species_value_, pi_value_product_, 
                        i_product_, (*exponentsProduct)[i_product_]
                      );
     }

    double q_i_element = (k_f_i_[grid_point_] * pi_value_reactant_)
                       - (k_r_i_[grid_point_] * pi_value_product_ );

    q_dot_i_[grid_point_] = q_i_element;

    if (verbose_)
        printf("q_i_element [%i] = %g \n",grid_point_,q_i_element);
    
}

// ----------------------------------------------------------------------
double ChemistryReactionSingle::integrateVolume_q_i(int particleID, double &particleVolume)
{
    //Integrate the volumetric reaction rate using spherical shells
    //centered around each grid point
    particleVolume = 0.0;
    double integral = 0.0;
    double ra(0.0), ri(0.0);
    double volume(0.0), deltaVolume(0.0);
    double dx = particleData().pRadius(particleID)/((particleMesh().nGridPoints())-1.0);
    for(int i=0; i<particleMesh().nGridPoints(); i++)
    {
        ra  = (double)i + 0.5; ra = min(ra,(double)(particleMesh().nGridPoints()-1));
        ra *= dx;
        ri  = (double)i - 0.5; ri = max(ri,0.0);
        ri *= dx;
        
        deltaVolume  = (ra*ra*ra-ri*ri*ri)*4.188790205; //4.188790205=4/3*pi
        particleVolume += deltaVolume;
        integral       += deltaVolume*q_dot_i_[i];
    } 
    
    return integral;
}

// ----------------------------------------------------------------------
void ChemistryReactionSingle::empty_vectors()
{
    delta_H_k_0_reactant_.clear();
    delta_S_k_0_reactant_.clear();
    delta_H_k_0_product_.clear();
    delta_S_k_0_product_.clear();
    reaction_order_reactant_.clear();
    reaction_order_product_.clear();
    original_IDs_product_.clear();
    original_IDs_reactant_.clear(); 
}


// ----------------------------------------------------------------------
ChemistryReactionSingle::~ChemistryReactionSingle()
{

}


