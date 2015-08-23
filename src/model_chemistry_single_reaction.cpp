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


#include "model_chemistry_single_reaction.h"

using namespace PASCAL_NS;

ModelChemistrySingleReaction::ModelChemistrySingleReaction(ParScale *ptr, char *name):
ModelChemistry(ptr,name),
reaction_(NULL),
grainModel_(NULL),
isIsoThermal_(true)
{

    ptr_ = ptr;

    input().openJsonFile("settings", "verbose", "verbose", global_properties_ );
    if(global_properties_["ModelChemistrySingleReaction"].isNull())
        error().throw_error_one(FLERR,"ERROR: could not read verbose settings for this class. \n",
                                "ModelChemistrySingleReaction");
    verbose_   = global_properties_["ModelChemistrySingleReaction"].toBool();    


}
// *************************************************************************************
ModelChemistrySingleReaction::~ModelChemistrySingleReaction()
{
    delete reaction_;
    delete grainModel_;
}

// *************************************************************************************
void ModelChemistrySingleReaction::init(int narg, char const* const* arg, int eqnType, int modelChemistryID)
{
    
    ModelChemistry::init(narg,arg,eqnType,modelChemistryID);
    read_chemisty_JSON();
   
    //Error Checks
    if( (heatEqnID_<0) && !isIsoThermal_)
        error().throw_error_one(FLERR,"ERROR: no heat equation found, but you like to run a non-isothermal simulations. This is impossible. \n");

    if( (heatEqnID_>=0) && isIsoThermal_)
        error().throw_error_one(FLERR,"ERROR: heat equation found, but you like to run an isothermal simulations. This is not meaningful. Deactivate the heat equation! \n");

}

// *************************************************************************************
void ModelChemistrySingleReaction::begin_of_step()
{

    int nrHeatEqn = modelEqnContainer().nrHeatEqns();

    //TODO: Must reset the original IDs for now, check in future to avoid this
    reaction_->resetOriginalIDs();
    reaction_->returnOriginalIDsreactant(original_IDs_reactant_);

    //Attention: first model Eqn is used to decide how to handle heat of reaction
    int     modelApproachFirst = modelEqnContainer().modelEqn
                                 (original_IDs_reactant_[0])->modelingApproach;


    for(int particleID=0; particleID<particleData().nbody(); particleID++)
    {
        int frontGridPoint     = 0; //for SHRINKINGCORE 
        int frontConcentration = 0; //for SHRINKINGCORE 
        
        if(verbose_) 
            printf("Now particle %i with heat eqn id %i \n", particleID, heatEqnID_); 


        if(isIsoThermal_)
        {
            double sourceSpecIConst = reaction_->arrheniusRate(temperature_);
            
            if (verbose_)
                printf("SINGLE ISOTHERM REACT: raction rate is = %g \n",sourceSpecIConst);

            for (int grid_point = 0; grid_point < particleMesh().nGridPoints(); grid_point++)
            {
              double source = sourceSpecIConst;
              if(modelApproachFirst==SHRINKINGCORE)
              {
                    //No action
              }
              else if(modelApproachFirst==CONTINUUM)
              {
                //Retrieve the local species concentration        
                particleData().retrieveIntraData(original_IDs_reactant_, 
                                                 particleID, grid_point, 
                                                 species_concentrations_);
                
                source *= reaction_->reactionProduct( species_concentrations_,
                                                      species_order_)
                        * grainModel_->f(species_concentrations_[grainModel_->solidID()]); 

                if (verbose_)
                {
                    printf("...retrieved concentrations. grainModel->f: %g \n",
                            grainModel_->f(species_concentrations_[grainModel_->solidID()]) 
                          );
                    for(int iC=0;iC<species_concentrations_.size();iC++)
                         printf("%g  \n", species_concentrations_[iC] );
                }
              }
              else
                   error().throw_error_one(FLERR,"ERROR: reaction handling for this modelingApproach not implemented \n");
              
              for (uint idSpec=0; idSpec < species_stoich_.size(); idSpec++)
              {
                    double sourceSpecI   = source;
                    double jacobianSpecI = 0.0;
                    int    speciesDataId = original_IDs_reactant_[idSpec];
                    sourceSpecI *= species_stoich_[idSpec];

                    if(modelApproachFirst==SHRINKINGCORE)
                    {
                         //no action
                    }
                    else if(modelApproachFirst==CONTINUUM)
                    {
                         sourceSpecI *= molar_mass_[idSpec];
                         //Jacobian (derivative of grain model IS NOT considered here, since no trivial analytical approach)
                         jacobianSpecI =  sourceSpecI 
                                         *  species_order_[idSpec]
                                         / (species_concentrations_[idSpec] + 1e-64);
                    }

                    if (verbose_)
                      printf("SINGLE ISOTHERM REACT: source = %g / jac: %g (stoich: %g) for Particle %i chem species ID %i, modelApproach: %i and grid Point %i \n",
                                sourceSpecI, 
                                
                                species_stoich_[idSpec],
                                particleID, 
                                speciesDataId, modelApproachFirst, 
                                grid_point
                              );

                    double linearizedSpeciesD = sourceSpecI - jacobianSpecI * species_concentrations_[idSpec];  
                    particleData().saveChemistryParticleData(
                                   speciesDataId,
                                   particleID,
                                   linearizedSpeciesD, //must save d of linear approximation!
                                   grid_point);
                    particleData().saveChemistryJacParticleData(
                                   speciesDataId,
                                   particleID,
                                   jacobianSpecI,
                                   grid_point);
                
              }
            } //end grid_point loop
        }

        else  //non-isothermal Case
        {
            particleData().setParticleIDPointer(heatEqnID_,particleID);	
    	    particleData().returnIntraData(tempIntraDataHeat_);
            reaction_->SetActualTemp(tempIntraDataHeat_);

            for (int grid_point = particleMesh().nGridPoints()-1; grid_point >= 0 ; grid_point--) 
            //MUST loop from back to front, to avoid problems with frontGridPoint used in SHRINKINGCORE
            {
                reaction_->calculate_k_f_i(grid_point);
                reaction_->returnKfi(k_f_i_);
                reaction_->calculate_q_i(grid_point, particleID);

                //retrieve concentrations since needed in Jacobian calculation
                particleData().retrieveIntraData(original_IDs_reactant_, 
                                                 particleID, grid_point,
                                                 species_concentrations_
                                                 );

                //apply grain model if active
                if(grainModel_->active())
                {

                    reaction_->correct_q_i(
                                grainModel_->f(species_concentrations_[grainModel_->solidID()]),
                                grid_point, particleID
                                          );
                    if (verbose_)
                    {
                        printf("...retrieved concentrations. grainModel->f: %g \n",
                                grainModel_->f(species_concentrations_[grainModel_->solidID()]) 
                              );
                        for(int iC=0;iC<species_concentrations_.size();iC++)
                             printf("%g  \n", species_concentrations_[iC] );
                    }
                }
                
                reaction_->returnQDoti(q_i_dot_);

                // 1 - Fill source terms for each involved species
                for (uint idSpec=0; idSpec < species_stoich_.size(); idSpec++)
                {
                    double sourceSpecI   = 0.0;
                    double jacobianSpecI = 0.0;
                    int    speciesDataId = original_IDs_reactant_[idSpec];
                    int    modelApproach = modelEqnContainer().modelEqn(speciesDataId)->modelingApproach;
                    if(modelApproach==SHRINKINGCORE)
                    {
                        //Determine grid point next to front
                        double rCore = 0.0;
                        rCore = particleData().retrieveIntraData(speciesDataId,particleID,0);
                        frontGridPoint =max(0, 
                                               min( particleMesh().nGridPoints()-1,
                                                (int)floor(
                                                          (double)(particleMesh().nGridPoints())
                                                        * rCore
                                                        )
                                               )
                                           );
                                           
                        //TODO: compute concentration on front with fluid concentrations
                        //and intra-partile resistances in SHRINKINGCORE model
                        frontConcentration = 0.0;
                                //this is the fluid concentration:
                                //particleData().retrieveIntraData(speciesDataId,particleID,XXX);
                        
                        if(grid_point==0) //just set for 0-th grid point --> radius of SCM
                        {
                            sourceSpecI = species_stoich_[idSpec]
                                        * k_f_i_[frontGridPoint];
                        }
                               
                        if (verbose_)
                            printf("...variables @ front of SCM: Rate = %g for Particle %i species data id: %i, front grid Point %i, rCore: %g \n",
                                sourceSpecI, particleID, 
                                speciesDataId,  
                                frontGridPoint,rCore
                                  );
                    }
                    else if(modelApproach==CONTINUUM)
                    {
                        sourceSpecI = species_stoich_[idSpec]
                               * molar_mass_[idSpec]
                               * q_i_dot_[grid_point];   

                        //Jacobian (derivative of grain model IS NOT considered here, since no trivial analytical approach)
                        jacobianSpecI =  sourceSpecI 
                                      *  species_order_[idSpec]
                                      / (species_concentrations_[idSpec] + 1e-64);
                 
                    }
                    else
                        error().throw_error_one(FLERR,"ERROR: reaction handling for this modelingApproach not implemented \n");



                    
                    if (verbose_)
                        printf("SINGLE REACT: Rate = %g / jac: %g for Particle %i chem species ID %i, data id: %i, modelApproach: %i and grid Point %i \n",
                                sourceSpecI, jacobianSpecI,
                                particleID, 
                                idSpec, speciesDataId, modelApproach, 
                                grid_point
                              );

                    //Save to data array 
                    double linearizedSpeciesD = sourceSpecI - jacobianSpecI * species_concentrations_[idSpec];  
                    particleData().saveChemistryParticleData(
                                   speciesDataId,
                                   particleID,
                                   linearizedSpeciesD, //must save d of linear approximation!
                                   grid_point);
                    particleData().saveChemistryJacParticleData(
                                   speciesDataId,
                                   particleID,
                                   jacobianSpecI,
                                   grid_point);
                } //End species loop
                
                // 2 - Fill source terms for heat equations
                {
                    double heatSource   = 0.0;
                    double jacobianHeat = 0.0;
                    if(modelApproachFirst==SHRINKINGCORE)
                    {
                        if(grid_point==frontGridPoint) //Just set source at the front!
                        {
                            //TODO: check frontConcentration normalize appropriately
                            //in order to return correct volumetric heat source
                            heatSource =-k_f_i_[frontGridPoint]
                                       * frontConcentration
                                       * deltaH_r;
                        }
                    }
                    else if(modelApproachFirst==CONTINUUM)
                    {
                        heatSource =-q_i_dot_[grid_point]
                                   * deltaH_r;
                        jacobianHeat = heatSource 
                                     * reaction_->computeArrheniusRateDerivativePrefactor(grid_point);
                    }
                    else
                        error().throw_error_one(FLERR,"ERROR: reaction handling for this modelingApproach not implemented \n");

                    

                    if (verbose_)
                        printf("SINGLE REACT: heat source = %g for Particle %i, modelApproachFirst: %i. \n",
                                heatSource, particleID, 
                                modelApproachFirst
                              );  

                    //Save to data array                     
                    double linearizedHeatD = heatSource - jacobianHeat * tempIntraDataHeat_[grid_point];                         
                    particleData().saveChemistryParticleData(
                                   heatEqnID_,
                                   particleID,
                                   linearizedHeatD, //must save d of linear approximation!
                                   grid_point);
                    particleData().saveChemistryJacParticleData(
                                   heatEqnID_,
                                   particleID,
                                   jacobianHeat,
                                   grid_point);
                }
            } //end loop over grid points
            
            //Compute growth
            if(growthMolesFormed_!=0.0)
            {
                double particleVolume(0.0);
                double molarGain = growthMolesFormed_*
                                   reaction_->integrateVolume_q_i(particleID,particleVolume);
                double newParticleVolume = particleVolume
                                         + molarGain
                                         * growthProductMolarWeight_
                                         / growthProductDensity_;
                                         
                //Update particle properties and grow particles
                particleData().normalizeInternalFields(particleID, particleVolume/newParticleVolume);
                particleData().pRadius(particleID) =
                            pow( newParticleVolume
                                *0.238732415 // 0.238732415 = 3/(4pi)
                                ,
                                0.333333333333333);
                if (verbose_)
                    printf("SINGLE REACT GROWTH: molarGain/massGain/volumeGain = %g/%g/%g, old pVolume: %g, new pRadius: %g \n",
                                molarGain,
                                molarGain*growthProductMolarWeight_,
                                molarGain*growthProductMolarWeight_/growthProductDensity_,
                                particleVolume, particleData().pRadius(particleID)
                          );
            }
        }// end non-isothermal
    } //end particle loop

}

// *************************************************************************************
void ModelChemistrySingleReaction::read_chemisty_JSON()
{
    //init new chemistry reaction
    reaction_   = new ChemistryReactionSingle(ptr_);   
    grainModel_ = new ChemistryGrainModel(ptr_); 

    QJsonObject myJSONDict = readQJsonObject("chemistry_single_reaction", "reactants");
    QJsonArray  myReactants = myJSONDict["names"].toArray();
    if(myReactants.size()==0)
        error().throw_error_one(FLERR,"ERROR: names not found in 'species_properties' file. \n");

    QJsonArray::iterator i;
    deltaH_r = 0.0;
    for (i = myReactants.begin(); i != myReactants.end(); ++i)
    {
        //Names
        species_names_.push_back((*i).toString().toStdString());
        reaction_->SetSpeciesNamesReactant(species_names_.back());
    
        //Molar Mass
        double myMass;  
        read_species_properties(species_names_.back().c_str(),"molarMass", &myMass);
        molar_mass_.push_back(myMass);

        //Reaction Order & Stoichiometry
        myJSONDict     = readQJsonObject("chemistry_single_reaction", species_names_.back().c_str());
        if(myJSONDict["reactionOrder"].isNull() || myJSONDict["stoichiometry"].isNull())
        error().throw_error_one(FLERR,"ERROR for a species: 'reactionOrder' or 'stoichiometry' not found in file. \n",
                              species_names_.back().c_str(),
                              "settings/chemistry_single_reaction.json");

        species_order_.push_back(myJSONDict["reactionOrder"].toDouble());
        species_stoich_.push_back(myJSONDict["stoichiometry"].toDouble());
        reaction_->SetSpeciesStoichReactant(species_stoich_.back());
                
        //Reaction enthalyp
        double myEnthalpy;
        read_species_properties(species_names_.back().c_str(),"H_i", &myEnthalpy);
        deltaH_r += myEnthalpy*species_stoich_.back();
    }

    species_concentrations_.assign(species_order_.size(),0.0); //set concentrations to zero
    reaction_->SetReactionOrderReactant(species_order_);

    read_chemistry_single_react_json_file("isIsoThermal",  &isIsoThermal_, true);
    read_chemistry_single_react_json_file("cMinimum",      &cMinimum_,     false);
    read_chemistry_single_react_json_file("Arrhenius_A",   &arrhenius_A_,  true);
    read_chemistry_single_react_json_file("Arrhenius_beta",&arrhenius_beta_, true);
    read_chemistry_single_react_json_file("Arrhenius_E_A", &arrhenius_E_A_, true);
    
    read_chemistry_single_react_json_file("growthMolesFormed", &growthMolesFormed_, true);
    if(growthMolesFormed_!=0.0)
    {
        read_chemistry_single_react_json_file("growthProductMolarWeight", &growthProductMolarWeight_, true);
        read_chemistry_single_react_json_file("growthProductDensity", &growthProductDensity_, true);
    }

    temperature_ = -1;
    if(isIsoThermal_)
        read_chemistry_single_react_json_file("temperature",  &temperature_, true);

    //Fill the reaction 
    reaction_->setArrheniusConstants(arrhenius_A_,arrhenius_beta_,arrhenius_E_A_);
    reaction_->setCMinimum(cMinimum_);
    reaction_->SetElementaryReaction(false);
    
    //Init grain model
    grainModel_->init();
        
    //Report Status
    reaction_->printStatus();
    grainModel_->printStatus();
}
// *************************************************************************************
