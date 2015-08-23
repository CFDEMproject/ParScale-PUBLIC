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

#include "model_phasechange_equilibrium_sharp.h"
#include "particle_data.h"
#include "control.h"
#include "simulation_state.h"

//#define verbose_ FALSE
using namespace PASCAL_NS;
#define MIN_DENSE_PHASE_TOEVAPORATE  1e-16
#define MIN_DILUTE_PHASE_TOEVAPORATE 1e-16

ModelPhaseChangeEquilibrium_sharp::ModelPhaseChangeEquilibrium_sharp(ParScale *ptr, char *name):
ModelPhaseChange(ptr,name),
update_phase_fraction_(false)

{    
    ptr_ = ptr;

    input().openJsonFile("settings", "verbose", "verbose", global_properties_);
    if(global_properties_["ModelPhaseChangeEquilibriumSharp"].isNull())
        error().throw_error_one(FLERR,"ERROR: could not read verbose settings for this class. \n",
                                "ModelPhaseChangeEquilibriumSharp");

    verbose_   = global_properties_["ModelPhaseChangeEquilibriumSharp"].toBool();    

    

}
// *************************************************************************************
ModelPhaseChangeEquilibrium_sharp::~ModelPhaseChangeEquilibrium_sharp()
{
}

// *************************************************************************************
void ModelPhaseChangeEquilibrium_sharp::init(int narg, char const* const* arg, int eqnType, int modelPhaseChangeID)
{
    
    ModelPhaseChange::init(narg,arg,eqnType,modelPhaseChangeID);

    for (int i_ = 0;i_ < modelEqnContainer().nrEqns(); i_++)
    {
        if(modelEqnContainer().modelEqn(i_)->updatePhaseFraction)
            update_phase_fraction_ = modelEqnContainer().modelEqn(i_)->updatePhaseFraction;
            
        if(verbose_)
            printf("Model eqn id for updating phase fraction is %i \n",i_);
    }

    if(!update_phase_fraction_)
        error().throw_error_one(FLERR,"For Phase Change model equilibrium sharp you have to have at least one species for which you update the phase fraction",
                                "ModelPhaseChangeChangeEquilibriumSharp");
        
}

// *************************************************************************************
void ModelPhaseChangeEquilibrium_sharp::begin_of_step()
{
    ModelPhaseChange::begin_of_step();

    readEquilibriumSharpJSON();

    for(int particleID=0; particleID<particleData().nbody(); particleID++)
    {
        particleData().setParticleIDPointerPhaseFraction(particleID);

        if(verbose_) 
            printf("setting evaporation rate for particle %i with heat eqn id %i \n", particleID, heatEqnID_); 

        if(!isIsoThermal_)
        {
            particleData().setParticleIDPointer(heatEqnID_,particleID);	
            particleData().returnIntraData(tempIntraDataHeat_);
        }

        //fill in the source
        for (int grid_point = 0; grid_point < particleMesh().nGridPoints(); grid_point++)
        {

                double rateDensePhase    = 0;
                double rateDilutePhase   = 0;
                double jacobiDensePhase  = 0;
                double jacobiDilutePhase = 0;

                //Retrieve the local species concentration        
                particleData().retrieveIntraData(speciesParticleDataID_, 
                                                 particleID, grid_point, 
                                                 species_concentrations_);

                /*if(phaseID_[0] == 0) //if dense is solid
                {
                    double temp_phaseFracGas = 0.0;
                    double temp_phaseFracLiq = 0.0;

                    particleData().returnPhaseFractionDataGridpoint(temp_phaseFracGas,temp_phaseFracLiq,grid_point);
                    //printf("Phase Fraction 1 = %g, Phase Fraction 2 = %g \n",temp_phaseFracGas,temp_phaseFracLiq); 

                    phaseFractions_[0]=1.0-temp_phaseFracGas-temp_phaseFracLiq;

                    if(phaseID_[1] == 1)                        //dilute is gas 
                        phaseFractions_[1] = temp_phaseFracGas;
                    else                                        //dilute must be liquid
                        phaseFractions_[1] = temp_phaseFracLiq;
                }*/

                if(phaseID_[1] == 0) //if dilute is solid - for equilibrium sharp it has to be!
                {
                    double temp_phaseFracGas = 0.0;
                    double temp_phaseFracLiq = 0.0;

                    particleData().returnPhaseFractionDataGridpoint(temp_phaseFracGas,temp_phaseFracLiq,grid_point);           
                    //printf("Phase Fraction 1 = %g, Phase Fraction 2 = %g \n",temp_phaseFracGas,temp_phaseFracLiq); 

                    phaseFractions_[1]=1.0-temp_phaseFracGas-temp_phaseFracLiq;     

                    if(phaseID_[0] == 1 || phaseID_[0] == 0)    //dilute is gas or solid 
                        error().throw_error_one(FLERR,"For Phase Change model equilibrium sharp your dense phase needs to be a liquid species but it is gasous or solid.",
                                "ModelPhaseChangeChangeEquilibriumSharp");
                    else                                        //dilute must be liquid
                        phaseFractions_[0] = temp_phaseFracLiq;
                }
                else
                {
                    error().throw_error_one(FLERR,"For Phase Change model equilibrium sharp your dilute phase needs to be solid species!",
                                "ModelPhaseChangeChangeEquilibriumSharp");
                    /*particleData().retrievePhaseFractionData(phaseID_, 
                                                     particleID, grid_point,
                                                     phaseFractions_);*/
                }

                if(verbose_) 
                    printf("Phase Fraction 1 = %g, Phase Fraction 2 = %g \n",phaseFractions_[0],phaseFractions_[1]);

                double densephaseFactor  = 0;
                double dilutephaseFactor = 0;
                if(phaseFractions_[0]>MIN_DENSE_PHASE_TOEVAPORATE)
                    densephaseFactor = 1.0;

                if(phaseFractions_[1]>MIN_DILUTE_PHASE_TOEVAPORATE) //boiling or have some gas
                    dilutephaseFactor = 1.0;


                if(verbose_)
                    printf("sat_concentration = %g \n",sat_concentration_liq_);
   
                delta_concentration_ = species_concentrations_[0] - sat_concentration_liq_;
                
                if(verbose_)
                    printf("delta_concentration = %g,species_concentrations_[1] = %g \n",delta_concentration_,species_concentrations_[0]);

                double time_step_ = control().simulationState().deltaT();

                if(delta_concentration_ >= 0.0 && species_concentrations_[1] <= sat_concentration_solid_)
                {
                    double nano_rate_ = 1/((sat_concentration_solid_-species_concentrations_[1])/(species_concentrations_[1]));
                    //printf("Nano rate for gridpoint %i = %g ,sat_concentration_solid_ = %g, species_concentrations_[1]= %g \n",grid_point,nano_rate_,sat_concentration_solid_,species_concentrations_[1]);
                    rateDilutePhase = (delta_concentration_/time_step_)/1.0;    //Species is added to dilute phase
                                                                                //No multiplication by phase fraction since dilute phase is solid phase [kmol/m³_tot]
                    rateDensePhase  = -1.0 * phaseFractions_[0] * (delta_concentration_/time_step_); //multiply with -1 since dense phase is falling out
                    //rateDensePhase  = -1.0 * nano_rate_ * phaseFractions_[0] * (delta_concentration_/time_step_); //multiply with -1 since dense phase is falling out
                    
                }
                else //species concentration is smaller than saturation concentration
                {
                    rateDilutePhase = 0.0;
                    rateDensePhase  = 0.0;
                }

                particleData().savePhaseChangeParticleData(   phaseID_[0], particleID, rateDensePhase,   grid_point);
                particleData().savePhaseChangeParticleDataJac(phaseID_[0], particleID, jacobiDensePhase, grid_point);

                particleData().savePhaseChangeParticleData(   phaseID_[1], particleID, rateDilutePhase,  grid_point);
                particleData().savePhaseChangeParticleDataJac(phaseID_[1], particleID, jacobiDilutePhase,grid_point);

                if(verbose_)
                    printf("Grid point = %i, Rate dense phase = %g, rate dilute phase = %g, timestep = %g \n",grid_point,rateDensePhase,rateDilutePhase,time_step_); 

        } //End gridpoint loop

    if(verbose_)
        printf("\n");
    } //End particle loop
    
}

// *************************************************************************************
void ModelPhaseChangeEquilibrium_sharp::readEquilibriumSharpJSON()
{
    QJsonObject  myJSONDict = readQJsonObject( name(), "SaturationConcentrationLiquid");   
    if(myJSONDict["Concentration"].isNull())
        error().throw_error_one(FLERR,"For Phase Change model equilibrium sharp you have to set a saturation concentration for the liquid dense species",
                                "ModelPhaseChangeChangeEquilibriumSharp");    
    sat_concentration_liq_ = myJSONDict["Concentration"].toDouble();


    myJSONDict = readQJsonObject( name(), "SaturationConcentrationSolid");   
    if(myJSONDict["Concentration"].isNull())
        error().throw_error_one(FLERR,"For Phase Change model equilibrium sharp you have to set a saturation concentration for the solid dilute species",
                                "ModelPhaseChangeChangeEquilibriumSharp"); 
    sat_concentration_solid_ = myJSONDict["Concentration"].toDouble();
}
