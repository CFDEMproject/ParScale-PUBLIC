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

#include "model_phasechange_evaporation.h"
#include "particle_data.h"
#include "control.h"
#include "model_eqn_container.h"

//#define verbose_ FALSE
using namespace PASCAL_NS;
#define MIN_DENSE_PHASE_TOEVAPORATE  1e-16
#define MIN_DILUTE_PHASE_TOEVAPORATE 1e-16

ModelPhaseChangeEvaporation::ModelPhaseChangeEvaporation(ParScale *ptr, char *name):
ModelPhaseChange(ptr,name)
{
    pressureAmbient_ = 1e5;
    evaporationTemp_ = 273;

    ptr_ = ptr;

    input().openJsonFile("settings", "verbose", "verbose", global_properties_);
    if(global_properties_["ModelPhaseChangeEvaporation"].isNull())
        error().throw_error_one(FLERR,"ERROR: could not read verbose settings for this class. \n",
                                "ModelPhaseChangeEvaporation");

    verbose_   = global_properties_["ModelPhaseChangeEvaporation"].toBool();    

}
// *************************************************************************************
ModelPhaseChangeEvaporation::~ModelPhaseChangeEvaporation()
{
}

// *************************************************************************************
void ModelPhaseChangeEvaporation::init(int narg, char const* const* arg, int eqnType, int modelPhaseChangeID)
{
    
    ModelPhaseChange::init(narg,arg,eqnType,modelPhaseChangeID);
    readEvaporationJSON();
}

// *************************************************************************************
void ModelPhaseChangeEvaporation::begin_of_step()
{
    ModelPhaseChange::begin_of_step();
    double pSat = 0;

    if( control().simulationState().deltaT() > (100/evaporationRateConstant_))
        error().throw_error_one(FLERR,"ERROR: evaporationRateConstant is too large for this time step (i.e., rate is > 0.1/deltaT). Reduce time step or evaporationRateConstant. \n",
                                "ModelPhaseChangeEvaporation");

    double phaseFractionMinimumConvection = modelEqnContainer().modelEqn(speciesModelEqnID_[1])->phaseFractionMinimumConvection;
    double convectionBoundMinStefanFactor = modelEqnContainer().modelEqn(speciesModelEqnID_[1])->convectionBoundMinStefanFactor;

    for(int particleID=0; particleID<particleData().nbody(); particleID++)
    {
       
        if(verbose_) 
            printf("setting evaporation rate for particle %i with heat eqn id %i \n", particleID, heatEqnID_); 


        double currentTemp = evaporationTemp_;
        pSat               = (this->*vaporPressure)(currentTemp);

        if(!isIsoThermal_)
        {
            particleData().setParticleIDPointer(heatEqnID_,particleID);	
            particleData().returnIntraData(tempIntraDataHeat_);
        }

        //fill in the source
        for (int grid_point = 0; grid_point < particleMesh().nGridPoints(); grid_point++)
        {

                if(!isIsoThermal_)
                {   
                    currentTemp = tempIntraDataHeat_[grid_point];
                    pSat        = (this->*vaporPressure)(currentTemp);
                }
 
                double rateDensePhase    = 0;
                double rateDilutePhase   = 0;
                double jacobiDensePhase  = 0;
                double jacobiDilutePhase = 0;

                //Retrieve the local species concentration        
                particleData().retrieveIntraData(speciesParticleDataID_, 
                                                 particleID, grid_point, 
                                                 species_concentrations_);

                if(phaseID_[0] == 0 || phaseID_[1] == 0) //if dense is solid
                {
                    error().throw_error_one(FLERR,"No phase change of solid phase possilbe if phase changemodel is evaporation \n", "ModelPhaseChangeEvaporation");
                }

                particleData().retrievePhaseFractionData(phaseID_, 
                                                 particleID, grid_point,
                                                 phaseFractions_);
                double densephaseFactor  = 0;
                double dilutephaseFactor = 0;
                if(phaseFractions_[0]>MIN_DENSE_PHASE_TOEVAPORATE)
                    densephaseFactor = 1.0;

                if(pSat>pressureAmbient_ || phaseFractions_[1]>MIN_DILUTE_PHASE_TOEVAPORATE) //boiling or have some gas
                    dilutephaseFactor = 1.0;


                //Must save LINEARIZED source term due to phase change    
                jacobiDilutePhase= -1*evaporationRateConstant_ // phaseFractions_[1]
                                   *densephaseFactor*dilutephaseFactor;
                rateDilutePhase  = -1*jacobiDilutePhase * pSat * INVERSE_UNIVERSAL_GAS_CONSTANT / currentTemp;

                rateDensePhase   = rateDilutePhase 
                                 - evaporationRateConstant_ * species_concentrations_[1]  // phaseFractions_[1]
                                  *densephaseFactor*dilutephaseFactor;
                rateDensePhase  *= -1 / phaseFractions_[1]; //multiply with -1/phaseFraction since dense phase is evaporating
                                                            //and dilutephase rate is normalized with phasefraction

                //Jacobi that would yield zero evaporation for dense at zero phase fraction
                jacobiDensePhase  = -1*rateDensePhase/(phaseFractions_[0]+MIN_DENSE_PHASE_TOEVAPORATE);
                jacobiDensePhase *= jacobiForDensePhase_;

                particleData().savePhaseChangeParticleData(   phaseID_[0], particleID, rateDensePhase,   grid_point);
                particleData().savePhaseChangeParticleDataJac(phaseID_[0], particleID, jacobiDensePhase, grid_point);

                particleData().savePhaseChangeParticleData(   phaseID_[1], particleID, rateDilutePhase,  grid_point);
                particleData().savePhaseChangeParticleDataJac(phaseID_[1], particleID, jacobiDilutePhase,grid_point);

                //Save the old temperature
                particleData().savePhaseChangeParticleDataJacLastSlot(particleID, currentTemp, grid_point);

                if(verbose_) 
                    printf("species_concentrations_: %g/%g, phaseID_: %d/%d, pSat: %g, rateDense/Dilute: %g/%g, jacobi dense/dilute: %g/%g, oldTemperature: %g. \n", 
                            species_concentrations_[0], species_concentrations_[1],
                            phaseID_[0], phaseID_[1],
                            pSat, 
                            rateDensePhase,   rateDilutePhase,
                            jacobiDensePhase, jacobiDilutePhase, currentTemp
                          ); 
        }
    }
    
}

// *************************************************************************************
void ModelPhaseChangeEvaporation::readEvaporationJSON()
{

    QJsonObject  myJSONDict = readQJsonObject( name(), "vaporPressureModel");
    QString      myModel    = myJSONDict["name"].toString();
    QJsonArray   myParams   = myJSONDict["parameters"].toArray();
    if(myJSONDict["evaporationRateConstant"].isNull())
        error().throw_error_one(FLERR,"ERROR: could not read 'evaporationRateConstant' from JSON file for this class. \n",
                                "ModelPhaseChangeEvaporation");
    evaporationRateConstant_   = myJSONDict["evaporationRateConstant"].toDouble();
    jacobiForDensePhase_       = myJSONDict["jacobiForDensePhase"].toDouble();

    if(!(myJSONDict["evaporationTemp"].isNull()))
        evaporationTemp_    = myJSONDict["evaporationTemp"].toDouble();

    if(!(myJSONDict["pressureAmbient"].isNull()))
        pressureAmbient_    = myJSONDict["pressureAmbient"].toDouble();

    particleData().referencePressure = pressureAmbient_;

    if(strcmp(qPrintable(myModel),"antoine") == 0 )
    {

        vaporPressure=&ModelPhaseChangeEvaporation::vaporPressureAntoine;
        evaporationDerivativeConcTemp=&ModelPhaseChangeEvaporation::derivativeConcTempAntoine;
        QJsonArray::iterator i;
        for (i = myParams.begin(); i != myParams.end(); ++i)
        {
            vaporPressureParameters_.push_back((*i).toDouble());
        }
        printf("ModelPhaseChangeEvaporation: have read model '%s', with %d parameters. \n \n",
               qPrintable(myModel), vaporPressureParameters_.size());

        if(vaporPressureParameters_.size()!=3)
            error().throw_error_one(FLERR,"ERROR: You have not provided exactly 3 parameters for the model 'antoine'.");
    }
    else
         error().throw_error_one(FLERR,"ERROR: please specify the name of an implemented vaporPressureModel.");


    myJSONDict = readQJsonObject( name(), "enthalpyModel");
    myModel    = myJSONDict["name"].toString();

    if(strcmp(qPrintable(myModel),"constant") == 0 )
    {
        evaporationDeltaH=&ModelPhaseChangeEvaporation::evaporationDeltaHConstant;
        evaporationDeltaHParameters_.push_back(myJSONDict["value"].toDouble());
    }
    else
         error().throw_error_one(FLERR,"ERROR: please specify the name of an implemented enthalpyModel.");



}

// *************************************************************************************
double ModelPhaseChangeEvaporation::vaporPressureAntoine(double temperatureK) const
{
    //Parameters must be for pressure to be in mmHg and temperature in Kelvin
    //Result is then in Pascal

    return 133.32236 * pow(10.0,                        //convert from mmHg to Pa
                              vaporPressureParameters_[0] 
                         -    vaporPressureParameters_[1] 
                           / (vaporPressureParameters_[2]+temperatureK)
                          );
}

// *************************************************************************************
double ModelPhaseChangeEvaporation::derivativeConcTempAntoine(double temperatureK) const
{
    //Parameters must be for pressure to be in mmHg and temperature in Kelvin
    //Result is then in Pascal
    //printf("Temperature = %g \n",temperatureK);
    double pSat = vaporPressureAntoine(temperatureK);
    //printf("Saturation Pressure = %g \n",pSat);
    double inverseTemp = 1/temperatureK;

    return  pSat * INVERSE_UNIVERSAL_GAS_CONSTANT * inverseTemp
          * (
                 2.30258509299 //ln(10) 
                -inverseTemp
            );
}
// *************************************************************************************
double ModelPhaseChangeEvaporation::evaporationDeltaHConstant(double empty) const
{
    return evaporationDeltaHParameters_[0];
}
// *************************************************************************************
