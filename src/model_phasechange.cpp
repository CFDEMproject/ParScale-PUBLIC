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


#include "model_phasechange.h"
#include "model_eqn_container.h"

using namespace PASCAL_NS;

ModelPhaseChange::ModelPhaseChange(ParScale *ptr, char *name):
ModelBase(ptr, name),
verbose_(false),
isActive_(true),
isSet_(false),
isIsoThermal_(true),
heatEqnID_(-1),
tempIntraDataHeat_(NULL),
tempIntraDataSpecies_(NULL)
{
    create<double>(tempIntraDataHeat_,    particleMesh().nGridPoints() );
    create<double>(tempIntraDataSpecies_, particleMesh().nGridPoints() );
}
// *************************************************************************************
ModelPhaseChange::~ModelPhaseChange()
{
    destroy<double>(tempIntraDataHeat_);
    destroy<double>(tempIntraDataSpecies_);
}

// *************************************************************************************
void ModelPhaseChange::init(int narg, char const* const* arg, int eqnType, int modelPhaseChangeID)
{
    //Check if heatEqn is found and can be set
    if(modelEqnContainer().nrHeatEqns() > 0)
        heatEqnID_ = 0; //take the first heat equation for the chemistry!
    else
    {
        printf("********ModelPhaseChange::WARNING: you do not have a heat equation. \n");
        printf("********Can only run in iso-thermal mode for ALL involved particles! \n\n");
    }

    readPhaseChangeBasicsJSON();
}

// *************************************************************************************
void ModelPhaseChange::readPhaseChangeBasicsJSON()
{

    QJsonObject myJSONDict = readQJsonObject( name(), "involvedPhases");
    myPhases_  = myJSONDict["phaseNames"].toArray();
    mySpecies_ = myJSONDict["speciesNames"].toArray();
    if(myPhases_.size()!=2 )
        error().throw_error_one(FLERR,"ERROR: phaseNames not found in file, or number of phaseNames is not exactly 2. \n");
    if(mySpecies_.size()!=2)
        error().throw_error_one(FLERR,"ERROR: speciesNames not found in file, or number of speciesNames is not exactly 2. \n");

}

// *************************************************************************************
void ModelPhaseChange::setPhaseChangeBasics()
{
    if(isSet_) return;

    for(int iPhase=0;iPhase<2; iPhase++) //loop involved phases
    {
        if( strcmp(qPrintable(myPhases_[iPhase].toString()),"gas") == 0 )
            phaseID_.push_back(GAS);
        else if( strcmp(qPrintable(myPhases_[iPhase].toString()),"liquid") == 0 )
            phaseID_.push_back(LIQUID);
        else if( strcmp(qPrintable(myPhases_[iPhase].toString()),"solid") == 0 )
            phaseID_.push_back(SOLID);
        else
            error().throw_error_one(FLERR,"ERROR: could not find one of the 'phaseNames' provided. \n");


        int eqnId = modelEqnContainer().lookupEqn( mySpecies_[iPhase].toString().toStdString().c_str() );
        if (eqnId<0)
            error().throw_error_one(FLERR,"ERROR: could not find the equn that was specified in 'speciesNames'. \n");

        speciesParticleDataID_.push_back( modelEqnContainer().modelEqn(eqnId)->particleDataID() );
        speciesModelEqnID_.push_back( eqnId );
        phaseFractions_.push_back(0.0);
        species_concentrations_.push_back(0.0);
    }

    if(phaseID_.size()!=2)
        error().throw_error_one(FLERR,"ERROR: phaseID_.size()!=2. \n");
    if(speciesParticleDataID_.size()!=2 )
        error().throw_error_one(FLERR,"ERROR: speciesParticleDataID_.size()!=2. \n");
    if(speciesModelEqnID_.size()!=2)
        error().throw_error_one(FLERR,"ERROR: speciesModelEqnID_.size()!=2. \n");

    printf("ModelPhaseChange::readPhaseChangeBasicsJSON successfully identified %d speciesParticleDataIDs: %d/%d, speciesModelEqnlID: %d/%d and %d phaseIDs: %d/%d. \n \n",
           (int)speciesParticleDataID_.size(),
           speciesParticleDataID_[0], speciesParticleDataID_[1],
           speciesModelEqnID_[0], speciesModelEqnID_[1],
           (int)phaseID_.size(),phaseID_[0],phaseID_[1]
          );

    if(modelEqnContainer().nrHeatEqns()>0) //if heat equation, perform non-isothermal calculation
        isIsoThermal_ = false;

    isSet_ = true;

}

// *************************************************************************************
void ModelPhaseChange::postParticleDataBuild()
{
    setPhaseChangeBasics(); //must do here since particleData does not know before of phase information
}

// *************************************************************************************
void ModelPhaseChange::begin_of_step()
{
    if(!isSet_)
        error().throw_error_one(FLERR,"ERROR: cannot do 'begin_of_step()' since ModelPhaseChange is not set. \n");
}
