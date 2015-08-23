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


#include "coupling_model_json.h"
#include "input.h"
#include "output.h"
#include "particle_data.h"
#include "control.h"

using namespace PASCAL_NS;

/* ----------------------------------------------------------------------
   ModelBase Constructor
------------------------------------------------------------------------- */

CouplingModelJSON::CouplingModelJSON(ParScale *ptr, const char *_name)
:
CouplingModel(ptr,_name)
{

}


/* --------------------------------------------------------------------- */
void CouplingModelJSON::read()
{
    CouplingModel::read();
    printf("read in json coulpling \n");
    input().openJsonFile("settings", scope(), "setBCs", setBCs_);

    //Loop through setBCs and report
    QStringList setBCStrings = setBCs_.keys();
    QJsonObject::const_iterator i;
    int iCurr = 0;
    for (i = setBCs_.begin(); i != setBCs_.end(); ++i)
    {
        settingBCName_.push_back(setBCStrings[iCurr]);
        QJsonObject currSetting = (*i).toObject();

        settingBCEqn_.push_back(currSetting["dataField"].toString());
        settingBCType_.push_back(currSetting["type"].toString());
        QJsonArray  currTimes   = currSetting["times"].toArray();
        settingBCTime_.push_back(currTimes);
        QJsonArray  currValues  = currSetting["values"].toArray();
        settingBCValue_.push_back(currValues);

        //Settings for where to apply   
        applyAt_ = -1;
        if( !currSetting["applyAt"].isNull() )
        {
            if( strcmp(currSetting["applyAt"].toString().toStdString().c_str(),"fluid")==0 )
                applyAt_ = 9999999;
            else if( strcmp(currSetting["applyAt"].toString().toStdString().c_str(),"surface")==0 )
                applyAt_ = -2;
            else
                applyAt_ = currSetting["applyAt"].toInt();
        }
        settingBCApplyAt_.push_back(applyAt_); 

        //check for correct length of values
        if(!(currValues.size()==currTimes.size()))
            error().throw_error_one(FLERR,"Not enough/too much values in file settings/coupling*.json! \n");

        if(verbose_)
        {
            printf( "CouplingModelJSON: reading current setBC: %s ", settingBCName_[iCurr].toStdString().c_str() );
            printf( " for dataField '%s' with type '%s' and applyAt setting at: %d \n",
                       settingBCEqn_[settingBCEqn_.size()-1].toStdString().c_str(),
                       settingBCType_[settingBCType_.size()-1].toStdString().c_str(),
                       settingBCApplyAt_[settingBCType_.size()-1]
                  );
            for(int j = 0; j < settingBCTime_[settingBCTime_.size()-1].size(); j++)
            {

                printf("CouplingModelJSON: found time/value pairs %g / %g \n",
                        settingBCTime_ [settingBCTime_.size()-1 ][j].toDouble(),
                        settingBCValue_[settingBCValue_.size()-1][j].toDouble()
                      );
            }
        }
        iCurr++;
    }

}

/* ----------------------------------------------------------------------
   pull/push functions
------------------------------------------------------------------------- */
bool CouplingModelJSON::fill_container_from_coupling(class ContainerBase &container) const
{
    if(pullNext_ > control().simulationState().timeStep() )
        return false;

    const ParticleDataContainerProperties& containerProps = container.prop();

    //Check correct container
    //const char *id               = container.prop().id();
    const char *scope            = container.prop().scope();
    const bool  element_property = container.prop().element_property();

    //Loop all settings for BCs
    for(unsigned int iSetting=0; iSetting < settingBCName_.size(); iSetting++)
    {
        //if scope equals the dataField and it is an element_property, actually do the fill
        if( (strcmp(settingBCEqn_[iSetting].toStdString().c_str(),scope) == 0) && element_property )
        {
                if(verbose_)
                    printf("...processing container with id '%s', needsPull: %d, needsRead: %d double: %d, int: %d. \n",
                        containerProps.id(),
                        containerProps.needsPull(),
                        containerProps.needsRead(),
                        container.isDoubleData(), container.isIntData()
                    );

                //Set intra-particle position to apply setting
                applyAt_ = settingBCApplyAt_[iSetting];
                if( applyAt_ >= container.lenVecUsed())
                    applyAt_  = container.lenVecUsed()-1; //apply in the fluid

                if( settingBCApplyAt_[iSetting] < -1)
                    applyAt_  = container.lenVecUsed()+settingBCApplyAt_[iSetting];  //apply at a distance from the surface

                double currTime  = control().simulationState().time()
                                 + control().simulationState().deltaT(); //Use next time!
                double currValue = MathExtraPascal::interpolateLinearly(settingBCTime_[iSetting],
                                                                        settingBCValue_[iSetting],
                                                                        currTime);

                double *  ptr     = (double*)      container.begin_slow_dirty();   //ptr to scalar per-particle data
                double ***ptr3    = (double***)    container.begin_slow_dirty();   //ptr to vectorial per-particle data
 
                //This type will apply the same BC to all particles
                if(strcmp("applyToAllParticles",settingBCType_[iSetting].toStdString().c_str()) == 0)
                {
                    if(container.lenVec()>1)
                        for (int iP = 0; iP < particleData().nbody(); ++iP)
                            insertSingleValueInArray(currValue, container, iP, ptr3 );  
                    else              
                        for (int iP = 0; iP < particleData().nbody(); ++iP)
                            insertSingleValue(currValue, container, iP, ptr );  
                }
                else
                    error().throw_error_one(FLERR,"This type is not implemented! See src code for implemented types. \n");

        }

    }

    hasPulled_ = true;
    return true;
}

/* -------------------------------------------------------------------- */
bool CouplingModelJSON::dump_container_to_coupling(class ContainerBase &container) const
{
    //TODO: Implement filling of the container
    return false;
}


/* -------------------------------------------------------------------- */
void CouplingModelJSON::insertSingleValue(  double currValue, 
                                            class  ContainerBase &container, 
                                            int    iP,
                                            double *  ptr
                                         ) const
{
    if(iP>=container.size())
        error().throw_error_one(FLERR,"particleId>=container.size(). Something went wrong.\n");

    ptr[iP] = currValue;

    if(verbose_)
        printf("...inserted single value (%g) in a scalar container.  \n \n", 
                currValue);
}

/* -------------------------------------------------------------------- */
void CouplingModelJSON::insertSingleValueInArray(  double currValue, 
                                            class  ContainerBase &container, 
                                            int    iP,
                                            double *** ptr3
                                         ) const
{
    if(iP>=container.size())
        error().throw_error_one(FLERR,"particleId>=container.size(). Something went wrong.\n");

    ptr3[iP][0][applyAt_] = currValue; //Insert at a certain position
    if(verbose_)
        printf("...inserted single value (%g) at position %d of an array with length %d.  \n \n", 
                currValue, applyAt_, container.lenVec());

}
