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

#include "chemistry_grainmodel.h"
#include "input.h"

using namespace PASCAL_NS;

/* ----------------------------------------------------------------------
   Constructor / Destructor
------------------------------------------------------------------------- */

ChemistryGrainModel::ChemistryGrainModel(ParScale *ptr)
:
   ParScaleBaseAccessible(ptr),
   cSolidInit_(1.0)
{
    active_     = false;
    solidID_    = 0;
    psi0_       = -1;
    alpha_      = -1;
}


// ----------------------------------------------------------------------
ChemistryGrainModel::~ChemistryGrainModel()
{
}

// ----------------------------------------------------------------------
void ChemistryGrainModel::init()
{
    //Read verbose file
    input().openJsonFile("settings", "verbose", "verbose", global_properties_ );
    if(global_properties_["chemistry_grainmodel"].isNull())
        error().throw_error_one(FLERR,"ERROR: could not read verbose settings for class. \n",
                                 "chemistry_grainmodel");
    verbose_   = global_properties_["chemistry_grainmodel"].toBool();


    //Set the model type and corresponding parameters
    input().openJsonFile("settings", "chemistry_grainmodel", "modelParameters", grainmodel_properties_);
    if(grainmodel_properties_["type"].isNull())
        error().throw_error_one(FLERR,"ERROR: could not read model/type from settings/chemistry_grainmodel.json. \n");

    //General model settings
    active_ = true;
    if(  grainmodel_properties_["cSolidInit"].isNull()
      && strcmp(qPrintable(grainmodel_properties_["type"].toString()),"none")
      )
        error().throw_error_one(FLERR,"ERROR: could not read cSolidInit from settings/chemistry_grainmodel.json. \n");
    cSolidInit_ = grainmodel_properties_["cSolidInit"].toDouble();

    solidID_ = 0;
    if(  grainmodel_properties_["solidID"].isNull()
      && strcmp(qPrintable(grainmodel_properties_["type"].toString()),"none")
      )
        error().throw_error_one(FLERR,"ERROR: could not read solidID from settings/chemistry_grainmodel.json. \n");
    solidID_     = grainmodel_properties_["solidID"].toInt();

    //Model Types & Params
    psi0_  = -1;
    alpha_ = -1;
    if(strcmp(qPrintable(grainmodel_properties_["type"].toString()),"none") == 0 )
    {
        active_ = false;
        FX = &ChemistryGrainModel::fxNone;
    }
    else if(strcmp(qPrintable(grainmodel_properties_["type"].toString()),"volumetric") == 0 )
    {
        FX = &ChemistryGrainModel::fxVolumetric;
    }
    else if(strcmp(qPrintable(grainmodel_properties_["type"].toString()),"grainClassical") == 0 )
    {
        FX = &ChemistryGrainModel::fxGrainClassical;
    }
    else if(strcmp(qPrintable(grainmodel_properties_["type"].toString()),"randomPore") == 0 )
    {
        FX = &ChemistryGrainModel::fxRandomPore;
        if(grainmodel_properties_["psi0"].isNull())
            error().throw_error_one(FLERR,"ERROR: could not read psi0 from settings/  chemistry_grainmodel.json. This parameter is required for the randomPore model! \n");
        psi0_ = grainmodel_properties_["psi0"].toDouble();
    }
    else if(strcmp(qPrintable(grainmodel_properties_["type"].toString()),"powerLaw") == 0 )
    {
        FX = &ChemistryGrainModel::fxPowerLaw;
        if(grainmodel_properties_["alpha"].isNull())
            error().throw_error_one(FLERR,"ERROR: could not read alpha from settings/  chemistry_grainmodel.json. This parameter is required for the powerLaw model! \n");
        alpha_ = grainmodel_properties_["alpha"].toDouble();
    }
    else
        error().throw_error_one(FLERR,"ERROR: could not find correct type of grain model. \n");

}


// ----------------------------------------------------------------------
void ChemistryGrainModel::printStatus() const
{
        printf(" *** ChemistryGrainModel - Status ***\n\n");
        printf("  type      : %s \n", qPrintable(grainmodel_properties_["type"].toString()));
        printf("  cSolidInit: %g \n", cSolidInit_);
        printf("  solidID   : %d \n", solidID_);
        if(strcmp(qPrintable(grainmodel_properties_["type"].toString()),"powerLaw") == 0 )
            printf("  alpha     : %g\n",
                 alpha_);
        if(strcmp(qPrintable(grainmodel_properties_["type"].toString()),"randomPore") == 0 )
            printf("  psi0      : %g\n",
                    psi0_);
        printf("\n *** ChemistryGrainModel - Status ***\n\n");
}
