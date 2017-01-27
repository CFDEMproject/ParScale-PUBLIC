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


#include "stdlib.h"
#include "memory_ns.h"
#include "model_chemistry.h"
#include "model_eqn_container.h"
#include "control.h"
#include "simulation_state.h"
#include "integrator_simple.h"
#include "integrator_cvode.h"
#include "model_eqn.h"
#include "model_container.h"

using namespace PASCAL_NS;
using namespace PASCAL_MEMORY_NS;

#define FEQUENZFACTOR         1
#define LARGENUMBER           1e32
/* ----------------------------------------------------------------------
   ModelChemistry Constructor
------------------------------------------------------------------------- */

ModelChemistry::ModelChemistry(ParScale *ptr,  char *name) :
       ModelBase(ptr, name),
       eqnReactionActive_(true),
       tempIntraDataHeat_(NULL),
       tempIntraDataSpecies_(NULL),
       kSurface_arrhenius_(0.0)
{
  create<double>(tempIntraDataHeat_,                   particleMesh().nGridPoints() );
  create<double>(tempIntraDataSpecies_,                particleMesh().nGridPoints() );
}

ModelChemistry::~ModelChemistry()
{
    destroy<double>(tempIntraDataHeat_);
    destroy<double>(tempIntraDataSpecies_);
}
/* ----------------------------------------------------------------------
   MemberFunctions
------------------------------------------------------------------------- */
void ModelChemistry::init(int narg, char const* const* arg, int eqnType, int modelChemistryID)
{
    //Check if heatEqn is found and can be set
    heatEqnID_ = -1;
    if(modelEqnContainer().nrHeatEqns() > 0)
        heatEqnID_ = 0; //take the first heat equation for the chemistry!
    else
    {
        printf("********ModelChemistry::WARNING: you do not have a heat equation. \n");
        printf("********Can only run in iso-thermal mode for ALL involved particles! \n\n");
    }

}

// ----------------------------------------------------------------------
void ModelChemistry::begin_of_step()
{

}
