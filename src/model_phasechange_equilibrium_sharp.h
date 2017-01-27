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
/*
    Contributing author:    Thomas Forgber, TU Graz
    Phase change model for sharp phase changes at saturation concentration
    (i.e., liquid <-> solid) phenomena

\*-----------------------------------------------------------------------------------*/

#ifdef MODEL_PHASECHANGE_CLASS

ModelPhaseChangeStyle(Equilibrium_sharp, ModelPhaseChangeEquilibrium_sharp)

#else

#ifndef PASC_PHASECHANGE_EQUILIBRIUM_SHARP_H
#define PASC_PHASECHANGE_EQUILIBRIUM_SHARP_H

#include "model_phasechange_container.h"
#include "input.h"
#include "simulation_state.h"
#include "model_eqn_container.h"
#include "model_eqn.h"

namespace PASCAL_NS
{

#define INVERSE_UNIVERSAL_GAS_CONSTANT 0.0001202723625

class ModelPhaseChangeEquilibrium_sharp : public ModelPhaseChange
{
    public:

     ModelPhaseChangeEquilibrium_sharp(ParScale *ptr, char *name);
     ~ModelPhaseChangeEquilibrium_sharp();

     void begin_of_step();
     void init(int narg, char const* const* arg, int eqnType, int modelPhaseChangeID);

    private:

     double         evaporationRateConstant_;   //Constant to control evaporation rate
     double         jacobiForDensePhase_;       //Jacobi for dense phase

     double         sat_concentration_liq_;          //Concentraiton of liquid phase at saturation pressure
     double         sat_concentration_solid_;        //Concentraiton of solid phase at saturation pressure

     double         delta_concentration_;

     bool           update_phase_fraction_;     // Find out if phase fraction is updated
     int            update_fraction_id;         // Model eqn id for which fraction is solved

     ParScale       *ptr_;

     QJsonObject    global_properties_;

     //Functions
     void           readEquilibriumSharpJSON();

    };

} //end PASCAL_NS

#endif

#endif
