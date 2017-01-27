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
    Base Class for phase change models (e.g., evaporation, solidification,
    precipitation...) between phase A (more dense, e.g., liquid)
    and phase B (less dense, e.g., gas).

\*-----------------------------------------------------------------------------------*/

#ifndef PASC_PHASECHANGE_H
#define PASC_PHASECHANGE_H

#include "memory_ns.h"
#include "particle_mesh.h"
#include "model_base.h"
#include "control.h"
#include "simulation_state.h"

using namespace PASCAL_MEMORY_NS;

namespace PASCAL_NS
{

class ModelPhaseChange : public ModelBase
{
    public:

     ModelPhaseChange(ParScale *ptr, char *name);
     ~ModelPhaseChange();

     void postParticleDataBuild();

     virtual void begin_of_step();
     virtual void init(int narg, char const* const* arg, int eqnType, int modelPhaseChangeID);

     //Access
     const vector<int>*   phaseID() const               { return &phaseID_; };
     const vector<int>*   speciesParticleDataID() const { return &speciesParticleDataID_; };
     const vector<int>*   speciesModelEqnlID()    const { return &speciesModelEqnID_; };

     virtual void         phaseChangeDeltaH(double, double & d) const {};
     virtual void         derivativeConcTemp(double, double & d) const {};

    protected:

     ParScale      *ptr_;

     bool           verbose_;

     bool           isActive_;

     bool           isSet_;

     bool           isIsoThermal_;

     QJsonArray     myPhases_;
     QJsonArray     mySpecies_;

     vector<int>    phaseID_;                //id of phase,              0: more dense; 1: less dense
     vector<int>    speciesParticleDataID_;  //id of species in phase;   0: more dense; 1: less dense
     vector<int>    speciesModelEqnID_;  //id of species in phase;   0: more dense; 1: less dense

     vector<double> phaseFractions_;         //local phase fractions
     vector<double> species_concentrations_; //concentration of species

     int            heatEqnID_;         //id of heat eqn.

     double*        tempIntraDataSpecies_;
     double*        tempIntraDataHeat_;

     void           readPhaseChangeBasicsJSON();
     void           setPhaseChangeBasics();

    };

} //end PASCAL_NS

#endif
