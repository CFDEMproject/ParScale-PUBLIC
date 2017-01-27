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

#ifndef PASC_CONTROL_H
#define PASC_CONTROL_H

#include "stdio.h"
#include "driver.h"
#include "simulation_state.h"
#include "pascal_base.h"
#include "pascal_base_interface.h"

namespace PASCAL_NS
{

class Control : public ParScaleBase, public ParScaleBaseInterface
{
    public:

      Control(ParScale *ptr);
      ~Control();

      void init();

      void parse_command(int narg,char const* const* arg);

      SimulationState& simulationState()  const {return simulationState_;}

      //Access function
      bool doImplicitIteration()
      {
            currImplicitIterations_+=1;
            if(currImplicitIterations_> numberImplicitIterations_)
            {
                currImplicitIterations_ = 0;
                return false;
            }
            else
                return true;
      };

    private:

      Driver& driver_;

      SimulationState& simulationState_;

      int numberImplicitIterations_;

      int currImplicitIterations_;
};

} //end PASCAL_NS

#endif
