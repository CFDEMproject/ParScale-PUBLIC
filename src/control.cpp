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

#include "control.h"
#include "stdlib.h"
#include "input.h"
#include "coupling.h"
#include "output.h"

using namespace PASCAL_NS;

/* ----------------------------------------------------------------------
   Control Constructor
------------------------------------------------------------------------- */

Control::Control(ParScale *ptr) : ParScaleBase(ptr),
    driver_(* new Driver(ptr)),
    simulationState_(* new SimulationState(ptr)),
    numberImplicitIterations_(0),
    currImplicitIterations_(0)
{

}

Control::~Control()
{
    delete &driver_;
    delete &simulationState_;
}

/* ----------------------------------------------------------------------
    Initialization
------------------------------------------------------------------------- */

void Control::init()
{


}

/* ----------------------------------------------------------------------
    Parse a command
------------------------------------------------------------------------- */

void Control::parse_command(int narg,char const* const* arg)
{
    double endTime(0.0);
    double timeStep(1.0);

  // parse optional args
  int iarg = 0;

  while (iarg < narg)
  {
       if (strcmp(arg[iarg],"run")==0)
       {

            iarg++;

            if(iarg + 1 > narg)
                error().throw_error_all(FLERR,"Not enough arguments for 'control run' command\n");
            endTime = atof(arg[iarg++]);

            if(endTime <= 0.)
                error().throw_error_all(FLERR,"endTime > 0 required\n");

            bool do_init = true;

            if(iarg + 2 <= narg)
            {
                if(strcmp(arg[iarg++],"init"))
                    error().throw_error_all(FLERR,"Expecting keyword 'init' for 'control run' command\n");

                printf("arg[%d]: %s \n", iarg,arg[iarg]);
                if(0 == strcmp(arg[iarg],"yes"))
                    do_init = true;
                else if(0 == strcmp(arg[iarg],"no"))
                    do_init = false;
                else
                    error().throw_error_all(FLERR,"Expecting 'yes' or 'no' after 'init' for 'control run' command\n");
            }

            if(coupling().external_code_in_control() && do_init)
//                error().throw_error_all(FLERR,"Can not use 'control run' command with a coupling model, "
//                                              "since the coupling model fully controls time-stepping\n");

            //check if time step has been set
            if(simulationState().deltaT()>0.0f)
            {
            }
            else
            {
                output().write_screen_one("Control: forcing to set time step to 1.0");
                simulationState().setDeltaT(1.0);
            }

            // Run the simulation
            // this also inits the simulation
            output().write_screen_one("Control is firing up ParScale.");
            driver_.run(endTime,do_init);
            output().write_screen_one("ParScale run done.");

            iarg++;
      }
      else if (strcmp(arg[iarg],"timeStep")==0)
      {
          timeStep=atof(arg[iarg+1]);
          simulationState().setDeltaT(timeStep);
          iarg++;
      }
      else if (strcmp(arg[iarg],"outputTimeStep")==0)
      {
          timeStep=atof(arg[iarg+1]);
          simulationState().setOutputTimeStep(timeStep);
          iarg++;
      }
      else  iarg++;
  };

}
