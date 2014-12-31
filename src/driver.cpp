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

#include "error.h"
#include "driver.h"
#include "input.h"
#include "output.h"
#include "coupling.h"
#include "control.h"
#include "particle_data.h"
#include "simulation_state.h"
#include "model_container.h"
#include "model_eqn_container.h"
#include "model_chemistry_container.h"
#include "comm.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

using namespace PASCAL_NS;

/* ----------------------------------------------------------------------
   Driver Constructor
------------------------------------------------------------------------- */

Driver::Driver(ParScale *ptr)
:
 ParScaleBaseAccessible(ptr),
 minTimeStep_(1e-16)
{

}

/* ----------------------------------------------------------------------
    Run application for a given time
------------------------------------------------------------------------- */

void Driver::run(double _time,bool _do_init)
{
    double time_left = _time;
    double dt;

    if(_do_init)
    {   
        // init will read all properties from file
        // init also pulls data from LIGGGHTS via Coupling::read()
        // which is local on each processor, thus no need to bcast
        pascal_ptr()->init();
        if(!particleData().runtimeContainersExist())
        error().throw_error_one(FLERR,"internal error: particleData is not allocated. Cannot proceed. \n");
    }
    else
    {
        // pull new domain boundaries and proc info, then exchange atoms
        comm().pull();
        // pull atom data via coupling
        coupling().pull();
    }

    while(time_left > 0.)
    {
        printf("time : %g \n", control().simulationState().time());
        //printf("timestep : %g \n", control().simulationState().deltaT());

        dt = control().simulationState().deltaT();
        if(dt > time_left)
        {
            dt = time_left;
            control().simulationState().setDeltaT(dt);
        }

        if(dt > minTimeStep_)
            run_one_timestep();
        else
            printf("Driver:: supressing time step because it is smaller than %g.\n", minTimeStep_);

        time_left -= dt;

        control().simulationState().progress_time(dt);

        if(control().simulationState().writeContainers())
            particleData().write();

        printf("\n");
    }
   //printf("endtime : %g \n", control().simulationState().time());

   //TODO: finalize output to files
}

/* ----------------------------------------------------------------------
    Run application for a timestep
------------------------------------------------------------------------- */

void Driver::run_one_timestep()
{
    // TODO: pass non-const reference to modifyable data as function arguments
    coupling().pull_tags();

    comm().exchange();

    coupling().pull();

    modelEqnContainer().setupParticle();

//    printf("Advancing models for all particles ... \n");
      do
      {
        //Begin Section
        modelContainer().begin_of_step();

        particleData().resetChemicalSourceTerms();
        modelChemistryContainer().begin_of_step();

        modelEqnContainer().begin_of_step();

        //Middle Section     
        modelContainer().pre_middle_of_step();
        modelEqnContainer().pre_middle_of_step();

        //Post Section
        modelContainer().post_middle_of_step();
        modelEqnContainer().post_middle_of_step();

        //Finalize
        modelContainer().end_of_step();
      }

      while ( control().doImplicitIteration() );

    //printf("Computing properties and saving data for all particles ... \n");
    modelEqnContainer().computeParticleProps();

    coupling().push();
    
}
