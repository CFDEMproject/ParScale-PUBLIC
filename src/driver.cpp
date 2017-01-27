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
#include "model_phasechange_container.h"
#include "model_chemistry_container.h"
#include "comm.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

using namespace PASCAL_NS;

#define VERBOSE_MAIN false
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
        char msgstr[500];
        sprintf(msgstr,": initalizig ParScale... \n");
        output().write_screen_all(msgstr);

        pascal_ptr()->init();

        sprintf(msgstr,": ParScale initalized! \n");
        output().write_screen_all(msgstr);

        if(!particleData().runtimeContainersExist())
            error().throw_error_one(FLERR,"internal error: particleData is not allocated. Cannot proceed. \n");

        modelPhaseChangeContainer().postParticleDataBuild();
    }
    else
    {
        // pull new domain boundaries and proc info, then exchange atoms
        comm().pull();
    }

    while(time_left > 0.)
    {
        if(comm().is_proc_0())
        {
#if VERBOSE_MAIN
            printf("ParScale time : %g \n", control().simulationState().time());
#endif
        }

        dt = control().simulationState().deltaT();
        if(dt > time_left)
        {
            dt = time_left;
            control().simulationState().setDeltaT(dt);
        }

        if(dt > minTimeStep_)
            run_one_timestep();
        else
        {
#if VERBOSE_MAIN
            printf("Driver:: supressing time step because it is smaller than %g.\n", minTimeStep_);
#endif
        }

        time_left -= dt;

        control().simulationState().progress_time(dt);

        if(control().simulationState().writeContainers())
            particleData().write();

    }

   //TODO: finalize output to files
}

/* ----------------------------------------------------------------------
    Run application for a timestep
------------------------------------------------------------------------- */

void Driver::run_one_timestep()
{

    // TODO: pass non-const reference to modifyable data as function arguments
    coupling().pull_tags(); //TODO-NOT implemented yet

    comm().exchange();

    coupling().pull();

    modelEqnContainer().setupParticle();

//     printf("[%d/%d]: Advancing models ... \n",
//            comm().me(),comm().nprocs()
//           );
      do
      {
        //Begin Section - 1 Setup Models
        modelContainer().begin_of_step();

        particleData().resetPhaseChangeSourceTerms();
        modelPhaseChangeContainer().begin_of_step();

        particleData().resetChemicalSourceTerms();
        modelChemistryContainer().begin_of_step();

        particleData().computeConvection();

        //Begin Section - 2 Evaluate Equations
        modelEqnContainer().begin_of_step();

        //Middle Section
        modelContainer().pre_middle_of_step();
        modelEqnContainer().pre_middle_of_step();

        //Post Section
        modelContainer().post_middle_of_step();
        modelEqnContainer().post_middle_of_step();

        //Finalize
        modelContainer().end_of_step();
        modelEqnContainer().end_of_step();
      }

      while ( control().doImplicitIteration() );

    //printf("Computing properties and saving data for all particles ... \n");
    modelEqnContainer().computeParticleProps();

    coupling().push();

}
