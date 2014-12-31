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

#include "pascal.h"
#include "input.h"
#include "output.h"
#include "control.h"
#include "comm.h"
#include "coupling.h"
#include "error.h"
#include "particle_mesh.h"
#include "particle_data.h"
#include "fluid_data.h"
#include "model_container.h"
#include "model_eqn_container.h"
#include "model_chemistry_container.h"
#include "model_base.h"
#include "pascal_base_interface_vector.h"
#include "chemistry_reader.h"
#include "coupling_model.h"

using namespace PASCAL_NS;

/* ----------------------------------------------------------------------
   Central ParScale Class
    - allocates memory, error, input, etc.
------------------------------------------------------------------------- */

ParScale::ParScale(int narg, char **arg, MPI_Comm communicator,void* caller) :
  callingProgram_(caller),
  paScalBaseInterfaceVector_( new ParScaleBaseInterfaceVector),  
  comm_     ( paScalBaseInterfaceVector_->add<Comm>(new Comm(this,communicator),"comm")),
  input_    ( paScalBaseInterfaceVector_->add<Input>(new Input(this),"input")),
  output_   ( paScalBaseInterfaceVector_->add<Output>(new Output(this),"output")),
  control_  ( paScalBaseInterfaceVector_->add<Control>(new Control(this),"control")),
  error_    ( paScalBaseInterfaceVector_->add<Error>(new Error(this),"error")),
  particleMesh_ ( paScalBaseInterfaceVector_->add<ParticleMesh>(new ParticleMesh(this),"particle_mesh")),
  modelContainer_( paScalBaseInterfaceVector_->add<ModelContainer>(new ModelContainer(this),"model")),
  modelEqnContainer_( paScalBaseInterfaceVector_->add<ModelEqnContainer>(new ModelEqnContainer(this),"modelEqn")),
  modelChemistryContainer_( paScalBaseInterfaceVector_->add<ModelChemistryContainer>(new ModelChemistryContainer(this),"modelchemistry")),
  // must come before particleData_ so coupling_.read() is executed before particleData_.read()
  // and coupling.pull() is executed before particleData_.pull(
  coupling_ ( paScalBaseInterfaceVector_->add<Coupling>(new Coupling(this),"coupling")),
  particleData_ ( paScalBaseInterfaceVector_->add<ParticleData>(new ParticleData(this),"particle_data")),
  fluidData_ ( paScalBaseInterfaceVector_->add<FluidData>(new FluidData(this),"fluid_data"))

{
}

/* ----------------------------------------------------------------------
   shutdown ParScale
   delete top-level classes
------------------------------------------------------------------------- */
ParScale::~ParScale()
{
    delete paScalBaseInterfaceVector_;

    delete input_;
    delete output_;
    delete control_;
    delete coupling_;
    delete error_;
    delete particleMesh_;
    delete particleData_;
    delete fluidData_;
    delete modelContainer_;
    delete modelEqnContainer_;
    delete modelChemistryContainer_;
    //delete couplingModel_;
    delete comm_;          // delete this one last
   // delete modelBase_;
}

/* **********************************************************************
       MEMBER FUNCTIONS
********************************************************************** */
/* ----------------------------------------------------------------------
   initialize top-level classes
   do not initialize Timer class, other classes like Run() do that explicitly
------------------------------------------------------------------------- */
void ParScale::init()
{

    printf("\n*******************************************\n");
    printf("\nCreating a Pascal object \n");
    printf("Will initialize individual components of Pascal now... \n");

    paScalBaseInterfaceVector_->init();

    printf("\n*******************************************\n\n");

}

/* ----------------------------------------------------------------------
   Set input & directory
------------------------------------------------------------------------- */
void ParScale::set_dir(const char * dirname)
{
    input_->set_runDirectory(dirname);
}

// * * * * * * * * * * * * * * * * * * * * * * * *
void ParScale::set_input(const char * filename)
{
    input_->set_input_script(filename);
}

/* ----------------------------------------------------------------------
   Process input
------------------------------------------------------------------------- */
void ParScale::input()
{
   input_->process_input_script();
}

/* ----------------------------------------------------------------------
   Command execution
------------------------------------------------------------------------- */
void ParScale::runCommand(const char * command)
{
   input_->process_input(command);
}

/* ----------------------------------------------------------------------
   help message for command line options and styles present in executable
------------------------------------------------------------------------- */
void ParScale::help()
{

  printf( //(screen,
          "\nCommand line options:\n\n"
          "-in filename                : read input from file, not stdin NOTIMPLEMENTED (-i)\n"
          "-help                       :  help message NOTIMPLEMENTED (-h)\n"
          "-log none/filename          : where to send log output NOTIMPLEMENTED (-l)\n"
          "-screen none/filename       : where to send screen output NOTIMPLEMENTED (-sc)\n"
          "-other command line Options          : NOTIMPLEMENTED\n\n");

}

