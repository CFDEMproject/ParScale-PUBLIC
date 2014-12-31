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

#ifndef PASC_PASCAL_H
#define PASC_PASCAL_H

#include "stdio.h"
#include "mpi.h"

namespace PASCAL_NS
{

class ParScale
{
 friend class ParScaleBase;

 public:

  ParScale(int narg, char **args, MPI_Comm communicator, void* caller);
  ~ParScale();

  void init();

 private:

  void help();


  //Pointer to the calling program (if used as a library). 
  //Used to call functions of calling program if needed in a sub-model
  void* callingProgram_;

  // class to register all members of this class
  class ParScaleBaseInterfaceVector *paScalBaseInterfaceVector_;

  // code framework

  class Comm     *comm_;       // inter-processor communication
  class Input    *input_;      // input to simulation - input script processing
  class Output   *output_;     // simulation output
  class Control  *control_;    // simulation control
  class Error    *error_;      // error handling

  // data and physical/chemical models
  class ParticleMesh      *particleMesh_;     //default particle mesh
  class ModelContainer    *modelContainer_;   // container for physical/chemical models
  class ModelEqnContainer *modelEqnContainer_; // container for (discretized) equations
  class ModelChemistryContainer *modelChemistryContainer_; //Container for chemistry models
  class ChemistryReader *chemistryReader_;
  class Coupling          *coupling_;       // coupling to other codes
                                            // must come before particleData_ so coupling_.read()
                                            // is executed before particleData_.read()
  class CouplingModel     *couplingModel_;
  class ParticleData      *particleData_;     // storage for particle data
  class FluidData         *fluidData_;        // storage for fluid data

  public:

    //Top level functionalities
    void  set_dir(const char *);      // set the input
    void  set_input(const char *);    // set the input
    void  input();                    // process all input
    void  runCommand(const char *);   // run a single command
};

} //end PASCAL_NS

#endif

