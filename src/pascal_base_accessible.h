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

#ifndef PASC_BASE_ACCESSIBLE_H
#define PASC_BASE_ACCESSIBLE_H

#include "pascal_base.h"

namespace PASCAL_NS
{

class ParScaleBaseAccessible : public ParScaleBase
{
  public:

    ParScaleBaseAccessible(ParScale *ptr) : ParScaleBase(ptr) {}

    virtual ~ParScaleBaseAccessible() {}

  protected:

    inline ParScale*         pascal_ptr()     {return &pscl_;}

    inline Comm&           comm()           {return *comm_;}
    inline Input&          input()          {return *input_;}
    inline Output&         output()         {return *output_;}
    inline Control&        control()        {return *control_;}
    inline Coupling&       coupling()       {return *coupling_;}
    inline Error&          error()          {return *error_;}
    inline ParticleMesh&   particleMesh()   {return *particleMesh_;}
    inline ParticleData&   particleData()   {return *particleData_;}
    inline FluidData&      fluidData()      {return *fluidData_;}
    inline ModelContainer& modelContainer() {return *modelContainer_;}
    inline ModelEqnContainer& modelEqnContainer() {return *modelEqnContainer_;}
    inline ModelPhaseChangeContainer& modelPhaseChangeContainer() {return *modelPhaseChangeContainer_;}
    inline ModelChemistryContainer& modelChemistryContainer() {return *modelChemistryContainer_;}
    inline CouplingModel&       couplingModel()       {return *couplingModel_;}
    inline ChemistryReaderCHEMKIN& chemistryReaderCHEMKIN() {return *chemistryReaderCHEMKIN_;}


    inline void*            callingProgram()  { return callingProgram_;};

  private:

};

} //end PASCAL_NS

#endif
