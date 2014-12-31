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
#include "model_eqn.h"
#include "model_container.h"
#include "integrator_cvode.h"

using namespace PASCAL_NS;

/* ----------------------------------------------------------------------
   ModelEqn Constructor
------------------------------------------------------------------------- */

ModelEqn::ModelEqn(ParScale *ptr,  char *name) :
          ModelBase(ptr, name),
          integrator_(NULL),
          tempIntraData_(NULL),
          tempAvData_(0.0),
          tempPartFlux_(0.0),
          diffusivity(NULL),
          thermal_solid_conductivity(NULL),
          thermal_gas_conductivity(NULL),
          transfer_coeff(NULL),
          porosity(NULL),
          capacity_solid(NULL),
          capacity_gas(NULL),
          density_solid(NULL),
          density_gas(NULL),
          tortuosity(NULL),
          flux(NULL),
          eqnSolveThis_(true),
          eqnID_(-1),
          eqnType_(-1),
          particleDataID_(-1),
          haveParticleDataID(false),
          modelingApproach(CONTINUUM),
          nGridPointsUsed_(0)
{


}

/* ----------------------------------------------------------------------
   ModelEqn Destructor
------------------------------------------------------------------------- */

ModelEqn::~ModelEqn()
{

    delete integrator_;
    delete tempIntraData_;
    // TODO destroy all members of models_
}

/* ----------------------------------------------------------------------
   ModelEqn Member Functions
------------------------------------------------------------------------- */

void ModelEqn::init(int narg, char const* const* arg, int eqnType, int modelEqnID)
{
  eqnType_ = eqnType;
  eqnID_   = modelEqnID;
    
  //Allocate mem
  tempIntraData_   = create<double>(tempIntraData_, nGridPointsUsed_); 

  //Check models for physical properties one might need
  //TODO:throw errors,effective values possible - check for that in model_eqn_1D_spherical
  for(int iModel=0; iModel < modelContainer().modelCount(); iModel++)
  {
      //A - these models belong to all types of eqns.
      if(strstr(modelContainer().model(iModel)->name(), "porosity") != NULL)
      {
              printf("found porosity in model %s \n", modelContainer().model(iModel)->name());
              porosity=modelContainer().model(iModel);
              printf("with value: %g \n", porosity->value());
      }
    
      //B - these models are specific for a certain type of eqn.
      if(strstr(modelContainer().model(iModel)->name(), name()) != NULL)
      {
		 if(strstr(modelContainer().model(iModel)->name(), "Diffusivity") != NULL)
          {
			  printf("found diffusivity in model %s \n", modelContainer().model(iModel)->name());
              diffusivity=modelContainer().model(iModel);
              printf("with value: %g \n", diffusivity->value());
          }

          if(strstr(modelContainer().model(iModel)->name(), "Tortuosity") != NULL)
          {
			  printf("found tortuosity in model %s \n", modelContainer().model(iModel)->name());
              tortuosity=modelContainer().model(iModel);
              printf("with value: %g \n", tortuosity->value());
          }

          if(strstr(modelContainer().model(iModel)->name(), "Capacity_solid") != NULL)
          {
			 
              printf("found solid capacity in model %s \n", modelContainer().model(iModel)->name());
              capacity_solid=modelContainer().model(iModel);
              printf("with value: %g \n", capacity_solid->value());
          }

          if(strstr(modelContainer().model(iModel)->name(), "Capacity_gas") != NULL)
          {
              printf("found gas capacity in model %s \n", modelContainer().model(iModel)->name());
              capacity_gas=modelContainer().model(iModel);
              printf("with value: %g \n", capacity_gas->value());
          }

		  if(strstr(modelContainer().model(iModel)->name(), "ThermalConductivity_solid") != NULL)
          {
              printf("found solid thermal conductivity in model %s \n", modelContainer().model(iModel)->name());
              thermal_solid_conductivity=modelContainer().model(iModel);
              printf("with value: %g \n", thermal_solid_conductivity->value());
          }

           if(strstr(modelContainer().model(iModel)->name(), "ThermalConductivity_gas") != NULL)
          {
              printf("found gas thermal conductivity in model %s \n", modelContainer().model(iModel)->name());
              thermal_gas_conductivity=modelContainer().model(iModel);
              printf("with value: %g \n", thermal_gas_conductivity->value());
          }

		  if(strstr(modelContainer().model(iModel)->name(), "TransferCoeff") != NULL)
          {
              printf("found heat/mass transfer coefficient in model %s \n", modelContainer().model(iModel)->name());
              transfer_coeff=modelContainer().model(iModel);
              printf("with value: %g \n", transfer_coeff->value());
          }

          if(strstr(modelContainer().model(iModel)->name(), "Density_solid") != NULL)
          {
              printf("found solid desity in model %s \n", modelContainer().model(iModel)->name());
              density_solid=modelContainer().model(iModel);
              printf("with value: %g \n", density_solid->value());
          }

            if(strstr(modelContainer().model(iModel)->name(), "Density_gas") != NULL)
          {
              printf("found gas desity in model %s \n", modelContainer().model(iModel)->name());
              density_gas=modelContainer().model(iModel);
              printf("with value: %g \n", density_gas->value());
          }

          if(strstr(modelContainer().model(iModel)->name(), "flux") != NULL)
          {
              printf("found flux in model %s \n", modelContainer().model(iModel)->name());
              flux=modelContainer().model(iModel);
              printf("with value: %g \n", flux->value());
          }

		   if(strstr(modelContainer().model(iModel)->name(), "global_properties") != NULL)
          {
			  printf("found global_properties in model %s \n", modelContainer().model(iModel)->name());
              global_properties=modelContainer().model(iModel);
              //printf("with value: %g \n", diffusivity->value());
          }
      } //end have found model eqn with appropriate namespace


  } //end loop over model

  // parse boundary conditions and id of storage
  int iarg = 0;
  ICvalue=-999;
  BC[0]=BC[1]=-1;
  BCvalue[0]=BCvalue[1]=-999;
  BCvalueExtra[0]=BCvalueExtra[1]=-999;



  char BCname[4];


  while (iarg < narg)
  {

    for(int iBdry=0; iBdry<2;iBdry++)
    {
      sprintf(BCname,"BC%d",iBdry);
       if (strcmp(arg[iarg],BCname)==0)
       {
            BC[iBdry]=atoi(arg[iarg+1]);
            iarg+=1;

            printf("Boundary condition %d defined for %s. \n",
                    BC[iBdry],
                    BCname
                  );
       }
    }
    iarg++;
  };

  if( ((modelingApproach==CONTINUUM) && BC[0]==-1) || BC[1]==-1)
  {
    printf("ERROR: please specify one of the following BCs:\n%s(%d)\n%s(%d)\n%s(%d)\n",
            "NEUMANN",    NEUMANN,
            "DIRICHLET",  DIRICHLET,
            "CONVECTIVE", CONVECTIVE);

    error().throw_error_one(FLERR,"Boundary conditions not specified.\n");
  }


  //init the integrator (must happen after BCs are defined!)
  //TODO: allow a more flexible selection of the integrator!
  //for now, cvode is the standard integrator
  if(!integrator_)
      integrator_ = new IntegratorCvode(pascal_ptr(), nGridPointsUsed_);

  integrator_->init(0.0, *this);

  printf("*********ModelEqn initialized*****************\n\n");

}


// ----------------------------------------------------------------------
void ModelEqn::setupParticle()
{
    if(!haveParticleDataID)
        error().throw_error_one(FLERR,"You must specify a 'particleDataID' to let this model know where to save/request data.\n");

    return;
}


// ----------------------------------------------------------------------
void ModelEqn::computeParticleProps()
{

    if(!haveParticleDataID)
	error().throw_error_one(FLERR,"You must specify a 'particleDataID' to let this model know where to save/request data.\n");	

	//compute the averages
	computeParticleAverages(); 

	//compute flux at surface
	computeSurfaceFluxes();  

    return;
}
