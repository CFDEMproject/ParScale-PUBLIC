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
#include "integrator_simple.h"
#include "integrator_cvode.h"

using namespace PASCAL_NS;

/* ----------------------------------------------------------------------
   ModelEqn Constructor
------------------------------------------------------------------------- */

ModelEqn::ModelEqn(ParScale *ptr,  char *name) :
          ModelBase(ptr, name),
          integrator_(NULL),
          tempIntraData_(NULL),
          tempPhaseDataGas_(NULL),
          tempPhaseDataLiquid_(NULL),
          tempAvData_(0.0),
          tempPartFlux_(0.0),
          diffusivity(NULL),
          thermal_solid_conductivity(NULL),
          thermal_gas_conductivity(NULL),
          transfer_coeff(NULL),
          phaseFraction(NULL),
          capacity_solid(NULL),
          capacity_gas(NULL),
          capacity_liquid(NULL),
          density_solid(NULL),
          density_gas(NULL),
          tortuosity(NULL),
          eqnSolveThis_(true),
          eqnID_(-1),
          eqnType_(-1),
          particleDataID_(-1),
          haveParticleDataID(false),
          modelingApproach(CONTINUUM),
          nGridPointsUsed_(0),
          inPhase_(-1),
          surface_area(NULL),
          segment_vol(NULL),
          heatEqnID_(-1),
          molar_mass(NULL),
          permeability(NULL),
          viscosity(NULL),
          surface_tension(NULL),
          film_flow(NULL),
          pore_radius(NULL),
          particleID(-1),
          IsoThermal_(FALSE),
          diffu_eff_(NULL)
{
    phaseFractionMinimum = SMALLNUMBER;
    phaseFractionMinimumConvection = SMALLNUMBER;
    convectionBoundMinStefanFactor = 1e-2;
    solveMe             = true;
    updatePhaseFraction = false;
    averagePhaseFraction= false;
    solveConvectiveFlux = false;
    writeDebugContainers= false;
    normalizeDuringGrowth = false;

    verbose_ = false;

}

/* ----------------------------------------------------------------------
   ModelEqn Destructor
------------------------------------------------------------------------- */

ModelEqn::~ModelEqn()
{
    delete integrator_;
    destroy<double>(tempIntraData_);
    destroy<double>(tempPhaseDataGas_);
    destroy<double>(tempPhaseDataLiquid_);
}

/* ----------------------------------------------------------------------
   ModelEqn Member Functions
------------------------------------------------------------------------- */
void ModelEqn::end_of_step()
{
    //Loops particles to ensure positive species concentration
    if(eqnType_==SPECIES)
    {
       for(particleID=0; particleID<particleData().nbody(); particleID++)
       {
        if(updatePhaseFraction)
        {
           particleData().setParticleIDPointerPhaseFraction(particleID);
           particleData().returnPhaseFractionData(tempPhaseDataGas_, tempPhaseDataLiquid_);
           for(int i=0; i<particleMesh().nGridPoints(); i++)
           {
              tempPhaseDataGas_[i]    = max(0.0, min(1.0, tempPhaseDataGas_[i]));       //bound between zero and unity
              tempPhaseDataLiquid_[i] = max(0.0, min(1.0, tempPhaseDataLiquid_[i]));    //bound concentrations to zero
           }
//           tempPhaseDataGas_[particleMesh().nGridPoints()]    = tempPhaseDataGas_[particleMesh().nGridPoints()-1];      //copy last element
//           tempPhaseDataLiquid_[particleMesh().nGridPoints()] = tempPhaseDataLiquid_[particleMesh().nGridPoints()-1];   //copy last element
           particleData().savePhaseFractionData(particleID, tempPhaseDataGas_, tempPhaseDataLiquid_);
        }
        else
        {
          //Pull out data and correct
          particleData().setParticleIDPointer(particleDataID_,particleID);	
          particleData().returnIntraData(tempIntraData_); 
    
          for(int i=0; i<particleMesh().nGridPoints(); i++)
              tempIntraData_[i] = max(0.0, tempIntraData_[i]); //bound concentrations to zero
         
          particleData().saveIntraParticleData(particleDataID_, particleID, tempIntraData_);
        }
       }
    }
}

// *****************************************************************
void ModelEqn::init(int narg, char const* const* arg, int eqnType, int modelEqnID)
{
  eqnType_ = eqnType;
  eqnID_   = modelEqnID;
    
  if(modelEqnContainer().nrHeatEqns() > 0)
  {
    heatEqnID_ = 0; //take the first heat equation for the temperature dependency
  }  
  else
  { 
    read_chemistry_single_react_json_file("isIsoThermal",  &IsoThermal_, true);
  }  
  //TODO: make this more pretty - runs only if chemistry eqn is there
   
  //Error Checks
    if( (heatEqnID_<0) && !IsoThermal_)
        error().throw_error_one(FLERR,"ERROR: no heat equation found, but you like to run a non-isothermal simulations. This is impossible. \n");

    if( (heatEqnID_>=0) && IsoThermal_)
        error().throw_error_one(FLERR,"ERROR: heat equation found, but you like to run an isothermal simulations. This is not meaningful. Deactivate the heat equation! \n");


  //Allocate mem
  tempIntraData_       = create<double>(tempIntraData_, nGridPointsUsed_); 
  tempPhaseDataGas_    = create<double>(tempPhaseDataGas_, nGridPointsUsed_); 
  tempPhaseDataLiquid_ = create<double>(tempPhaseDataLiquid_, nGridPointsUsed_); 

  //Check models for physical properties one might need
  //TODO:throw errors,effective values possible - check for that in model_eqn_1D_spherical
  for(int iModel=0; iModel < modelContainer().modelCount(); iModel++)
  {
      //A - these models belong to all types of eqns.

      if(strstr(modelContainer().model(iModel)->name(), "PhaseFraction") != NULL)
      {
              printf("found phaseFraction in model %s \n", modelContainer().model(iModel)->name());
              phaseFraction=modelContainer().model(iModel);
              printf("with value: %g \n", phaseFraction->value());
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

          if(strstr(modelContainer().model(iModel)->name(), "Capacity_liquid") != NULL)
          {
              printf("found liquid capacity in model %s \n", modelContainer().model(iModel)->name());
              capacity_liquid=modelContainer().model(iModel);
              printf("with value: %g \n", capacity_liquid->value());
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

           if(strstr(modelContainer().model(iModel)->name(), "ThermalConductivity_liquid") != NULL)
          {
              printf("found liquid thermal conductivity in model %s \n", modelContainer().model(iModel)->name());
              thermal_liquid_conductivity=modelContainer().model(iModel);
              printf("with value: %g \n", thermal_liquid_conductivity->value());
          }

		  if(strstr(modelContainer().model(iModel)->name(), "TransferCoeff") != NULL)
          {
              printf("found heat/mass transfer coefficient in model %s \n", modelContainer().model(iModel)->name());
              transfer_coeff=modelContainer().model(iModel);
              printf("with value: %g \n", transfer_coeff->value());
          }

          if(strstr(modelContainer().model(iModel)->name(), "Density_solid") != NULL)
          {
              printf("found solid density in model %s \n", modelContainer().model(iModel)->name());
              density_solid=modelContainer().model(iModel);
              printf("with value: %g \n", density_solid->value());
          }

          if(strstr(modelContainer().model(iModel)->name(), "Density_gas") != NULL)
          {
              printf("found gas density in model %s \n", modelContainer().model(iModel)->name());
              density_gas=modelContainer().model(iModel);
              printf("with value: %g \n", density_gas->value());
          }

          if(strstr(modelContainer().model(iModel)->name(), "Density_liquid") != NULL)
          {
              printf("found liquid density in model %s \n", modelContainer().model(iModel)->name());
              density_liquid=modelContainer().model(iModel);
              printf("with value: %g \n", density_liquid->value());
          }

		  if(strstr(modelContainer().model(iModel)->name(), "global_properties") != NULL)
          {
			  printf("found global_properties in model %s \n", modelContainer().model(iModel)->name());
              global_properties=modelContainer().model(iModel);
              //printf("with value: %g \n", diffusivity->value());
          }

          if(strstr(modelContainer().model(iModel)->name(), "Pore_Radius") != NULL)
          {
              printf("found Pore_Radius in model %s \n", modelContainer().model(iModel)->name());
              pore_radius=modelContainer().model(iModel);
              printf("with value: %g \n", pore_radius->value());
          }

          if(strstr(modelContainer().model(iModel)->name(), "Molar_Mass") != NULL)
          {
              printf("found Molar Mass in model %s \n", modelContainer().model(iModel)->name());
              molar_mass=modelContainer().model(iModel);
              printf("with value: %g \n", molar_mass->value());
          }

          if(strstr(modelContainer().model(iModel)->name(), "Permeability") != NULL)
          {
              printf("found permeability in model %s \n", modelContainer().model(iModel)->name());
              permeability=modelContainer().model(iModel);
              printf("with value: %g \n", permeability->value());
          }

          if(strstr(modelContainer().model(iModel)->name(), "Viscosity") != NULL)
          {
              printf("found viscosity in model %s \n", modelContainer().model(iModel)->name());
              viscosity=modelContainer().model(iModel);
              printf("with value: %g \n", viscosity->value());
          }

          if(strstr(modelContainer().model(iModel)->name(), "Surface_tension") != NULL)
          {
              printf("found surface_tension in model %s \n", modelContainer().model(iModel)->name());
              surface_tension=modelContainer().model(iModel);
              printf("with value: %g \n", surface_tension->value());
          }

          if(strstr(modelContainer().model(iModel)->name(), "Film_flow") != NULL)
          {
              printf("found film_flow in model %s \n", modelContainer().model(iModel)->name());
              film_flow=modelContainer().model(iModel);
              printf("with value: %g \n", film_flow->value());
          }


      } //end have found model eqn with appropriate namespace


  } //end loop over model

  // parse boundary conditions and id of storage
  int iarg = 0;
  ICvalue=-999;
  BC[0]=BC[1]=-1;
  BCvalue[0]=BCvalue[1]=-999;

  char BCname[10];
  char SpeciesType[10];

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
    };

    if(eqnType_==SPECIES)
    {
        sprintf(SpeciesType,"solid");
        if (strcmp(arg[iarg],SpeciesType)==0)
            inPhase_ = SOLID;

        sprintf(SpeciesType,"gas");
        if (strcmp(arg[iarg],SpeciesType)==0)
            inPhase_ = GAS;

        sprintf(SpeciesType,"liquid");
        if (strcmp(arg[iarg],SpeciesType)==0)
            inPhase_ = LIQUID;

        if (strcmp(arg[iarg],"updatePhaseFraction")==0)
            updatePhaseFraction = true;

        if (strcmp(arg[iarg],"averagePhaseFraction")==0)
            averagePhaseFraction = true;
    };

    if (strcmp(arg[iarg],"phaseFractionMinimum")==0)
        phaseFractionMinimum = atof(arg[iarg+1]);

    if (strcmp(arg[iarg],"phaseFractionMinimumConvection")==0)
        phaseFractionMinimumConvection = atof(arg[iarg+1]);

    if (strcmp(arg[iarg],"convectionBoundMinStefanFactor")==0)
        convectionBoundMinStefanFactor = atof(arg[iarg+1]);

    if (strcmp(arg[iarg],"inactive")==0)
         solveMe = false;
 
    if (strcmp(arg[iarg],"solveConvectiveFlux")==0)
         solveConvectiveFlux = true;
            
    if (strcmp(arg[iarg],"writeDebugContainers")==0)
         writeDebugContainers = true;

    if (strcmp(arg[iarg],"normalizeDuringGrowth")==0)
         normalizeDuringGrowth = true;

    if (strcmp(arg[iarg],"verbose")==0)
         verbose_ = true;

    iarg++;
  };

  // Error Checks 
  if(updatePhaseFraction && eqnType_!=SPECIES)
    error().throw_error_one(FLERR,"In case 'updatePhaseFraction' = true, you must have a species eqn.\n");

  if( ((modelingApproach==CONTINUUM) && BC[0]==-1) || BC[1]==-1)
  {
    printf("ERROR: please specify one of the following BCs:\n%s(%d)\n%s(%d)\n%s(%d)\n",
            "NEUMANN",    NEUMANN,
            "DIRICHLET",  DIRICHLET,
            "CONVECTIVE", CONVECTIVE);

    error().throw_error_one(FLERR,"Boundary conditions not specified.\n");
  }

  if((modelingApproach==CONTINUUM) && inPhase_==-1 &&eqnType_==SPECIES)
  {
    printf("ERROR: please specify the type (solid, gas or fluid) of your species equation!\n");
    error().throw_error_one(FLERR,"Species type not specified.\n");
  }

  if(pore_radius!=NULL && molar_mass==NULL)
  {  
    printf("Since you specefied a pore radius and a molar mass, Knudsen Diffusion is activated\n");
  }

  if(pore_radius==NULL && molar_mass!=NULL)
  {  
     error().throw_error_one(FLERR,"You specified a molar mass in one of you species Eqn. In order to activate Knudsen diffusion you must also set a pore radius!\n");
  }

  if(pore_radius!=NULL && molar_mass==NULL)
  {  
     error().throw_error_one(FLERR,"You specified a pore radius in one of you species Eqn. In order to activate Knudsen diffusion you must also set a molar mass!\n");
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
