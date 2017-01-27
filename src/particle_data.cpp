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

#include "particle_mesh.h"
#include "particle_data.h"
#include "model_eqn_container.h"
#include "model_phasechange_container.h"
#include "model_chemistry_container.h"
#include "input.h"
#include "comm.h"
#include "output.h"
#include "mpi_pascal.h"
#include "model_eqn.h"

using namespace PASCAL_NS;

/* ----------------------------------------------------------------------
   ParticleData Constructor
------------------------------------------------------------------------- */

ParticleData::ParticleData(ParScale *ptr) :
   ParScaleBase(ptr),
   particle_data_tracker_ (*new CustomValueTracker(ptr,this)),
   id_ (*particle_data_tracker_.addElementProperty< ContainerScalar<int> >
            (
                "id", "comm_exchange_borders", "frame_general", "restart_yes", "coupling_pull",
                "read_no", "output_yes"
            )
        ),
   radius_(*particle_data_tracker_.addElementProperty< ContainerScalar<double> >
             (
                 "radius","comm_exchange_borders","frame_general","restart_yes","coupling_pull_push",
                 "read_yes", "output_no"
             )
          ),
   runtimeContainersExist_(false)
{

    nBodyPerProcess_ = NULL;
    nBodyInProcessesBelowMe_ = 0;
    verbose_    = false;
    hasPulled_  = false;

    haveGasPhase=false;
    haveLiquidPhase=false;
    numberOfPhases = 1; //by default we assume a single solid phase
    eqnIdFirstLiquid = -1;

    referencePressure = 1e5;

    haveConvectiveFluxGasPhase    = false;
    haveConvectiveFluxLiquidPhase = false;

    writeDebugGasPhase    = false;
    writeDebugLiquidPhase = false;

}

ParticleData::~ParticleData()
{
    delete &particle_data_tracker_;
    destroy<int>(nBodyPerProcess_);
    destroy<double>(tmpSourceVariable_);
}

/////////////////////////////////////////////////////////////////////////////
                          // MEMBER functions
/////////////////////////////////////////////////////////////////////////////

/* ----------------------------------------------------------------------
   Allocate all necessary data in containers
------------------------------------------------------------------------- */

void ParticleData::allocate()
{
    addRunTimeContainers();
    data().allocate();


   //Setting GLOBAL IDs to -1 (i.e., indicate it is undefined)
   int _id_length = id_.size();
   for(int iter=0; iter < _id_length; iter++)
        id_.begin()[iter] = -1;

   return;
}

/* ----------------------------------------------------------------------
 Called via coupling
------------------------------------------------------------------------- */
void ParticleData::grow_nbody(int _nbody,int _nbody_all)
{
    data().grow(_nbody+1);
    data().set_n_body(_nbody,_nbody_all);
}

/* ----------------------------------------------------------------------
   Read all necessary global and particle properties so all MPI procs have them
------------------------------------------------------------------------- */

void ParticleData::read()
{
    OperationProperties op(OPERATION_READ,false,false,false);
    if(verbose_)
        output().write_screen_all("ParticleData is now reading your input...");
    data().read(op);
}

/* ----------------------------------------------------------------------
   Pull properties from LIGGGHTS
   called after Comm::pull() and Coupling::pull()
------------------------------------------------------------------------- */

void ParticleData::pull()
{
    //Read particle information from coupling
    OperationProperties op(OPERATION_PULL,false,false,false);
    if(verbose_)
        output().write_screen_all("ParticleData attempts to pull ...");

    data().pull(op);

    hasPulled_ = true;

}

/* ----------------------------------------------------------------------
   Write all necessary global properties so all MPI procs have them
------------------------------------------------------------------------- */

void ParticleData::write()
{
    //Write particle information to file
    OperationProperties op(OPERATION_OUTPUT,false,false,false);
    if(verbose_)
        output().write_screen_all("ParticleData attempts to write some output ...");
    data().write(op);
}

/* ----------------------------------------------------------------------
   Push properties to LIGGGHTS
------------------------------------------------------------------------- */

void ParticleData::push()
{
    //Send particle information to coupling
    OperationProperties op(OPERATION_PUSH,false,false,false);
    if(verbose_)
        output().write_screen_all("ParticleData attempts to push ...");

    data().push(op);
}

/* ----------------------------------------------------------------------
   BCast global properties so all MPI procs have them
------------------------------------------------------------------------- */

void ParticleData::scatter()
{
    // is done in Comm::bcast(), not here
}

/* ----------------------------------------------------------------------
   Tasks important for parallelization of the simulation
------------------------------------------------------------------------- */

void ParticleData::parallelize()
{
    // nothing to do here

//    if(comm().nprocs() > 1)
//        error().throw_error_one(FLERR,"TODO: check if all per-particle data in correct subdomain");
}

/* ----------------------------------------------------------------------
   Intitialization phase - done before tun
------------------------------------------------------------------------- */

void ParticleData::parse_command(int narg,char const* const* arg)
{

  int iarg = 0;

  while (iarg < narg)
  {
      if (strcmp(arg[iarg],"number_particles")==0)
      {
            iarg++;

            if(narg - iarg < 1)
                error().throw_error_all(FLERR,"Not enough arguments for command: ",arg[iarg-1]);

            int nbody = atoi(arg[iarg]);
            data().set_n_body(nbody,nbody);

            iarg ++;
      }
      else if (strcmp(arg[iarg],"verbose")==0)
      {
            iarg++;

            if(narg - iarg < 1)
                error().throw_error_all(FLERR,"Not enough arguments for command: ",arg[iarg-1]);

            verbose_ = bool(atoi(arg[iarg]));
            iarg ++;
      }
      else
        error().throw_error_all(FLERR,"Unknown command: ",arg[iarg]);
  }
}

/* ----------------------------------------------------------------------
   Intitialization phase - done before run
------------------------------------------------------------------------- */

void ParticleData::init()
{

    //set IDs to value in case it is not pulled from coupling
    countBodiesOnMachine();
    resetIds();

    return;
}

/* ----------------------------------------------------------------------
   Reset the global particle IDs
------------------------------------------------------------------------- */
void ParticleData::countBodiesOnMachine() const
{

    char msgstr[500];
    sprintf(msgstr,"Counting bodies on the machine. using communicator: %d. \n",
            comm().nprocs()
           );
    if(verbose_)
       output().write_screen_all(msgstr);

    //Compute global particle counts
    if(nBodyPerProcess_)
        destroy<int>(nBodyPerProcess_);
    create<int>(nBodyPerProcess_,comm().nprocs());
    nBodyInProcessesBelowMe_ = 0;

    if(comm().nprocs()>1)
    {
    for(int iProc=0; iProc<comm().nprocs(); iProc++)
    {
        int currNBody=0;
        if(iProc==comm().me())
            currNBody = nbody();

//        MPI_Sum_Scalar(currNBody,comm().world()); //do an allreduce TODO-DEFECTIVE

        nBodyPerProcess_[iProc] = currNBody;
        if( iProc < comm().me() )
            nBodyInProcessesBelowMe_ += currNBody;
    }
    }
    else
        nBodyPerProcess_[0] = nbody();

   if(verbose_)
       printf("**me: %d /%d : nBodyPerProcess computed to be: %d, bodies below me: %d \n",
              comm().me(), comm().nprocs(),
              nBodyPerProcess_[comm().me()], nBodyInProcessesBelowMe_
              );

}


/* ----------------------------------------------------------------------
   Reset the global particle IDs
------------------------------------------------------------------------- */
void ParticleData::resetIds() const
{
    if(verbose_)
       printf("Resetting ids for process %d. \n",comm().me());

    for(int iP=0;iP < nbody(); iP++)
    {
        id_.begin()[iP] = iP+nBodyInProcessesBelowMe_+1; //by default, IDs start with 1
    }
}

/* ----------------------------------------------------------------------
   Reset the global particle IDs
------------------------------------------------------------------------- */
void ParticleData::scanPhaseInformation() const
{
    printf("ParticleData: scanning %d ModelEqns in modelEqnContainer() for phase information. \n",
           modelEqnContainer().nrEqns());

    int myPhase;
    int eqn_type;

    for(int iEqn=0; iEqn < modelEqnContainer().nrEqns(); iEqn++)
    {
        myPhase    = modelEqnContainer().modelEqn(iEqn)->return_myPhase();
        eqn_type   = modelEqnContainer().modelEqn(iEqn)->return_eqn_type();
        bool writeDebug = modelEqnContainer().modelEqn(iEqn)->writeDebugContainers;

        bool solveConvective = modelEqnContainer().modelEqn(iEqn)->solveConvectiveFlux;

        if(solveConvective && myPhase==GAS && eqn_type == 1)
            haveConvectiveFluxGasPhase = true;

        if(solveConvective && myPhase==LIQUID && eqn_type == 1)
        {
            haveConvectiveFluxLiquidPhase = true;
            if(eqnIdFirstLiquid<0)
                eqnIdFirstLiquid = iEqn;
        }

        if(writeDebug && myPhase==GAS && eqn_type == 1)
            writeDebugGasPhase    = true;

        if(writeDebug && myPhase==LIQUID && eqn_type == 1)
            writeDebugLiquidPhase = true;

        printf("Yor are looking at a %i model eqn, if species is in phase %i!\n", eqn_type, myPhase);

        //count phase fractions to be added
        if(eqn_type == 1)
        {
            if(myPhase==GAS && !haveGasPhase)
            {
                haveGasPhase=true;
                numberOfPhases++;
            }

            if(myPhase==LIQUID && !haveLiquidPhase)
            {
                haveLiquidPhase=true;
                numberOfPhases++;
            }
        }

    }

    runtimeContainersExist_ = false;
}


/* ----------------------------------------------------------------------
   add containers that need information to be read at run-time (script)
------------------------------------------------------------------------- */
void ParticleData::addRunTimeContainers() const
{

    //Scan all the phase information
    scanPhaseInformation();

    char outputSettingDebug [50];

    printf("ParticleData: will allocate memory for %d ModelEqns in modelEqnContainer(). \n", modelEqnContainer().nrEqns());
    for(int iEqn=0; iEqn < modelEqnContainer().nrEqns(); iEqn++)
    {
        modelEqnContainer().modelEqn(iEqn)->setParticleDataID(iEqn);
        printf("...allocating mem for modelEqn '%s' with particleDataId %d \n",
               modelEqnContainer().modelEqn(iEqn)->name(),
               modelEqnContainer().modelEqn(iEqn)->particleDataID()
              );
        char avName [50],fluxName [50],transCoeffName [50], chemName [50], chemJacobiName [50];

        if(modelEqnContainer().modelEqn(iEqn)->averagePhaseFraction)
            sprintf(avName,         "%sPhaseAv", modelEqnContainer().modelEqn(iEqn)->name());
        else
            sprintf(avName,         "%sAv", modelEqnContainer().modelEqn(iEqn)->name());
            
        sprintf(fluxName,       "%sFlux", modelEqnContainer().modelEqn(iEqn)->name());
        sprintf(transCoeffName, "%sTransCoeff", modelEqnContainer().modelEqn(iEqn)->name());
        sprintf(chemName,       "%sChem",  modelEqnContainer().modelEqn(iEqn)->name());
        sprintf(chemJacobiName, "%sChemJac", modelEqnContainer().modelEqn(iEqn)->name());

        char readSetting [50], outputSetting [50], outputSettingDebug [50];
        char couplingSettingIntra [50], couplingSettingAv [50], couplingSettingFlux [50], couplingSettingCoeff [50];
        sprintf(outputSetting,         "output_yes");

        sprintf(outputSettingDebug,    "output_no");
        if(modelEqnContainer().modelEqn(iEqn)->writeDebugContainers)
            sprintf(outputSettingDebug,    "output_yes");

        sprintf(readSetting,           "read_yes");
        sprintf(couplingSettingIntra,  "coupling_push_min_max_pull_reset");
        sprintf(couplingSettingAv,     "coupling_push");
        sprintf(couplingSettingFlux,   "coupling_pull_push");
        sprintf(couplingSettingCoeff,  "coupling_pull");
        if(!modelEqnContainer().modelEqn(iEqn)->solveMe)
        {
              sprintf(outputSetting, "output_no");
              sprintf(readSetting,   "read_yes");
              sprintf(couplingSettingIntra,  "coupling_none");
              sprintf(couplingSettingAv,     "coupling_none");
              sprintf(couplingSettingFlux,   "coupling_none");
              sprintf(couplingSettingCoeff,  "coupling_none");
        }


        intraPartMem_.push_back
        (
            particle_data_tracker_.addElementProperty< ContainerCvode<double,N_HISTORY_DEF,N_INTRA_GRIDPOINTS_DEF> >
                 (
                     modelEqnContainer().modelEqn(iEqn)->name(),
                     "comm_exchange_borders",
                     "frame_general","restart_yes",couplingSettingIntra,
                     readSetting, outputSetting
                 )
        );

        chemistryMem_.push_back
        (
            particle_data_tracker_.addElementProperty< ContainerChemistry<double,1,N_INTRA_GRIDPOINTS_DEF> >
                 (
                     chemName,
                     "comm_exchange_borders",
                     "frame_general","restart_yes","coupling_none",
                     "read_no", outputSettingDebug
                 )
        );

        chemistryMemJacobi_.push_back
        (
            particle_data_tracker_.addElementProperty< ContainerChemistry<double,1,N_INTRA_GRIDPOINTS_DEF> >
                 (
                     chemJacobiName,
                     "comm_exchange_borders",
                     "frame_general","restart_yes","coupling_none",
                     "read_no", outputSettingDebug
                 )
        );
        intraPartAv_.push_back
        (
             particle_data_tracker_.addElementProperty< ContainerScalar<double> >
                (
                    avName,
                    "comm_exchange_borders",
                    "frame_general","restart_yes",couplingSettingAv,
                    "read_no", outputSetting
                )
        );

        intraPartFlux_.push_back
        (
             particle_data_tracker_.addElementProperty< ContainerScalar<double> >
                (
                    fluxName,
                    "comm_exchange_borders",
                    "frame_general","restart_yes",couplingSettingFlux,
                    "read_no", outputSettingDebug
                )
        );

        intraPartTransCoeff_.push_back
        (
             particle_data_tracker_.addElementProperty< ContainerScalar<double> >
                (
                    transCoeffName,
                    "comm_exchange_borders",
                    "frame_general","restart_yes",couplingSettingCoeff,
                    "read_no", outputSettingDebug
                )
        );
    }

    //generate phaseFraction containers
    if(haveGasPhase)
    {
        phaseList.push_back(GAS);

        sprintf(outputSettingDebug,    "output_no");
        if(writeDebugGasPhase)
            sprintf(outputSettingDebug,    "output_yes");

        printf("creating phase fraction mem for gas species equation now \n");
        intraPhaseFractionMem_.push_back
        (
            particle_data_tracker_.addElementProperty< ContainerCvode<double,1,N_INTRA_GRIDPOINTS_DEF> >
            (
                "gasPhaseFraction",
                "comm_exchange_borders",
                "frame_general","restart_yes","coupling_push_min_max_pull_reset",
                "read_yes", "output_yes"
            )
        );

        printf("creating phase change rate mem / Jacobi for gas species equation now \n");
        intraPhaseChangeRateMem_.push_back
        (
            particle_data_tracker_.addElementProperty< ContainerCvode<double,1,N_INTRA_GRIDPOINTS_DEF> >
            (
                "gasPhaseChangeRate",
                "comm_exchange_borders",
                "frame_general","restart_yes","coupling_none",
                "read_no", outputSettingDebug
            )
        );
        intraPhaseChangeRateMemJacobi_.push_back
        (
            particle_data_tracker_.addElementProperty< ContainerCvode<double,1,N_INTRA_GRIDPOINTS_DEF> >
            (
                "gasPhaseChangeRateJac",
                "comm_exchange_borders",
                "frame_general","restart_yes","coupling_none",
                "read_no", outputSettingDebug
            )
        );
        intraPhaseChangeRateVolumetric_.push_back
        (
            particle_data_tracker_.addElementProperty< ContainerCvode<double,1,N_INTRA_GRIDPOINTS_DEF> >
            (
                "gasPhaseChangeRateVol",
                "comm_exchange_borders",
                "frame_general","restart_yes","coupling_none",
                "read_no", outputSettingDebug
            )
        );

        if(haveConvectiveFluxGasPhase || haveConvectiveFluxLiquidPhase)
        intraConvectiveFluxMem_.push_back
        (
            particle_data_tracker_.addElementProperty< ContainerCvode<double,1,N_INTRA_GRIDPOINTS_DEF> >
            (
                "gasConvectiveFlux",
                "comm_exchange_borders",
                "frame_general","restart_yes","coupling_none",
                "read_no", outputSettingDebug
            )
        );
    }
    if(haveLiquidPhase)
    {
        phaseList.push_back(LIQUID);
        sprintf(outputSettingDebug,    "output_no");
        if(writeDebugLiquidPhase)
            sprintf(outputSettingDebug,    "output_yes");

        printf("creating phase fraction mem for liquid species equation now. \n");
        intraPhaseFractionMem_.push_back
        (
            particle_data_tracker_.addElementProperty< ContainerCvode<double,1,N_INTRA_GRIDPOINTS_DEF> >
            (
                "liquidPhaseFraction",
                "comm_exchange_borders",
                "frame_general","restart_yes","coupling_push_min_max_pull_reset",
                "read_yes", "output_yes"
            )
        );

        printf("creating phase change rate mem / Jacobi for liquid species equation now \n");
        intraPhaseChangeRateMem_.push_back
        (
            particle_data_tracker_.addElementProperty< ContainerCvode<double,1,N_INTRA_GRIDPOINTS_DEF> >
            (
                "liquidPhaseChangeRate",
                "comm_exchange_borders",
                "frame_general","restart_yes","coupling_none",
                "read_no", outputSettingDebug
            )
        );
        intraPhaseChangeRateMemJacobi_.push_back
        (
            particle_data_tracker_.addElementProperty< ContainerCvode<double,1,N_INTRA_GRIDPOINTS_DEF> >
            (
                "liquidPhaseChangeRateJac",
                "comm_exchange_borders",
                "frame_general","restart_yes","coupling_none",
                "read_no", outputSettingDebug
            )
        );
        intraPhaseChangeRateVolumetric_.push_back
        (
            particle_data_tracker_.addElementProperty< ContainerCvode<double,1,N_INTRA_GRIDPOINTS_DEF> >
            (
                "liquidPhaseChangeRateVol",
                "comm_exchange_borders",
                "frame_general","restart_yes","coupling_none",
                "read_no", outputSettingDebug
            )
        );

        if( haveConvectiveFluxGasPhase || haveConvectiveFluxLiquidPhase )
        intraConvectiveFluxMem_.push_back
        (
            particle_data_tracker_.addElementProperty< ContainerCvode<double,1,N_INTRA_GRIDPOINTS_DEF> >
            (
                "liquidConvectiveFlux",
                "comm_exchange_borders",
                "frame_general","restart_yes","coupling_none",
                "read_no", outputSettingDebug
            )
        );
    }
    if(numberOfPhases>1) //only needed if there are at least 2 phases
    {
        sprintf(outputSettingDebug,    "output_no");
        if(writeDebugLiquidPhase || writeDebugGasPhase)
            sprintf(outputSettingDebug,    "output_yes");

      //push back the solid phase
      intraPhaseChangeRateMem_.push_back
      (
        particle_data_tracker_.addElementProperty< ContainerCvode<double,1,N_INTRA_GRIDPOINTS_DEF> >
        (
                "solidPhaseChangeRate",
                "comm_exchange_borders",
                "frame_general","restart_yes","coupling_none",
                "read_no", outputSettingDebug
        )
      );
      intraPhaseChangeRateMemJacobi_.push_back
      (
        particle_data_tracker_.addElementProperty< ContainerCvode<double,1,N_INTRA_GRIDPOINTS_DEF> >
        (
                "solidPhaseChangeRateJac",
                "comm_exchange_borders",
                "frame_general","restart_yes","coupling_none",
                "read_no", outputSettingDebug
        )
      );

      intraPhaseChangeRateVolumetric_.push_back
      (
            particle_data_tracker_.addElementProperty< ContainerCvode<double,1,N_INTRA_GRIDPOINTS_DEF> >
            (
                "solidPhaseChangeRateVol",
                "comm_exchange_borders",
                "frame_general","restart_yes","coupling_none",
                "read_no", outputSettingDebug
            )
      );

      //A container for previous-time temperature needed for Evaluating the Jacobi of phase-change models
      //this must be the last container in the memory for the Jacobi!
      intraPhaseChangeRateMemJacobi_.push_back
      (
        particle_data_tracker_.addElementProperty< ContainerCvode<double,1,N_INTRA_GRIDPOINTS_DEF> >
        (
                "temperaturePrevious",
                "comm_exchange_borders",
                "frame_general","restart_no","coupling_none",
                "read_no", "output_no"
        )
      );

    }
    phaseList.push_back(SOLID);

    generatePhaseIDMap();

    //Empty data pointers to current state (to be filled later)
    datapointer_.resize(modelEqnContainer().nrEqns());
    datapointerJac_.resize(modelEqnContainer().nrEqns());
    datapointerPhaseFraction_.resize(numberOfPhases-1);

    // Set containers to correct vector length (only for intra-particle containers)
    for(int iEqn = 0; iEqn < modelEqnContainer().nrEqns(); iEqn++)
    {
        intraPartMem_[iEqn]->setLenVecUsed( modelEqnContainer().modelEqn(iEqn)->nGridPointsUsed()  );
        chemistryMem_[iEqn]->setLenVecUsed( modelEqnContainer().modelEqn(iEqn)->nGridPointsUsed()  );
        chemistryMemJacobi_[iEqn]->setLenVecUsed( modelEqnContainer().modelEqn(iEqn)->nGridPointsUsed()  );

        if(intraPartMem_[iEqn]->isFull())
            error().throw_error_one(FLERR,
                 "You overfilled the 'intraPartMem' container. Increase the value for N_INTRA_GRIDPOINTS_DEF");
        if(chemistryMem_[iEqn]->isFull())
            error().throw_error_one(FLERR,
                 "You overfilled the 'chemistryMem' container. Increase the value for N_INTRA_GRIDPOINTS_DEF");
        if(chemistryMemJacobi_[iEqn]->isFull())
            error().throw_error_one(FLERR,
                 "You overfilled the 'chemistryMemJacobi' container. Increase the value for N_INTRA_GRIDPOINTS_DEF");
    }

    //Set phase fraction container
    for(uint iPhase = 0; iPhase < intraPhaseFractionMem_.size(); iPhase++) //Only save phase fraction for fluid phases!
    {
        intraPhaseFractionMem_[iPhase]->setLenVecUsed(modelEqnContainer().modelEqn(0)->nGridPointsUsed() );
        if(intraPhaseFractionMem_[iPhase]->isFull())
            error().throw_error_one(FLERR, "You overfilled the 'intraPhaseFractionMem_' container. Increase the value for N_INTRA_GRIDPOINTS_DEF");

        if( haveConvectiveFluxGasPhase || haveConvectiveFluxLiquidPhase )
        {
          intraConvectiveFluxMem_[iPhase]->setLenVecUsed(modelEqnContainer().modelEqn(0)->nGridPointsUsed() );

          if(intraConvectiveFluxMem_[iPhase]->isFull())
            error().throw_error_one(FLERR, "You overfilled the 'intraConvectiveFluxMem_' container. Increase the value for N_INTRA_GRIDPOINTS_DEF");
        }

    }
    for(uint iPhase = 0; iPhase < intraPhaseChangeRateMem_.size(); iPhase++)
    {
        intraPhaseChangeRateMem_[iPhase]->setLenVecUsed(modelEqnContainer().modelEqn(0)->nGridPointsUsed() );
        if(intraPhaseChangeRateMem_[iPhase]->isFull())
            error().throw_error_one(FLERR, "You overfilled the 'intraPhaseChangeRateMem_' container. Increase the value for N_INTRA_GRIDPOINTS_DEF");
    }
    for(uint iPhase = 0; iPhase < intraPhaseChangeRateMemJacobi_.size(); iPhase++)
    {
        intraPhaseChangeRateMemJacobi_[iPhase]->setLenVecUsed(modelEqnContainer().modelEqn(0)->nGridPointsUsed() );
        if(intraPhaseChangeRateMemJacobi_[iPhase]->isFull())
            error().throw_error_one(FLERR, "You overfilled the 'intraPhaseFractionMemJacobi_' container. Increase the value for N_INTRA_GRIDPOINTS_DEF");
    }
    for(uint iPhase = 0; iPhase < intraPhaseChangeRateVolumetric_.size(); iPhase++)
    {
        intraPhaseChangeRateVolumetric_[iPhase]->setLenVecUsed(modelEqnContainer().modelEqn(0)->nGridPointsUsed() );
        if(intraPhaseChangeRateVolumetric_[iPhase]->isFull())
            error().throw_error_one(FLERR, "You overfilled the 'intraPhaseChangeRateVolumetric_' container. Increase the value for N_INTRA_GRIDPOINTS_DEF");
    }

    tmpSourceVariable_ = create<double>(tmpSourceVariable_,
                                        modelEqnContainer().modelEqn(0)->nGridPointsUsed()
                                       );

    //set to true since container handles everything
    runtimeContainersExist_ = true;

    output().write_screen_one("ParticleData allocation completed. \n");
}

/* ----------------------------------------------------------------------
   Set Data pointer and acess Intra-Particle data
------------------------------------------------------------------------- */
void ParticleData::setParticleIDPointer(int particleDataID_, int particleID)
{
    current_particleDataID_ = particleDataID_;
    current_particleID_     = particleID;

    if(!runtimeContainersExist_)
        error().throw_error_one(FLERR,"Internal error: requesting intra particle data but containers not existing");

    ptr3 = (double***)  intraPartMem_[particleDataID_]->begin_slow_dirty();   //ptr to vectorial per-particle data
    datapointer_[particleDataID_] = &ptr3[particleID][0][0]; //TODO: put into separate function
//    printf("pointer for particleDataID %i and particle ID %i with length %d set\n ",
//             particleDataID_,
//             particleID,
//             intraPartMem_[particleDataID_]->lenVecUsed()
//          );
}

/* ----------------------------------------------------------------------
   Set Data pointer and acess Intra-Particle data
------------------------------------------------------------------------- */
void ParticleData::setParticleIDPointerPhaseFraction(int _particleID)
{
    if(!runtimeContainersExist_)
        error().throw_error_one(FLERR,"Internal error: requesting intra particle data but containers not existing");

    int phaseId=0;
    if(haveGasPhase)
    {
        ptr3GasPhase_ = (double***)  intraPhaseFractionMem_[phaseId]->begin_slow_dirty();     //ptr to vectorial per-particle data
        datapointerPhaseFraction_[phaseId] = &ptr3GasPhase_[_particleID][0][0];
        phaseId++;
    }
    if(haveLiquidPhase)
    {
        ptr3LiquidPhase_ = (double***)  intraPhaseFractionMem_[phaseId]->begin_slow_dirty();  //ptr to vectorial per-particle data
        datapointerPhaseFraction_[phaseId] = &ptr3LiquidPhase_[_particleID][0][0];
        phaseId++;
    }
}

/* ---------------------------------------------------------------------- */
void ParticleData::resetParticleIDPointer(int particleDataID_)
{
    if(!runtimeContainersExist_)
        error().throw_error_one(FLERR,"Internal error: requesting intra particle data but containers not existing");

    current_particleDataID_ = particleDataID_;

    ptr3 = (double***)  intraPartMem_[current_particleDataID_]->begin_slow_dirty();   //ptr to vectorial per-particle data
    datapointer_[current_particleDataID_] = &ptr3[current_particleDataID_][0][0];
}

/* ---------------------------------------------------------------------- */
void ParticleData::returnDataPoint(double & data, int _j)
{
    data = ptr3[current_particleID_][0][_j]; //get data
    return;
}

/* ---------------------------------------------------------------------- */
void ParticleData::setDataPoint(double data, int _j)
{
    ptr3[current_particleID_][0][_j] = data; //set data
    return;
}

/* ---------------------------------------------------------------------- */
void ParticleData::retrieveIntraData(vector<int> _dataIDs, int _particleID, int _gridPoint, vector<double>& output)
{
    static double *** myPtr;
    for(uint i = 0; i<_dataIDs.size();i++)
    {
        myPtr = (double***)  intraPartMem_[_dataIDs[i]]->begin_slow_dirty();
        output[i] = myPtr[_particleID][0][_gridPoint];
    }
}

/* ---------------------------------------------------------------------- */
double ParticleData::retrieveIntraData(int _dataID, int _particleID, int _gridPoint)
{
    static double *** myPtr;
    myPtr = (double***)  intraPartMem_[_dataID]->begin_slow_dirty();
    return myPtr[_particleID][0][_gridPoint];
}

/* ---------------------------------------------------------------------- */
void ParticleData::retrievePhaseFractionData(vector<int> _dataIDs, int _particleID, int _gridPoint, vector<double>& output)
{
    static double *** myPtr;
    for(uint i = 0; i<_dataIDs.size();i++)
    {
        myPtr = (double***)  intraPhaseFractionMem_[phaseIDMap[_dataIDs[i]]]->begin_slow_dirty();
        output[i] = myPtr[_particleID][0][_gridPoint];
    }
}

/* ---------------------------------------------------------------------- */
double ParticleData::retrievePhaseFractionData(int _dataID, int _particleID, int _gridPoint)
{
    static double *** myPtr;
    myPtr = (double***)  intraPhaseFractionMem_[phaseIDMap[_dataID]]->begin_slow_dirty();
    return myPtr[_particleID][0][_gridPoint];
}

/* ---------------------------------------------------------------------- */
void ParticleData::returnchemistryDataPoint(int _particleDataID, int particleID, int _gridPoint, double & data)
{
    current_particleDataID_ = _particleDataID;
    current_particleID_     = particleID;

    if(!runtimeContainersExist_)
        error().throw_error_one(FLERR,"Internal error: requesting intra particle data, but containers not existing");

    ptr3 = (double***)  chemistryMem_[current_particleDataID_]->begin_slow_dirty();   //ptr to vectorial per-particle data
    datapointer_[current_particleDataID_] = &ptr3[current_particleID_][0][0]; //TODO: put into separate function
    data = ptr3[current_particleID_][0][_gridPoint]; //get data

    return;
}


/* ---------------------------------------------------------------------- */
void ParticleData::returnchemistryJacDataPoint(int _particleDataID, int particleID, int _gridPoint, double & data)
{
    current_particleDataID_ = _particleDataID;
    current_particleID_     = particleID;

    if(!runtimeContainersExist_)
        error().throw_error_one(FLERR,"Internal error: requesting intra particle data, but containers not existing");

    ptr3Jac = (double***)  chemistryMemJacobi_[current_particleDataID_]->begin_slow_dirty();   //ptr to vectorial per-particle data
    datapointerJac_[current_particleDataID_] = &ptr3Jac[current_particleID_][0][0]; //TODO: put into separate function
    data = ptr3Jac[current_particleID_][0][_gridPoint]; //get data

    return;
}

/* ---------------------------------------------------------------------- */
void ParticleData::returnIntraData(double * data)
{

    for(int j=0;j < intraPartMem_[current_particleDataID_]->lenVecUsed();j++) //Must only loop over used IDs!
    {
        data[j] = ptr3[current_particleID_][0][j]; //get data
    }
    return;
}

/* ---------------------------------------------------------------------- */
void ParticleData::returnPhaseFractionData(double * phaseFracGas, double * phaseFracLiq)
{

    for(int j=0;j < intraPartMem_[0]->lenVecUsed();j++) //Must only loop over used IDs!
    {
       if(haveGasPhase)
           phaseFracGas[j] = ptr3GasPhase_[current_particleID_][0][j]; //get data
       else
           phaseFracGas[j] = 0;

       if(haveLiquidPhase)
           phaseFracLiq[j] = ptr3LiquidPhase_[current_particleID_][0][j]; //get data
       else
           phaseFracLiq[j] = 0;
    }

    return;
}

/* ---------------------------------------------------------------------- */
void ParticleData::returnPhaseFractionDataGridpoint(double &phaseFracGas, double &phaseFracLiq, int grid_point)
{
    if(haveGasPhase)
           phaseFracGas = ptr3GasPhase_[current_particleID_][0][grid_point]; //get data
       else
           phaseFracGas = 0;

    if(haveLiquidPhase)
           phaseFracLiq = ptr3LiquidPhase_[current_particleID_][0][grid_point]; //get data
       else
           phaseFracLiq = 0;

    return;
}
/* ---------------------------------------------------------------------- */
void ParticleData::returnPhaseChangeRateDataPoint(int _phaseId, int particleID, int _gridPoint, double & data)
{
    //input: phaseId is enueration refering to the phase (e.g., SOLID=0, GAS=1, LIQUID=2 )
    if(!runtimeContainersExist_)
        error().throw_error_one(FLERR,"Internal error: requesting intraPhaseChangeRate data, but containers not existing");
    if(phaseIDMap[_phaseId]<0)
        error().throw_error_one(FLERR,"ParticleData: phaseIDMap not set for the requested phase");

    double *** aTempPtr = (double***)  intraPhaseChangeRateMem_[phaseIDMap[_phaseId]]->begin_slow_dirty();   //ptr to vectorial per-particle data
    data = aTempPtr[particleID][0][_gridPoint]; //get data

    return;
}

/* ---------------------------------------------------------------------- */
void ParticleData::returnPhaseChangeRateJacDataPoint(int _phaseId, int particleID, int _gridPoint, double & data)
{
    if(!runtimeContainersExist_)
        error().throw_error_one(FLERR,"Internal error: requesting intraPhaseChangeRateJac data, but containers not existing");
    if(phaseIDMap[_phaseId]<0)
        error().throw_error_one(FLERR,"ParticleData: phaseIDMap not set for the requested phase");

    double *** aTempPtr = (double***)  intraPhaseChangeRateMemJacobi_[phaseIDMap[_phaseId]]->begin_slow_dirty();   //ptr to vectorial per-particle data
    data = aTempPtr[particleID][0][_gridPoint]; //get data

    return;
}


/* ---------------------------------------------------------------------- */
void ParticleData::returnPhaseChangeRateJacLastSlotDataPoint(int particleID, int _gridPoint, double & data)
{
    if(!runtimeContainersExist_)
        error().throw_error_one(FLERR,"Internal error: requesting intraPhaseChangeRateJac data, but containers not existing");

    double *** aTempPtr = (double***)  intraPhaseChangeRateMemJacobi_[intraPhaseChangeRateMemJacobi_.size()-1]->begin_slow_dirty();   //ptr to last lost
    data = aTempPtr[particleID][0][_gridPoint]; //get data

    return;
}

/* ---------------------------------------------------------------------- */
void ParticleData::returnIntraFlux(int particleID, double &flux)
{
    flux = intraPartFlux_[current_particleDataID_]->get(particleID);
    return;
}


/* ---------------------------------------------------------------------- */
void ParticleData::returnIntraTransCoeff(int particleID, double &transCoeff)
{
    transCoeff = intraPartTransCoeff_[current_particleDataID_]->get(particleID);
    return;
}

/* ----------------------------------------------------------------------
   Save Intra-Particle data
------------------------------------------------------------------------- */
void ParticleData::saveIntraParticleData(int particleDataID_, int particleID, double * data)
{
    if(!runtimeContainersExist_)
        error().throw_error_one(FLERR,"Internal error: trying to save intra particle data but containers not existing");

    //ave intra-particle data
    double *** ptr3 = (double***)  intraPartMem_[particleDataID_]->begin_slow_dirty();   //ptr to vectorial per-particle data

    for(int j=0;j < intraPartMem_[particleDataID_]->lenVecUsed();j++)
    {
           ptr3[particleID][0][j] = data[j]; //save data
    }

    return;
}

/* ---------------------------------------------------------------------- */
void ParticleData::saveChemistryParticleData(int _particleDataID, int particleID, double &data, int gridPoint)
{
    if(!runtimeContainersExist_)
        error().throw_error_one(FLERR,"Internal error: trying to save intra particle data but containers not existing");

    //ave intra-particle data
    double *** aPtr3 = (double***) chemistryMem_[_particleDataID]->begin_slow_dirty();   //ptr to vectorial per-particle data
    aPtr3[particleID][0][gridPoint] += data;      //save data

    return;
}

/* ----------------------------------------------------------------------
   Save phase fraction data
------------------------------------------------------------------------- */
void ParticleData::savePhaseFractionData(int particleID, double * dataGas, double * dataLiquid)
{
    if(!runtimeContainersExist_)
        error().throw_error_one(FLERR,"Internal error: trying to save intra particle data but containers not existing");

    int phaseId=0;
    if(haveGasPhase)
    {
        ptr3GasPhase_ = (double***)  intraPhaseFractionMem_[phaseId]->begin_slow_dirty();     //ptr to vectorial per-particle data
        for(int j=0;j < intraPhaseFractionMem_[phaseId]->lenVecUsed();j++)
            ptr3GasPhase_[particleID][0][j] = dataGas[j]; //save data

        phaseId++;
    }
    if(haveLiquidPhase)
    {
        ptr3LiquidPhase_ = (double***)  intraPhaseFractionMem_[phaseId]->begin_slow_dirty();  //ptr to vectorial per-particle data
        for(int j=0;j < intraPhaseFractionMem_[phaseId]->lenVecUsed();j++)
            ptr3LiquidPhase_[particleID][0][j] = dataLiquid[j]; //save data

        phaseId++;
    }

    return;
}

/* ----------------------------------------------------------------------
   Update time derivates of  phase fraction
------------------------------------------------------------------------- */
void ParticleData::phaseFractionChangeRateStart(int particleID, double * dataGas, double * dataLiquid)
{
    if(!runtimeContainersExist_)
        error().throw_error_one(FLERR,"Internal error: trying to save intra particle data but containers not existing");

    int phaseId=0;
    if(haveGasPhase)
    {
        ptr3GasPhase_ = (double***)  intraPhaseChangeRateVolumetric_[phaseId]->begin_slow_dirty();     //ptr to vectorial per-particle data
        for(int j=0;j < intraPhaseChangeRateVolumetric_[phaseId]->lenVecUsed();j++)
            ptr3GasPhase_[particleID][0][j] = -1*dataGas[j]; //add old

        phaseId++;
    }
    if(haveLiquidPhase)
    {
        ptr3LiquidPhase_ = (double***)  intraPhaseChangeRateVolumetric_[phaseId]->begin_slow_dirty();  //ptr to vectorial per-particle data
        for(int j=0;j < intraPhaseChangeRateVolumetric_[phaseId]->lenVecUsed();j++)
            ptr3LiquidPhase_[particleID][0][j]  = -1*dataLiquid[j]; //add old

        phaseId++;
    }

    return;
}

/* ---------------------------------------------------------------------- */
void ParticleData::phaseFractionChangeRateEnd(int particleID, double * dataGas, double * dataLiquid, double deltaT)
{
    if(!runtimeContainersExist_)
        error().throw_error_one(FLERR,"Internal error: trying to save intra particle data but containers not existing");

    int phaseId=0;
    if(haveGasPhase)
    {
        ptr3GasPhase_ = (double***)  intraPhaseChangeRateVolumetric_[phaseId]->begin_slow_dirty();     //ptr to vectorial per-particle data
        for(int j=0;j < intraPhaseChangeRateVolumetric_[phaseId]->lenVecUsed();j++)
        {
            ptr3GasPhase_[particleID][0][j] += dataGas[j]; //add new
            ptr3GasPhase_[particleID][0][j] /= deltaT;
        }

        phaseId++;
    }
    if(haveLiquidPhase)
    {
        ptr3LiquidPhase_ = (double***)  intraPhaseChangeRateVolumetric_[phaseId]->begin_slow_dirty();  //ptr to vectorial per-particle data
        for(int j=0;j < intraPhaseChangeRateVolumetric_[phaseId]->lenVecUsed();j++)
        {
            ptr3LiquidPhase_[particleID][0][j] += dataLiquid[j]; //add new
            ptr3LiquidPhase_[particleID][0][j] /= deltaT;
        }

        phaseId++;
    }

    return;
}

/* ---------------------------------------------------------------------- */
void ParticleData::savePhaseChangeParticleData(int _phaseID, int particleID, double &data, int gridPoint)
{
    if(!runtimeContainersExist_)
        error().throw_error_one(FLERR,"Internal error: trying to save phase change data but containers not existing");

    //ave intra-particle data
    double *** aPtr3 = (double***) intraPhaseChangeRateMem_[phaseIDMap[_phaseID]]->begin_slow_dirty();   //ptr to vectorial per-particle data
    aPtr3[particleID][0][gridPoint] += data;      //save data

    return;
}

/* ---------------------------------------------------------------------- */
void ParticleData::savePhaseChangeParticleDataJac(int _phaseID, int particleID, double &data, int gridPoint)
{
    if(!runtimeContainersExist_)
        error().throw_error_one(FLERR,"Internal error: trying to save phase change data Jacobi but containers not existing");

    //ave intra-particle data
    double *** aPtr3 = (double***) intraPhaseChangeRateMemJacobi_[phaseIDMap[_phaseID]]->begin_slow_dirty();   //ptr to vectorial per-particle data
    aPtr3[particleID][0][gridPoint] += data;      //save data

    return;
}


/* ---------------------------------------------------------------------- */
void ParticleData::savePhaseChangeParticleDataJacLastSlot(int particleID, double &data, int gridPoint)
{
    if(!runtimeContainersExist_)
        error().throw_error_one(FLERR,"Internal error: trying to save phase change data Jacobi but containers not existing");

    //ave intra-particle data
    double *** aPtr3 = (double***) intraPhaseChangeRateMemJacobi_[intraPhaseChangeRateMemJacobi_.size()-1]->begin_slow_dirty();   //ptr to vectorial per-particle data
    aPtr3[particleID][0][gridPoint] += data;      //save data

    return;
}

/* ---------------------------------------------------------------------- */
void ParticleData::saveChemistryJacParticleData(int _particleDataID, int particleID, double &data, int gridPoint)
{
    if(!runtimeContainersExist_)
        error().throw_error_one(FLERR,"Internal error: trying to save intra particle data but containers not existing");

    //ave intra-particle data
    double *** aPtr3 = (double***) chemistryMemJacobi_[_particleDataID]->begin_slow_dirty();   //ptr to vectorial per-particle data
    aPtr3[particleID][0][gridPoint] += data;      //save data

    return;
}

/* ---------------------------------------------------------------------- */
void ParticleData::saveIntraParticleAv(int particleDataID_, int particleID, double dataAv)
{
    //save averages
    intraPartAv_[particleDataID_]->set(particleID,dataAv);

    return;
}

/* ---------------------------------------------------------------------- */
void ParticleData::saveIntraParticleFlux(int particleDataID_, int particleID, double dataflux)
{
    //save fluxes
    intraPartFlux_[particleDataID_]->set(particleID,dataflux);

    return;
}


/* ---------------------------------------------------------------------- */
void ParticleData::saveIntraParticleTransCoeff(int particleDataID_, int particleID, double transCoeff)
{
    //save fluxes
    intraPartTransCoeff_[particleDataID_]->set(particleID,transCoeff);

    return;
}

/* ----------------------------------------------------------------------
   Clean phase change data
------------------------------------------------------------------------- */

void ParticleData::resetPhaseChangeSourceTerms()
{
    if(!runtimeContainersExist_)    return;
    for (uint m = 0; m < intraPhaseChangeRateMem_.size() ;m++) //particle Data ID
    {
        //printf("cleaning for %i model eqn \n", m);
        double *** aPtr3    = (double***) intraPhaseChangeRateMem_[m]->begin_slow_dirty();

        for (int i=0;i < nbody();i++)  //particle ID
            for(int j=0;j < intraPhaseChangeRateMem_[m]->lenVecUsed();j++)
                   aPtr3[i][0][j]    = 0; //reset to zero
    }

    for (uint m = 0; m < intraPhaseChangeRateMemJacobi_.size() ;m++) //particle Data ID
    {
        //printf("cleaning for %i model eqn \n", m);
        double *** bPtr3    = (double***) intraPhaseChangeRateMemJacobi_[m]->begin_slow_dirty();

        for (int i=0;i < nbody();i++)  //particle ID
            for(int j=0;j < intraPhaseChangeRateMemJacobi_[m]->lenVecUsed();j++)
                   bPtr3[i][0][j]    = 0; //reset to zero
    }

    for (uint m = 0; m < intraPhaseChangeRateVolumetric_.size() ;m++) //particle Data ID
    {
        //printf("cleaning for %i model eqn \n", m);
        double *** aPtr3    = (double***) intraPhaseChangeRateVolumetric_[m]->begin_slow_dirty();

        for (int i=0;i < nbody();i++)  //particle ID
            for(int j=0;j < intraPhaseChangeRateVolumetric_[m]->lenVecUsed();j++)
                   aPtr3[i][0][j]    = 0; //reset to zero
    }


    if(!runtimeContainersExist_)
        error().throw_error_one(FLERR,"Internal error: trying to clean phase change data but containers not existing");

    return;
}

/* ----------------------------------------------------------------------
   Clean chemistry data
------------------------------------------------------------------------- */

void ParticleData::resetChemicalSourceTerms()
{
    if(!runtimeContainersExist_)    return;

    for (int m=0; m<modelEqnContainer().nrEqns(); m++) // no need to pull out particle Data ID, since we do not care about order
    {
        //printf("cleaning for %i model eqn \n", m);
        double *** aPtr3    = (double***)   chemistryMem_[m]->begin_slow_dirty();
        double *** aPtrJac3 = (double***)   chemistryMemJacobi_[m]->begin_slow_dirty();

        for (int i=0;i < nbody();i++)  //particle ID
        {
            //printf("cleaning for %i particle \n", i);
            for(int j=0;j < chemistryMem_[m]->lenVecUsed();j++)
            {
                   aPtr3[i][0][j]    = 0; //reset to zero
                   aPtrJac3[i][0][j] = 0; //reset to zero
            }
        }
    }

    if(!runtimeContainersExist_)
        error().throw_error_one(FLERR,"Internal error: trying to clean chemical particle data but containers not existing");

    return;
}

/* ----------------------------------------------------------------------
   Update the convective Fluxes
------------------------------------------------------------------------- */
void ParticleData::generatePhaseIDMap() const
{
    phaseIDMap.assign(4,-1); //init negative, i.e., don't have phase

    printf("ParticleData::generatePhaseIDMap: phaseList: \n");
    for(uint iPhaseList=0; iPhaseList<phaseList.size(); iPhaseList++)
        printf("phase[%d]: %d \n", iPhaseList, phaseList[iPhaseList]);

    for(int iPhase=0; iPhase<4; iPhase++)
    {
        //search phase list and enter id of each phase
        for(uint iPhaseList=0; iPhaseList<phaseList.size(); iPhaseList++)
            if( phaseList[iPhaseList]==iPhase )
            {
                if(phaseIDMap[iPhase] != -1) //was set before, not possible
                    error().throw_error_one(FLERR,"ParticleData: Cannot generatePhaseIDMap since a phase was inserted twice.");
                phaseIDMap[iPhase]=iPhaseList;
            }
    }
    printf("ParticleData::generatePhaseIDMap: phaseIDMap: \n");
    for(int iPhase=0; iPhase<4; iPhase++)
        printf("phaseIDMap[%d]: %d \n", iPhase, phaseIDMap[iPhase]);
}

/* ----------------------------------------------------------------------
   Update the convective Fluxes
------------------------------------------------------------------------- */
void ParticleData::computeConvection()
{
    if(!runtimeContainersExist_)    return;
    if( !(haveConvectiveFluxGasPhase || haveConvectiveFluxLiquidPhase) ) return;

    for(uint iPhase = 0; iPhase < intraPhaseFractionMem_.size(); iPhase++) //Only save phase fraction for fluid phases!
    {
        int containerPhase = phaseList[iPhase];

        //loop over phase change models to identify involved species, in order to compute evaporation rate
        int speciesModelEqnlID=-1;
        const  vector<int>* phaseID;
        const  vector<int>* specModelEqn;
        for(int iPhaseChangeModel=0;iPhaseChangeModel<modelPhaseChangeContainer().nrPhaseChangeEqns(); iPhaseChangeModel++)
        {
            //double pDot    = 0.0;
            //double pDotJac = 0.0;
            phaseID = modelPhaseChangeContainer().modelPhaseChangeEqn(iPhaseChangeModel)->phaseID();
            specModelEqn = modelPhaseChangeContainer().modelPhaseChangeEqn(iPhaseChangeModel)->speciesModelEqnlID();

            //Check dense phase
            if( (*phaseID)[0]==containerPhase )
            {
              speciesModelEqnlID    = (*specModelEqn)[0];
              continue;
            }
            else if( (*phaseID)[1]==containerPhase )  //Check dilute phase
            {
              speciesModelEqnlID    = (*specModelEqn)[1];
              continue;
            }
        }

        //evaluate the phase flux for this phase
        modelEqnContainer().modelEqn(speciesModelEqnlID)->evaluatePhaseFlux();

    }

    return;
}

/* ----------------------------------------------------------------------
   Normalize the internal (species) fields
------------------------------------------------------------------------- */
void ParticleData::normalizeInternalFields(int _particleID, double _factor)
{
        for(int iEqn=0; iEqn < modelEqnContainer().nrEqns(); iEqn++)
        {
            if( !modelEqnContainer().modelEqn(iEqn)->normalizeDuringGrowth )
                continue;

            double*** aPtr = (double***)  intraPartMem_[iEqn]->begin_slow_dirty();
            for(int i=0; i<modelEqnContainer().modelEqn(iEqn)->nGridPointsUsed()-1; i++)
            {
                aPtr[_particleID][0][i] *= _factor;
            }
        }
}
