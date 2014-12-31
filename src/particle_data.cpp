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
#include "model_chemistry_container.h"
#include "input.h"
#include "comm.h"
#include "output.h"
#include "mpi_pascal.h"

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
    verbose_ = false;

}

ParticleData::~ParticleData()
{
 //TODO: clean data arrays if necessary
//   delete &particle_data_tracker_;
//   delete &radius_;
//   delete &intraPartAv0_;
//   delete &intraPartMem_[0];
//   delete &chemistryMem_[0];
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
        output().write_screen_one("\nParticleData is now reading your input...");
    data().read(op);
}

/* ----------------------------------------------------------------------
   Pull properties from LIGGGHTS
   called after Comm::pull() and Coupling::pull()
------------------------------------------------------------------------- */

void ParticleData::pull()
{
    //Read particle information from file
    OperationProperties op(OPERATION_PULL,false,false,false);
    if(verbose_)
        output().write_screen_one("\nParticleData attempts to pull ...");

    data().pull(op);
}

/* ----------------------------------------------------------------------
   Write all necessary global properties so all MPI procs have them
------------------------------------------------------------------------- */

void ParticleData::write()
{
    //Write particle information to file
    OperationProperties op(OPERATION_OUTPUT,false,false,false);
    if(verbose_)
        output().write_screen_one("ParticleData attempts to writing some output ...");
    data().write(op);
}

/* ----------------------------------------------------------------------
   Push properties to LIGGGHTS
------------------------------------------------------------------------- */

void ParticleData::push()
{
    //Read particle information from file
    OperationProperties op(OPERATION_PUSH,false,false,false);
    if(verbose_)
        output().write_screen_one("\nParticleData attempts to push ...");
        
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

    if(comm().nprocs() > 1)
        error().throw_error_one(FLERR,"TODO: check if all per-particle data in correct subdomain");
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
    if(verbose_)
       printf("Counting bodies on the machine. \n");
    
    //Compute global particle counts
    create<int>(nBodyPerProcess_,comm().nprocs());
    nBodyInProcessesBelowMe_ = 0;
    
    if(comm().nprocs()>1)
    {
    for(int iProc=0; iProc<comm().nprocs(); iProc++)
    { 
        int currNBody=0;
        if(iProc==comm().me())
            currNBody = nbody();
            
        MPI_Sum_Scalar(currNBody,comm().world()); //do an allreduce
        
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
   add containers that need information to be read at run-time (script)
------------------------------------------------------------------------- */
void ParticleData::addRunTimeContainers() const
{

    printf("ParticleData: will allocate memory for %d ModelEqns in modelEqnContainer(). \n", modelEqnContainer().nrEqns());
    for(int iEqn=0; iEqn < modelEqnContainer().nrEqns(); iEqn++)
    {

        modelEqnContainer().modelEqn(iEqn)->setParticleDataID(iEqn);
        printf("...allocating mem for modelEqn '%s' with particleDataId %d \n",
               modelEqnContainer().modelEqn(iEqn)->name(),
               modelEqnContainer().modelEqn(iEqn)->particleDataID()
              );
        char avName [50],fluxName [50], chemName [50]; int nString(0);
        nString = sprintf(avName,   "%sAv", modelEqnContainer().modelEqn(iEqn)->name());
        nString = sprintf(fluxName, "%sFlux", modelEqnContainer().modelEqn(iEqn)->name());
        nString = sprintf(chemName, "%sChem", modelEqnContainer().modelEqn(iEqn)->name());

        intraPartMem_.push_back
        (
            particle_data_tracker_.addElementProperty< ContainerCvode<double,N_HISTORY_DEF,N_INTRA_GRIDPOINTS_DEF> >
                 (
                     modelEqnContainer().modelEqn(iEqn)->name(),
                     "comm_exchange_borders",
                     "frame_general","restart_yes","coupling_push_min_max",
                     "read_yes", "output_yes"
                 )
        );

        chemistryMem_.push_back
        (
            particle_data_tracker_.addElementProperty< ContainerChemistry<double,N_HISTORY_DEF,N_INTRA_GRIDPOINTS_DEF> >
                 (
                     chemName,
                     "comm_exchange_borders",
                     "frame_general","restart_yes","coupling_none",
                     "read_no", "output_yes"
                 )
        );

        intraPartAv_.push_back
        (
             particle_data_tracker_.addElementProperty< ContainerScalar<double> >
                (
                    avName,
                    "comm_exchange_borders",
                    "frame_general","restart_yes","coupling_push",
                    "read_no", "output_yes"
                )
        );

        intraPartFlux_.push_back
        (
             particle_data_tracker_.addElementProperty< ContainerScalar<double> >
                (
                    fluxName,
                    "comm_exchange_borders",
                    "frame_general","restart_yes","coupling_pull_push",
                    "read_no", "output_yes"
                )
        );

    }

    //Empty data pointers to current state (to be filled later)
    datapointer_.resize(modelEqnContainer().nrEqns());

    // Set containers to correct vector length (only for intra-particle containers)
    for(int iEqn = 0; iEqn < modelEqnContainer().nrEqns(); iEqn++)
    {
        intraPartMem_[iEqn]->setLenVecUsed( modelEqnContainer().modelEqn(iEqn)->nGridPointsUsed()  );
        chemistryMem_[iEqn]->setLenVecUsed( modelEqnContainer().modelEqn(iEqn)->nGridPointsUsed()  );

        if(intraPartMem_[iEqn]->isFull())
            error().throw_error_one(FLERR,
                 "You overfilled the container. Increase the value for N_INTRA_GRIDPOINTS_DEF");

        if(chemistryMem_[iEqn]->isFull())
            error().throw_error_one(FLERR,
                 "You overfilled the container. Increase the value for N_INTRA_GRIDPOINTS_DEF");
    }

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
//	printf("pointer for particleDataID %i and particle ID %i with length %d set\n ",
//	         particleDataID_,
//	         particleID,
//	         intraPartMem_[particleDataID_]->lenVecUsed()
//	      );
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
vector<double> ParticleData::retrieveIntraData(vector<int> _dataIDs, int _particleID, int _gridPoint)
{
    vector<double> output;
    for(int i = 0; i<_dataIDs.size();i++)
    {
        double *** myPtr = (double***)  intraPartMem_[_dataIDs[i]]->begin_slow_dirty();
        double     data = myPtr[_particleID][0][_gridPoint];
        output.push_back(data);
    }
    
    return output;
}

/* ---------------------------------------------------------------------- */
double ParticleData::retrieveIntraData(int _dataID, int _particleID, int _gridPoint)
{
    double *** myPtr = (double***)  intraPartMem_[_dataID]->begin_slow_dirty();
    return myPtr[_particleID][0][_gridPoint];
}

/* ---------------------------------------------------------------------- */
void ParticleData::returnchemistryDataPoint(int particleDataID_, int particleID, int gridPoint, double & data)
{
    current_particleDataID_ = particleDataID_;
	current_particleID_     = particleID;

	if(!runtimeContainersExist_)
        error().throw_error_one(FLERR,"Internal error: requesting intra particle data but containers not existing");

    ptr3 = (double***)  chemistryMem_[particleDataID_]->begin_slow_dirty();   //ptr to vectorial per-particle data
    datapointer_[particleDataID_] = &ptr3[current_particleID_][0][0]; //TODO: put into separate function
    data = ptr3[current_particleID_][0][gridPoint]; //get data
    
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
void ParticleData::returnIntraFlux(int particleID, double &flux)
{
    flux = intraPartFlux_[current_particleDataID_]->get(particleID);
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
void ParticleData::saveChemicalParticleData(int particleDataID_, int particleID, double data, int gridPoint)
{
    if(!runtimeContainersExist_)
        error().throw_error_one(FLERR,"Internal error: trying to save intra particle data but containers not existing");

    //ave intra-particle data
    double *** ptr3 = (double***) chemistryMem_[particleDataID_]->begin_slow_dirty();   //ptr to vectorial per-particle data
    ptr3[particleID][0][gridPoint] += data;      //save data
    
//    printf("PARTICLE DATA: Data (%g) for Particle %i species ID %i and grid Point %i \n",data, particleID, particleDataID_, gridPoint);
    //printf("\n");
    return;
}

/* ---------------------------------------------------------------------- */
void ParticleData::saveIntraParticleAv(int particleDataID_, int particleID, double dataAv)
{
    //save averages
    intraPartAv_[particleDataID_]->set(particleID,dataAv);

    return;
}

void ParticleData::saveIntraParticleFlux(int particleDataID_, int particleID, double dataflux)
{
    //save fluxes
    intraPartFlux_[particleDataID_]->set(particleID,dataflux);

    return;
}

/* ----------------------------------------------------------------------
   Clean chemistry data
------------------------------------------------------------------------- */

void ParticleData::resetChemicalSourceTerms()
{
    for (int m=0;m < modelEqnContainer().nrEqns() ;m++) //particle Data ID
    {
        //printf("cleaning for %i model eqn \n", m); 
        double *** ptr3 = (double***)   chemistryMem_[m]->begin_slow_dirty();   //ptr to vectorial per-particle data
        for (int i=0;i < nbody();i++)  //particle ID
        {
            //printf("cleaning for %i particle \n", i); 
            for(int j=0;j < chemistryMem_[m]->lenVecUsed();j++)
            {
                   ptr3[i][0][j] = 0; //reset to zero
            }
        }
    }

    if(!runtimeContainersExist_)
        error().throw_error_one(FLERR,"Internal error: trying to clean chemical particle data but containers not existing");

    return;
}

