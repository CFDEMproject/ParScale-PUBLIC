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
/*-----------------------------------------------------------------------------------
Description
	This class holds all scalar (e.g., local temperature) and vectorial (e.g.,
	particle velocity) information.
-----------------------------------------------------------------------------------*/

#ifndef PASC_PARTICLE_DATA_H
#define PASC_PARTICLE_DATA_H

#define N_HISTORY_DEF           1
#define N_INTRA_GRIDPOINTS_DEF  30

#include "stdio.h"
#include "container.h"
#include "custom_value_tracker.h"
#include "pascal_base.h"
#include "pascal_base_interface.h"
#include <vector>


//#include "particle_mesh.h"

using std::vector;

namespace PASCAL_NS
{

class ParticleData : public ParScaleBase, public ParScaleBaseInterface
{
    public:

      ParticleData(ParScale *ptr);
      ~ParticleData();

      void allocate();

      void read();

      void write();

      void scatter();

      void parallelize();

      void init();

      void push();
      
      void pull();

      void parse_command(int narg,char const* const* arg);

      bool runtimeContainersExist() const
      {return runtimeContainersExist_;} ;

      CustomValueTracker& data() const
      { return particle_data_tracker_; }

      int nbody() const
      { return data().nbody(); }

      int nbody_all() const
      { return data().nbody_all(); }

      void grow_nbody(int _nbody,int _nbody_all);

      void requestIntraParticleData(int particleDataID_, int particleID, double * data);

	  void setParticleIDPointer(int particleDataID_, int particleID);

 	  void returnIntraData(double * data);
 	  
 	  void returnIntraFlux(int particleID, double &flux);
 	  
 	  inline int returnId(int _particleID) const {return id_.get(_particleID);};
 	  
 	  void resetIds() const;
 	  
 	  void countBodiesOnMachine() const;

      void resetParticleIDPointer(int particleDataID_);

      void saveIntraParticleData(int particleDataID_, int particleID, double * data);

      void saveIntraParticleAv(int particleDataID_, int particleID, double dataAv);

      void saveIntraParticleFlux(int particleDataID_, int particleID, double dataflux);
    
      void saveChemicalParticleData(int particleDataID_, int particleID, double data, int gridPoint);

      void returnDataPoint(double & data, int j);
      
      void setDataPoint(double data, int j);

      void returnchemistryDataPoint(int particleDataID_, int particleID, int gridPoint, double & data);

      vector<double> retrieveIntraData(vector<int> dataIDs_, int particleID, int gridPoint);
      
      double retrieveIntraData(int dataID, int particleID, int gridPoint);

      void resetChemicalSourceTerms();

      //inline functions to access scalar and vector particle properties
      inline double& pRadius(int i)
      { return radius_(i); }

      double *** ptr3;
	  int current_particleID_ ;
	  int current_particleDataID_;

    private:
      // add containers which depend on settings which are read during run-time
      void addRunTimeContainers() const;

      void addParticles(int particleCount, int particleCountGlobal) const;

      // ************ DATA STORED FOR EACH PARTICLE *****************
      // class holding all fields
      // basic properties (see below) register here
      // models may also register additional properties
      CustomValueTracker &particle_data_tracker_;

      // A - Basic data for each particle, one for each particle! Scalar or vector quantity
      ContainerScalar<int>        &id_;        // pull
      ContainerScalar<double>     &radius_;    // pull - push

      mutable vector<double*> datapointer_; //pointers to current particle data, for each particleDataID_ !

      // B - Intra-Particle Data for each particle, one array for each particle and model_eqn!
      mutable vector< ContainerCvode<double,N_HISTORY_DEF,N_INTRA_GRIDPOINTS_DEF>* > intraPartMem_;

      mutable vector< ContainerChemistry<double,N_HISTORY_DEF,N_INTRA_GRIDPOINTS_DEF>* > chemistryMem_;

      // C - Important intra-particle data, always scalar
      mutable vector< ContainerScalar<double>* > intraPartAv_;    //Volume-average particle information
      mutable vector< ContainerScalar<double>* > intraPartFlux_;  //Fluxes at the boundary of the particle


      // ************ DATA STORED FOR EACH PARTICLE *****************
      mutable int* nBodyPerProcess_;            //number of particles per process AT THE TIME OF ALLOC!
      mutable int  nBodyInProcessesBelowMe_;    //number of particles per process AT THE TIME OF ALLOC!
      mutable bool runtimeContainersExist_;

      bool verbose_;
      // global-local lookup
      //int mapTagMax_;
      //int *mapArray_;
};

} //end PASCAL_NS

#endif
