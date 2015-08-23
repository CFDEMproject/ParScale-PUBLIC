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
#define N_INTRA_GRIDPOINTS_DEF  65

#include "stdio.h"
#include "container.h"
#include "custom_value_tracker.h"
#include "pascal_base.h"
#include "pascal_base_interface.h"
#include <vector>
#include "model_base.h"

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

      bool hasPulled() const    
      {return hasPulled_;  }

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
      void setParticleIDPointerPhaseFraction( int particleID);

 	  void returnIntraData(double * data);

      void returnPhaseFractionData(double * phaseFracGas, double * phaseFracLiq);

      void returnPhaseFractionDataGridpoint(double &phaseFracGas, double &phaseFracLiq, int grid_point);
 	  
 	  void returnIntraFlux(int particleID, double &flux);
 	  
  	  void returnIntraTransCoeff(int particleID, double &transCoeff);
 	  
 	  inline int returnId(int _particleID) const {return id_.get(_particleID);};
 	  
 	  void resetIds() const;

      void scanPhaseInformation() const;
 	  
 	  void countBodiesOnMachine() const;

      void resetParticleIDPointer(int particleDataID_);

      void saveIntraParticleData(int particleDataID_, int particleID, double * data);
      void saveIntraParticleAv(int particleDataID_, int particleID, double dataAv);
      void saveIntraParticleFlux(int particleDataID_, int particleID, double dataflux);
      void saveIntraParticleTransCoeff(int particleDataID_, int particleID, double transCoeff);
      void saveChemistryParticleData(int particleDataID_, int particleID, double &data, int gridPoint);
      void saveChemistryJacParticleData(int particleDataID_, int particleID, double &data, int gridPoint);

      void savePhaseFractionData(int particleID, double * dataGas, double * dataLiquid);
      void savePhaseChangeParticleData(   int _phaseID, int particleID, double &data, int gridPoint);
      void savePhaseChangeParticleDataJac(int _phaseID, int particleID, double &data, int gridPoint);
      void savePhaseChangeParticleDataJacLastSlot(int particleID, double &data, int gridPoint);
      void phaseFractionChangeRateEnd(int particleID, double * dataGas, double * dataLiquid, double deltaT);
      void phaseFractionChangeRateStart(int particleID, double * dataGas, double * dataLiquid);

      void returnDataPoint(double & data, int j);
      
      void setDataPoint(double data, int j);

      void returnchemistryDataPoint(int particleDataID_, int particleID, int gridPoint, double & data);
      void returnchemistryJacDataPoint(int particleDataID_, int particleID, int gridPoint, double & data);

      void   retrieveIntraData(vector<int> dataIDs_, int particleID, int gridPoint, vector<double>& output);
      
      double retrieveIntraData(int dataID, int particleID, int gridPoint);

      void   retrievePhaseFractionData(vector<int> _dataIDs, int _particleID, int _gridPoint, vector<double>& output);
      double retrievePhaseFractionData(int _dataID, int _particleID, int _gridPoint);
      void   returnPhaseChangeRateDataPoint(int _phaseId, int particleID, int gridPoint, double & data);
      void   returnPhaseChangeRateJacDataPoint(int _phaseId, int particleID, int _gridPoint, double & data);
      void   returnPhaseChangeRateJacLastSlotDataPoint(int particleID, int _gridPoint, double & data);

      //Access to key containers
      double *** accessIntraPartMem(int _particleDataId) const
      {  return (double***)  intraPartMem_[_particleDataId]->begin_slow_dirty(); }

      double *** accessChemistryMem(int _particleDataId) const
      {  return (double***)  chemistryMem_[_particleDataId]->begin_slow_dirty(); }

      double *** accessChemistryMemJac(int _particleDataId) const
      {  return (double***)  chemistryMemJacobi_[_particleDataId]->begin_slow_dirty(); }

      double *** accessPhaseChangeRate(int _phaseId) const
      {  return (double***)  intraPhaseChangeRateMem_[phaseIDMap[_phaseId]]->begin_slow_dirty(); }

      double *** accessPhaseChangeRateJac(int _phaseId) const
      {  return (double***)  intraPhaseChangeRateMemJacobi_[phaseIDMap[_phaseId]]->begin_slow_dirty(); }

      double *** accessPhaseChangeRateVolumetric(int _phaseId) const
      {  return (double***)  intraPhaseChangeRateVolumetric_[phaseIDMap[_phaseId]]->begin_slow_dirty(); }

      double *** accessIntraConvectiveFlux(int _phaseId) const
      {  return (double***)  intraConvectiveFluxMem_[phaseIDMap[_phaseId]]->begin_slow_dirty(); }

      double *** accessPhaseFractionMem(int _phaseId) const
      {  return (double***)  intraPhaseFractionMem_[phaseIDMap[_phaseId]]->begin_slow_dirty(); }

      void   resetPhaseChangeSourceTerms();
      void   resetChemicalSourceTerms();

      void   computeConvection();

      void   normalizeInternalFields(int _particleID, double _factor);

      //inline functions to access scalar and vector particle properties
      inline double& pRadius(int i)
      { return radius_(i); }

      inline double& pRadiusConst(int i) const
      { return radius_(i); }

      double *** ptr3;
      double *** ptr3Jac;
	  int current_particleID_ ;
	  int current_particleDataID_;

      double *** ptr3GasPhase_;
      double *** ptr3LiquidPhase_;

      mutable bool haveGasPhase;
      mutable bool haveLiquidPhase;
      mutable int  numberOfPhases;
      mutable int  eqnIdFirstLiquid;          //id of first liquid-phase equation (=solvent), useful for many liquid flux calculations
      mutable vector<int>  phaseList;         //list of phases
      mutable vector<int>  phaseIDMap;        //list of ids for each phase, 0...SOLID, 1...GAS, 2 ... LIQUID, 3 ... NONE, 
                                              //this is the inverse map of phaseList
      mutable double referencePressure;       //a reference pressure. useful for many flux calculations. to be set by 
    
      mutable bool haveConvectiveFluxGasPhase;
      mutable bool haveConvectiveFluxLiquidPhase;

      mutable bool writeDebugGasPhase;
      mutable bool writeDebugLiquidPhase;

      int sizePhaseChangeRateMemJacobi() const {return intraPhaseChangeRateMemJacobi_.size();};

    private:
      // add containers which depend on settings which are read during run-time
      void addRunTimeContainers() const;

      void addParticles(int particleCount, int particleCountGlobal) const;

      void generatePhaseIDMap() const;

      // ************ DATA STORED FOR EACH PARTICLE *****************
      // class holding all fields
      // basic properties (see below) register here
      // models may also register additional properties
      CustomValueTracker &particle_data_tracker_;

      // A - Basic data for each particle, one for each particle! Scalar or vector quantity
      ContainerScalar<int>        &id_;        // pull
      ContainerScalar<double>     &radius_;    // pull - push

      mutable vector<double*> datapointer_; //pointers to current particle data, for each particleDataID_ !
      mutable vector<double*> datapointerJac_; //pointers to current particle data, for each particleDataID_ !
      mutable vector<double*> datapointerPhaseFraction_; //pointers to current particle data, for each phase!

      // B - Intra-Particle Data for each particle, one array for each particle and model_eqn!
      mutable vector< ContainerCvode<double,N_HISTORY_DEF,N_INTRA_GRIDPOINTS_DEF>* > intraPartMem_;

      mutable vector< ContainerChemistry<double,1,N_INTRA_GRIDPOINTS_DEF>* > chemistryMem_;
      mutable vector< ContainerChemistry<double,1,N_INTRA_GRIDPOINTS_DEF>* > chemistryMemJacobi_;

      // C - Important intra-particle data, always scalar
      mutable vector< ContainerScalar<double>* > intraPartAv_;    //Volume-average particle information
      mutable vector< ContainerScalar<double>* > intraPartFlux_;  //Fluxes at the boundary of the particle, only EXPLICIT part on pull, on EXPLICIT+IMPLICIT on push
      mutable vector< ContainerScalar<double>* > intraPartTransCoeff_;  //Transfer coefficient at the boundary of the particle, for IMPLICIT flux calculation, only pull
    
      mutable vector< ContainerCvode<double,1,N_INTRA_GRIDPOINTS_DEF>* > intraPhaseFractionMem_;
      mutable vector< ContainerCvode<double,1,N_INTRA_GRIDPOINTS_DEF>* > intraPhaseChangeRateMem_;
      mutable vector< ContainerCvode<double,1,N_INTRA_GRIDPOINTS_DEF>* > intraPhaseChangeRateMemJacobi_;
      mutable vector< ContainerCvode<double,1,N_INTRA_GRIDPOINTS_DEF>* > intraPhaseChangeRateVolumetric_;

      //Convective fluxes
      mutable vector< ContainerCvode<double,1,N_INTRA_GRIDPOINTS_DEF>* > intraConvectiveFluxMem_;
      
      // ************ DATA STORED FOR EACH PARTICLE *****************
      mutable double * tmpSourceVariable_;
      mutable int* nBodyPerProcess_;            //number of particles per process AT THE TIME OF ALLOC!
      mutable int  nBodyInProcessesBelowMe_;    //number of particles per process AT THE TIME OF ALLOC!
      mutable bool runtimeContainersExist_;

      bool verbose_;
      bool hasPulled_;
      
      // global-local lookup
      //int mapTagMax_;
      //int *mapArray_;
};

} //end PASCAL_NS

#endif
