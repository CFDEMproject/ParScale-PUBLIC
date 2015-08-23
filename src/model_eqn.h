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
	This class is the base class for all physical sub-models (e.g., for properties),
	as well as overall particle model equations (e.g., for transient heat conduction
	inside the particle).
-----------------------------------------------------------------------------------*/

#ifndef PASC_MODEL_EQN_H
#define PASC_MODEL_EQN_H
#define VLARGENUMBER 1e64
#define LARGENUMBER 1e32
#define SMALLNUMBER 1e-32
#define VSMALLNUMBER 1e-64

#include "error.h"
#include "memory_ns.h"
#include "model_base.h"
#include "particle_data.h"
#include "integrator.h"
#include "particle_mesh.h"
#include <sundials/sundials_direct.h>
#include "coupling.h"
#include "model_chemistry_container.h"
#include "model_eqn_container.h"
//#include "model_base.h"


using namespace PASCAL_MEMORY_NS;

namespace PASCAL_NS
{

#define INVERSE_UNIVERSAL_GAS_CONSTANT 0.0001202723625

//numbering of boundary condition types
enum{ NEUMANN,		//0
      DIRICHLET,	//1
      CONVECTIVE	//2
    };

//numbering of modeling approaches, i.e., variable modelingApproach
enum{CONTINUUM, SHRINKINGCORE};

enum{SOLID, GAS, LIQUID, NONE};

class ModelEqn : public ModelBase
{
    friend class ModelEqn1DSpherical;
    friend class ModelEqn1DCartesian;
    friend class ModelEqnShrinkingCore;

    public:

      ModelEqn(ParScale *ptr, char *name);
      virtual ~ModelEqn();

      virtual void init(int narg, char const* const* arg, int eqnType, int modelEqnID);

      virtual void setupParticle();
      virtual void computeParticleProps();
      virtual void computeParticleAverages() = 0;
      virtual void computeSurfaceFluxes() = 0;

      virtual void begin_of_step() {};
      
      virtual void pre_middle_of_step() {};
      
      virtual void post_middle_of_step() {};
      
      virtual void end_of_step();

      /*
         * Evaluate the right-hand-side function. Called by the integrator.
         * @param[in]  t time.
         * @param[in]  udata solution vector, length neq()
         * @param[out] dudata rate of change of solution vector, length neq()
         * @param[in]  p spare parameter vector
      */
      virtual void eval(double t, double* udata, double* dudata, double* p){};
      virtual void returnJac(long int N, long int mu, long int ml,
                   double t, double* udata, double* fudata,
                   DlsMat J, double* p,
                   double* tmp1data, double* tmp2data, double* tmp3data){};

      virtual void  evaluatePhaseFlux()       const {};

      //Access functions
      Integrator&   integrator()    { return *integrator_; }

      void          setParticleDataID(int ID) const {   particleDataID_ = ID; haveParticleDataID=true; };
      int           particleDataID()          const {   return particleDataID_; };
      
      int           nGridPointsUsed()         const {   return nGridPointsUsed_; };
      int           modelingApproach;

      int           return_myPhase()          const {    return inPhase_;         };
      int           return_eqn_type()         const {    return eqnType_;         };

      double        phaseFractionMinimum;   //smallest phase fraction to solve equation

      double        phaseFractionMinimumConvection;   //smallest phase fraction to add convection term

      double        convectionBoundMinStefanFactor;   //smallest factor to account for Stefan effect

      bool          solveMe;                //main switch to solve equation

      bool          updatePhaseFraction;    //switch to update the phase fraction of the phase the species is in, not the concentration
                                            //this is useful for solvents

      bool          averagePhaseFraction;   //average will be that of phase fraction (not that of species)

      bool          solveConvectiveFlux;    //switch to activate convective flux calculation
      
      bool          writeDebugContainers;   //switch to activate dumping of secondary containers
      
      bool          normalizeDuringGrowth;  //switch to activate normalization of concentration during growth

     private:
      Integrator     * integrator_;         //ptr - the (time) integrator that will be used to solve this equation

      //Initial & Boundary conditions
      double    ICvalue;            //IC value
      int       BC[2];              //BC type
      double    BCvalue[2];         //BC value
      double    rMAX;                              //radius of the current particle
      double    environmentU;         			   //enviroment temperature/concentration 
      double    environmentFlux;                   //flux to environment, EXPLICIT or EXPLICIT+IMPLICIT 
                                                   //flux is SURFACE-AREA SPECIFIC !!
      double    environmentTransCoeff;             //transfer coefficient to environment, only for fixed calculation
                                                      //(i.e., heat/species flux across particle interface)   
      int       heatEqnID_;                         //ID refering to the heat equation
      bool      IsoThermal_;                        //Determine whether species is isothermal or not
      double    isoTemp_;                           //Temperature if species system is isothermal

      //Temporary array for intra-particle data 
      double *  tempIntraData_;
      double *  tempPhaseDataGas_;
      double *  tempPhaseDataLiquid_;

      double    tempAvData_;
      double    tempPartFlux_;
      double    diffu_eff_;

      double    surface_area;            //surface area of particle [m²]
      double    segment_vol;             //volume of segment, depending on number of grid points and radius of particle [m³]

      int       particleID;

      bool      verbose_;

      //Pointers to key physical models
      ModelBase*    diffusivity;
      ModelBase*    phaseFraction;
      ModelBase*    capacity_solid;
      ModelBase*    capacity_gas;
      ModelBase*    capacity_liquid;
      ModelBase*	thermal_solid_conductivity; 
      ModelBase*	thermal_gas_conductivity;
      ModelBase*	thermal_liquid_conductivity;
      ModelBase* 	transfer_coeff;
      ModelBase*	density_solid;
      ModelBase* 	density_gas;
      ModelBase* 	density_liquid;
	  ModelBase*	global_properties;
      ModelBase*	tortuosity;
      ModelBase*    pore_radius;
      ModelBase*    molar_mass;
      ModelBase*    permeability;
      ModelBase*    viscosity;
      ModelBase*    surface_tension;
      ModelBase*    film_flow;

      mutable int   inPhase_;          //Species type - enum - 0-Solid, 1-Gas, 2-Fluid, 3-None
      bool          eqnSolveThis_;     //switch to activate or deactive solution of this equation
      int           eqnID_;            //id of the eqn
      int           eqnType_;          //integer to identify if heat, species or other transport eqn is solved
      mutable int   particleDataID_;   //id of the object in the particleData class that saves the information
      mutable bool  haveParticleDataID;
      mutable int   nGridPointsUsed_;   //number of grid points used to define state of the model
      
};

} //end PASCAL_NS

#endif
