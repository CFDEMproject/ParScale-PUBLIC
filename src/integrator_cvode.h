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


#ifndef PASC_INTEGRATOR_CVODE_H
#define PASC_INTEGRATOR_CVODE_H

#include "integrator.h"
#include "model_eqn.h"
#include <cvode/cvode.h>             		/* prototypes for CVODE fcts., consts. */
#include <cvode/cvode_band.h>        		/* prototype for CVBand */
#include <nvector/nvector_serial.h>  		/* serial N_Vector types, fcts., macros */
#include <sundials/sundials_band.h>  		/* definitions of type DlsMat and macros */
#include <sundials/sundials_types.h> 		/* definition of type realtype */
#include <sundials/sundials_math.h>  		/* definition of ABS and EXP */
#include "particle_data.h"


namespace PASCAL_NS
{

class modelData;


class IntegratorCvode : public Integrator
{

    
    public:

    IntegratorCvode(ParScale *ptr, int nGrid);
	
	~IntegratorCvode();
     
    //int returnParticleID() {return particleID;}; 

    private:

 	FILE *fp;			    //pointer for writing data
	N_Vector u;             //solution-vector 
    double dx; 		        //data to be saved in an instance of particleMesh 
    void *cvode_mem;
    realtype reltol;        //relative tolerance for iteration
  	realtype abstol;		//absolute tolerance for iteration
  	realtype t;			    //actual time, doesnt matter
	realtype T0;			//initial time, 		
	realtype T1;			//first output time 		
	realtype umax;			//max norm of u (N_VMaxNorm(u))	
    long int nst;			//number of integration steps 
	int flag;			    //error flag
	int MX;				    //mesh dimensions (number intervals) 
	int NEQ;			    //number of equations (NEQ=MX)
	int h;
	void *flagvalue;
	char *funcname;
	int opt;

    modelData* m_data_;
     
  	virtual void init(double t0, ModelEqn& m_eqn);

   	virtual void integrate_begin(const char* stateType, int nGridPointsUsed_, int particleDataID_);
  	virtual void integrate_middle() {};
   	virtual void integrate_end() {};	

	/* Private Helper*/ 
	void PrintHeader(realtype reltol, realtype abstol, realtype umax);
	void PrintOutput(realtype t, realtype umax, long int nst);
	void PrintFinalStats(void *cvode_mem);
    void N_VPrint_Serial(N_Vector u);
    void Reinitialize(double T0, modelData* m_data_);
	/*Private function to check function return values */
	int check_flag(void *flagvalue, char *funcname, int opt);

	protected:

};

} //end PASCAL_NS

#endif
