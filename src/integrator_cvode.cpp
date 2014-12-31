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



#include "integrator_cvode.h"
#include <stdio.h>
#include "control.h"
#include "output.h"
#include <sys/stat.h>

namespace PASCAL_NS
{

class modelData
{
public:

   modelData(ModelEqn* m_eqn)
   {
      m_func = m_eqn;
   }
    virtual ~modelData() {}
//    vector_fp m_pars; //TODO: add parameters of the model if necessary
      ModelEqn* m_func;
};
}

extern "C" 
{
    /*
     *  Functions called by cvodes to evaluate ydot given y.  The cvode
     *  integrator allows passing in a void* pointer to access
     *  external data. This pointer can be cast to a pointer ("d") to an instance
     *  of class "modelEqn". The equations to be integrated should be
     *  specified by deriving a class from class "modelEqn" that evaluates 
     *  the desired equations.
     *  @Radl, IPPT
  */
    int cvodes_rhs(realtype t, N_Vector u, N_Vector udot, void* f_model)
    {
        realtype *udata, *dudata;
        udata = NV_DATA_S(u);
        dudata = NV_DATA_S(udot);
        
        //get access to ModelEqn object that called the integrator
        PASCAL_NS::modelData* m_data = (PASCAL_NS::modelData*)f_model;
        PASCAL_NS::ModelEqn*  m_eqn = m_data->m_func;
   	    
        m_eqn->eval(t, udata, dudata, NULL);
      
	    return 0;
    }

    int cvodes_jac(long int N, long int mu, long int ml,
                   realtype t, N_Vector u, N_Vector fu, 
                   DlsMat J, void *f_model,
                   N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
    {
        realtype *udata, *tmp1data, *tmp2data, *tmp3data, *fudata;
        udata = NV_DATA_S(u);
        fudata = NV_DATA_S(fu);
        tmp1data = NV_DATA_S(tmp1);
        tmp2data = NV_DATA_S(tmp2);
        tmp3data = NV_DATA_S(tmp3);
       
        //get access to ModelEqn object that called the integrator
        PASCAL_NS::modelData* m_data = (PASCAL_NS::modelData*)f_model;
        PASCAL_NS::ModelEqn*  m_eqn = m_data->m_func;
        
        m_eqn->returnJac(N, mu, ml,
                   t, udata, fudata, 
                   J, NULL,
                   tmp1data, tmp2data, tmp3data);
        
  	return 0;
    }

} //end extern "C"


using namespace PASCAL_NS;

/* ----------------------------------------------------------------------
   ModelBase Constructor
------------------------------------------------------------------------- */

IntegratorCvode::IntegratorCvode(ParScale *ptr, int nGrid) : 
    Integrator(ptr),
    cvode_mem(NULL),			
    u(NULL),		
    MX(nGrid),
    NEQ(nGrid),
    T0(control().simulationState().time()),
    t(0.0),
    T1(0.0),
    abstol(1.0e-5),  
    reltol(0.0),
    m_data_(0),
	flag(0)
{
  //check the particle mesh
   if(MX<1)
      output().write_screen_one("WARNING: Your integrator has no grid points! \n"); 
}

/////////////////////////////////////////////////////////////////////////////
                             // Destructor 
/////////////////////////////////////////////////////////////////////////////

IntegratorCvode::~IntegratorCvode()
{
    fclose(fp);
    delete m_data_;
    N_VDestroy_Serial(u);   				// Free the u vector 
    CVodeFree(&cvode_mem);  				// Free the integrator memory 
    delete  tempIntraData_;
}

/////////////////////////////////////////////////////////////////////////////
                          // MEMBER functions 
/////////////////////////////////////////////////////////////////////////////
void IntegratorCvode::init(double T0, ModelEqn& m_eqn)
{

   if (u)
   {
   N_VDestroy_Serial(u); 			// Free the u vector 
   }  

   if (cvode_mem)
   {
   CVodeFree(&cvode_mem);  			// Free the integrator memory 
   } 	
			
   delete m_data_;				    // pass a pointer to func in m_data_
   m_data_ = new modelData(&m_eqn); 		// m_eqn.nparams()); //TODO: could hand over data if necessary

   						//TODO: hand over data from the ModelEqn to set solution vectors
 
   u = N_VNew_Serial(NEQ+1);  			// Create and allocate a serial vector, (NEQ+1) needed because of behavior of CVODE on boundary points!!
     
   if(check_flag((void*)u, "N_VNew_Serial", 0)) return; //TODO: throw error
   //if(check_flag((void *)m_eqn, "malloc", 2)) return; //TODO: throw error
      
   //Call CVodeCreate to create the solver memory and specify the 
   //Backward Differentiation Formula and the use of a Newton iteration 
   cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
   if(check_flag((void *)cvode_mem, "CVodeCreate", 0)) return;
           
   // Call CVodeInit to initialize the integrator memory and specify the
   // user's right hand side function in u'=f(t,u), the inital time T0, and
   // the initial dependent variable vector u. 
   flag = CVodeInit(cvode_mem, cvodes_rhs, T0, u);
   if(check_flag(&flag, "CVodeInit", 1)) return;

   // Call CVodeSStolerances to specify the scalar relative tolerance
   // and scalar absolute tolerance 
   flag = CVodeSStolerances(cvode_mem, reltol, abstol);
   if (check_flag(&flag, "CVodeSStolerances", 1)) return; 
   
   // Set the pointer to user-defined data 
   flag = CVodeSetUserData(cvode_mem, m_data_);
   if(check_flag(&flag, "CVodeSetUserData", 1)) return; 
  
   // Call CVBand to specify the CVBAND band linear solver 
   flag = CVBand(cvode_mem, NEQ, 1, 1);
   if(check_flag(&flag, "CVBand", 1)) return; 	

   // Set the user-supplied Jacobian routine Jac 
   flag = CVDlsSetBandJacFn(cvode_mem, cvodes_jac);
   if(check_flag(&flag, "CVDlsSetBandJacFn", 1)) return; 

   printf("IntegratorCvode successfully initialized with %d grid points \n", MX );
   
   return;
}

///////////////////////////////////////////////////////////////////////////

void IntegratorCvode::integrate_begin(const char* stateType, int nGridPointsUsed_, int particleDataID_)
{
   tempIntraData_   = create<double>(tempIntraData_, nGridPointsUsed_); 

   double r;
   int h;
   
   for(particleID=0; particleID<particleData().nbody(); particleID++)
   {
       //printf("number of particles in integrator = %i, actual particle ID \n",particleData().nbody(),);
       particleData().setParticleIDPointer(particleDataID_,particleID);	
       particleData().returnIntraData(tempIntraData_);

      // printf(" requesting data with id %d for particle %d on this CPU with %d particles \n", particleDataID_, particleID, particleData().nbody());

       if(strcmp(stateType,"udata") == 0)
       {
          //set the current udata
          realtype *udata;
          udata = NV_DATA_S(u);

          for (int h=1; h <= MX; h++) 
          {
            udata[h] = tempIntraData_[h-1];
            //printf("Udata[%i] in integrator for particle Data ID %i before integration = %g at t = %g sec\n", h,particleDataID_,udata[h],control().simulationState().time());
          }

       //Reinitialize with CVodeReInit
       Reinitialize(control().simulationState().time(), m_data_);
       }

       else if(strcmp(stateType,"radius") == 0)
       {
           error().throw_error_one(FLERR,"IntegratorCvode does not manage the radius in this one.\n");
       }
		  
	   umax = N_VMaxNorm(u);
	   T0 = control().simulationState().time();
	   T1 = T0+control().simulationState().deltaT();

	   realtype *udata;
	   udata = NV_DATA_S(u);
		
	   flag = CVode(cvode_mem, T1, u, &t, CV_NORMAL); 				
	   if(check_flag(&flag, "CVode", 1)) return;


		double radialDist(0.0), deltaR(0.0);
		deltaR = 1.0/(MX-1);
		  
	   flag = CVodeGetNumSteps(cvode_mem, &nst);
	   check_flag(&flag, "CVodeGetNumSteps", 1);
      
        for (int h=1; h <= MX; h++) 
        {
            //printf(" udata[%i] in integrator after integration = %g \n", h,udata[h]);
            tempIntraData_[h-1]= udata[h];
        }

        particleData().saveIntraParticleData(particleDataID_, particleID, tempIntraData_);
   }
}

///////////////////////////////////////////////////////////////////////////

void IntegratorCvode::Reinitialize(double T0, modelData* m_data_)
{
    int flag;
    flag = CVodeReInit(cvode_mem, T0, u);
}


/////////////////////////////////////////////////////////////////////////////////////

void IntegratorCvode::PrintHeader(realtype reltol, realtype abstol, realtype umax)
{
   printf("Mesh dimensions = %d\n", MX);
   printf("Total system size = %d\n", NEQ);
 #if defined(SUNDIALS_EXTENDED_PRECISION)
   printf("Tolerance parameters: reltol = %Lg   abstol = %Lg\n\n",
         reltol, abstol);
   printf("At t = %Lg      max.norm(u) =%14.6Le \n", T0, umax);
 #elif defined(SUNDIALS_DOUBLE_PRECISION)
   printf("Tolerance parameters: reltol = %lg   abstol = %lg\n\n",
         reltol, abstol);
   printf("At t = %lg      max.norm(u) =%14.6le \n", T0, umax);
 #else
   printf("Tolerance parameters: reltol = %g   abstol = %g\n\n", reltol, abstol);
   printf("At t = %g      max.norm(u) =%14.6e \n", T0, umax);
 #endif

  return;
}

///////////////////////////////////////////////////////////////////////////
	
void IntegratorCvode::PrintOutput(realtype t, realtype umax, long int nst)
{
 #if defined(SUNDIALS_EXTENDED_PRECISION)
   printf("At t = %4.4Lf   max.norm(u) =%14.6Le   nst = %4ld\n", t, umax, nst);
 #elif defined(SUNDIALS_DOUBLE_PRECISION)
   printf("At t = %4.4f   max.norm(u) =%14.6le   nst = %4ld\n", t, umax, nst);
 #else
   printf("At t = %4.4f   max.norm(u) =%14.6e   nst = %4ld\n", t, umax, nst);
 #endif
   return;
}

///////////////////////////////////////////////////////////////////////////

void IntegratorCvode::PrintFinalStats(void *cvode_mem)
{
   int flag;
   long int nst, nfe, nsetups, netf, nni, ncfn, nje, nfeLS;

   flag = CVodeGetNumSteps(cvode_mem, &nst);
   check_flag(&flag, "CVodeGetNumSteps", 1);
   flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
   check_flag(&flag, "CVodeGetNumRhsEvals", 1);
   flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
   check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
   flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
   check_flag(&flag, "CVodeGetNumErrTestFails", 1);
   flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
   check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
   flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
   check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

   flag = CVDlsGetNumJacEvals(cvode_mem, &nje);
   check_flag(&flag, "CVDlsGetNumJacEvals", 1);
   flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);
   check_flag(&flag, "CVDlsGetNumRhsEvals", 1);

   printf("\nFinal Statistics:\n");
   printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
	 nst, nfe, nsetups, nfeLS, nje);
   printf("nni = %-6ld ncfn = %-6ld netf = %ld\n \n",
	 nni, ncfn, netf);

   return;
}

///////////////////////////////////////////////////////////////////////////

int IntegratorCvode::check_flag(void *flagvalue, char *funcname, int opt)
{
   int *errflag;

   // Check if SUNDIALS function returned NULL pointer - no memory allocated 
   if (opt == 0 && flagvalue == NULL) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
     return(1); }

  // Check if flag < 0 
   else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
              funcname, *errflag);
      return(1); }}

   // Check if function returned NULL pointer - no memory allocated 
   else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

   return(0);
}


