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
    fp(NULL),
    u(NULL),		
    MX(nGrid),
    NEQ(nGrid),
    T0(control().simulationState().time()),
    t(0.0),
    T1(0.0),
    mxsteps_(500),
    maxord_(5),
    deltaTInit_(1e-16),
    deltaTMin_(1e-16),
    deltaTMax_(1e3),
    maxNonLinearIterations_(3),
    reltol(0.0e-5),
    abstol(1.0e-5),  
    m_data_(0),
	flag(0),
    umax(0.0)
{
  //check the particle mesh
   if(MX<1)
      output().write_screen_one("WARNING: Your integrator has no grid points! \n");

  linearSolver_ = CVDIAG;

}

/////////////////////////////////////////////////////////////////////////////
                             // Destructor 
/////////////////////////////////////////////////////////////////////////////

IntegratorCvode::~IntegratorCvode()
{
    if(fp)
        fclose(fp);
    N_VDestroy_Serial(u);   				// Free the u vector 
    CVodeFree(&cvode_mem);  				// Free the integrator memory 

    delete  m_data_;
    destroy<double>(tempIntraData_);
    destroy<double>(tempPhaseDataGas_);
    destroy<double>(tempPhaseDataLiquid_);
    destroy<double>(tempPhaseDataSolid_);
}

/////////////////////////////////////////////////////////////////////////////
                          // MEMBER functions 
/////////////////////////////////////////////////////////////////////////////
void IntegratorCvode::init(double T0, ModelEqn& m_eqn)
{
   //Initialize Integrator
   Integrator::init(T0, m_eqn);

   //Set main integrator
   if(parameters_["absTol"].isNull())
    
   //Set all relevant integration parameters
   if(integrator_["type"].isNull())
     error().throw_error_one(FLERR,"'type' was not specified in settings/integrator.json/integrator.\n");
   if(strcmp(qPrintable(integrator_["type"].toString()),"CVODE") != 0 )
     error().throw_error_one(FLERR,"'type' incorrectly specified (i.e., not 'CVODE') in settings/integrator.json/integrator.\n");

   if(integrator_["linearSolver"].isNull())
     error().throw_error_one(FLERR,"'linearSolver' was not specified in settings/integrator.json/integrator.\n");

   if(strcmp(qPrintable(integrator_["linearSolver"].toString()),"CVDiag") == 0 )
        linearSolver_ = CVDIAG;
   else if(strcmp(qPrintable(integrator_["linearSolver"].toString()),"CVDense") == 0 )
        linearSolver_ = CVDENSE;
   else if(strcmp(qPrintable(integrator_["linearSolver"].toString()),"CVBand") == 0 )
        linearSolver_ = CVBAND;
   else if(strcmp(qPrintable(integrator_["linearSolver"].toString()),"CVSpgmr") == 0 )
        linearSolver_ = CVSPGMR;
   else if(strcmp(qPrintable(integrator_["linearSolver"].toString()),"CVSpbcg") == 0 )
        linearSolver_ = CVSPBCG;
   else if(strcmp(qPrintable(integrator_["linearSolver"].toString()),"CVSptfqmr") == 0 )
        linearSolver_ = CVSPTFQMR;
   else 
        error().throw_error_one(FLERR,"'linearSolver' specified in settings/integrator.json/integrator/linearSolver was not recognized.\n");

   //settings for this integrator
   if(parameters_["relTol"].isNull())
     printf("WARNING: IntegratorCvode will use default relTol, since not found in settings/integrator.json/parameters. \n");
   else
       reltol=(parameters_["relTol"].toDouble());  

   if(parameters_["mxsteps"].isNull())
     printf("WARNING: IntegratorCvode will use default mxsteps, since not found in settings/integrator.json/parameters. \n");
   else
       mxsteps_=(int)(parameters_["mxsteps"].toDouble());    

   if(parameters_["maxord"].isNull())
     printf("WARNING: IntegratorCvode will use default maxord, since not found in settings/integrator.json/parameters. \n");
   else
    maxord_=(int)(parameters_["maxord"].toDouble());   

   if(parameters_["deltaTInit"].isNull())
     printf("WARNING: IntegratorCvode will use default deltaTInit, since not found in settings/integrator.json/parameters. \n");
   else
    deltaTInit_=(parameters_["deltaTInit"].toDouble());  

   if(parameters_["deltaTMin"].isNull())
     printf("WARNING: IntegratorCvode will use default deltaTMin, since not found in settings/integrator.json/parameters. \n");
   else
   deltaTMin_=(parameters_["deltaTMin"].toDouble());  

   if(parameters_["deltaTMax"].isNull())
     printf("WARNING: IntegratorCvode will use default deltaTMax, since not found in settings/integrator.json/parameters. \n");
   else
    deltaTMax_=(parameters_["deltaTMax"].toDouble());  

   if(parameters_["maxNonLinearIterations"].isNull())
     printf("WARNING: IntegratorCvode will use default maxNonLinearIterations, since not found in settings/integrator.json/parameters. \n");
   else
    maxNonLinearIterations_=(int)(parameters_["maxNonLinearIterations"].toDouble()); 


   //Allocate mem
   tempIntraData_      = create<double>(tempIntraData_,       NEQ+1); 
   tempPhaseDataGas_   = create<double>(tempPhaseDataGas_,    NEQ+1); 
   tempPhaseDataLiquid_= create<double>(tempPhaseDataLiquid_, NEQ+1); 
   tempPhaseDataSolid_ = create<double>(tempPhaseDataSolid_,  NEQ+1); 

   if (u)
       N_VDestroy_Serial(u); 			// Free the u vector 

   if (cvode_mem)
       CVodeFree(&cvode_mem);  			// Free the integrator memory 


   delete m_data_;				            // pass a pointer to func in m_data_
   m_data_ = new modelData(&m_eqn); 		// m_eqn.nparams()); //TODO: could hand over data if necessary

   						//TODO: hand over data from the ModelEqn to set solution vectors
  
   u = N_VNew_Serial(NEQ+1);  			// Create and allocate a serial vector, (NEQ+1) needed because of behavior of CVODE on boundary points!!
     
   if(check_flag((void*)u, "N_VNew_Serial", 0)) return; //TODO: throw error
      
   //Init u
   for(uint iVecId=0; iVecId<(NEQ+1); iVecId++)
       NV_Ith_S(u,iVecId) = 0.0;

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
   // and scalar absolute tolerance, and max steps 
   flag = CVodeSStolerances(cvode_mem, reltol, abstol);
   if (check_flag(&flag, "CVodeSStolerances", 1)) return;

   flag = CVodeSetMaxNumSteps(cvode_mem, mxsteps_);
   if (check_flag(&flag, "CVodeSetMaxNumSteps", 1)) return;

   flag = CVodeSetMaxOrd(cvode_mem, maxord_);
   if (check_flag(&flag, "CVodeSetMaxOrd", 1)) return;
 
   flag = CVodeSetInitStep(cvode_mem, deltaTInit_);
   if (check_flag(&flag, "CVodeSetInitStep", 1)) return;

   flag = CVodeSetMinStep(cvode_mem, deltaTMin_);
   if (check_flag(&flag, "CVodeSetMinStep", 1)) return;

   flag = CVodeSetMaxStep(cvode_mem, deltaTMax_);
   if (check_flag(&flag, "CVodeSetMaxStep", 1)) return;

   flag = CVodeSetMaxNonlinIters(cvode_mem, maxNonLinearIterations_);
   if (check_flag(&flag, " CVodeSetMaxNonlinIters(", 1)) return;

   // Set the pointer to user-defined data 
   flag = CVodeSetUserData(cvode_mem, m_data_);
   if(check_flag(&flag, "CVodeSetUserData", 1)) return; 

   //Select the linear solver
   if(linearSolver_==CVDIAG)
   {
       flag = CVDiag(cvode_mem);
       if(check_flag(&flag, "CVDiag", 1)) return; 	
       printf("IntegratorCvode will use 'CVDiag' as linear solver.  \n");
   }
   else if(linearSolver_==CVDENSE)
   {
       flag = CVDense(cvode_mem, NEQ);
       if(check_flag(&flag, "CVDense", 1)) return; 	

       // Set the user-supplied Jacobian routine Jac, TODO: currently not done! 
       printf("IntegratorCvode will use 'CVDense' as linear solver.  \n");
   }
   else if(linearSolver_==CVBAND)
   {
       flag = CVBand(cvode_mem, NEQ, 1, 1);
       if(check_flag(&flag, "CVBand", 1)) return; 	

       // Set the user-supplied Jacobian routine Jac 
       flag = CVDlsSetBandJacFn(cvode_mem, cvodes_jac);
       if(check_flag(&flag, "CVDlsSetBandJacFn", 1)) return; 
       printf("IntegratorCvode will use 'CVBand' as linear solver.  \n");
   }
   else if(linearSolver_==CVSPGMR)
   {
       flag = CVSpgmr(cvode_mem, PREC_NONE, 0);
       if(check_flag(&flag, "CVSpgmr", 1)) return; 	
       printf("IntegratorCvode will use 'CVSpgmr' as linear solver.  \n");
   }
   else if(linearSolver_==CVSPBCG)
   {
       flag = CVSpbcg(cvode_mem, PREC_NONE, 0);
       if(check_flag(&flag, "CVSpbcg", 1)) return; 	
       printf("IntegratorCvode will use 'CVSpbcg' as linear solver.  \n");
   }
   else if(linearSolver_==CVSPTFQMR)
   {
       flag = CVSptfqmr(cvode_mem, PREC_NONE, 0);
       if(check_flag(&flag, "CVSptfqmr", 1)) return; 	
       printf("IntegratorCvode will use 'CVSptfqmr' as linear solver.  \n");
   }
   else
       error().throw_error_one(FLERR,"linearSolver_ incorrectly specified in the CVODE integrator. Check: settings/integrator.json/integrator/linearSolver.\n");  

   printf("IntegratorCvode successfully initialized with %d grid points \n", MX );
   
   return;
}

///////////////////////////////////////////////////////////////////////////

void IntegratorCvode::integrate_begin(const char* stateType, int nGridPointsUsed, int dataID, bool updatePhaseFraction)
{
   double* intraDataPtr;
   //printf("number of particles in integrator = %i, actual particle ID \n",particleData().nbody());
   for(particleID=0; particleID<particleData().nbody(); particleID++)
   {
       if(updatePhaseFraction)
       {
           particleData().setParticleIDPointerPhaseFraction(particleID);
           particleData().returnPhaseFractionData(tempPhaseDataGas_, tempPhaseDataLiquid_);
           particleData().phaseFractionChangeRateStart(particleID,tempPhaseDataGas_, tempPhaseDataLiquid_);
           for(int iGrid = 0; iGrid<(nGridPointsUsed-1); iGrid++) //no need to save outer grid point
           {
             tempPhaseDataSolid_[iGrid] = 1.
                                        - tempPhaseDataGas_[iGrid]
                                        - tempPhaseDataLiquid_[iGrid];
             //printf("tempPhaseDataSolid_ [%i] = %g, tempPhaseDataGas_ = %g, tempPhaseDataLiquid_[iGrid] = %g \n",iGrid,tempPhaseDataSolid_[iGrid],tempPhaseDataGas_[iGrid],tempPhaseDataLiquid_[iGrid]);
           }
           if(dataID==GAS)  
             intraDataPtr = tempPhaseDataGas_;    //just set pointer
           else if(dataID==LIQUID) 
             intraDataPtr = tempPhaseDataLiquid_; //just set pointer
           else if(dataID==SOLID)
           {
               error().throw_error_one(FLERR,"IntegratorCvode detected that you like to update the solid phase fraction. This is not allowed. \n");
           }
           else
               error().throw_error_one(FLERR,"IntegratorCvode was handed over an invalid phase id.\n");

       }
       else
       {
           particleData().setParticleIDPointer(dataID,particleID);	
           particleData().returnIntraData(tempIntraData_);
           intraDataPtr = tempIntraData_;
       }

	   realtype *udata;
       if(strcmp(stateType,"udata") == 0)
       {
          //set the current udata
          udata = NV_DATA_S(u);

          for (int h=1; h <= MX; h++) 
          {
            //printf("udata [%i] = %g \n",h,udata[h]);
            udata[h] = intraDataPtr[h-1];
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

	   udata = NV_DATA_S(u);
		
	   flag = CVode(cvode_mem, T1, u, &t, CV_NORMAL); 				
	   if(check_flag(&flag, "CVode", 1)) 
       {
            printf("Problem when integrating for dataID %d.\n", dataID);
            error().throw_error_one(FLERR,"IntegratorCvode could not finish command 'CVode' without a flag. Must abort.\n");
       }


	   double radialDist(0.0), deltaR(0.0);
	   deltaR = 1.0/(MX-1);
		  
	   flag = CVodeGetNumSteps(cvode_mem, &nst);
	   check_flag(&flag, "CVodeGetNumSteps", 1);
      
       for (int h=1; h <= MX; h++) 
            intraDataPtr[h-1] = udata[h];

       if(updatePhaseFraction)
       {
           if(dataID==GAS)   //have written into tempPhaseDataGas_
            for(int iGrid = 0; iGrid<(nGridPointsUsed-1); iGrid++)
                tempPhaseDataLiquid_[iGrid] = 1.0 
                                            - tempPhaseDataSolid_[iGrid]
                                            - tempPhaseDataGas_[iGrid];
           else if(dataID==LIQUID) //have written into tempPhaseDataLiquid_
            for(int iGrid = 0; iGrid<(nGridPointsUsed-1); iGrid++)
                tempPhaseDataGas_[iGrid] = 1.0 
                                         - tempPhaseDataSolid_[iGrid]
                                         - tempPhaseDataLiquid_[iGrid];

           particleData().savePhaseFractionData(particleID, tempPhaseDataGas_, tempPhaseDataLiquid_);
           particleData().phaseFractionChangeRateEnd(particleID, tempPhaseDataGas_, tempPhaseDataLiquid_, control().simulationState().deltaT());
       }
       else
          particleData().saveIntraParticleData(dataID, particleID, intraDataPtr);



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
   check_flag(&flag, "CVDlfpsGetNumRhsEvals", 1);

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


