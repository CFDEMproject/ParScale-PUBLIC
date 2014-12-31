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



#include "model_eqn_shrinkingCore.h"
#include "output.h"
#include "math_extra_pascal.h"

using namespace PASCAL_NS;

#define SMALLESTRADIUS 1e-64
#define LARGENUMBER    1e32
/* ----------------------------------------------------------------------
   ModelExample Constructor
------------------------------------------------------------------------- */

ModelEqnShrinkingCore::ModelEqnShrinkingCore(ParScale *ptr, char *name) :
      ModelEqn(ptr, name)
{
}

//------------------------------------
void ModelEqnShrinkingCore::init(int narg, char const* const* arg, int eqnType, int modelEqnID)
{
    nGridPointsUsed_ = 2; // allocate one interior (used for the relative position of the reaction front) and one exterior grid point
    modelingApproach = SHRINKINGCORE;
    ModelEqn::init(narg,arg,eqnType,modelEqnID);

}
//------------------------------------
void ModelEqnShrinkingCore::begin_of_step()
{
    //updateProperties();             //this will update model-internal properties.
    integrator().integrate_begin("udata", nGridPointsUsed_,particleDataID_); //this will start the integrator
}

//------------------------------------
void ModelEqnShrinkingCore::eval(double t, double* udata, double* dudata, double* p)
{

    if(modelChemistryContainer().nrChemistryEqns()!=1)
        error().throw_error_one(FLERR,"ERROR: Shrinking core more requires exactly one chemistry eqn! \n");

    particleID = integrator().returnParticleID();
    updateProperties(particleID);             //this will update model-internal properties.

    //update chemistry
    particleData().returnchemistryDataPoint(particleDataID_, particleID, 0, kSurface_);

    //Returns the time derivative of the DIMENSIONLESS core radius u!
   	h          = 1;
    int hFluid = 2; //index of fluid variable

    double Apart = rMAX*rMAX;      //4*pi missing on purpose
    double rCore = max(SMALLESTRADIUS, udata[h]) * rMAX;
    double ACore = rCore*rCore;
    const  double invCapacity = 1.0     //4*pi missing on purpose
                            /(
                                  ACore
                                * rho_solid
                                * rMAX      //normalization factor, 
                                           //since we return dimensionless change of rCore=u*rMAX
                             );

     //fluid variable must not be changed
     dudata[hFluid] =0.0;

     if(BC[1]==NEUMANN)
        dudata[h] = -( 
                        flux->value() 
                       *Apart
                     ) 
                     * invCapacity;
                     
    double resistance = 0.0; //multiplied by 4*pi!
    if( BC[1]==DIRICHLET || BC[1]==CONVECTIVE )
    {
         // A - Resistance due to pore diffusion
         resistance += (1.0/rCore - 1.0/rMAX) 
                     /  diffu_eff_;

         // B - Resistance due to diffusion in the boundary layer
         if(BC[1]==CONVECTIVE)
             resistance += 1.0 
                        / (
                             max(1e-64, transfer_coeff->value() )
                            *Apart
                          );
                          
         // C - Add chemistry
         resistance += 1.0 
                        / (
                             max(1e-64, -kSurface_)
                            *ACore
                          );

         // D - Compute time derivative of rCore
         dudata[h] = -environmentU / resistance * invCapacity;
    }

#if 0
    printf("udata: %g %g, dudata: %g rCore: %g, rMAX: %g, resistance: %g, diffu_eff_: %g, nrChe: %d, kSurface: %g\n",
           udata[1],udata[2],
           dudata[1],
           rCore,
           rMAX,
           resistance,
           diffu_eff_,
           modelChemistryContainer().nrChemistryEqns(),
           kSurface_
           );
#endif
}

/////////////////////////////////////////////////////////////////////////////////
void ModelEqnShrinkingCore::returnJac(long int N, long int mu, long int ml,
                   double t, double* udata, double* fudata,
                   DlsMat J, double* p,
                   double* tmp1data, double* tmp2data, double* tmp3data)
{

    //set Jac coefficients at centre of sphere
    j=0;
    col_j = BAND_COL(J,j);
    
    if(BC[1]==NEUMANN)
        BAND_COL_ELEM(col_j,j,j) = 0.0;

    double ks = LARGENUMBER;
    double kg = LARGENUMBER;
    double u  = udata[1];
    double Nval  = 0.0;
    if( BC[1]==DIRICHLET || BC[1]==CONVECTIVE )
    {
         if(modelChemistryContainer().nrChemistryEqns()!=0)
            ks = max( 1e-64, -kSurface_);
          
         if(BC[1]==CONVECTIVE)
            kg = max( 1e-64, transfer_coeff->value() );
            
         Nval =  u*u * diffu_eff_ * ks
              + (1.0 - u)*u*rMAX*kg*ks
              + kg * diffu_eff_;
          
         BAND_COL_ELEM(col_j,j,j) = environmentU 
                                  * kg * diffu_eff_ * ks
                                  / (rMAX * rho_solid)
                                  * (
                                         2.0 * u  * diffu_eff_ * ks
                                      + (1.0 + u) * rMAX * kg * ks
                                    )
                                  / (Nval*Nval);
    }

    /*printf("BAND_COL_ELEM(col_j,j,j): %g, kg: %g, ks: %g\n",
           BAND_COL_ELEM(col_j,j,j),
           kg, ks
          );*/
           
}

/////////////////////////////////////////////////////////////////////////////////

void ModelEqnShrinkingCore::computeParticleAverages()
{
    //This is the concentration AT THE REACTING CORE, i.e., 
    //this concentration can be used to calculate the reaction rate
    
    for(int particleID=0; particleID<particleData().nbody(); particleID++)
    {
        particleData().setParticleIDPointer(particleDataID_,particleID);	
        particleData().returnIntraData(tempIntraData_);

        double rc = tempIntraData_[0]; //dim-less radius of reacting core
        rMAX     = particleData().pRadius(particleID);

        double divBy = 1.0;
        if(BC[1]==CONVECTIVE)
        {
                divBy +=  -kSurface_
                        * rc * rc
                        / max( 1e-64, transfer_coeff->value() );
        }
        divBy += -kSurface_
               * (rMAX * rc) * (1.0 - rc)
               /  max( 1e-64, diffu_eff_ );

        tempAvData_ = tempIntraData_[1]
                    / divBy;


        particleData().saveIntraParticleAv(particleDataID_, particleID, tempAvData_);    
    }
}

/////////////////////////////////////////////////////////////////////////////////

//calculate heat and species flux over the surface of each particle [J/(s*m^2)]!!
void ModelEqnShrinkingCore::computeSurfaceFluxes()
{
    for(int particleID=0; particleID<particleData().nbody(); particleID++)
    {

        //TODO: Compute flux and the gas concentration at rCore
        tempPartFlux_ = 0.0;
        {
            if(BC[1]==DIRICHLET)
            {
                tempPartFlux_ = 0.0;
            }
            else if(BC[1]==NEUMANN)
            {
                tempPartFlux_ = environmentFlux;
            }
            else if(BC[1]==CONVECTIVE)
            {
                tempPartFlux_ =  alpha * (tempIntraData_[0] - environmentU);
               // printf("heat flux in ModelEqnShrinkingCore::computeSurfaceFluxes() = %g \n", tempPartFlux_);
            }
        }
        particleData().saveIntraParticleFlux(particleDataID_, particleID, tempPartFlux_);
    }
}

/////////////////////////////////////////////////////////////////////////////////
void ModelEqnShrinkingCore::updateProperties(int particleID)
{
    for(int particleID=0; particleID<particleData().nbody(); particleID++)
    {
        //set grid coefficients and BC condition
        rMAX = particleData().pRadius(particleID);

        particleData().setParticleIDPointer(particleDataID_,particleID);	
        particleData().returnIntraData(tempIntraData_);

        BCvalue[1] = tempIntraData_[0]; //set environment variable, re-use intra-particle container for this purpose
        environmentU=BCvalue[1];
        //TODO: all this will not work in case we have locally-defined (i.e., intra-particle)
        //      variable properties!
        //set properties
        if (porosity!=NULL)
        {
            if(eqnType_==HEAT)
            {
                //get thermal conductivity from in.fil,calculate effective thermal conductivity
                lambda_solid=thermal_solid_conductivity->value();
                lambda_gas=thermal_gas_conductivity->value();
                lambda_eff=(1.0-(porosity->value()*porosity->value()))*lambda_solid
                      +(porosity->value()*porosity->value())*lambda_gas; //TODO: check this eqn.

                //get heat capacity from in.fil,calculate effective heat capacity
                c_p_solid=capacity_solid->value();
                c_p_gas=capacity_gas->value();
                c_p_eff=(1.0-porosity->value())*c_p_solid
                        +porosity->value()*c_p_gas;

                //get density from in.fil,calculate effective density
                rho_gas=density_gas->value();
                rho_eff=(1.0-porosity->value())*rho_solid
                        +porosity->value()*rho_gas;

                diffu_eff_=lambda_eff/(c_p_eff*rho_eff);
            }
            else if(eqnType_==SPECIES)
            {
                if (diffusivity==NULL)
                    error().throw_error_one(FLERR,"ERROR: You must specify a diffusivity! \n");

            }
                diffu_eff_ = diffusivity->value();
            rho_solid  = density_solid->value();
        }
        else
        {
            if(eqnType_==HEAT)
            {
                lambda_eff=thermal_solid_conductivity->value();
		        c_p_eff=capacity_solid->value();
                rho_eff=density_solid->value();
                diffu_eff_=lambda_eff/(c_p_eff*rho_eff);

            }
            else if(eqnType_==SPECIES)
            {
                if (diffusivity==NULL)
                    error().throw_error_one(FLERR,"ERROR: You must specify a diffusivity! \n");

                diffu_eff_ = diffusivity->value();
                rho_solid  = density_solid->value();
	         }
        }

        //set boundary properties
        if(BC[1]==CONVECTIVE)
        {
            if(transfer_coeff!=NULL)
		    {
                alpha=transfer_coeff->value();
		    }
		    else
            {
                printf("WARNING: you have not specified a transfer coefficient. Assuming 0.\n");
                alpha=0.;
            }
        }
        else
        {
            alpha=0.0;
        }

        if(BC[1]==NEUMANN)
        {
            if(flux!=NULL)
                environmentFlux=flux->value();
            else
            {
                printf("WARNING: you have not specified a flux. Assuming 0.\n");
                environmentFlux=0.;
            }
        }
        else
        {
            environmentFlux=0.0;
        }

        //printf("using diffu_eff_: %g, alpha: %g, biot_num: %g.\n",
        //       diffu_eff_, alpha, biot_num);
    }
    return;
}

