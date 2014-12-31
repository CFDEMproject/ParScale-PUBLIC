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



#include "model_eqn_1D_spherical.h"
#include "output.h"
#include "math_extra_pascal.h"

using namespace PASCAL_NS;

/* ----------------------------------------------------------------------
   ModelExample Constructor
------------------------------------------------------------------------- */

ModelEqn1DSpherical::ModelEqn1DSpherical(ParScale *ptr, char *name) :
      ModelEqn(ptr, name),
      boundZero_(true) //TODO: read from json, currently not used
{

  //check the particle mesh
   if(particleMesh().nGridPoints()<1)
   output().write_screen_one("WARNING: ModelEqn1DSpherical:your particle mesh has no grid points! \n");

}

//------------------------------------
void ModelEqn1DSpherical::init(int narg, char const* const* arg, int eqnType, int modelEqnID)
{
    //Allocate temporary memory to exchange particle data with particleData class
    nGridPointsUsed_ = particleMesh().nGridPoints() + 1; // allocate one interior and one exterior grid point
    printf("nGridPointsUsed_ = %i \n",nGridPointsUsed_);
    modelingApproach = CONTINUUM;
    ModelEqn::init(narg,arg,eqnType,modelEqnID);
}

//------------------------------------
void ModelEqn1DSpherical::begin_of_step()
{
    //printf("nGridPointsUsed_ = %i, particleDataID_=%i \n",nGridPointsUsed_,particleDataID_);
    integrator().integrate_begin("udata", nGridPointsUsed_,particleDataID_);
}

//------------------------------------
void ModelEqn1DSpherical::eval(double t, double* udata, double* dudata, double* p)
{
    particleID = integrator().returnParticleID();
    //printf("particle ID = %i\n",particleID);
    updateProperties(particleID);             //this will update model-internal properties.
    //printf("Coupling model name = %s \n",coupling().couplingModel().name());
    if (strcmp(coupling().couplingModel().name(),"liggghts") == 0)// && environmentFlux != 0)
    {
        //printf("Before receiving flux from LIGGGHTS\n");
        particleData().returnIntraFlux(particleID,environmentFlux);
    }

    // Loop over all grid points, from inside to outside
   	for (h=1;  h <= particleMesh().nGridPoints(); h++)
    {
          //specify BC condition at surface
          if (h == (particleMesh().nGridPoints()))
          {
            if(BC[1]==DIRICHLET)
            {
                udata[h] = BCvalue[1];
                dudata[h] = 0.0;
            }
            if(BC[1]==NEUMANN)
            {
                dudata[h] = (
                                (  2.0
                                 *(
                                       udata[h-1]
                                     - ((environmentFlux*dx)/diffu_eff_)
                                     - udata[h]
                                  )
                                )
                                /(dx*dx)
                             )
                            -(
                                 (2*environmentFlux)
                                /(diffu_eff_*rMAX)
                             );
            }
            if(BC[1]==CONVECTIVE)
            {
                //TODO:check formular here!!!!
                dudata[h] =  2.0*diffu_eff_/(dx*dx)
                           *(
                               - 1.0*biot_num*(udata[h]-environmentU)
                               + udata[h-1]-udata[h]
                            )
                           +diffu_eff_/(rMAX*dx)
                           *(
                               - 2.0*biot_num*(udata[h]
                               - environmentU)
                            );
               // printf("dudata before it gets LIGGGTHS flux at outer grid point= %g \n",dudata[h]);
               // printf("diffu_eff_ = %g, biot_num = %g, udata[h] = %g, environmentU = %g, udata[h-1] = %g, rMAX = %g, dx = %g \n",diffu_eff_,biot_num,udata[h],environmentU,udata[h-1],rMAX,dx);  
                //Flux due to possible Coupling to LIGGGHTS 
                if (strcmp(coupling().couplingModel().name(),"liggghts") == 0)
                { 
                    dudata[h] +=
                        (environmentFlux/(rho_eff*c_p_eff));
                }
            }
          }

          //specify BC condition in middle of sphere
          else if (h == 1)
          {
            dudata[h] =  diffu_eff_
                        *6.0
                        *(
                             ( 1.0/(dx*dx) )
                           * ( udata[h+1]-udata[h] )
                         );
          }

          //apply numerical scheme (see docu) for all points between middle and surface
          else
          {
              x_coeff_1st_dev = dx*(h-1);
              coeff_1st_dev   = (2.0)/(2.0*x_coeff_1st_dev*dx);
              dudata[h] = diffu_eff_*
                        (
                            (  coeff_2nd_dev*(udata[h-1]
                             - 2.0*udata[h]
                             + udata[h+1])
                            )
                           +(  coeff_1st_dev*(udata[h+1]
                             - udata[h-1])
                            )
                        );
             //printf("diffu_eff_ = %g, biot_num = %g, udata[h] = %g, environmentU = %g, udata[h-1] = %g, rMAX = %g, dx = %g \n",diffu_eff_,biot_num,udata[h],environmentU,udata[h-1],rMAX,dx);  
            //printf("dudata before it gets LIGGGTHS flux at all grid points= %g \n",dudata[h]);
          }

         //add chemistry
         double b_ = 0.0;
         if(modelChemistryContainer().nrChemistryEqns()!=0)
         {
            particleData().returnchemistryDataPoint(particleDataID_, particleID, h-1, b_);
            //printf("returned value for gridpoint %i = %g \n",h-1,b_);
          
            if(eqnType_==SPECIES)
            {
                //TODO: check to get correct dimensions!
                dudata[h] += b_;
            }
            else if (eqnType_==HEAT)
            {
                //TODO: check to get correct dimensions!
                dudata[h] += b_
                          / (c_p_eff*rho_eff);
            } 
            
         }

#if 0
    printf("udata[%d]: %g, dudata: %g, b: %g.\n",
           h,
           udata[h],
           dudata[h],
           b_
           );
#endif
    } //End loop grid points
    //printf("udata eval [surfcae]=%g \n",udata[10]);
}

/////////////////////////////////////////////////////////////////////////////////
void ModelEqn1DSpherical::returnJac(long int N, long int mu, long int ml,
                   double t, double* udata, double* fudata,
                   DlsMat J, double* p,
                   double* tmp1data, double* tmp2data, double* tmp3data)
{
    particleID = integrator().returnParticleID();

    //set Jac coefficients at centre of sphere
    j=0;
    col_j = BAND_COL(J,j);
    BAND_COL_ELEM(col_j,j,j)   = -6.0*diffu_eff_*coeff_2nd_dev;
    BAND_COL_ELEM(col_j,j+1,j) = 6.0*diffu_eff_*coeff_2nd_dev;

    j=particleMesh().nGridPoints()-1;
    if(BC[1]==DIRICHLET)
    {
        BAND_COL_ELEM(col_j,j,j) = 0.0;
    }
    else if(BC[1]==NEUMANN)
    {
        BAND_COL_ELEM(col_j,j,j) = -(2.0/dx*dx);
        BAND_COL_ELEM(col_j,j-1,j) = (2.0/dx*dx);
    }
    else if(BC[1]==CONVECTIVE)
    {
        //printf("jac Convective!\n");
        BAND_COL_ELEM(col_j,j,j) = ((-2.0*diffu_eff_)/dx)*(biot_num/rMAX + biot_num/dx - 1.0/dx);        
        BAND_COL_ELEM(col_j,j-1,j) = (2.0*diffu_eff_)/(dx*dx);
   }

    for (h=2; h < particleMesh().nGridPoints(); h++)
    {
        x_coeff_1st_dev=dx*(h-1);
        coeff_1st_dev = (2.0)/(2.0*x_coeff_1st_dev*dx);
        j = h-1;					            //h runs from 1 to MX
                                    		    //j (Jacobian index) runs from 0 to MX-1
        col_j = BAND_COL(J,j);				    //BAND_COL BAND_COL_ELEM supposed to be faster - could be replaced by
        BAND_COL_ELEM(col_j,j,j) = -2.0*diffu_eff_*(coeff_2nd_dev); 	//diagonal elements of Jac
      	BAND_COL_ELEM(col_j,j-1,j) = diffu_eff_*(coeff_2nd_dev-coeff_1st_dev); //below the diagonal
       	BAND_COL_ELEM(col_j,j+1,j) = diffu_eff_*(coeff_2nd_dev+coeff_1st_dev); //above the diagonal
    }


    //Add chemical term
    if(modelChemistryContainer().nrChemistryEqns()!=0)
    {
      for(h=1; h <= particleMesh().nGridPoints(); h++)
      {

         double b_ = 0.0;
         particleData().returnchemistryDataPoint(particleDataID_, particleID, h-1, b_);
         if(eqnType_==SPECIES)
         {
                //TODO: check to get correct value
         }
         else if (eqnType_==HEAT)
         {
                //TODO: check to get correct value
         } 
        //TODO:
        //j = h-1;
        //col_j = BAND_COL(J,j);
        //printf("chemical source term jac: %g \n", chemistry_->returnSourceJac(eqnType_, eqnID_)[j]);
        //BAND_COL_ELEM(col_j,j,j) += chemistry_->returnSourceJac(eqnType_, eqnID_)[j];
      }
    }
}

/////////////////////////////////////////////////////////////////////////////////

void ModelEqn1DSpherical::computeParticleAverages()
{
    for(int particleID=0; particleID<particleData().nbody(); particleID++)
    {
        particleData().setParticleIDPointer(particleDataID_,particleID);	
        particleData().returnIntraData(tempIntraData_);

        tempAvData_ = 0.0;
        double ra(0.0), ri(0.0);
        double volume(0.0), deltaVolume(0.0);

        //integrate over spherical shells, volume and radii are normalized with dx
        for(int i=0; i<particleMesh().nGridPoints(); i++)
        {
            ra = (double)i + 0.5; ra = min(ra,(double)(particleMesh().nGridPoints()-1));
            ri = (double)i - 0.5; ri = max(ri,0.0);
            deltaVolume  = ra*ra*ra-ri*ri*ri;
            volume      += deltaVolume;
            tempAvData_ += deltaVolume*tempIntraData_[i];
    //        printf("i: %d, ri: %g ra: %g, volume: %g, tempAvData: %g. \n",
    //                i, ri, ra, volume, tempAvData_);
        }
        tempAvData_ /= volume;
        
        particleData().saveIntraParticleAv(particleDataID_, particleID, tempAvData_);
    }
}

/////////////////////////////////////////////////////////////////////////////////

//calculate heat and species flux over the surface of each particle [J/(s*m^2)]!!

void ModelEqn1DSpherical::computeSurfaceFluxes()
{
    for(int particleID=0; particleID<particleData().nbody(); particleID++)
    {
        tempPartFlux_ = 0.0;
        particleData().setParticleIDPointer(particleDataID_,particleID);	
        particleData().returnIntraData(tempIntraData_);

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
            tempPartFlux_ =  alpha * (tempIntraData_[particleMesh().nGridPoints()-1] - environmentU);
            //printf("heat flux in ModelEqn1DSpherical::computeSurfaceFluxes() = %g \n", tempPartFlux_);
        }
        particleData().saveIntraParticleFlux(particleDataID_, particleID, tempPartFlux_);
    }
}

/////////////////////////////////////////////////////////////////////////////////

void ModelEqn1DSpherical::updateProperties(int particleID)
{
        //set grid coefficients and BC condition
        rMAX=particleData().pRadius(particleID);
#if 0
        printf("rMAX = %g \n",rMAX);
        printf("particleDataID = %i, particleID = %i  of %i particles \n", particleDataID_, particleID,particleData().nbody());
#endif

        particleData().setParticleIDPointer(particleDataID_,particleID);	
        particleData().returnIntraData(tempIntraData_); 
     
        BCvalue[1] = tempIntraData_[particleMesh().nGridPoints()]; //set environment variable, re-use intra-particle container for this purpose
        //printf("BC value = %g \n",BCvalue[1]); 
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
                rho_solid=density_solid->value();
                rho_gas=density_gas->value();
                rho_eff=(1.0-porosity->value())*rho_solid
                        +porosity->value()*rho_gas;

                diffu_eff_=lambda_eff/(c_p_eff*rho_eff);
            }
            else if(eqnType_==SPECIES)
            {
                if (diffusivity==NULL)
                    error().throw_error_one(FLERR,"ERROR: You must specify a diffusivity! \n");


                diffusivity_= diffusivity->value();
                tortuosity_ = tortuosity->value();

                if (tortuosity_==NULL)
                    error().throw_error_one(FLERR,"ERROR: You must specify a tortuosity if you want to calculate a species equation within a pourous particle! \n");

                diffu_eff_ = diffusivity_*(porosity->value()/tortuosity_);
            }
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
                 
	         }
        }

        //set boundary properties
        dx = rMAX/((particleMesh().nGridPoints())-1.0);

        //coeffcient results off 2nd deviation O(dx^2)
        coeff_2nd_dev = (1.0)/(dx*dx);

        if(BC[1]==CONVECTIVE)
        {
            if(transfer_coeff!=NULL)
		    {
                alpha=transfer_coeff->value();
			    //printf("alpha = %g \n", alpha);
		    }
		    else
            {
                printf("WARNING: you have not specified a transfer coefficient. Assuming 0.\n");
                alpha=0.;
            }

            //Biot number
            if(eqnType_==HEAT)
                biot_num = (alpha*dx)/lambda_eff;   //TODO: check this
            else if(eqnType_==SPECIES)
                biot_num = (alpha*dx)/diffu_eff_;   //TODO: check this
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
        /*else
        {
            environmentFlux=0.0;
        }*/

        //printf("using diffu_eff_: %g, alpha: %g, biot_num: %g.\n",
        //       diffu_eff_, alpha, biot_num);
        return;
    
}


