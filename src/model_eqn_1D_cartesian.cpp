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



#include "model_eqn_1D_cartesian.h"

using namespace PASCAL_NS;

/* ----------------------------------------------------------------------
   ModelExample Constructor
------------------------------------------------------------------------- */

ModelEqn1DCartesian::ModelEqn1DCartesian(ParScale *ptr, char *name) :
      ModelEqn(ptr, name)
{
        printf("...this is a ModelEqn1DCartesian. \n");
        debug_ = true;
}

//------------------------------------
void ModelEqn1DCartesian::init(int narg, char const* const* arg, int eqnType, int modelEqnID)
{
    //Allocate temporary memory to exchange particle data with particleData class
    nGridPointsUsed_ = particleMesh().nGridPoints() + 1; // allocate one interior and one exterior grid point
    modelingApproach = CONTINUUM;    
    ModelEqn::init(narg,arg,eqnType,modelEqnID);
}

//------------------------------------
void ModelEqn1DCartesian::begin_of_step()
{
        if(debug_)
            printf("begin_of_step: advancing modelEqn %s ... \n", name());

       integrator().integrate_begin("udata", nGridPointsUsed_,particleDataID_);
}

//////////////////////////////////////////////////////////////////////////////////////////

void ModelEqn1DCartesian::eval(double t, double* udata, double* dudata, double* p)
{
	 double uh, ult, urt, hordc, horac, hdiff, hadv;
        if(debug_)
            printf("ModelEqn1DCartesian::eva  ... \n");
	// Loop over all grid points.
        for (int h=0;  h < particleMesh().nGridPoints(); h++)
        {
            //dudata[h]=udata[h];

            // Extract u at x_i and two neighboring points
            uh = udata[h];
            ult = (h == 1)  ? 0 : udata[h-1];
            urt = (h == particleMesh().nGridPoints()) ? 0 : udata[h+1];

            // Set diffusion and advection terms and load into udot
            hdiff = hordc*(ult - 2.0*uh + urt);	        //CDS
            hadv = horac*(urt - ult);			//CDS
            dudata[h] = hdiff + hadv;
            }
}

//////////////////////////////////////////////////////////////////////////////////////////

void ModelEqn1DCartesian::computeParticleAverages()
{
    for(int particleID=0; particleID<particleData().nbody(); particleID++)
    {
        tempAvData_ = 0.0;
        double xa(0.0), xi(0.0);
        double volume(0.0), deltaVolume(0.0);

        //integrate over 1D layer, volume and position (xi, xa) are normalized with dx
        for(int i=0; i<particleMesh().nGridPoints(); i++)
        {
            xa = (double)i + 0.5; xa = min(xa,(double)(particleMesh().nGridPoints()-1));
            xi = (double)i - 0.5; xi = max(xi,0.0);
            deltaVolume  = xa-xi;
            volume      += deltaVolume;
            tempAvData_ += deltaVolume*tempIntraData_[i];
    //        printf("i: %d, xi: %g xa: %g, volume: %g, tempAvData: %g. \n",
    //                i, xi, xa, volume, tempAvData_);
        }
        tempAvData_ /= volume;

        particleData().saveIntraParticleAv(particleDataID_, particleID, tempAvData_);
   }
    
}

//////////////////////////////////////////////////////////////////////////////////////////

void ModelEqn1DCartesian::computeSurfaceFluxes()
{
    for(int particleID=0; particleID<particleData().nbody(); particleID++)
    {
        tempPartFlux_ = 0.0;
        particleData().saveIntraParticleFlux(particleDataID_, particleID, tempPartFlux_);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////

void ModelEqn1DCartesian::updateProperties()
{
    for(int particleID=0; particleID<particleData().nbody(); particleID++)
    {
        //set grid coefficients and BC condition
        rMAX=particleData().pRadius(particleID);

        particleData().setParticleIDPointer(particleDataID_,particleID);	
        particleData().returnIntraData(tempIntraData_);

        BCvalue[1] = tempIntraData_[particleMesh().nGridPoints()]; //set environment variable, re-use intra-particle container for this purpose
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

                diffu_eff_ = diffusivity->value();
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
        else
        {
            environmentFlux=0.0;
        }

        //printf("using diffu_eff_: %g, alpha: %g, biot_num: %g.\n",
        //       diffu_eff_, alpha, biot_num);
    }
    return;
}
