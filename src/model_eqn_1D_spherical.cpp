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


#include "model_phasechange_container.h"
#include "model_eqn_1D_spherical.h"
#include "output.h"
#include "comm.h"
#include "math_extra_pascal.h"

using namespace PASCAL_NS;
#define VERBOSE_1DSPHERICAL TRUE    //developer to hard-code here for debug
#define APPLYCONVECTIVEFLUX TRUE   //developer to hard-code here for testing
/* ----------------------------------------------------------------------
   ModelExample Constructor
------------------------------------------------------------------------- */

ModelEqn1DSpherical::ModelEqn1DSpherical(ParScale *ptr, char *name) :
      ModelEqn(ptr, name),
      boundZero_(true), //TODO: read from json, currently not used
      tempDiff_eff_(NULL),
      tempConvectionSpeed_eff_(NULL),
      tempTempData_(NULL),
      tempPhaseFractGas_(NULL),
      tempPhaseFractLiq_(NULL),
      tempPhaseFractSolid_(NULL),
      current_modelEqn_ID_(0)
{

  //check the particle mesh
   if(particleMesh().nGridPoints()<1)
       output().write_screen_one("WARNING: ModelEqn1DSpherical:your particle mesh has no grid points! \n");

}
//------------------------------------
ModelEqn1DSpherical::~ModelEqn1DSpherical()
{
    destroy<double>(tempDiff_eff_);
    destroy<double>(tempConvectionSpeed_eff_);
    destroy<double>(tempTempData_);
    destroy<double>(tempPhaseFractGas_);
    destroy<double>(tempPhaseFractLiq_);
    destroy<double>(tempPhaseFractSolid_);
}

//------------------------------------
void ModelEqn1DSpherical::init(int narg, char const* const* arg, int eqnType, int modelEqnID)
{
    //current_modelEqn_ID_ = modelEqnID;
    //Allocate temporary memory to exchange particle data with particleData class
    nGridPointsUsed_ = particleMesh().nGridPoints() + 1; // allocate one interior and one exterior grid point
    //printf("nGridPointsUsed_ = %i \n",nGridPointsUsed_);
    modelingApproach = CONTINUUM;
    ModelEqn::init(narg,arg,eqnType,modelEqnID);
    tempDiff_eff_        = create<double>(tempDiff_eff_, nGridPointsUsed_-1);
    tempConvectionSpeed_eff_ = create<double>(tempConvectionSpeed_eff_, nGridPointsUsed_-1);
    tempTempData_        = create<double>(tempTempData_, nGridPointsUsed_-1);
    tempPhaseFractGas_   = create<double>(tempPhaseFractGas_, nGridPointsUsed_-1);
    tempPhaseFractLiq_   = create<double>(tempPhaseFractLiq_, nGridPointsUsed_-1);
    tempPhaseFractSolid_ = create<double>(tempPhaseFractSolid_, nGridPointsUsed_-1);

    //Set pointers to sub-functions. TODO: put into separate class for more flexibility
    phaseFluxGas    = &ModelEqn1DSpherical::phaseFluxIntegrate;
    phaseFluxLiquid = &ModelEqn1DSpherical::phaseFluxCapillarity;

    if(IsoThermal_)
        read_chemistry_single_react_json_file("temperature",  &isoTemp_, true);
}

//------------------------------------
void ModelEqn1DSpherical::begin_of_step()
{
    if(!solveMe)    return;
    if(updatePhaseFraction)
        integrator().integrate_begin("udata", nGridPointsUsed_, inPhase_,        true);
    else
        integrator().integrate_begin("udata", nGridPointsUsed_, particleDataID_, false);

}

//------------------------------------
void ModelEqn1DSpherical::eval(double t, double* udata, double* dudata, double* p)
{
    if(verbose_)
    {
          if (eqnType_==HEAT)
            printf("Eval heat here!\n");
          if (eqnType_==SPECIES && inPhase_==SOLID)
            printf("Eval a SOLID here!\n");

          if (eqnType_==SPECIES && inPhase_==LIQUID)
          {
            printf("Eval a LIQUID here!\n");
            if (updatePhaseFraction)
                printf("Update Phase Fraction activated \n");
            else if (!updatePhaseFraction)
                printf("No Update Phase Fraction \n");
           }

          if (eqnType_==SPECIES && inPhase_==GAS)
            printf("Eval a GAS here!\n");
    }
    double tOld;
    double *** phaseChangeVol(NULL);
    if(modelPhaseChangeContainer().nrPhaseChangeEqns()>0)
        phaseChangeVol = ParScaleBase::particleData().accessPhaseChangeRateVolumetric(inPhase_);

    particleID = integrator().returnParticleID();
    //printf("particle ID = %i\n",particleID);
    updateProperties(particleID);       //this will update model-internal properties.

    //Reset dudata in order to store data
       for (h=1;  h <= particleMesh().nGridPoints(); h++)
        dudata[h]=0;
    // Loop over all grid points, from inside to outside
       for (h=1;  h <= particleMesh().nGridPoints(); h++)
    {
          int       iGrid = h-1; //position in the data grid

          //0 - Term due to change of volume of the phase
          if(modelPhaseChangeContainer().nrPhaseChangeEqns()>0)
          {
            if (inPhase_==LIQUID)
              dudata[h] -= (phaseChangeVol[particleID][0][iGrid]*udata[h]) / tempPhaseFractLiq_[iGrid];

            if (inPhase_ == GAS)
              dudata[h] -= (phaseChangeVol[particleID][0][iGrid]*udata[h]) / tempPhaseFractGas_[iGrid];
          }

          //specify BC condition at surface
          if (h == particleMesh().nGridPoints() )
          {
            if(BC[1]==DIRICHLET)
                udata[h]  = BCvalue[1];

            if(BC[1]==NEUMANN)
            {
                //1 - Diffusive terms
                dudata[h]+= tempDiff_eff_[h-1]*((
                                (  2.0
                                 *(
                                       udata[h-1]
                                     - (((-environmentFlux/surface_area)*dx_)/lambda_eff_)
                                     - udata[h]
                                  )
                                )
                                /(dx_*dx_)
                             )
                            -(
                                 (2.0*(-environmentFlux/surface_area))        //environmentFlux [W/m²], positive = heating
                                /(lambda_eff_*rMAX)
                             ));

                  //2 - Convective terms

                  currentRadius_   = dx_*(h-1);
#if APPLYCONVECTIVEFLUX
                  if(solveConvectiveFlux &&
                     control().simulationState().timeStep()>1)
                  {
                    double convContrib
                          = -2.0 / currentRadius_
                          * (tempConvectionSpeed_eff_[iGrid-1]+tempConvectionSpeed_eff_[iGrid])*0.5
                          * udata[h]
                          - 1.0 / dx_
                           * (
                                   tempConvectionSpeed_eff_[iGrid]   * udata[h]     //upwind
                                  -tempConvectionSpeed_eff_[iGrid-1] * udata[h-1]   //upwind
                             );
                      if (inPhase_==LIQUID && tempPhaseFractLiq_[iGrid]>phaseFractionMinimumConvection)
                          dudata[h] += convContrib / tempPhaseFractLiq_[iGrid];
                      if (inPhase_ == GAS && tempPhaseFractGas_[iGrid]>phaseFractionMinimumConvection)
                          dudata[h] += convContrib / tempPhaseFractGas_[iGrid];
                  }
#endif

                  if(verbose_)
                      printf("dudata from NEUMANN BC= %g, environmentFlux: %g \n",dudata[h],-environmentFlux);

            }
            if(BC[1]==CONVECTIVE)
            {
                dudata[h]+=  2.0*tempDiff_eff_[h-1]/(dx_*dx_)
                           *(
                               - biot_num_*(udata[h]-environmentU)
                               + udata[h-1]-udata[h]
                            )
                           +tempDiff_eff_[h-1]/(rMAX*dx_)
                           *(
                               - 2.0*biot_num_*(udata[h]
                               - environmentU)
                            );

                  currentRadius_   = dx_*(h-1);
#if APPLYCONVECTIVEFLUX
                  if(solveConvectiveFlux &&
                     control().simulationState().timeStep()>1)
                  {
                    double convContrib
                          = -2.0 / currentRadius_
                          * (tempConvectionSpeed_eff_[iGrid-1]+tempConvectionSpeed_eff_[iGrid])*0.5
                          * udata[h]
                          - 1.0 / dx_
                           * (
                                   tempConvectionSpeed_eff_[iGrid]   * udata[h]     //upwind
                                  -tempConvectionSpeed_eff_[iGrid-1] * udata[h-1]   //upwind
                             );
                        if (inPhase_==LIQUID && tempPhaseFractLiq_[iGrid]>phaseFractionMinimumConvection)
                            dudata[h] += convContrib / tempPhaseFractLiq_[iGrid];
                        if (inPhase_ == GAS && tempPhaseFractGas_[iGrid]>phaseFractionMinimumConvection)
                            dudata[h] += convContrib / tempPhaseFractGas_[iGrid];
                  }
#endif
                //Flux due to possible Coupling to LIGGGHTS- unit of environmentFlux has to be [W/m³]
                if (strcmp(coupling().couplingModel().name(),"liggghts") == 0 && eqnType_==HEAT)
                {
                    if(verbose_)
                    {
                      printf("dudata before it gets LIGGGTHS heat flux at outer grid point= %g \n",dudata[h]);
                      printf("transfer_coeff: %g, environmentU = %g. \n",
                              transfer_coeff->value(), environmentU);
                    }

                    double ra(0.0), ri(0.0);
                    ra = (h-1)*dx_;
                    ri = (h-2)*dx_;
                    segment_vol = 4.18879020479*(ra*ra*ra-ri*ri*ri); //volume of the shell, 4/3*pi*(ra^3-ri^3)
                    dudata[h] +=
                        ((environmentFlux/segment_vol)/(rho_eff_*c_p_eff_));
                }

                if(verbose_)
                {
                  printf("dudata for convective BC = %g \n",dudata[h]);
                  printf("Biot number = %g\n",biot_num_);
                  printf("transfer_coeff: %g, environmentU = %g. \n",
                          transfer_coeff->value(), environmentU);
                }
            }
          }

          //specify BC condition in middle of sphere
          else if (h == 1)
          {
            dudata[h]+=  tempDiff_eff_[h-1]
                        *6.0
                        *(
                             ( 1.0/(dx_*dx_) )
                           * ( udata[h+1]-udata[h] )
                         );
          }

          //apply numerical scheme (see docu) for all points between center and surface
          else
          {
              currentRadius_   = dx_*(h-1);
              coeff_1st_dev_   = 1.0/(currentRadius_*dx_);

              //1 - Diffusive terms
              dudata[h]+= tempDiff_eff_[iGrid]*
                        (
                            (  coeff_2nd_dev_*(udata[h-1]
                             - 2.0*udata[h]
                             + udata[h+1])
                            )
                           +(  coeff_1st_dev_*(udata[h+1]
                             - udata[h-1])
                            )
                        );
              //printf("dudata[%i] = %g, tempDiff_eff_[%i] = %g, coeff_2nd_dev_ = %g, coeff_1st_dev_ = %g \n",h,dudata[h], iGrid, tempDiff_eff_[iGrid],coeff_2nd_dev_,coeff_1st_dev_);

#if APPLYCONVECTIVEFLUX
              //2 - Convective terms (only apply if non-zero, handle fronts correctly)
              double convContrib = 0;
              if(solveConvectiveFlux &&
                 control().simulationState().timeStep()>1)
              {
               if(  (inPhase_==LIQUID && tempPhaseFractLiq_[iGrid]>phaseFractionMinimumConvection)
                  ||(inPhase_==GAS    && tempPhaseFractGas_[iGrid]>phaseFractionMinimumConvection)
                 )
               {

                 bool upwindFront = false;
                 if(tempPhaseFractGas_[iGrid-1]<phaseFractionMinimumConvection)
                    upwindFront = true;
                 //bool downwindFront = false;
                 //if(tempPhaseFractGas_[iGrid+1]<phaseFractionMinimumConvection)
                 //   downwindFront = true;

                 convContrib = -2.0 / currentRadius_
                             * (tempConvectionSpeed_eff_[iGrid-1] + tempConvectionSpeed_eff_[iGrid]) *0.5
                             * udata[h];

                 if(upwindFront) //treat upwind front
                 {
                  double convSpeedUpwind = 2.0*tempConvectionSpeed_eff_[iGrid]-tempConvectionSpeed_eff_[iGrid+1]; //extrapolate
                  convContrib +=
                              - 1.0 / dx_
                              * (
                                   tempConvectionSpeed_eff_[iGrid]   * (udata[h]  +udata[h+1])*0.5
                                  -convSpeedUpwind                   * (udata[h-1]+udata[h]  )*0.5
                                );
                 }
                 else
                  convContrib +=
                              - 1.0 / dx_
                              * (
                                   tempConvectionSpeed_eff_[iGrid]   * (udata[h]  +udata[h+1])*0.5
                                  -tempConvectionSpeed_eff_[iGrid-1] * (udata[h-1]+udata[h]  )*0.5
                                );

                if (inPhase_==LIQUID)
                    dudata[h] += convContrib / tempPhaseFractLiq_[iGrid];

                if (inPhase_ == GAS)
                    dudata[h] += convContrib / tempPhaseFractGas_[iGrid];
               }

              }
#endif

         } //end loop over inner points

         for(int iPhaseChangeModel=0;iPhaseChangeModel<modelPhaseChangeContainer().nrPhaseChangeEqns(); iPhaseChangeModel++)
         {
            double pDot    = 0.0;
            double pDotJac = 0.0;
            double derivativeConcTemp = 0;
            const vector<int>* phaseID = modelPhaseChangeContainer().modelPhaseChangeEqn(iPhaseChangeModel)->phaseID(); //this is the enueration refering to the phase (e.g., SOLID=0, GAS=1, LIQUID=2 )
            const vector<int>* dataIds = modelPhaseChangeContainer().modelPhaseChangeEqn(iPhaseChangeModel)->speciesParticleDataID();

            if( (*dataIds)[0]==particleDataID_ )     //Dense phase
            {
                 particleData().returnPhaseChangeRateDataPoint(    (*phaseID)[0], particleID, h-1, pDot);
                 particleData().returnPhaseChangeRateJacDataPoint((*phaseID)[0], particleID, h-1, pDotJac);
            }
            else if( (*dataIds)[1]==particleDataID_ ) //Dilute phase
            {
                 particleData().returnPhaseChangeRateDataPoint(    (*phaseID)[1], particleID, h-1, pDot);
                 particleData().returnPhaseChangeRateJacDataPoint((*phaseID)[1], particleID, h-1, pDotJac);
            }
            else if( eqnType_==HEAT )
            {
                 double pDotDilute    = 0;
                 double pDotJacDilute = 0;
                 double diluteConc    = 0;
                 particleData().returnPhaseChangeRateDataPoint(    (*phaseID)[1], particleID, h-1, pDotDilute);
                 particleData().returnPhaseChangeRateJacDataPoint((*phaseID)[1], particleID, h-1, pDotJacDilute);
                 diluteConc = particleData().retrieveIntraData(   (*dataIds)[1], particleID, h-1);
                 pDot = pDotDilute + pDotJacDilute*diluteConc;

                 //old temperature
                 particleData().returnPhaseChangeRateJacLastSlotDataPoint(particleID, h-1, tOld);

                 //retrieve derivative of saturation concentration
                 modelPhaseChangeContainer().modelPhaseChangeEqn(iPhaseChangeModel)->derivativeConcTemp(udata[h], derivativeConcTemp);
                 pDotJac = -1*pDotJacDilute*derivativeConcTemp;

                 //retrieve heat of evaporation
                 double deltaH = 1.0;
                 modelPhaseChangeContainer().modelPhaseChangeEqn(iPhaseChangeModel)->phaseChangeDeltaH(udata[h],deltaH);
                 dudata[h] += ( pDot + (udata[h]-tOld) * pDotJac )
                            * -1 * deltaH
                            / (c_p_eff_*rho_eff_);
            }
            else
                continue;

            if(eqnType_==SPECIES)
            {
                if(updatePhaseFraction ) //assume concentration is constant, and only vary the phase fraction
                {
                    if(control().simulationState().timeStep()>1) //cannot update phasefraction at first time step
                    {
                        double cLocal = particleData().retrieveIntraData(particleDataID_, particleID, h-1);
                        dudata[h] = (pDot + udata[h]*pDotJac)
                                  / cLocal;
                    }
                    else
                        dudata[h] = 0.0;
                }

                else
                {
                    if (inPhase_==SOLID )
                    {
                        if (udata[h] >= 0.0)
                            dudata[h] += pDot + udata[h]*pDotJac;
                        else
                        {
                            dudata[h] = 0.0; //bound to zero in case of negative concentration!
                            udata[h] = 0.0;
                        }
                    }

                    if (inPhase_==LIQUID)
                    {
                        if (udata[h] >= 0.0)
                        {
                            if (tempPhaseFractLiq_[h-1] < phaseFractionMinimum)
                                tempPhaseFractLiq_[h-1] = phaseFractionMinimum;

                            dudata[h] += (pDot + udata[h]*pDotJac ) / tempPhaseFractLiq_[h-1];
                        }
                        else
                        {
                            dudata[h] = 0; //bound to zero!
                            udata[h] = 0.0;
                        }
                    }

                    if (inPhase_ == GAS)
                    {
                        if (udata[h] >= 0)
                        {
                            if (tempPhaseFractGas_[h-1] < phaseFractionMinimum)
                                tempPhaseFractGas_[h-1] = phaseFractionMinimum;

                            dudata[h] += (pDot + udata[h]*pDotJac) / tempPhaseFractGas_[h-1];
                        }
                        else
                        {
                            dudata[h] = 0; //bound to zero!
                            udata[h] = 0.0;
                        }
                    }
                }
            }

                if(verbose_)
                {
                printf("iPhaseChangeModel: %d, phaseId: %d/%d, dataIds: %d/%d, particleDataID_: %d, pDot: %g, pDotJac: %g, tOld: %g, derivativeConcTemp: %g, phaseFractionMinimum: %g. \n",
                       iPhaseChangeModel,
                       (*phaseID)[0],(*phaseID)[1],(*dataIds)[0],(*dataIds)[1],
                       particleDataID_,
                       pDot, pDotJac,
                       tOld, derivativeConcTemp, phaseFractionMinimum
                      );
                }
         }

         //add chemistry
         double b_    = 0.0;
         double bJac_ = 0.0;
         if(modelChemistryContainer().nrChemistryEqns()!=0)
         {
            particleData().returnchemistryDataPoint(   particleDataID_, particleID, h-1, b_);
            particleData().returnchemistryJacDataPoint(particleDataID_, particleID, h-1, bJac_);

            if(eqnType_==SPECIES)
            {
                if (inPhase_ == SOLID)
                {
                    dudata[h] += b_ + bJac_ * udata[h]; //Linearized form of the chemistry source term! b_ is NOT the current source!

                    if (udata[h] < 0)
                    {
                        udata[h] = 0; //bound to zero in case of negative concentration!
                        if(dudata[h]<0)
                            dudata[h] = 0.0;
                    }
                }
                if (inPhase_ == LIQUID)
                {
                    if (tempPhaseFractLiq_[h-1] < phaseFractionMinimum)
                        tempPhaseFractLiq_[h-1] = phaseFractionMinimum;

                    dudata[h] += ( b_ + bJac_ * udata[h] ) / tempPhaseFractLiq_[h-1];

                    if (udata[h] < 0)
                    {
                        udata[h] = 0; //bound to zero in case of negative concentration!
                        if(dudata[h]<0)
                            dudata[h] = 0.0;
                    }
                }
                if (inPhase_ == GAS)
                {
                    if (tempPhaseFractGas_[h-1] < phaseFractionMinimum)
                        tempPhaseFractGas_[h-1] = phaseFractionMinimum;

                    dudata[h] += ( b_ + bJac_ * udata[h] ) / tempPhaseFractGas_[h-1];

                    if (udata[h] < 0)
                    {
                        udata[h] = 0; //bound to zero in case of negative concentration!
                        if(dudata[h]<0)
                            dudata[h] = 0.0;
                    }
                }
            }
            else if (eqnType_==HEAT)
            {
                dudata[h] += ( b_ + bJac_ * udata[h] )
                          / (c_p_eff_*rho_eff_);
            }
         }

    if(verbose_)
    {
        printf("eval final: udata[%d/%d]: %g, dudata: %g, b: %g, bJac: %g, s_chem: %g.\n",
               h, particleMesh().nGridPoints(),
               udata[h],
               dudata[h],
               b_, bJac_,  b_ + bJac_ * udata[h]
               );
    }
    } //End loop grid points
    //printf("\n");
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
    BAND_COL_ELEM(col_j,j,j)   = -6.0*tempDiff_eff_[j]*coeff_2nd_dev_;
    BAND_COL_ELEM(col_j,j+1,j) =  6.0*tempDiff_eff_[j]*coeff_2nd_dev_;

    j=particleMesh().nGridPoints()-1;
    if(BC[1]==DIRICHLET)
    {
        BAND_COL_ELEM(col_j,j,j) = 0.0;
    }
    else if(BC[1]==NEUMANN)
    {
        BAND_COL_ELEM(col_j,j,j)   = -(2.0/(dx_*dx_));
        BAND_COL_ELEM(col_j,j-1,j) =   (2.0/(dx_*dx_));
    }
    else if(BC[1]==CONVECTIVE)
    {
        //printf("jac Convective!\n");
        BAND_COL_ELEM(col_j,j,j) = ((-2.0*tempDiff_eff_[j])/dx_)*(biot_num_/rMAX + biot_num_/dx_ - 1.0/dx_);
        BAND_COL_ELEM(col_j,j-1,j) = (2.0*tempDiff_eff_[j])/(dx_*dx_);
        //printf("BAND_COL_ELEM(col_j,j,j) = %g , BAND_COL_ELEM(col_j,j-1,j) = %g \n", BAND_COL_ELEM(col_j,j,j), BAND_COL_ELEM(col_j,j,j));
   }

    for (h=2; h < particleMesh().nGridPoints(); h++)
    {
        currentRadius_=dx_*(h-1);
        coeff_1st_dev_ = (2.0)/(2.0*currentRadius_*dx_);
        j = h-1;                                //runs from 1 to MX-2
                                                //j (Jacobian index) runs from 0 to MX-1
        col_j = BAND_COL(J,j);                    //BAND_COL BAND_COL_ELEM supposed to be faster - could be replaced by
        BAND_COL_ELEM(col_j,j,j) = -2.0*tempDiff_eff_[j]*(coeff_2nd_dev_);     //diagonal elements of Jac
          BAND_COL_ELEM(col_j,j-1,j) = tempDiff_eff_[j]*(coeff_2nd_dev_-coeff_1st_dev_); //below the diagonal
           BAND_COL_ELEM(col_j,j+1,j) = tempDiff_eff_[j]*(coeff_2nd_dev_+coeff_1st_dev_); //above the diagonal
    }


    //Add chemical term - explicit source term so derivation = 0
    if(modelChemistryContainer().nrChemistryEqns()!=0)
    {
      for(h=1; h <= particleMesh().nGridPoints(); h++)
      {
         double bJac_ = 0.0;
         j = h-1;    //runs from 0 to MX-1
         particleData().returnchemistryJacDataPoint(particleDataID_, particleID, j, bJac_);
         if(eqnType_==SPECIES)
         {
            if (inPhase_==SOLID && udata[h] >= 0)
            {
                col_j                     = BAND_COL(J,j);
                BAND_COL_ELEM(col_j,j,j) += bJac_;
            }

            if (inPhase_==GAS && udata[h] >= 0)
            {
                col_j                     = BAND_COL(J,j);
                if(tempPhaseFractGas_[j] < phaseFractionMinimum)
                    tempPhaseFractGas_[j] = phaseFractionMinimum;

                BAND_COL_ELEM(col_j,j,j) += bJac_  / tempPhaseFractGas_[j];
            }
         }
         else if (eqnType_==HEAT)
         {
            col_j                     = BAND_COL(J,j);
            BAND_COL_ELEM(col_j,j,j) +=  bJac_ / (c_p_eff_*rho_eff_);
         }
      }
    }
}

/////////////////////////////////////////////////////////////////////////////////
void ModelEqn1DSpherical::computeParticleAverages()
{
    //printf("Try to get access to particle data eqn_type %i and phase %i \n",eqnType_,inPhase_);
    //eqnType_==SPECIES && inPhase_==SOLID

    double *** myPhaseFrac = NULL;
    if( inPhase_== GAS || inPhase_== LIQUID)
        myPhaseFrac = ParScaleBase::particleData().accessPhaseFractionMem(inPhase_);
    else if(inPhase_== LIQUID && averagePhaseFraction)
        error().throw_error_one(FLERR,"You want to average the phase fraction for a solid, this is not implemented yet since it is not directy saved in the particleData object.\n");
    else
         if(averagePhaseFraction)
            error().throw_error_one(FLERR,"You want to average the phase fraction, but it is not set for this modelEqn.\n");

    double sum_avg_ = 0.0;        // Variable to sum up all avg intra particle properties
    double avg_ = 0.0;            // Average value for integral variable

    for(int particleID=0; particleID<particleData().nbody(); particleID++)
    {
        particleData().setParticleIDPointer(particleDataID_,particleID);
        particleData().returnIntraData(tempIntraData_);

        tempAvData_ = 0.0;
        double ra(0.0), ri(0.0);
        double volume(0.0), deltaVolume(0.0);

        //integrate over spherical shells, volume and radii are normalized with dx_
        for(int i=0; i<particleMesh().nGridPoints(); i++)
        {
            ra = (double)i + 0.5; ra = min(ra,(double)(particleMesh().nGridPoints()-1));
            ri = (double)i - 0.5; ri = max(ri,0.0);
            deltaVolume  = ra*ra*ra-ri*ri*ri;
            volume      += deltaVolume;
            if(    averagePhaseFraction )
                tempAvData_ += deltaVolume*myPhaseFrac[particleID][0][i];
            else
                tempAvData_ += deltaVolume*tempIntraData_[i];
        }
        tempAvData_ /= volume;
        sum_avg_ +=  tempAvData_;
        particleData().saveIntraParticleAv(particleDataID_, particleID, tempAvData_);
    }
    //current_modelEqn_ID_ = modelEqnID;
    int n_body_ = particleData().nbody();
    avg_ = sum_avg_/n_body_;
    //printf("modelEqnContainer().nrEqns() = %i \n",modelEqnContainer().nrEqns());
    if(writeVolAvgProp)
    {
        //printf("modelEqnContainer().modelEqn(eqnID_)->name() = %s, avg_ = %g, time = %g \n",modelEqnContainer().modelEqn(eqnType_)->name(),avg_,control().simulationState().time());
        output().write_integral_value(name_,"vol_avg",control().simulationState().time(),avg_); 
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
            if (strcmp(coupling().couplingModel().name(),"liggghts") == 0) 
                particleData().returnIntraFlux(particleID,environmentFlux);
            else
                environmentFlux = BCvalue[1];

            tempPartFlux_ = environmentFlux;
        }
        else if(BC[1]==CONVECTIVE)
        {
            if( strcmp(coupling().couplingModel().name(),"liggghts") == 0 )
                 particleData().returnIntraTransCoeff(particleID,environmentTransCoeff);
            else if(transfer_coeff!=NULL)
                 environmentTransCoeff=transfer_coeff->value();
            else
                 environmentTransCoeff=0.0;
            double particleRadius = particleData().pRadius(particleID);

            //This is total flux! Positive if heating
            tempPartFlux_ = particleRadius * particleRadius * 12.5663706144 //4*pi = 12.566
                          * environmentTransCoeff
                          * (
                                  tempIntraData_[particleMesh().nGridPoints()]      //ambient fluid temp/conc
                                - tempIntraData_[particleMesh().nGridPoints()-1]    //particle surface temp/conc
                            );
            //printf("heat flux in ModelEqn1DSpherical::computeSurfaceFluxes() = %g \n", tempPartFlux_);
        }
        particleData().saveIntraParticleFlux(particleDataID_, particleID, tempPartFlux_);
    }
}

/////////////////////////////////////////////////////////////////////////////////
void ModelEqn1DSpherical::updateProperties(int particleID)
{
        //Clean temporary diffusivity / convection array
        for(int i=0; i<(nGridPointsUsed_-1); i++)
        {
            tempDiff_eff_[i] = LARGENUMBER;
            tempConvectionSpeed_eff_[i] = 0;
        }

        //set grid coefficients and BC condition
        rMAX=particleData().pRadius(particleID);

        surface_area = rMAX*rMAX* 12.5663706144;         //4*pi*r²
        particleData().setParticleIDPointer(particleDataID_,particleID);
        particleData().returnIntraData(tempIntraData_);

        if (phaseFraction!=NULL)
        {
            for(int i=0; i<(nGridPointsUsed_-1); i++)
            {
                tempPhaseFractGas_[i] = phaseFraction->value();
                tempPhaseFractLiq_[i] = 0.0;
                tempPhaseFractSolid_[i] = 1.0 - phaseFraction->value();
            }
        }
        else
        {
            particleData().setParticleIDPointerPhaseFraction(particleID);
            particleData().returnPhaseFractionData(tempPhaseFractGas_,tempPhaseFractLiq_);
            for(int i=0; i<(nGridPointsUsed_-1); i++)
                tempPhaseFractSolid_[i] = 1.0 - tempPhaseFractGas_[i] - tempPhaseFractLiq_[i];

            if(verbose_)
             for(int i=0; i<(nGridPointsUsed_-1); i++)
                printf("Phase Fraction(s) [%d] = %g/%g \n", i, tempPhaseFractGas_[i], tempPhaseFractLiq_[i]);

        }


        //pull out all environment properties
        BCvalue[1] = tempIntraData_[particleMesh().nGridPoints()]; //set environment variable, re-use intra-particle container for this purpose
        environmentU=BCvalue[1];

        //Set Neumann BC via flux
        if(BC[1]==NEUMANN)
        {
            if (strcmp(coupling().couplingModel().name(),"liggghts") == 0)// && environmentFlux != 0)
                particleData().returnIntraFlux(particleID,environmentFlux);
            else
                environmentFlux = BCvalue[1];

            if( environmentFlux==0.0 )
                if(verbose_)
                    printf("WARNING: you have not specified a flux. Assuming 0.\n");
        }


        if(eqnType_==HEAT)
        {
            if(particleData().haveGasPhase)
            {
              if (thermal_gas_conductivity==NULL)
                 error().throw_error_one(FLERR,"ERROR: You must specify thermal_gas_conductivity for your heat equation\n");
              if (capacity_gas==NULL)
                 error().throw_error_one(FLERR,"ERROR: You must specify capacity_gas for your heat equation\n");
              if (density_gas==NULL)
                 error().throw_error_one(FLERR,"ERROR: You must specify density_gas for your heat equation\n");
            }

            if(particleData().haveLiquidPhase)
            {
              if (thermal_liquid_conductivity==NULL)
                 error().throw_error_one(FLERR,"ERROR: You must specify thermal_liquid_conductivity for your heat equation\n");
              if (capacity_liquid==NULL)
                 error().throw_error_one(FLERR,"ERROR: You must specify capacity_liquid for your heat equation\n");
              if (density_liquid==NULL)
                 error().throw_error_one(FLERR,"ERROR: You must specify density_liquid for your heat equation\n");
            }

            for(int i=0; i<(nGridPointsUsed_-1); i++)
            {
                lambda_eff_ = thermal_solid_conductivity->value() * tempPhaseFractSolid_[i];
                c_p_eff_    = capacity_solid->value()             * tempPhaseFractSolid_[i];
                rho_eff_    = density_solid->value()              * tempPhaseFractSolid_[i];

                if(particleData().haveGasPhase)
                {
                    lambda_eff_ += thermal_gas_conductivity->value() * tempPhaseFractGas_[i];
                    c_p_eff_    += capacity_gas->value()             * tempPhaseFractGas_[i];
                    rho_eff_    += density_gas->value()              * tempPhaseFractGas_[i];
                }

                if(particleData().haveLiquidPhase)
                {
                    lambda_eff_ += thermal_liquid_conductivity->value() * tempPhaseFractLiq_[i];
                    c_p_eff_    += capacity_liquid->value()             * tempPhaseFractLiq_[i];
                    rho_eff_    += capacity_liquid->value()             * tempPhaseFractLiq_[i];
                }
                tempDiff_eff_[i]=lambda_eff_/(c_p_eff_*rho_eff_);

                if(solveConvectiveFlux) //save for later usage
                    tempConvectionSpeed_eff_[i] = 1.0
                                                / (c_p_eff_*rho_eff_);
            }
        }

        if(eqnType_==SPECIES)
        {
            if (diffusivity==NULL)
                    error().throw_error_one(FLERR,"ERROR: You must specify a diffusivity for every species equation! \n");
            diffusivity_= diffusivity->value();

            if( (inPhase_==GAS || inPhase_==LIQUID) && tortuosity==NULL)
                error().throw_error_one(FLERR,"ERROR: You must specify a species tortuosity in porous particles when looking at a fluid phase! \n");

            if (inPhase_ == GAS)
            {
                tortuosity_ =  tortuosity->value();

                for(int i=0; i<(nGridPointsUsed_-1); i++)
                    tempDiff_eff_[i]=diffusivity_/tortuosity_;

                if(pore_radius!=NULL && molar_mass!=NULL)
                {
                    if(!IsoThermal_)
                    {
                         particleData().setParticleIDPointer(heatEqnID_,particleID);
                         particleData().returnIntraData(tempTempData_);
                    }
                    else
                         for(int i=0; i<(nGridPointsUsed_-1); i++)
                             tempTempData_[i] = isoTemp_;

                    //Calculate Knudsen Diffusion Coeff in Medium, store in temporary array
                    for(int i=0; i<(nGridPointsUsed_-1); i++)
                    {
                        if (tempPhaseFractGas_[i] < phaseFractionMinimum)
                            tempPhaseFractGas_[i] = phaseFractionMinimum;

                           double dKnudsen = 0.66666666666666667*pore_radius->value()
                                           * sqrt((8.0/3.14159654)
                                                  *(8.3144621*tempTempData_[i]/molar_mass->value())
                                                 )
                                           *  tempPhaseFractGas_[i] / tortuosity->value();

                            tempDiff_eff_[i] = 1.0
                                             / ( 1.0/dKnudsen + 1.0/tempDiff_eff_[i] );

                        }
                  }//end for Knudsen diffusion

                  //Reset setParticleIDPointer to original position
                  particleData().setParticleIDPointer(particleDataID_,particleID);
            }
            else if (inPhase_ == LIQUID)
            {
                tortuosity_ =  tortuosity->value();
                for(int i=0; i<(nGridPointsUsed_-1); i++)
                    tempDiff_eff_[i]= diffusivity->value() / tortuosity_;
            }
            else //must be SOLID
            {
                for(int i=0; i<(nGridPointsUsed_-1); i++)
                    tempDiff_eff_[i]=diffusivity->value()*tempPhaseFractSolid_[i];
            }
        }

        //Copy molecular diffusion in array, Knudsen diffusion is 0 (array is empty)
        //if(pore_radius->value()==NULL && molar_mass->value()==NULL && eqnType_ != SPECIES && inPhase_ != GAS)
        tempDiff_eff_[nGridPointsUsed_-1] = 0.0; //last point --> outside, so set zero
        for(int i=0; i<=(nGridPointsUsed_-1); i++)
        {
            if(verbose_)
                printf("Temp Diffu [%i] = %g \n", i,tempDiff_eff_[i]);

        }

        //set boundary properties
        dx_ = rMAX/((particleMesh().nGridPoints())-1.0);

        //coeffcient results off 2nd deviation O(dx_^2)
        coeff_2nd_dev_ = (1.0)/(dx_*dx_);

        if(BC[1]==CONVECTIVE)
        {
            if( strcmp(coupling().couplingModel().name(),"liggghts") == 0 )
            {
                particleData().returnIntraTransCoeff(particleID,environmentTransCoeff);
            }
            else if(transfer_coeff!=NULL)
            {
                environmentTransCoeff=transfer_coeff->value();
                //printf("environmentTransCoeff = %g \n", environmentTransCoeff);
            }
            else
            {
                printf("WARNING: you have not specified a transfer coefficient. Assuming 0.\n");
                environmentTransCoeff=0.;
            }

            if(environmentTransCoeff<0. || environmentTransCoeff>1.e20)
            {
                printf("environmentTransCoeff: %g. \n", environmentTransCoeff);
                error().throw_error_one(FLERR,"The transfer coefficient you specified is negative or very large. this is a problem. \n");
            }

            //Biot number
            if(eqnType_==HEAT)
            {
                biot_num_ = (environmentTransCoeff*dx_)/lambda_eff_;
            }
            else if(eqnType_==SPECIES)
            {
                  //must multiply here with phaseFraction, since diffu_eff_ has a different meaning for gas!
               if(inPhase_ == GAS)
                  biot_num_ = (environmentTransCoeff*dx_)
                           / (
                               ( tempDiff_eff_[nGridPointsUsed_-2] + VSMALLNUMBER )
                              *  tempPhaseFractGas_[nGridPointsUsed_-2]
                             );
               else
                  biot_num_ = environmentTransCoeff*dx_
                           / ( tempDiff_eff_[nGridPointsUsed_-2] + VSMALLNUMBER );
            }
        }


        // * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        //Set convective fluxes if needed
    if(solveConvectiveFlux)
    {
        double *** convectiveFluxGas = NULL;
        if(particleData().haveGasPhase)
            convectiveFluxGas = ParScaleBase::particleData().accessIntraConvectiveFlux(GAS);

        double *** convectiveFluxLiquid = NULL;
        if(particleData().haveLiquidPhase)
            convectiveFluxLiquid = ParScaleBase::particleData().accessIntraConvectiveFlux(LIQUID);

        if(eqnType_==HEAT)
        {
            for(int i=0; i<(nGridPointsUsed_-1); i++)
            {
                double totalFluxTimesCp = 0.0;
                if(particleData().haveGasPhase)
                    totalFluxTimesCp  += capacity_gas->value()    * convectiveFluxGas   [particleID][0][i];

                if(particleData().haveLiquidPhase)
                    totalFluxTimesCp  += capacity_liquid->value() * convectiveFluxLiquid[particleID][0][i];

                tempConvectionSpeed_eff_[i] *= totalFluxTimesCp; //must multiply here! c_p_eff_*rho_eff_ already in denominator

                 if(verbose_)
                    printf("tempConvectionSpeed_eff_[%i] = %g, totalFluxTimesCp: %g. \n",
                           i, tempConvectionSpeed_eff_[i], totalFluxTimesCp);
            }
        }
        if(eqnType_==SPECIES && inPhase_==GAS)
        {
            //obtain the local temperature
            if(!IsoThermal_)
            {
                 particleData().setParticleIDPointer(heatEqnID_,particleID);
                 particleData().returnIntraData(tempTempData_);
            }
            else
                 for(int i=0; i<(nGridPointsUsed_-1); i++)
                      tempTempData_[i] = isoTemp_;

            for(int i=0; i<(nGridPointsUsed_-1); i++)
            {
                double cTotalGas = particleData().referencePressure
                                 * INVERSE_UNIVERSAL_GAS_CONSTANT
                                 / tempTempData_[i];
                if (tempPhaseFractGas_[i] < phaseFractionMinimum)
                    tempPhaseFractGas_[i] = phaseFractionMinimum;

                tempConvectionSpeed_eff_[i] = convectiveFluxGas   [particleID][0][i]
                                            / ( cTotalGas );

                 if(verbose_)
                    printf("tempConvectionSpeed_eff_[%i] = %g, convectiveFluxGas: %g, cTotalGas: %g. \n",
                           i, tempConvectionSpeed_eff_[i], convectiveFluxGas[particleID][0][i], cTotalGas);
            }
        }
        if(eqnType_==SPECIES && inPhase_==LIQUID)
        {
            if(particleData().eqnIdFirstLiquid<0)
                error().throw_error_one(FLERR,"ERROR: particleData().eqnIdFirstLiquid<0! This is a problem! \n");

            double *** solventConcentration = particleData().accessIntraPartMem
                                              (
                                                particleData().eqnIdFirstLiquid
                                              );
            for(int i=0; i<(nGridPointsUsed_-1); i++)
            {
                if (tempPhaseFractGas_[i] < phaseFractionMinimum)
                    tempPhaseFractGas_[i] = phaseFractionMinimum;

                tempConvectionSpeed_eff_[i] = convectiveFluxLiquid[particleID][0][i]
                                            / ( solventConcentration[particleID][0][i]);
                if(verbose_)
                    printf("tempConvectionSpeed_eff_[%i] = %g, solventConcentration: %g. \n",
                           i, tempConvectionSpeed_eff_[i], solventConcentration[particleID][0][i]);
            }
        }

        tempConvectionSpeed_eff_[nGridPointsUsed_-1] = 0.0; //last point --> outside, so set zero
    } //end if(solveConvectiveFlux)
    return;
}


/////////////////////////////////////////////////////////////////////////////////
void ModelEqn1DSpherical::phaseFluxIntegrate()  const
//                                            bool haveGasPhase, bool haveLiquidPhase,
//                                             vector<double*> datapointerPhaseFraction,
//                                             double* totalSource, double* totalFlux, double &rParticle
//                                            ) const
{
    if(!solveConvectiveFlux) return;

    if(ParScaleBase::modelPhaseChangeContainer().nrPhaseChangeEqns()<1) return;
        //Phase flux model based on integrating the total mole balance equation
        //assumes uniform pressure distribution, and that the change of the phase fraction/concentration
        //is small

        double *** phaseFracMyPhase   = ParScaleBase::particleData().accessPhaseFractionMem(inPhase_);
        double *** phaseChangePDot    = ParScaleBase::particleData().accessPhaseChangeRate(inPhase_);
        double *** phaseChangePDotJac = ParScaleBase::particleData().accessPhaseChangeRateJac(inPhase_);
        double *** intraPartMem       = ParScaleBase::particleData().accessIntraPartMem(particleDataID_);
        double *** convectiveFlux     = ParScaleBase::particleData().accessIntraConvectiveFlux(inPhase_);
        double *** phaseChangeVol     = ParScaleBase::particleData().accessPhaseChangeRateVolumetric(inPhase_);

        double *** chemistrySource    = NULL;
        double *** chemistrySourceJac = NULL;

        if(ParScaleBase::modelChemistryContainer().nrChemistryEqns()!=0)
        {
            chemistrySource    = ParScaleBase::particleData().accessChemistryMem(particleDataID_);
            chemistrySourceJac = ParScaleBase::particleData().accessChemistryMemJac(particleDataID_);
        }

        for (int i=0; i<ParScaleBase::particleData().nbody(); i++)  //particle loop
        {
            double IntegralTerm = 0.0;
            double deltaR       = ParScaleBase::particleData().pRadiusConst(i)/(nGridPointsUsed_-2);

            //Do calculation @ the center of the particle
            double totalSource       = 0.0;
            convectiveFlux[i][0][0]    = 0;

            for(int jGrid=1; jGrid<(nGridPointsUsed_-1); jGrid++) //start with 2nd grid point
            {
                double radiusJPlus05   = ((double)jGrid+0.5)*deltaR;
                double radiusJMinus05  = radiusJPlus05-deltaR;

                totalSource = phaseChangePDot[i][0][jGrid]
                            + phaseChangePDotJac[i][0][jGrid]
                            * intraPartMem[i][0][jGrid] ;
                totalSource /= max(1e-16,phaseFracMyPhase[i][0][jGrid]); //Re-scale phase change term: now specific to phase volume

                double cTotalGas = ParScaleBase::particleData().referencePressure
                                 * INVERSE_UNIVERSAL_GAS_CONSTANT
                                 / tempTempData_[jGrid];

#if APPLYCONVECTIVEFLUX
                //scale source in case computed in highly saturated region (e.g., a dense wet core)
                //(no convection modeled in this region for stability reasons)
                if(phaseFracMyPhase[i][0][jGrid]<phaseFractionMinimumConvection)
                    totalSource /= max(convectionBoundMinStefanFactor,
                                       1.0 - intraPartMem[i][0][jGrid]/cTotalGas); //avoid division by small number
#endif
                if(ParScaleBase::modelChemistryContainer().nrChemistryEqns()!=0)
                {
                    totalSource += chemistrySource   [i][0][jGrid]
                                 + chemistrySourceJac[i][0][jGrid]
                                 * intraPartMem      [i][0][jGrid];
                }

                //Correct change due to generation of new phase volume
                totalSource -= phaseChangeVol[i][0][jGrid] * cTotalGas; //this is accurate as long as ddt(c_g) is small

                IntegralTerm += 0.3333333333333333 // 1/3; do not multiply with 4*PI, since we divide later only by r^2
                              *(
                                   radiusJPlus05 *radiusJPlus05 *radiusJPlus05
                                 - radiusJMinus05*radiusJMinus05*radiusJMinus05
                               )
                              * totalSource ;

                //Re-scale phase change term: now superficial flux @ upper face
                convectiveFlux[i][0][jGrid] = IntegralTerm / radiusJPlus05 / radiusJPlus05
                                            * (phaseFracMyPhase[i][0][jGrid]+phaseFracMyPhase[i][0][jGrid+1])*0.5;
                                            //Re-scale phase change term: now superficial flux @ upper face

                if(verbose_)
                    printf("phaseFluxIntegrate: jGrid=%d, deltaR: %g radiusJMinus05: %g, IntegralTerm: %g, conc: %g, source/flux: %g/%g. \n",
                        jGrid, deltaR, radiusJMinus05,
                        IntegralTerm*12.56637, //multiply with 4*Pi to display correct value
                        intraPartMem[i][0][jGrid],
                        totalSource,
                        convectiveFlux[i][0][jGrid]);
            }

        } // end particle loop

    return;
}

/////////////////////////////////////////////////////////////////////////////////
void ModelEqn1DSpherical::phaseFluxCapillarity()   const
{
    //Phase flux model based on a capillary pressure distribution
    //only useful for liquid phase. Currently, assumes uniform pressure distribution
    //in the gas phase

    if(solveConvectiveFlux)
    {
        if( !ParScaleBase::particleData().haveGasPhase || !ParScaleBase::particleData().haveLiquidPhase)
            ParScaleBase::error().throw_error_one(FLERR,"phaseFluxCapillarity could not find gas and liquid phase. \n");

        if (permeability==NULL)
            ParScaleBase::error().throw_error_one(FLERR,"ERROR: You must specify permeability for your liquid-phase equation\n");
        if (viscosity==NULL)
            ParScaleBase::error().throw_error_one(FLERR,"ERROR: You must specify viscosity for your liquid-phase equation\n");
        if (viscosity->value() < SMALLNUMBER)
            ParScaleBase::error().throw_error_one(FLERR,"ERROR: You must specify a viscosity that is larger than SMALLNUMBER.\n");
        if (surface_tension==NULL)
            ParScaleBase::error().throw_error_one(FLERR,"ERROR: You must specify surface_tension for your liquid-phase equation\n");
        if (surface_tension->parameters().size()<3)
            ParScaleBase::error().throw_error_one(FLERR,"ERROR: You must specify at least 3 parameters in the surface_tension model for your liquid-phase equation.\n");
        if (film_flow==NULL)
            ParScaleBase::error().throw_error_one(FLERR,"ERROR: You must specify film_flow for your liquid-phase equation\n");
        if (film_flow->parameters().size()<2)
            ParScaleBase::error().throw_error_one(FLERR,"ERROR: You must specify at least 2 parameters in the film_flow model for your liquid-phase equation.\n");

        double deltaR       = 0.0;
        double *** phaseFracGas    = ParScaleBase::particleData().accessPhaseFractionMem(GAS);
        double *** phaseFracLiquid = ParScaleBase::particleData().accessPhaseFractionMem(LIQUID);
        double *** convectiveFlux  = ParScaleBase::particleData().accessIntraConvectiveFlux(inPhase_);
        double *** intraPartMem    = ParScaleBase::particleData().accessIntraPartMem(particleDataID_);

        for (int i=0; i<ParScaleBase::particleData().nbody(); i++)  //particle loop
        {
            deltaR                = ParScaleBase::particleData().pRadiusConst(i)/(nGridPointsUsed_-2);

            phaseFracLiquid[i][0][0] = phaseFracLiquid[i][0][1]; //set liquid phase fraction at center
            convectiveFlux[i][0][0]  = 0;
            for(int jGrid=1; jGrid<(nGridPointsUsed_-1); jGrid++) //start with 2nd grid point, and end with last
            {
                //We use centra differences here
                double currDeltaR(1e16);

                double liqPhaseFrac   = phaseFracLiquid      [i][0][jGrid-1];
                double pressureMinus1 = surface_tension->parameters()[0]
                                  * surface_tension->value()
                                  * pow( liqPhaseFrac+surface_tension->parameters()[2], surface_tension->parameters()[1] );

                if(jGrid==(nGridPointsUsed_-2)) //last grid point
                {
                    liqPhaseFrac          = phaseFracLiquid[i][0][jGrid];
                    currDeltaR            = deltaR;
                }
                else
                {
                    liqPhaseFrac          = phaseFracLiquid[i][0][jGrid+1];
                    currDeltaR            = deltaR*2.;
                }
                double pressurePlus1  = surface_tension->parameters()[0]
                                      * surface_tension->value()
                                      * pow( liqPhaseFrac+surface_tension->parameters()[2], surface_tension->parameters()[1] ) ;

                //Compute saturation and Keff at current grid point
                double saturation    =  phaseFracLiquid[i][0][jGrid]
                                     / ( phaseFracLiquid[i][0][jGrid]
                                        +phaseFracGas[i][0][jGrid]
                                        +SMALLNUMBER);
                double Keff          = saturation * saturation * saturation;

                double alphaFilm = 1.0;
                if(saturation<film_flow->parameters()[1])
                    alphaFilm = (saturation                 - film_flow->parameters()[0])
                              / (film_flow->parameters()[1] - film_flow->parameters()[0]);
                if(saturation<film_flow->parameters()[0])
                    alphaFilm = 0.0;

                convectiveFlux[i][0][jGrid] = alphaFilm
                                            * intraPartMem[i][0][jGrid]
                                            * permeability->value()
                                            * Keff
                                            / viscosity->value()
                                            * (pressureMinus1-pressurePlus1) //negative gradient!
                                            / currDeltaR ;

                if(verbose_)
                  printf("phaseFluxCapillarity: iPArticle=%d, jGrid=%d, gasPhaseFrac: %g, liquidPhaseFrac: %g, pressure-1/+1: %g/%g, Keff: %g, alphaFilm: %g, convectiveFlux: %g. \n",
                        i, jGrid,
                        phaseFracGas[i][0][jGrid], phaseFracLiquid[i][0][jGrid],
                        pressureMinus1, pressurePlus1, Keff,alphaFilm,
                        convectiveFlux[i][0][jGrid]);

            }

        } // end particle loop
    }//end if(solveConvectiveFlux)

    return;
}
