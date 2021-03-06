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


#ifdef MODEL_EQN_CLASS

ModelEqnStyle(1DSpherical, ModelEqn1DSpherical)

#else

#ifndef PASC_MODEL_1D_SPHERICAL_H
#define PASC_MODEL_1D_SPHERICAL_H

#include "model_base.h"
#include "model_eqn.h"
#include "model_container.h"
#include "error.h"
#include "pascal_base_accessible.h"
#include "coupling.h"
#include "particle_data.h"


namespace PASCAL_NS
{

class ModelEqn1DSpherical : public ModelEqn
{
    public:

      ModelEqn1DSpherical(ParScale *ptr, char *name);
      ~ModelEqn1DSpherical();

      void init(int narg, char const* const* arg, int eqnType, int modelEqnID);

      virtual void begin_of_step();

      virtual void eval(double t, double* udata, double* dudata, double* p);

      virtual void returnJac(long int N, long int mu, long int ml,
                   double t, double* udata, double* fudata,
                   DlsMat J, double* p,
                   double* tmp1data, double* tmp2data, double* tmp3data);

      void updateProperties(int particleID);
      void computeParticleAverages();
      void computeSurfaceFluxes();

      void  evaluatePhaseFlux() const
        {
            if(inPhase_==GAS)
                (this->*phaseFluxGas)();
            else if(inPhase_==LIQUID)
                (this->*phaseFluxLiquid)();
        }

        //TODO: set coeff correct related to spherical
        double dx_;                                    //dx: distance between grid points
        double coeff_2nd_dev_, coeff_1st_dev_;        //coefficients of first end second derivatives
        double currentRadius_;                        //actual radial position depending on h,MX
        int h;                                         //index of spatial position 1 ... MX
        int j;                                         //Jacobian matrix index 0...MX-1
        realtype *col_j;                            //jth collum of jacobian matrix

        double biot_num_;                            //Biot Number

        double lambda_eff_;                         //Effective conductivity

        double c_p_eff_;                            //Effective capacity

        double rho_eff_;                            //Effective density    

        double tortuosity_;                         //Tortousity
        double diffusivity_;                        //Binary Diffusivity
        double pore_radius_;                        //constant pore radius seen by species
        double molar_mass_;                         //constant molar_mass of spieces
        double porosity_;

        double * tempDiff_eff_;
        double * tempConvectionSpeed_eff_;          //for species: convective flux / total concentration / eps_phase
                                                    //for heat: Sum(heatCapacity*flux) / (effective density * effective heat capacity)
        double * tempTempData_;
        double * tempPhaseFractGas_;
        double * tempPhaseFractLiq_;
        double * tempPhaseFractSolid_;

    private:
        mutable int current_modelEqn_ID_;
        bool debug_;
        bool boundZero_; //TODO: read from json in order to swith on/off bounding

        //Ptrs to top-level phase flux model for gas and liquid phase
        void (ModelEqn1DSpherical::*phaseFluxGas)() const;
        void (ModelEqn1DSpherical::*phaseFluxLiquid)() const;

        //Individual phase flux models
        void  phaseFluxIntegrate() const;
        void  phaseFluxCapillarity() const;

};

} //end PASCAL_NS

#endif

#endif
