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

ModelEqnStyle(ShrinkingCore, ModelEqnShrinkingCore)

#else

#ifndef PASC_MODEL_SHRINKINGCORE_H
#define PASC_MODEL_SHRINKINGCORE_H

#include "model_base.h"
#include "model_eqn.h"
#include "model_container.h"
#include "error.h"
#include "pascal_base_accessible.h"


namespace PASCAL_NS
{

class ModelEqnShrinkingCore : public ModelEqn
{
    public:

      ModelEqnShrinkingCore(ParScale *ptr, char *name);

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

    private:
        bool debug_;

        int h;                                         //index of spatial position 1 (reaction front) ... 2 (fluid)
        int j;                                         //Jacobian matrix index 0...MX-1
        realtype *col_j;                            //jth collum of jacobian matrix

        double diffu_eff_;                          //effective diffusivity

        double lambda_solid_;                        //thermoconductivity solid,gas,effective
        double lambda_gas_;
        double lambda_eff_;

        double c_p_solid_;                           //heat capacity solid,gas,effective
        double c_p_gas_;
        double c_p_eff_;

        double rho_solid_;                           //density solid,gas,effective
        double rho_gas_;
        double rho_eff_;

        double kSurface_;

};

} //end PASCAL_NS

#endif

#endif
