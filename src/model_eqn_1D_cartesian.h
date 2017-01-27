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

ModelEqnStyle(1DCartesian, ModelEqn1DCartesian)

#else

#ifndef PASC_MODEL_1D_CARTESIAN_H
#define PASC_MODEL_1D_CARTESIAN_H

#include "model_eqn.h"

namespace PASCAL_NS
{

class ModelEqn1DCartesian : public ModelEqn
{
    public:

      ModelEqn1DCartesian(ParScale *ptr, char *name);

      void init(int narg, char const* const* arg, int eqnType, int modelEqnID);

      virtual void begin_of_step();

      virtual void eval(double t, double* udata, double* dudata, double* p);

      void updateProperties();
      void computeParticleAverages();
      void computeSurfaceFluxes();

        double dx;                                    //dx: distance between grid points
        double coeff_2nd_dev, coeff_1st_dev;        //coefficients of first end second derivatives
        double x_coeff_1st_dev;                        //actual radial position depending on h,MX
        int h;                                         //index of spatial position 1 ... MX
        int j;                                         //Jacobian matrix index 0...MX-1
        realtype *col_j;                            //jth collum of jacobian matrix

        double biot_num;                            //Biot Number
        double diffu_eff_;                          //effective diffusivity

        double lambda_solid;                        //thermoconductivity solid,gas,effective
        double lambda_gas;
        double lambda_eff;

        double c_p_solid;                           //heat capacity solid,gas,effective
        double c_p_gas;
        double c_p_eff;

        double rho_solid;                           //density solid,gas,effective
        double rho_gas;
        double rho_eff;

        const char* ptr_name;
        double * ptr_value;
        double steel_constant_;

    private:
        bool debug_;

};

} //end PASCAL_NS

#endif

#endif
