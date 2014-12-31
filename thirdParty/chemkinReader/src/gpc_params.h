#ifndef GPC_PARAMS_H
#define GPC_PARAMS_H

namespace Sprog
{
    const double PI = 3.1415926535897932384626433832795;
    const double kB	= 1.3806504e-23; // Boltzmann constant m^2 kg s^-2 K^-1
    const double EPSILON0 = 8.854187816e-12;
    const double Debye__ = 3.33564e-30;
    const double Angstroem__ = 1.0e-10;
}

namespace fastMath
{
    //! Faster pow functions for integer powers.
    inline double pow2(const double& x) {return x*x;}
    inline double pow3(const double& x) {return x*pow2(x);}
    inline double pow4(const double& x) {return pow2(x)*pow2(x);}
}

#endif
