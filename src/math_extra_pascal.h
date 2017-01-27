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

     Copyright (C): 2012 - 2014 DCS Computing GmbH (www.dcs-computing.com), Linz, Austria
                    2012 - 2014 Department of Particulate Flow Modelling, JKU Linz
                              (www.jku.at/pfm), Linz, Austria

   This file was originally part of LIGGGHTS (www.cfdem.com), and is now re-distributed
   under LGPL as part of ParScale with the permission of the copyright holders
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



#ifndef PSC_MATH_EXTRA_PASCAL_H
#define PSC_MATH_EXTRA_PASCAL_H

#include "math.h"
#include "stdio.h"
#include "string.h"
#include "ctype.h"
#include "qjson_includes.h"

namespace MathExtraPascal {

#define SMALLNUMBER 1e-32

  inline int min(int a,int b);
  inline int max(int a,int b);
  inline int abs(int a);
  inline double min(double a,double b);
  inline double max(double a,double b);
  inline double min(double a,double b,double c);
  inline double max(double a,double b,double c);
  inline double min(double a,double b,double c,double d);
  inline double min(double *input, int n,int &which);
  inline double max(double *input, int n,int &which);
  inline double abs(double a);

  inline void vec_quat_rotate(double *vec, double *quat, double *result);
  inline void vec_quat_rotate(double *vec, double *quat);
  inline void vec_quat_rotate(int *vec, double *quat) { UNUSED(vec); UNUSED(quat); }
  inline void vec_quat_rotate(bool *vec, double *quat) { UNUSED(vec); UNUSED(quat); }

  inline void quatquat(double *a, double *b, double *c);
  inline void quat_from_vec(const double *v, double *q);
  inline void vec_from_quat(const double *q, double *v);
  inline void qconjugate(double *q, double *qc);

  inline bool compDouble(double const a, double const b, double const prec);
  inline double fastPow(double base, double exponent);
  inline double interpolateLinearly(QJsonArray, QJsonArray, double);
};

/* ----------------------------------------------------------------------
   min max stuff
------------------------------------------------------------------------- */

  int MathExtraPascal::min(int a,int b) { if (a<b) return a; return b;}
  int MathExtraPascal::max(int a,int b) { if (a>b) return a; return b;}

  double MathExtraPascal::min(double a,double b) { if (a<b) return a; return b;}
  double MathExtraPascal::max(double a,double b) { if (a>b) return a; return b;}

  double MathExtraPascal::min(double a,double b,double c)
  {
      double ab = MathExtraPascal::min(a,b);
      if (ab<c) return ab;
      return c;
  }
  double MathExtraPascal::max(double a,double b,double c)
  {
      double ab = MathExtraPascal::max(a,b);
      if (ab<c) return c;
      return ab;
  }

  double MathExtraPascal::min(double a,double b,double c,double d)
  {
      double ab = MathExtraPascal::min(a,b);
      double cd = MathExtraPascal::min(c,d);
      if (ab<cd) return ab;
      return cd;
  }

  double MathExtraPascal::min(double *input, int n,int &which)
  {
      double min = input[0];
      which = 0;

      for(int i = 1; i < n; i++)
      {
          if(input[i] < min)
          {
              which = i;
              min = input[i];
          }
      }
      return min;
  }
  double MathExtraPascal::max(double *input, int n,int &which)
  {
      double max = input[0];
      which = 0;

      for(int i = 1; i < n; i++)
      {
          if(input[i] > max)
          {
              which = i;
              max = input[i];
          }
      }
      return max;
  }

  int MathExtraPascal::abs(int a) { if (a>0) return a; return -a;}
  double MathExtraPascal::abs(double a) { if (a>0.) return a; return -a;}

/*----------------------------------------------------------------------
   rotoate vector by quaternion
------------------------------------------------------------------------- */

void MathExtraPascal::vec_quat_rotate(double *vec, double *quat, double *result)
{
    double vecQ[4], resultQ[4], quatC[4], temp[4];

    // construct quaternion (0,vec)
    quat_from_vec(vec,vecQ);

    // conjugate initial quaternion
    qconjugate(quat,quatC);

    // rotate by quaternion multiplications
    quatquat(quat,vecQ,temp);
    quatquat(temp,quatC,resultQ);

    // return result
    vec_from_quat(resultQ,result);
}

/*----------------------------------------------------------------------
   rotoate vector by quaternion
------------------------------------------------------------------------- */

void MathExtraPascal::vec_quat_rotate(double *vec, double *quat)
{
    double vecQ[4], resultQ[4], quatC[4], temp[4], result[3];

    // construct quaternion (0,vec)
    quat_from_vec(vec,vecQ);

    // conjugate initial quaternion
    qconjugate(quat,quatC);

    // rotate by quaternion multiplications
    quatquat(quat,vecQ,temp);
    quatquat(temp,quatC,resultQ);

    // return result
    vec_from_quat(resultQ,result);
    vec[0] = result[0];
    vec[1] = result[1];
    vec[2] = result[2];
}

/* ----------------------------------------------------------------------
   quaternion-quaternion multiply: c = a*b
------------------------------------------------------------------------- */

void MathExtraPascal::quatquat(double *a, double *b, double *c)
{
  c[0] = a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3];
  c[1] = a[0]*b[1] + b[0]*a[1] + a[2]*b[3] - a[3]*b[2];
  c[2] = a[0]*b[2] + b[0]*a[2] + a[3]*b[1] - a[1]*b[3];
  c[3] = a[0]*b[3] + b[0]*a[3] + a[1]*b[2] - a[2]*b[1];
}

/* ----------------------------------------------------------------------
   construct quaternion4 from vector3
------------------------------------------------------------------------- */

void MathExtraPascal::quat_from_vec(const double *v, double *q)
{
  q[0] = 0.;
  q[1] = v[0];
  q[2] = v[1];
  q[3] = v[2];
}

/* ----------------------------------------------------------------------
   construct vector3 from quaternion4
------------------------------------------------------------------------- */

void MathExtraPascal::vec_from_quat(const double *q, double *v)
{
  v[0] = q[1];
  v[1] = q[2];
  v[2] = q[3];
}

/* ----------------------------------------------------------------------
   conjugate of a quaternion: qc = conjugate of q
   assume q is of unit length
------------------------------------------------------------------------- */

void MathExtraPascal::qconjugate(double *q, double *qc)
{
  qc[0] = q[0];
  qc[1] = -q[1];
  qc[2] = -q[2];
  qc[3] = -q[3];
}

/* -----------------------------------------------------------------------------
 * compare two doubles by using their integer representation
 * source: http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm
 -------------------------------------------------------------------------------*/

bool MathExtraPascal::compDouble(double const a, double const b, double const prec)
{
  if (a == b)
    return true;

  /*NP non-tested alternative
  double x = a-b;
  if(MathExtraLiggghts::abs(x) < prec)
    return true;
  return false;
  */
  if (b == 0)
    return a < prec && a > -prec;

  double x = (a-b);//b;

  return x < prec && x > -prec;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
inline double MathExtraPascal::fastPow(double base, double exponent)
{
    double result = 1.0;  int i_ = 0;
    if (floor(exponent) == exponent)
    {   while(i_ < floor(exponent))
        {
                result *= base;
                i_++;
        }
    }
    else
           result = pow(base, exponent);

    return result;

}
/* -----------------------------------------------------------------------------
 * Operations
 -------------------------------------------------------------------------------*/
inline double MathExtraPascal::interpolateLinearly(QJsonArray x, QJsonArray y, double xToInterpolate)
{
    //above last value, take the last value
    double result = y[x.count()-1].toDouble();

    for(int i= ( x.count()-1 ); 0 <= i; i--) //start from from last value
    {
        if( xToInterpolate < x[i].toDouble()  )
        {
            //below first value, take the first value
            if(i==0)
                result = y[0].toDouble();

            //otherwise, interpolate
            else
            {
                double deltaX = x[i].toDouble() - x[i-1].toDouble();
                if(deltaX<SMALLNUMBER) //avoid division by zero and extrapolation
                    result = y[i-1].toDouble();
                else
                    result = y[i-1].toDouble()
                          +  ( xToInterpolate - x[i-1].toDouble() )
                            /( deltaX )
                            *( y[i].toDouble() - y[i-1].toDouble() );

            }
        }
    }

    return result;
}



#endif
