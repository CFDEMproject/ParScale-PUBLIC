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

#ifndef LMP_VECTOR_PASCAL_H
#define LMP_VECTOR_PASCAL_H

#include "math.h"

namespace PASCAL_NS {

//================================================
//SOME VERY SIMPLE VECTOR OPERATIONS
//================================================

inline void vectorConstruct3D(double *v,double x, double y, double z)
{
  v[0] = x;
  v[1] = y;
  v[2] = z;
}

inline void vectorConstruct3D(int *v,int x, int y, int z)
{
  v[0] = x;
  v[1] = y;
  v[2] = z;
}

inline void vectorNormalize3D(double *v)
{
    double norm = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    double invnorm = (norm == 0.) ? 0. : 1./norm;
    v[0] *= invnorm;
    v[1] *= invnorm;
    v[2] *= invnorm;
}

inline double vectorMag3D(const double *v)
{
  return (  sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])  );
}

inline double vectorMag3DSquared(const double *v)
{
  return (  v[0]*v[0]+v[1]*v[1]+v[2]*v[2]  );
}

inline double vectorMag4D(const double *v)
{
  return (  sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]+v[3]*v[3])  );
}

inline double pointDistance(const double *point1, const double *point2)
{
  return
  (
     sqrt
     (
          (point1[0]-point2[0]) * (point1[0]-point2[0]) +
          (point1[1]-point2[1]) * (point1[1]-point2[1]) +
          (point1[2]-point2[2]) * (point1[2]-point2[2])
     )
  );
}

inline double vectorMag4DSquared(const double *v)
{
  return (  v[0]*v[0]+v[1]*v[1]+v[2]*v[2]+v[3]*v[3]  );
}

inline double vectorDot3D(const double *v1, const double *v2)
{
  return (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]);
}

inline double vectorDot2D(const double *v1, const double *v2)
{
  return (v1[0]*v2[0]+v1[1]*v2[1]);
}

inline void vectorCopy2D(const double *from, double *to)
{
  to[0]=from[0];
  to[1]=from[1];
}

inline void vectorCopy3D(const double *from, double *to)
{
  to[0]=from[0];
  to[1]=from[1];
  to[2]=from[2];
}

inline void vectorFlip3D(double *v)
{
  v[0]=-v[0];
  v[1]=-v[1];
  v[2]=-v[2];
}

inline void vectorCopyN(const double *from, double *to, int N)
{
    for(int i = 0; i < N; i++)
       to[i] = from[i];
}

inline void vectorCopy3D(const int *from, int *to)
{
  to[0]=from[0];
  to[1]=from[1];
  to[2]=from[2];
}

inline void vectorAbs3D(double *v)
{
    if(v[0] < 0) v[0] = -v[0];
    if(v[1] < 0) v[1] = -v[1];
    if(v[2] < 0) v[2] = -v[2];
}

inline double vectorMin3D(double *v)
{
    if(v[0] < v[1] && v[0] < v[2])
    return v[0];

    if(v[1] < v[2])
    return v[1];

    return v[2];
}

inline double vectorMax3D(double *v)
{
    if(v[0] > v[1] && v[0] > v[2])
    return v[0];

    if(v[1] > v[2])
    return v[1];

    return v[2];
}

inline void vectorComponentMin3D(double *v1,double *v2,double *min)
{
    if(v1[0] > v2[0])
        min[0] = v2[0];
    else
        min[0] = v1[0];

    if(v1[1] > v2[1])
        min[1] = v2[1];
    else
        min[1] = v1[1];

    if(v1[2] > v2[2])
        min[2] = v2[2];
    else
        min[2] = v1[2];
}

inline void vectorComponentMax3D(double *v1,double *v2,double *max)
{
    if(v1[0] > v2[0])
        max[0] = v1[0];
    else
        max[0] = v2[0];

    if(v1[1] > v2[1])
        max[1] = v1[1];
    else
        max[1] = v2[1];

    if(v1[2] > v2[2])
        max[2] = v1[2];
    else
        max[2] = v2[2];
}

inline void vectorCopy4D(const double *from, double *to)
{
  to[0]=from[0];
  to[1]=from[1];
  to[2]=from[2];
  to[3]=from[3];
}

inline void vectorScalarMultN(int n,double *v, double s)
{
    for(int i = 0; i < n; i++)
        v[i] = s*v[i];
}

inline void vectorScalarMult3D(double *v, double s)
{
  v[0]=s*v[0];
  v[1]=s*v[1];
  v[2]=s*v[2];
}

inline void vectorScalarMult3D(double *v, double s, double *result)
{
  result[0]=s*v[0];
  result[1]=s*v[1];
  result[2]=s*v[2];
}

inline void vectorScalarDiv3D(double *v, double s)
{
  double sinv = 1./s;
  v[0]=sinv*v[0];
  v[1]=sinv*v[1];
  v[2]=sinv*v[2];
}

inline void vectorScalarAdd3D(double *v, double s)
{
  v[0]+=s;
  v[1]+=s;
  v[2]+=s;
}

inline void vectorScalarSubtract3D(double *v, double s)
{
  v[0]-=s;
  v[1]-=s;
  v[2]-=s;
}

inline void vectorNegate3D(double *v, double *result)
{
  result[0]=-v[0];
  result[1]=-v[1];
  result[2]=-v[2];
}

inline void vectorNegate3D(double *v)
{
  v[0]=-v[0];
  v[1]=-v[1];
  v[2]=-v[2];
}

inline void vectorScalarDiv3D(double *v, double s, double *result)
{
  double sinv = 1./s;
  result[0]=sinv*v[0];
  result[1]=sinv*v[1];
  result[2]=sinv*v[2];
}

inline void vectorAdd3D(const double *v1, const double *v2, double *result)
{
  result[0]=v1[0]+v2[0];
  result[1]=v1[1]+v2[1];
  result[2]=v1[2]+v2[2];
}

inline void vectorAddN(double *v1, const double *v2, int n)
{
  for(int i = 0; i < n; i++)
    v1[i] += v2[i];
}

inline void vectorAddN(int *v1, const int *v2, int n)
{
  for(int i = 0; i < n; i++)
    v1[i] += v2[i];
}

inline void vectorAddMultiple3D(const double *v1, double v2factor, const double *v2, double *result)
{
  result[0]=v1[0]+v2factor*v2[0];
  result[1]=v1[1]+v2factor*v2[1];
  result[2]=v1[2]+v2factor*v2[2];
}

inline void vectorSubtract3D(const double *v1,const double *v2, double *result)
{
  result[0]=v1[0]-v2[0];
  result[1]=v1[1]-v2[1];
  result[2]=v1[2]-v2[2];
}

inline void vectorSubtract2D(const double *v1,const double *v2, double *result)
{
  result[0]=v1[0]-v2[0];
  result[1]=v1[1]-v2[1];
}

inline void vectorCross3D(const double *v1,const double *v2, double *result)
{
  result[0]=v1[1]*v2[2]-v1[2]*v2[1];
  result[1]=v1[2]*v2[0]-v1[0]*v2[2];
  result[2]=v1[0]*v2[1]-v1[1]*v2[0];
}

inline double vectorCrossMag3D(const double *v1,const double *v2)
{
  double res[3];
  res[0]=v1[1]*v2[2]-v1[2]*v2[1];
  res[1]=v1[2]*v2[0]-v1[0]*v2[2];
  res[2]=v1[0]*v2[1]-v1[1]*v2[0];
  return vectorMag3D(res);
}

inline void vectorZeroize2D(int *v)
{
  v[0]=0;
  v[1]=0;
}

inline void vectorZeroize2D(double *v)
{
  v[0]=0.;
  v[1]=0.;
}

inline void vectorZeroize3D(double *v)
{
  v[0]=0.;
  v[1]=0.;
  v[2]=0.;
}

inline void vectorZeroize3D(int *v)
{
  v[0]=0;
  v[1]=0;
  v[2]=0;
}

inline void vectorZeroize4D(double *v)
{
  v[0]=0.;
  v[1]=0.;
  v[2]=0.;
  v[3]=0.;
}

inline void vectorZeroizeN(double *v,int n)
{
  for(int i = 0; i < n; i++)
     v[i]=0.;
}

inline void vectorZeroizeN(int *v,int n)
{
  for(int i = 0; i < n; i++)
     v[i]=0;
}

inline void vectorInitialize3D(double *v,double init)
{
  v[0]=init;
  v[1]=init;
  v[2]=init;
}

inline void vectorInitializeN(int *v,int n,int init)
{
  for(int i = 0; i < n; i++)
     v[i]=init;
}

inline double vectorSumN(double *v,int n)
{
  double sum = 0.;
  for(int i = 0; i < n; i++)
     sum+=v[i];
  return sum;
}

inline void quatUnitize4D(double *q)
{
  q[0]=1.;
  q[1]=0.;
  q[2]=0.;
  q[3]=0.;
}

inline bool isUnitQuat4D(double *q)
{
    return
    (
        q[0] == 1. &&
        q[1] == 0. &&
        q[2] == 0. &&
        q[3] == 0.
    );
}

inline void normalize_bary(double *v)
{
  double mag = v[0]+v[1]+v[2];
  v[0]/=mag;
  v[1]/=mag;
  v[2]/=mag;
}

inline void vectorToBuf3D(double *vec,double *buf,int &m)
{
  buf[m++] = vec[0];
  buf[m++] = vec[1];
  buf[m++] = vec[2];
}

inline void bufToVector3D(double *vec,double *buf,int &m)
{
  vec[0] = buf[m++];
  vec[1] = buf[m++];
  vec[2] = buf[m++];
}

inline void vectorToBuf4D(double *vec,double *buf,int &m)
{
  buf[m++] = vec[0];
  buf[m++] = vec[1];
  buf[m++] = vec[2];
  buf[m++] = vec[3];
}

inline void bufToVector4D(double *vec,double *buf,int &m)
{
  vec[0] = buf[m++];
  vec[1] = buf[m++];
  vec[2] = buf[m++];
  vec[3] = buf[m++];
}

inline void printVec2D(FILE *out, const char *name, double *vec)
{
    fprintf(out," vector %s: %e %e\n",name,vec[0],vec[1]);
}

inline void printVec3D(FILE *out, const char *name, double *vec)
{
    fprintf(out," vector %s: %e %e %e\n",name,vec[0],vec[1],vec[2]);
}

inline void printVec3D(FILE *out, const char *name, int *vec)
{
    fprintf(out," vector %s: %d %d %d\n",name,vec[0],vec[1],vec[2]);
}

inline void printVec4D(FILE *out, const char *name, double *vec)
{
    fprintf(out," vector %s: %e %e %e %e\n",name,vec[0],vec[1],vec[2],vec[3]);
}

inline void printVecN(FILE *out, const char *name, double *vec, int n)
{
    fprintf(out," vector %s:",name);
    for(int i = 0; i < n; i++)
        fprintf(out,"%f ",vec[i]);
    fprintf(out,"\n");
}

inline void printVecN(FILE *out, const char *name, int *vec, int n)
{
    fprintf(out," vector %s:",name);
    for(int i = 0; i < n; i++)
        fprintf(out,"%d ",vec[i]);
    fprintf(out,"\n");
}

inline void printMat33(FILE *out, const char *name, double **mat)
{
    fprintf(out," matrix %s: %f %f %f\n",name,mat[0][0],mat[0][1],mat[0][2]);
    fprintf(out,"        %s: %f %f %f\n",name,mat[1][0],mat[1][1],mat[1][2]);
    fprintf(out,"        %s: %f %f %f\n",name,mat[2][0],mat[2][1],mat[2][2]);
}

}

#endif
