/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2012/03/07 21:41:19 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Vector constructors for the SUNDIALS Matlab interfaces.
 * -----------------------------------------------------------------
 */

#include <stdlib.h>
#include "nvm.h"

void PutData(N_Vector v, double *data, long int n)
{
  double *vdata;
  long int i;

  vdata = N_VGetArrayPointer(v);

  for(i=0;i<n;i++) vdata[i] = data[i];

  return;
}


void GetData(N_Vector v, double *data, long int n)
{
  double *vdata;
  long int i;

  vdata = N_VGetArrayPointer(v);

  for(i=0;i<n;i++) data[i] = vdata[i];

  return;
}
