/* Copyright (c) 2014 Carnegie Mellon University
 * Author: Chris Dellin <cdellin@gmail.com>
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 * - Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 * - Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the distribution.
 * - Neither the name of Carnegie Mellon University nor the names of its
 *   contributors may be used to endorse or promote products derived from
 *   this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE. 
 */

#include <cstdio>

#include <cstdlib>
#include <cstring>
#include <cblas.h>
#include <lapacke.h>
#include <pr_constraint/holonomic.h>
#define VEC_STACK_SIZE (16)
#define JAC_NUM_RADIUS (0.0001)

namespace pr_constraint
{

HolonomicConstraint::HolonomicConstraint()
{
}

HolonomicConstraint::~HolonomicConstraint()
{
}

int HolonomicConstraint::eval_dim(int * const n_config, int * const n_value)
{
   return ERR_NOT_IMPLEMENTED;
}

int HolonomicConstraint::eval(const double * const config, double * const value)
{
   return ERR_NOT_IMPLEMENTED;
}

int HolonomicConstraint::eval_jacobian(const double * const config, double * const jac)
{
   int i;
   int j;
   int err;
   int ret;
   int n_config;
   int n_value;
   double newconfig_stack[VEC_STACK_SIZE];
   double posvalue_stack[VEC_STACK_SIZE];
   double negvalue_stack[VEC_STACK_SIZE];
   double * newconfig;
   double * posvalue;
   double * negvalue;
   
   /* compute jacobian numerically */
   
   err = eval_dim(&n_config, &n_value);
   if (err)
      return err;
   
   newconfig = (n_config <= VEC_STACK_SIZE) ? newconfig_stack : new double[n_config];
   posvalue = (n_value <= VEC_STACK_SIZE) ? posvalue_stack : new double[n_value];
   negvalue = (n_value <= VEC_STACK_SIZE) ? negvalue_stack : new double[n_value];
   
   memcpy(newconfig, config, n_config*sizeof(double));
   for (j=0; j<n_config; j++)
   {
      newconfig[j] = config[j] + JAC_NUM_RADIUS;
      err = eval(newconfig, posvalue);
      if (err) { ret = err; goto error; }
      newconfig[j] = config[j] - JAC_NUM_RADIUS;
      err = eval(newconfig, negvalue);
      if (err) { ret = err; goto error; }
      newconfig[j] = config[j];
      for (i=0; i<n_value; i++)
         jac[i*n_config+j] = (posvalue[i] - negvalue[i]) / (2.0 * JAC_NUM_RADIUS);
   }
   
   ret = 0;
error:
   if (newconfig != newconfig_stack) delete newconfig;
   if (posvalue != posvalue_stack) delete posvalue;
   if (negvalue != negvalue_stack) delete negvalue;
   return ret;
}

int HolonomicConstraint::project(const double * const config, double * new_config)
{
   int i;
   int err;
   int ret;
   int n_config;
   int n_value;

   /* heap-allocated variables */
   double * heap;
#define VARS \
   X(value, n_value) \
   X(J, n_value*n_config) \
   X(JJT, n_value*n_value)
#define X(n,s) double * n;
   VARS
#undef X
   int * ipiv;
   
   /* do projection via newtons method */
   
   /* get sizes */
   err = eval_dim(&n_config, &n_value);
   if (err)
      return err;
   
   /* heap allocate */
   heap = (double *) malloc(sizeof(double) * (0
#define X(n,s) + (s)
      VARS
#undef X
      ) + n_value * sizeof(int));
   i = 0;
#define X(n,s) n = heap + i; i += (s);
   VARS
#undef X
   ipiv = (int *)(heap + i);
   
   
   memcpy(new_config, config, n_config*sizeof(double));
   
   
   /* one step */
   for (int iter=0; iter<3; iter++)
   {   
      err = eval(new_config, value);
      if (err) { ret = err; goto error; }
      
      err = eval_jacobian(new_config, J);
      if (err) { printf("error with eval_jacobian() ...\n"); ret = err; goto error; }
      
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
         n_value, n_value, n_config,
         1.0, J, n_config, J, n_config, 0.0, JJT, n_value);
      
      for (i=0; i<n_value; i++)
         JJT[i*n_value+i] += 1.0e-9;
      err = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n_value, 1, JJT, n_value, ipiv, value, 1);
      if (err) { printf("error %d with LAPACKE_dgesv() ...\n", err); err = 5; goto error; }
      
      cblas_dgemv(CblasRowMajor, CblasTrans,
         n_value, n_config,
         -1.0, J, n_config, value, 1, 1.0, new_config, 1);
   }
   
   ret = 0;
error:
   free(heap);
   return ret;
}

} /* namespace pr_constraint */
