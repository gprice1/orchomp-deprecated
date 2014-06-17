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

#include <cmath>
#include <cstring>
extern "C"
{
#include <libcd/kin.h>
#include <libcd/mat.h>
}
#include <pr_constraint/holonomic.h>
#include <pr_constraint_tsr/tsr.h>

namespace pr_constraint_tsr
{

TSRConstraint::TSRConstraint(double pose_0_w[7], double Bw[6][2], double pose_w_e[7])
{
   if (pose_0_w)
      memcpy(_pose_0_w, pose_0_w, sizeof(double [7]));
   else
      cd_kin_pose_identity(_pose_0_w);
   
   if (Bw)
      memcpy(_Bw, Bw, sizeof(double [6][2]));
   else
      cd_mat_set_zero(&_Bw[0][0], 6, 2);
   
   if (pose_w_e)
      memcpy(_pose_w_e, pose_w_e, sizeof(double [7]));
   else
      cd_kin_pose_identity(_pose_w_e);
   
   recalc();
}

TSRConstraint::~TSRConstraint()
{
}

void TSRConstraint::recalc()
{
   int i;
   /* calculate pose inverses */
   cd_kin_pose_invert(_pose_0_w, _pose_wprox_0);
   cd_kin_pose_invert(_pose_w_e, _pose_e_wdist);
   /* calculate dimensionalty of volume (Dmitry's _dimensionality) */
   _dim_volume = 0;
   for (i=0; i<6; i++)
      if (_Bw[i][0] != _Bw[i][1])
         _dim_volume++;
   /* calcuate dimensionalty of constraint */
   _dim = 0;
   for (i=0; i<6; i++)
   {
      if (-HUGE_VAL < _Bw[i][0] || _Bw[i][1] < HUGE_VAL)
      {
         _idim_Bw[_dim] = i;
         _dim++;
      }
   }
}

int TSRConstraint::eval_dim(int * const n_config, int * const n_value)
{
   if (n_config) *n_config = 7;
   if (n_value) *n_value = _dim;
   return 0;
}

int TSRConstraint::eval(const double * const config, double * const value)
{
   int i;
   int iBw;
   int ixyz;
   double pose_wprox_wdist[7];
   double xyzypr[6];
   cd_kin_pose_compose(_pose_wprox_0, config, pose_wprox_wdist);
   cd_kin_pose_compose(pose_wprox_wdist, _pose_e_wdist, pose_wprox_wdist);
   cd_kin_pose_to_xyzypr(pose_wprox_wdist, xyzypr);
   for (i=0; i<_dim; i++)
   {
      iBw = _idim_Bw[i];
      ixyz = (iBw < 3) ? iBw : 8 - iBw;
      if (xyzypr[ixyz] < _Bw[iBw][0])
         value[i] = xyzypr[ixyz] - _Bw[iBw][0];
      else if (_Bw[iBw][1] < xyzypr[ixyz])
         value[i] = xyzypr[ixyz] - _Bw[iBw][1];
      else
         value[i] = 0.0;
   }
   return 0;
}

} /* namespace pr_constraint_tsr */
