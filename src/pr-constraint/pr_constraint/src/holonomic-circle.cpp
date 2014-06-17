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
#include <pr_constraint/holonomic.h>
#include <pr_constraint/holonomic-circle.h>

namespace pr_constraint
{

HolonomicCircleConstraint::HolonomicCircleConstraint(double center[2], double radius)
{
   _center[0] = center[0];
   _center[1] = center[1];
   _radius = radius;
}

HolonomicCircleConstraint::~HolonomicCircleConstraint()
{
}

int HolonomicCircleConstraint::eval_dim(int * const n_config, int * const n_value)
{
   if (n_config) *n_config = 2;
   if (n_value) *n_value = 1;
   return 0;
}

int HolonomicCircleConstraint::eval(const double * const config, double * const value)
{
   double dx, dy;
   double dist;
   dx = config[0] - _center[0];
   dy = config[1] - _center[1];
   dist = sqrt(dx*dx + dy*dy);
   value[0] = dist - _radius;
   return 0;
}

} /* namespace pr_constraint */
