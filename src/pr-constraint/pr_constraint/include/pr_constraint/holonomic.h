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

namespace pr_constraint
{

enum Error
{
  SUCCESS = 0,
  ERR_NOT_IMPLEMENTED
};

class HolonomicConstraint
{
public:
   HolonomicConstraint();
   virtual ~HolonomicConstraint();
   
   /* the constraint as a function, with 0 the manifold itself;
    * implementing classes should implement eval_dim() and eval(),
    * and eval_jacobian() if it is analytically available
    */
   
   /* the dimensionality of the constraint eval function
    * config space is n_config
    * value is n_value (1 = scalar valued constraint)
    */
   virtual int eval_dim(int * const n_config, int * const n_value);
   
   /* a function that is 0-valued when the constraint is satisfied
    * input: config, a n_config-length configuration space coordinate
    * output: value, a n_value-length function value (n_value=1 is scalar)
    */
   virtual int eval(const double * const config, double * const value);
   
   /* the jacobian of the eval function
    * (will be numerically computed if necessary and not implemented)
    * jac is in row-major order
    * (a single row vector for scalar-valued functions)
    */
   virtual int eval_jacobian(const double * const config, double * const jac);
   
   /* hessian tensor? */
   
   
   /* for now, this uses newton's method */
   virtual int project(const double * const config, double * new_config);
};

} /* namespace pr_constraint */
