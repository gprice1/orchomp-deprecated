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
#include <string>
#include <openrave-core.h>
#include <pr_constraint/holonomic.h>
#include <pr_constraint_openrave/robot-link.h>

namespace pr_constraint_openrave
{

OpenraveRobotLinkConstraint::OpenraveRobotLinkConstraint(
   OpenRAVE::RobotBasePtr probot,
   std::string link_name, bool is_manip,
   pr_constraint::HolonomicConstraint * sub_constraint)
{
   _probot = probot;
   _orlink = probot->GetLink(link_name);
   _ortx_offset.identity();
   _ortx_offset.trans.z = 0.16;
   _sub_constraint = sub_constraint;
}

OpenraveRobotLinkConstraint::~OpenraveRobotLinkConstraint()
{
}

int OpenraveRobotLinkConstraint::eval_dim(int * const n_config, int * const n_value)
{
   int err;
   int sub_n_config;
   int sub_n_value;
   err = _sub_constraint->eval_dim(&sub_n_config, &sub_n_value);
   if (err) return err;
   if (sub_n_config != 7)
      return -5;
   *n_config = _probot->GetActiveDOF();
   *n_value = sub_n_value;
   return 0;
}

int OpenraveRobotLinkConstraint::eval(const double * const config, double * const value)
{
   int i;
   double pose[7];
   std::vector<OpenRAVE::dReal> dofvals(_probot->GetActiveDOF());
   for (i=0; i<_probot->GetActiveDOF(); i++)
      dofvals[i] = config[i];
   _probot->SetActiveDOFValues(dofvals);
   OpenRAVE::Transform ortx = _orlink->GetTransform() * _ortx_offset;
   pose[0] = ortx.trans.x;
   pose[1] = ortx.trans.y;
   pose[2] = ortx.trans.z;
   pose[3] = ortx.rot.y;
   pose[4] = ortx.rot.z;
   pose[5] = ortx.rot.w;
   pose[6] = ortx.rot.x;
   _sub_constraint->eval(pose, value);
   return 0;
}

} /* namespace pr_constraint_openrave */
