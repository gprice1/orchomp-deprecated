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
#include <ompl/base/spaces/RealVectorStateSpace.h>
#include <ompl/base/Planner.h>
#include <ompl/base/ProblemDefinition.h>
#include <ompl/base/ScopedState.h>
#include <ompl/base/SpaceInformation.h>
#include <ompl/geometric/planners/rrt/RRT.h>
#include <ompl/geometric/planners/rrt/RRTConnect.h>
#include <ompl/geometric/PathGeometric.h>
#include <pr_constraint/holonomic.h>
#include <pr_constraint/holonomic-circle.h>
#include <pr_constraint_ompl/geometric/planners/rrt/CRRT.h>

bool isStateValid(const ompl::base::State *state)
{
   /*printf("checking feasibility ...\n");*/
   return true;
}

/* see http://ompl.kavrakilab.org/geometricPlanningSE3.html
 * for ompl tutorial */
 
 /* run like
  * rosrun pr_constraint_tests test-ompl-circle
  * pdflatex -halt-on-error plot.tex
  */
int main()
{
   printf("started tets_ompl_circle!\n");
   
   /* state space */
   ompl::base::StateSpacePtr space(new ompl::base::RealVectorStateSpace(2));
   
   /* state space set bounds */
   {
      ompl::base::RealVectorBounds bounds(2);
      bounds.setLow(-2.0);
      bounds.setHigh(2.0);
      space->as<ompl::base::RealVectorStateSpace>()->setBounds(bounds);
   }
   
   /* space information container */
   ompl::base::SpaceInformationPtr si(new ompl::base::SpaceInformation(space));
   
   si->setStateValidityCheckingResolution(0.01);
   
   /* feasibility callback set on si container */
   si->setStateValidityChecker(boost::bind(&isStateValid, _1));
   
   /* create a problem from this si container */
   ompl::base::ProblemDefinitionPtr pdef(new ompl::base::ProblemDefinition(si));
   {
      ompl::base::ScopedState<ompl::base::RealVectorStateSpace> start(space);
      start[0] = -1.0;
      start[1] = 0.0;
      ompl::base::ScopedState<ompl::base::RealVectorStateSpace> goal(space);
      goal[0] = 1.0;
      goal[1] = 0.0;
      pdef->setStartAndGoalStates(start, goal);
   }
   
   /* create a planner from the si, set the problem definition */
#if 0
   ompl::base::PlannerPtr planner(new ompl::geometric::RRTConnect(si));
   planner->as<ompl::geometric::RRTConnect>()->setRange(0.1);
#endif
   
#if 0
   ompl::geometric::RRT * planner_rrt = new ompl::geometric::RRT(si);
   planner_rrt->setRange(0.1);
   ompl::base::PlannerPtr planner(planner_rrt);
#endif
   
   /* create the constraint */
   double cons_center[2] = {0.0, 0.0};
   pr_constraint::HolonomicCircleConstraint cons(cons_center, 1.0);
   
   pr_constraint_ompl::geometric::CRRT * planner_crrt = new pr_constraint_ompl::geometric::CRRT(si);
   planner_crrt->setRange(0.1);
   planner_crrt->setTrajectoryWideConstraint(&cons);
   ompl::base::PlannerPtr planner(planner_crrt);
   
   planner->setProblemDefinition(pdef);
   planner->setup();
   ompl::base::PlannerStatus solved = planner->solve(1.0); /* seconds */
   if (solved)
   {
      printf("path found!\n");
      FILE * fp;
      ompl::base::PathPtr path = pdef->getSolutionPath();
      ompl::geometric::PathGeometric * gpath = dynamic_cast<ompl::geometric::PathGeometric*>(path.get());
      fp = fopen("plot.tex","w");
      fprintf(fp,"\\documentclass{standalone}\n");
      fprintf(fp,"\\usepackage{tikz}\n");
      fprintf(fp,"\\usetikzlibrary{backgrounds,arrows,shapes,trees}\n");
      fprintf(fp,"\\begin{document}\n");
      fprintf(fp,"\\begin{tikzpicture}[x=1in,y=1in,show background rectangle]\n");
      fprintf(fp,"\\clip (-2,-2) rectangle (2,2); %% everything outside this is clipped out\n");
      for (int i=0; i<gpath->getStateCount()-1; i++)
      {
         ompl::base::ScopedState<ompl::base::RealVectorStateSpace> s1(space);
         ompl::base::ScopedState<ompl::base::RealVectorStateSpace> s2(space);
         s1 = gpath->getState(i);
         s2 = gpath->getState(i+1);
         fprintf(fp,"\\draw (%f,%f) -- (%f,%f);\n", s1[0],s1[1], s2[0],s2[1]);
      }
      fprintf(fp,"\\end{tikzpicture}\n");
      fprintf(fp,"\\end{document}\n");
      fclose(fp);
   }
   else
      printf("path not found in allotted time!\n");

   return 0;
}
