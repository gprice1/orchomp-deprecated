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
#include <boost/thread.hpp>
#include <openrave-core.h>
#include <openrave/utils.h>
#include <openrave/interface.h>
#include <openrave/planningutils.h>
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
#include <pr_constraint_tsr/tsr.h>
#include <pr_constraint_openrave/robot-link.h>
#include <pr_constraint_ompl/geometric/planners/rrt/CRRT.h>
#include <pr_constraint_ompl/geometric/planners/rrt/CRRTConnect.h>

OpenRAVE::EnvironmentBasePtr penv;
OpenRAVE::RobotBasePtr probot;

void viewermain()
{
  try
  {
    OpenRAVE::ViewerBasePtr viewer;
    viewer = OpenRAVE::RaveCreateViewer(penv, "qtcoin");
    penv->AddViewer(viewer);
    viewer->main(true);
  }
  catch (OpenRAVE::openrave_exception ex)
  {
    printf("Could not setup OpenRAVE viewer: %s\n", ex.what());
  }
}

bool isStateValid(const ompl::base::State *state)
{
   /*printf("checking feasibility ...\n");*/
   return true;
}

/* see http://ompl.kavrakilab.org/geometricPlanningSE3.html
 * for ompl tutorial */
int main()
{
   int i;
   
   printf("started test-ompl-or-tsr!\n");
   
   /* init openrave */
   OpenRAVE::RaveInitialize(true, OpenRAVE::Level_Info); /* plugins, level */
   penv = OpenRAVE::RaveCreateEnvironment();
   boost::thread viewerthread = boost::thread(viewermain);
   probot = penv->ReadRobotXMLFile("robots/barrettwam.robot.xml");
   if (!probot)
   {
      printf("Tried to open robot xml \"robots/barrettwam.robot.xml\", but loading failed!");
      return 2;
   }
   penv->Add(probot);
   
   {
      int d[] = { 0, 1, 2, 3, 4, 5, 6 };
      std::vector<int> dofindices(d, d+sizeof(d)/sizeof(d[0]));
      probot->SetActiveDOFs(dofindices);
   }
   
   printf("robot active dof: %d\n", probot->GetActiveDOF());
   
   /* state space */
   ompl::base::StateSpacePtr space(new ompl::base::RealVectorStateSpace(probot->GetActiveDOF()));
   
   /* state space set bounds */
   {
      std::vector<OpenRAVE::dReal> limits_lower;
      std::vector<OpenRAVE::dReal> limits_upper;
      probot->GetDOFLimits(limits_lower, limits_upper);
      ompl::base::RealVectorBounds bounds(probot->GetActiveDOF());
      for (i=0; i<probot->GetActiveDOF(); i++)
      {
         bounds.setLow(i, limits_lower[i]);
         bounds.setHigh(i, limits_upper[i]);
      }
      space->as<ompl::base::RealVectorStateSpace>()->setBounds(bounds);
   }
   
   /* space information container */
   ompl::base::SpaceInformationPtr si(new ompl::base::SpaceInformation(space));
   
   si->setStateValidityCheckingResolution(0.01);
   
   /* feasibility callback set on si container */
   si->setStateValidityChecker(boost::bind(&isStateValid, _1));
   
   /* construct the tsr constraint info */
   double tsr_pose_0_w[7]; /* filled from fk later */
   double tsr_Bw[6][2] = {{0.0, 0.0}, {-1.0, 1.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};
   double tsr_pose_w_e[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0};
   
   {
      OpenRAVE::ModuleBasePtr pikfast = RaveCreateModule(penv, "ikfast");
      penv->Add(pikfast,true,"");
      std::stringstream ssin,ssout;
      ssin << "LoadIKFastSolver " << probot->GetName() << " " << (int)OpenRAVE::IKP_Transform6D;
      if(!pikfast->SendCommand(ssout,ssin))
      {
         printf("failed to load iksolver\n");
         penv->Destroy();
         OpenRAVE::RaveDestroy();
         return 1;
      }
   }
   
   /* create a problem from this si container */
   ompl::base::ProblemDefinitionPtr pdef(new ompl::base::ProblemDefinition(si));
   {
      double d[] = { 0.0, 0.60, 0.0, 2.0, 0.0, 0.5415926536, 0.0 };
      std::vector<OpenRAVE::dReal> dofvals(d, d+sizeof(d)/sizeof(d[0]));
      probot->SetActiveDOFValues(dofvals);
      OpenRAVE::Transform ortx = probot->GetActiveManipulator()->GetEndEffectorTransform();
      tsr_pose_0_w[0] = ortx.trans.x;
      tsr_pose_0_w[1] = ortx.trans.y;
      tsr_pose_0_w[2] = ortx.trans.z;
      tsr_pose_0_w[3] = ortx.rot.y;
      tsr_pose_0_w[4] = ortx.rot.z;
      tsr_pose_0_w[5] = ortx.rot.w;
      tsr_pose_0_w[6] = ortx.rot.x;
      std::vector<OpenRAVE::dReal> vsolution;
      OpenRAVE::Transform ortx_start = ortx;
      ortx_start.trans.y += -0.3;
      if (!probot->GetActiveManipulator()->FindIKSolution(OpenRAVE::IkParameterization(ortx_start),vsolution,0))
      {
         printf("failed to load iksolver\n");
         penv->Destroy();
         OpenRAVE::RaveDestroy();
         return 1;
      }
      ompl::base::ScopedState<ompl::base::RealVectorStateSpace> start(space);
      for (i=0; i<probot->GetActiveDOF(); i++)
         start[i] = vsolution[i];
      OpenRAVE::Transform ortx_goal = ortx;
      ortx_goal.trans.y += 0.3;
      if (!probot->GetActiveManipulator()->FindIKSolution(OpenRAVE::IkParameterization(ortx_goal),vsolution,0))
      {
         printf("failed to load iksolver\n");
         penv->Destroy();
         OpenRAVE::RaveDestroy();
         return 1;
      }
      ompl::base::ScopedState<ompl::base::RealVectorStateSpace> goal(space);
      for (i=0; i<probot->GetActiveDOF(); i++)
         goal[i] = vsolution[i];
      pdef->setStartAndGoalStates(start, goal);
   }
   
   /* create constraints */
   pr_constraint_tsr::TSRConstraint cons_tsr(tsr_pose_0_w, tsr_Bw, tsr_pose_w_e);
   pr_constraint_openrave::OpenraveRobotLinkConstraint cons(probot, "wam7", false, &cons_tsr);

#if 0
   {
      int err;
      int n_config;
      int n_value;
      err = cons_tsr.eval_dim(&n_config, &n_value);
      if (err)
      {
         printf("failed to call cons.eval_dim()\n");
         penv->Destroy();
         OpenRAVE::RaveDestroy();
         return 1;
      }
      printf("n_config: %d, n_value: %d\n", n_config, n_value);
      std::vector<OpenRAVE::dReal> adofvals;
      probot->GetActiveDOFValues(adofvals);
      double config[7];
      double value[6];
      for (i=0; i<7; i++)
         config[i] = adofvals[i];
      err = cons.eval(config, value);
      printf("value: %f %f %f %f %f %f\n",
         value[0], value[1], value[2], value[3], value[4], value[5]);
   }
#endif
   
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
   
#if 0
   pr_constraint_ompl::geometric::CRRT * planner_derived = new pr_constraint_ompl::geometric::CRRT(si);
   planner_derived->setRange(0.1);
#endif

   pr_constraint_ompl::geometric::CRRTConnect * planner_derived = new pr_constraint_ompl::geometric::CRRTConnect(si);
   planner_derived->setRange(0.1);

   planner_derived->setTrajectoryWideConstraint(&cons);
   
   ompl::base::PlannerPtr planner(planner_derived);
   
   planner->setProblemDefinition(pdef);
   planner->setup();
   ompl::base::PlannerStatus solved = planner->solve(1.0); /* seconds */
   if (solved)
   {
      ompl::base::PathPtr path = pdef->getSolutionPath();
      ompl::geometric::PathGeometric * gpath = dynamic_cast<ompl::geometric::PathGeometric*>(path.get());
      OpenRAVE::TrajectoryBasePtr ptraj;
      
      printf("path found!\n");
      
      ptraj = OpenRAVE::RaveCreateTrajectory(penv);
      ptraj->Init(probot->GetActiveConfigurationSpecification());
      for (i=0; i<gpath->getStateCount(); i++)
      {
         ompl::base::ScopedState<ompl::base::RealVectorStateSpace> s(space);
         s = gpath->getState(i);
         std::vector<OpenRAVE::dReal> vec(&s[0], &s[0]+probot->GetActiveDOF());
         ptraj->Insert(i, vec);
      }
      OpenRAVE::planningutils::RetimeActiveDOFTrajectory(ptraj,probot,false,1.0,1.0,"","");
      
      while (1)
      {
         probot->GetController()->SetPath(ptraj);
         while (1)
         {
            if (probot->GetController()->IsDone())
               break;
            sleep(1);
         }
      }
   }
   else
      printf("path not found in allotted time!\n");

   sleep(9999);
   penv->Destroy();
   OpenRAVE::RaveDestroy();
   return 0;
}
