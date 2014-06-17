/*********************************************************************
* Software License Agreement (BSD License)
*
*  Copyright (c) 2008, Willow Garage, Inc.
*  All rights reserved.
*
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions
*  are met:
*
*   * Redistributions of source code must retain the above copyright
*     notice, this list of conditions and the following disclaimer.
*   * Redistributions in binary form must reproduce the above
*     copyright notice, this list of conditions and the following
*     disclaimer in the documentation and/or other materials provided
*     with the distribution.
*   * Neither the name of the Willow Garage nor the names of its
*     contributors may be used to endorse or promote products derived
*     from this software without specific prior written permission.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
*  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
*  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
*  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
*  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
*  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
*  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
*  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
*  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
*  POSSIBILITY OF SUCH DAMAGE.
*********************************************************************/

/* Author: Ioan Sucan */

/* Modifications from OMPL RRTConnect.cpp.
 * 
 * Copyright (c) 2014 Carnegie Mellon University
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

#include <pr_constraint/holonomic.h>
#include <pr_constraint_ompl/geometric/planners/rrt/CRRTConnect.h>
#include <ompl/base/goals/GoalSampleableRegion.h>
#include <ompl/base/spaces/RealVectorStateSpace.h>
#include <ompl/tools/config/SelfConfig.h>

pr_constraint_ompl::geometric::CRRTConnect::CRRTConnect(const ompl::base::SpaceInformationPtr &si) : ompl::base::Planner(si, "CRRTConnect")
{
    if (!dynamic_cast<ompl::base::RealVectorStateSpace*>(si->getStateSpace().get()))
      throw "for now, CRRTConnect only supports real vector state spaces!";
    
    specs_.recognizedGoal = ompl::base::GOAL_SAMPLEABLE_REGION;
    specs_.directed = true;

    maxDistance_ = 0.0;

    Planner::declareParam<double>("range", this, &CRRTConnect::setRange, &CRRTConnect::getRange, "0.:1.:10000.");
    connectionPoint_ = std::make_pair<ompl::base::State*, ompl::base::State*>(NULL, NULL);
    
    cons_ = NULL;
}

pr_constraint_ompl::geometric::CRRTConnect::~CRRTConnect(void)
{
    freeMemory();
}

void pr_constraint_ompl::geometric::CRRTConnect::setup(void)
{
    Planner::setup();
    ompl::tools::SelfConfig sc(si_, getName());
    sc.configurePlannerRange(maxDistance_);

    if (!tStart_)
        tStart_.reset(ompl::tools::SelfConfig::getDefaultNearestNeighbors<Motion*>(si_->getStateSpace()));
    if (!tGoal_)
        tGoal_.reset(ompl::tools::SelfConfig::getDefaultNearestNeighbors<Motion*>(si_->getStateSpace()));
    tStart_->setDistanceFunction(boost::bind(&CRRTConnect::distanceFunction, this, _1, _2));
    tGoal_->setDistanceFunction(boost::bind(&CRRTConnect::distanceFunction, this, _1, _2));
}

void pr_constraint_ompl::geometric::CRRTConnect::freeMemory(void)
{
    std::vector<Motion*> motions;

    if (tStart_)
    {
        tStart_->list(motions);
        for (unsigned int i = 0 ; i < motions.size() ; ++i)
        {
            if (motions[i]->state)
                si_->freeState(motions[i]->state);
            delete motions[i];
        }
    }

    if (tGoal_)
    {
        tGoal_->list(motions);
        for (unsigned int i = 0 ; i < motions.size() ; ++i)
        {
            if (motions[i]->state)
                si_->freeState(motions[i]->state);
            delete motions[i];
        }
    }
}

void pr_constraint_ompl::geometric::CRRTConnect::clear(void)
{
    Planner::clear();
    sampler_.reset();
    freeMemory();
    if (tStart_)
        tStart_->clear();
    if (tGoal_)
        tGoal_->clear();
    connectionPoint_ = std::make_pair<ompl::base::State*, ompl::base::State*>(NULL, NULL);
}

pr_constraint_ompl::geometric::CRRTConnect::GrowState pr_constraint_ompl::geometric::CRRTConnect::growTree(TreeData &tree, TreeGrowingInfo &tgi, Motion *rmotion)
{
    /* find closest state in the tree */
    Motion *nmotion = tree->nearest(rmotion);

    /* assume we can reach the state we go towards */
    bool reach = true;

    /* find state to add */
    ompl::base::State *dstate = rmotion->state;
    double d = si_->distance(nmotion->state, rmotion->state);
    if (d > maxDistance_)
    {
        si_->getStateSpace()->interpolate(nmotion->state, rmotion->state, maxDistance_ / d, tgi.xstate);
        dstate = tgi.xstate;
        reach = false;
    }
    
    /* do single-step constraint projection */
    if (cons_)
    {
        /* compute constraint value (~ disance to manifold)
         * (this is a total hack) */
        cons_->eval(
            dstate->as<ompl::base::RealVectorStateSpace::StateType>()->values,
            &cons_value_[0]);
        double norm2 = 0.0;
        for (unsigned int ui=0; ui<cons_value_.size(); ui++)
            norm2 += cons_value_[ui] * cons_value_[ui];
        if (norm2 > 0.0001 * 0.0001)
        {
            cons_->project(
               dstate->as<ompl::base::RealVectorStateSpace::StateType>()->values,
               tgi.pstate->as<ompl::base::RealVectorStateSpace::StateType>()->values);
            dstate = tgi.pstate;
            
            /* check distance between dstate and nmotion->state -- we could have gotten projected straight back! */
            double d_onmanifold = si_->distance(nmotion->state, dstate);
            if (d_onmanifold < 0.0001)
                return TRAPPED;
            
            reach = false;
        }
    }
    
    // if we are in the start tree, we just check the motion like we normally do;
    // if we are in the goal tree, we need to check the motion in reverse, but checkMotion() assumes the first state it receives as argument is valid,
    // so we check that one first
    bool validMotion = tgi.start ? si_->checkMotion(nmotion->state, dstate) : si_->getStateValidityChecker()->isValid(dstate) && si_->checkMotion(dstate, nmotion->state);

    if (validMotion)
    {
        /* create a motion */
        Motion *motion = new Motion(si_);
        si_->copyState(motion->state, dstate);
        motion->parent = nmotion;
        motion->root = nmotion->root;
        tgi.xmotion = motion;

        tree->add(motion);
        if (reach)
            return REACHED;
        else
            return ADVANCED;
    }
    else
        return TRAPPED;
}

ompl::base::PlannerStatus pr_constraint_ompl::geometric::CRRTConnect::solve(const ompl::base::PlannerTerminationCondition &ptc)
{
    checkValidity();
    ompl::base::GoalSampleableRegion *goal = dynamic_cast<ompl::base::GoalSampleableRegion*>(pdef_->getGoal().get());

    if (!goal)
    {
        OMPL_ERROR("%s: Unknown type of goal", getName().c_str());
        return ompl::base::PlannerStatus::UNRECOGNIZED_GOAL_TYPE;
    }

    while (const ompl::base::State *st = pis_.nextStart())
    {
        Motion *motion = new Motion(si_);
        si_->copyState(motion->state, st);
        motion->root = motion->state;
        tStart_->add(motion);
    }

    if (tStart_->size() == 0)
    {
        OMPL_ERROR("%s: Motion planning start tree could not be initialized!", getName().c_str());
        return ompl::base::PlannerStatus::INVALID_START;
    }

    if (!goal->couldSample())
    {
        OMPL_ERROR("%s: Insufficient states in sampleable goal region", getName().c_str());
        return ompl::base::PlannerStatus::INVALID_GOAL;
    }

    if (!sampler_)
        sampler_ = si_->allocStateSampler();

    OMPL_INFORM("%s: Starting with %d states", getName().c_str(), (int)(tStart_->size() + tGoal_->size()));

    TreeGrowingInfo tgi;
    tgi.xstate = si_->allocState();
    tgi.pstate = si_->allocState();

    Motion   *rmotion   = new Motion(si_);
    ompl::base::State *rstate = rmotion->state;
    bool startTree      = true;
    bool solved         = false;

    while (ptc == false)
    {
        TreeData &tree      = startTree ? tStart_ : tGoal_;
        tgi.start = startTree;
        startTree = !startTree;
        TreeData &otherTree = startTree ? tStart_ : tGoal_;

        if (tGoal_->size() == 0 || pis_.getSampledGoalsCount() < tGoal_->size() / 2)
        {
            const ompl::base::State *st = tGoal_->size() == 0 ? pis_.nextGoal(ptc) : pis_.nextGoal();
            if (st)
            {
                Motion* motion = new Motion(si_);
                si_->copyState(motion->state, st);
                motion->root = motion->state;
                tGoal_->add(motion);
            }

            if (tGoal_->size() == 0)
            {
                OMPL_ERROR("%s: Unable to sample any valid states for goal tree", getName().c_str());
                break;
            }
        }

        /* sample random state */
        sampler_->sampleUniform(rstate);

        GrowState gs = growTree(tree, tgi, rmotion);

        if (gs != TRAPPED)
        {
            /* remember which motion was just added */
            Motion *addedMotion = tgi.xmotion;

            /* attempt to connect trees */

            /* if reached, it means we used rstate directly, no need top copy again */
            if (gs != REACHED)
                si_->copyState(rstate, tgi.xstate);
            
            /* if constrained, then we always got to pstate */
            if (cons_)
                si_->copyState(rstate, tgi.pstate);

            GrowState gsc = ADVANCED;
            tgi.start = startTree;
            while (gsc == ADVANCED)
                gsc = growTree(otherTree, tgi, rmotion);

            Motion *startMotion = startTree ? tgi.xmotion : addedMotion;
            Motion *goalMotion  = startTree ? addedMotion : tgi.xmotion;

            /* if we connected the trees in a valid way (start and goal pair is valid)*/
            if (gsc == REACHED && goal->isStartGoalPairValid(startMotion->root, goalMotion->root))
            {
                // it must be the case that either the start tree or the goal tree has made some progress
                // so one of the parents is not NULL. We go one step 'back' to avoid having a duplicate state
                // on the solution path
                if (startMotion->parent)
                    startMotion = startMotion->parent;
                else
                    goalMotion = goalMotion->parent;

                connectionPoint_ = std::make_pair(startMotion->state, goalMotion->state);

                /* construct the solution path */
                Motion *solution = startMotion;
                std::vector<Motion*> mpath1;
                while (solution != NULL)
                {
                    mpath1.push_back(solution);
                    solution = solution->parent;
                }

                solution = goalMotion;
                std::vector<Motion*> mpath2;
                while (solution != NULL)
                {
                    mpath2.push_back(solution);
                    solution = solution->parent;
                }

                ompl::geometric::PathGeometric *path = new ompl::geometric::PathGeometric(si_);
                path->getStates().reserve(mpath1.size() + mpath2.size());
                for (int i = mpath1.size() - 1 ; i >= 0 ; --i)
                    path->append(mpath1[i]->state);
                for (unsigned int i = 0 ; i < mpath2.size() ; ++i)
                    path->append(mpath2[i]->state);

                pdef_->addSolutionPath(ompl::base::PathPtr(path), false, 0.0);
                solved = true;
                break;
            }
        }
    }

    si_->freeState(tgi.xstate);
    si_->freeState(tgi.pstate);
    si_->freeState(rstate);
    delete rmotion;

    OMPL_INFORM("%s: Created %u states (%u start + %u goal)", getName().c_str(), tStart_->size() + tGoal_->size(), tStart_->size(), tGoal_->size());

    return solved ? ompl::base::PlannerStatus::EXACT_SOLUTION : ompl::base::PlannerStatus::TIMEOUT;
}

void pr_constraint_ompl::geometric::CRRTConnect::getPlannerData(ompl::base::PlannerData &data) const
{
    Planner::getPlannerData(data);

    std::vector<Motion*> motions;
    if (tStart_)
        tStart_->list(motions);

    for (unsigned int i = 0 ; i < motions.size() ; ++i)
    {
        if (motions[i]->parent == NULL)
            data.addStartVertex(ompl::base::PlannerDataVertex(motions[i]->state, 1));
        else
        {
            data.addEdge(ompl::base::PlannerDataVertex(motions[i]->parent->state, 1),
                         ompl::base::PlannerDataVertex(motions[i]->state, 1));
        }
    }

    motions.clear();
    if (tGoal_)
        tGoal_->list(motions);

    for (unsigned int i = 0 ; i < motions.size() ; ++i)
    {
        if (motions[i]->parent == NULL)
            data.addGoalVertex(ompl::base::PlannerDataVertex(motions[i]->state, 2));
        else
        {
            // The edges in the goal tree are reversed to be consistent with start tree
            data.addEdge(ompl::base::PlannerDataVertex(motions[i]->state, 2),
                         ompl::base::PlannerDataVertex(motions[i]->parent->state, 2));
        }
    }

    // Add the edge connecting the two trees
    data.addEdge(data.vertexIndex(connectionPoint_.first), data.vertexIndex(connectionPoint_.second));
}
