/** \file orcdchomp_mod.h
 * \brief Interface to the orcdchomp module, an implementation of CHOMP
 *        using libcd.
 * \author Christopher Dellin
 * \edited by Temple Price for use with multigrid chomp
 * \date 2012
 */

/* (C) Copyright 2012-2013 Carnegie Mellon University */

/* This module (orcdchomp) is part of libcd.
 *
 * This module of libcd is free software: you can redistribute it
 * and/or modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This module of libcd is distributed in the hope that it will be
 * useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * A copy of the GNU General Public License is provided with libcd
 * (license-gpl.txt) and is also available at <http://www.gnu.org/licenses/>.
 */

/* requires:
 *  - openrave/openrave.h
 * */

#ifndef _ORCHOMP_MOD_H_
#define _ORCHOMP_MOD_H_


#include "chomp-multigrid/chomp/Chomp.h"
#include "orchomp_distancefield.h"
#include "orchomp_constraint.h"
#include "orchomp_kdata.h"
#include "orchomp_collision.h"

#include <openrave/openrave.h>
#include <openrave/planningutils.h>
#include <stack>

#include "utils/timer.h"


#define DEBUG_COUT 0
#define debugStream \
    if (DEBUG_COUT) {} \
    else std::cout



namespace orchomp
{
//this is a structure used to hold and initialize values
//  for an eventual call to chomp.
class ChompInfo {
public:
    // alpha : the gradient step
    // obstol: Relative error of the objective function - 
    //          once the error gets to a low enough point,
    //          reltive to the objective function, quit chomping.
    // t_total: the total amount of timesteps that
    //          N+1 timesteps will take            
    // gamma: how much the gradient term counts for.
    // epsilon: how far outside of obstacles costs start accruing.
    // epsilon_self: how far outside of the collision geometry
    // obs_factor :  The percentage that environmental collisions count
    //               towards total cost. Should be between 1 and 0.
    //               NOTE: obs_factor + obs_factor_self = 1.0
    //          self collision costs start accruing.
    // obs_factor_self :  The percentage that self collisions count
    //               towards total cost. Should be between 1 and 0.
    //               NOTE: obs_factor + obs_factor_self = 1.0
    // jointPadding : pads the joint limit values so that the constraints
    //                start before the limits are exceeded. 
    //                Must be a value between 0 and 1, however, it 
    //                should be very low. (between 0.01 and 0).
    double alpha, obstol, t_total, gamma, epsilon, epsilon_self, obs_factor,
           obs_factor_self, jointPadding, timeout_seconds;

    //n: the initial size of the trajectory,
    //n_max: the final size,
    //min_global_iter: the min # of global chomp interations,
    //max_global_iter: the max # of global chomp interations,
    //min_local_iter: the min # of local smoothing iterations
    //max_local_iter: the max # of local smoothing iterations
    size_t n, n_max, min_global_iter, max_global_iter,
                     min_local_iter, max_local_iter;

    //doGlobal/doLocal: whether or not global and/or local chomp should
    //                  be done.
    // doObserve : whether or not chomp should use an observer to 
    //             print out debug information
    // noFactory : if true, chomp will not be given a constraint factory,
    //             so no constraints will be enforced
    // noCollider : if true, chomp will not be given a gradient descender,
    //              so it will not descend out of collisions
    // noSelfCollision : the collider will not check for self-collision
    // noEnvironmentalCollision : if true, the collider will not check for
    //                           collisions between the robot and environ.
    bool doGlobal, doLocal, doObserve, noFactory, noCollider, 
         noSelfCollision, noEnvironmentalCollision;

    //a basic constructor to initialize values
    ChompInfo() :
        alpha(0.1), obstol(0.01), t_total(1.0), gamma(0.1),
        epsilon( 0.1 ), epsilon_self( 0.01 ), obs_factor( 0.7 ),
        obs_factor_self( 0.3 ), jointPadding( 0.001 ),
        timeout_seconds( -1.0),
        n(0), n_max(0),
        min_global_iter( 0 ), max_global_iter( size_t(-1) ), 
        min_local_iter( 0 ), max_local_iter( size_t(-1)), doGlobal( true ),
        doLocal( false), doObserve( false ), noFactory (false),
        noCollider( false ), noSelfCollision( false ),
        noEnvironmentalCollision( false )
        {}
};



/* the module itself */
class mod : public OpenRAVE::ModuleBase
{
public:

    //____________________PUBLIC MEMBER VARIABLES____________________//
    OpenRAVE::EnvironmentBasePtr environment; /* filled on module creation */
 
    //the trajectory, start, and endpoint.
    chomp::MatX trajectory, q0, q1;
    ORConstraintFactory * factory;
   
    //This holds basic info relating to an individual 
    //   run of chomp
    ChompInfo info;
  
    // This observes chomp for general debugging purposes
    chomp::DebugChompObserver observer;

    //these are useful for interfacing with the robot.
    OpenRAVE::RobotBasePtr robot;      // a pointer to the robot
    std::string robot_name;            // the name of the robot
    std::vector< int > active_indices; // the active indices of the robot
    size_t n_dof;                      // the degree of freedom of the bot.
 
    //holds the active manipulator of the robot. this is used for
    //  TSR constraints.
    OpenRAVE::RobotBase::ManipulatorPtr active_manip;
    

    //a vector containing the collision geometry
    std::vector< Sphere > active_spheres, inactive_spheres;
    
    //This holds information pertinent to the collision geometry, and
    //   assists chomp in computing collision gradients.
    SphereCollisionHelper * sphere_collider;
    chomp::ChompCollGradHelper * collisionHelper;

    //the upper and lower limits for the joints.
    std::vector< OpenRAVE::dReal > upperJointLimits, lowerJointLimits,
                        paddedUpperJointLimits, paddedLowerJointLimits;
                                    
    //This vector holds all of the sdf's.
    std::vector< DistanceField > sdfs;
    
    //this vector holds all of the TSR's 
    std::vector< ORTSRConstraint * > tsrs;
    
    //this is a timer for timing things.
    Timer timer;

    //this is a pointer to the chomp class that will pull most of the
    //   weight for the module.
    chomp::Chomp * chomper;
    
    //a pointer to an openrave trajectory, a call to gettraj, will fill
    //  this structure with the current chomp trajectory.
    OpenRAVE::TrajectoryBasePtr trajectory_ptr;


    //_______________________PUBLIC MEMBER FUNCTIONS___________________//

    ///////////////The following functions are in orchomp_mod.cpp

    //visualize the collision geometry 
    int viewspheres(std::ostream & sout, std::istream& sinput);

    //compute the distance field for use in collision detection, and
    //   descending the gradient out of collision
    int computedistancefield(std::ostream & sout, std::istream& sinput);
    
    //visualize a slice out of a signed distance field.
    int visualizeslice(std::ostream& sout, std::istream& sinput);

    // NOTE : Find out what this is supposed to do
    int addfield_fromobsarray(std::ostream & sout, std::istream& sinput);

    //create a chomp run, to prepare for running chomp.
    int create(std::ostream & sout, std::istream& sinput);

    //GO through one iteration of chomp
    int iterate(std::ostream & sout, std::istream& sinput);

    //Get the current trajectory
    int gettraj(std::ostream & sout, std::istream& sinput);

    //destroy the current chomp iteration.
    int destroy(std::ostream & sout, std::istream& sinput);
   
    //execute a construsted trajectory
    int execute(std::ostream & sout, std::istream& sinput);
    
    //execute a construsted trajectory
    int playback(std::ostream & sout, std::istream& sinput);
    
    //add a tsr to the factory
    int addtsr(std::ostream & sout, std::istream& sinput);
    //create a tsr, and add it to the factory. 
    //  This is identical to the above function,
    //  except that it uses a different method to parse the tsr.
    int createtsr(std::ostream & sout, std::istream& sinput);

    //creates a box in the environment to visualize a given TSR.
    int viewtsr(std::ostream & sout, std::istream& sinput);
    
    //Remove a constraint from the list of constraints 
    int removeconstraint(std::ostream & sout, std::istream& sinput);
    
    //view the collision geometry of the robot at a given robot
    //  configuration
    int viewspheresVec( const chomp::MatX & q,
                        const std::vector< OpenRAVE::dReal > & vec,
                        double time);
 
    //constructor that registers all of the commands to the openRave
    //   command line interface.
    mod(OpenRAVE::EnvironmentBasePtr penv);

    //Destructor
    virtual ~mod() {}
    void Destroy() { RAVELOG_INFO("module unloaded from environment\n"); }

    /* This is called on e.LoadProblem(m, 'command') */
    int main(const std::string& cmd)
       { RAVELOG_INFO("module init cmd: %s\n", cmd.c_str()); return 0; }

    //Get the collision geometry from the XML files
    void getSpheres();
  
  ////////////////////////////////////////////////////////////////////
  ////// these functions can be found in orchomp_mod_parse ///////////
  ////////////////////////////////////////////////////////////////////
  private: 
    // this is all helper code for the parsing.
    // The source code for these functions is in orchomp_mod_parse.cpp, 
    //      not the same file that contains many of the other functions.
    void parseCreate(std::ostream & sout, std::istream& sinput);
    void parseViewSpheres(std::ostream & sout, std::istream& sinput);
    void parseIterate(std::ostream & sout, std::istream& sinput);
    void parseGetTraj(std::ostream & sout, std::istream& sinput);
    void parseDestroy(std::ostream & sout, std::istream& sinput);
    void parseComputeDistanceField(std::ostream & sout,
                                   std::istream& sinput);
    void parseAddFieldFromObsArray(std::ostream & sout,
                                   std::istream& sinput);
    void parsePoint( std::istream & sinput, chomp::MatX & point);
    void parseExecute( std::ostream & sout , std::istream & sinput );
    void parseRobot( std::string & name );
    ORTSRConstraint * parseTSR( std::istream & sinput );

    
    ///////////////////////////////////////////////////////////////////
    /////The below functions are in the file orchomp_mod_utils.cpp/////
    ///////////////////////////////////////////////////////////////////
public:
    // A small helper function for creating a straight
    //  line trajectory between two endpoints:
    void createInitialTrajectory();
    //get the ith state in the trajectory matrix and turn it into an openrave
    //  state vector.
    void getIthStateAsVector( size_t i,
                std::vector< OpenRAVE::dReal > & state  );
    //get a state matrix, and turn it into an openrave state vector
    void getStateAsVector( const chomp::MatX & state,
                std::vector< OpenRAVE::dReal > & vec  );
    //get a random state that is within the robot's joint limits
    void getRandomState( chomp::MatX & vec );

    //turn an OpenRAVE::Vector to a chomp::MatX.
    void vectorToMat(const std::vector< OpenRAVE::dReal > & vec,
                             chomp::MatX & mat );

    //take a state matrix, and if any of the values exceed the joint limits,
    //  then clamp that value to the nearest joint limit
    void clampToLimits( chomp::MatX & state);
    //returns true if the given state matrix is within the joint limits
    bool isWithinLimits( const chomp::MatX & mat ) const;
    bool isWithinPaddedLimits( const chomp::MatX & mat ) const;
    //Take a state matrix, and set the robot's active DOF values
    void setActiveDOFValues( const chomp::MatX & qt );
    //get the ik for the currently active end effector for the transform
    void getIK( const OpenRAVE::Transform & xform, 
                     std::vector< OpenRAVE::dReal > & solution );



    //print out the trajectory 
    void coutTrajectory() const;
    //Checks to see if all of the points in the trajectory are within the
    //  joint limits. Print out the status of each point.
    void isTrajectoryWithinLimits() const;
    
    OpenRAVE::KinBodyPtr createBox( 
                                const OpenRAVE::Vector & pos,
                                const OpenRAVE::Vector & extents,
                                const OpenRAVE::Vector & color,
                                float transparency = 0);

    //Returns true if two joints are adjacent
    bool areAdjacent( int first, int second ) const ;
    

};



} /* namespace orchomp */
#endif 
