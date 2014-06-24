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
#include <openrave/openrave.h>
#include <openrave/planningutils.h>
#include <time.h>

extern "C" {
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "utils/os.h"
}

#define DEBUG_COUT 0
#define debugStream \
    if (DEBUG_COUT) {} \
    else std::cout

#define DEBUG_TIMING

#ifdef DEBUG_TIMING
#  define TIC() { struct timespec tic; struct timespec toc; clock_gettime(CLOCK_THREAD_CPUTIME_ID, &tic);
#  define TOC(tsptr) clock_gettime(CLOCK_THREAD_CPUTIME_ID, &toc); CD_OS_TIMESPEC_SUB(&toc, &tic); CD_OS_TIMESPEC_ADD(tsptr, &toc); }
#else
#  define TIC()
#  define TOC(tsptr)
#endif


namespace orchomp
{
class mod;

//this is a structure used to hold and initialize values
//  for an eventual call to chomp.
class ChompInfo {
public:
    // al: alpha
    // obstol: Relative error of the objective function - 
    //          once the error gets to a low enough point,
    //          reltive to the objective function, quit chomping.
    // t_total: the total amount of timesteps that
    //          N+1 timesteps will take            
    // gamma: how much the gradient term counts for.
    double alpha, obstol, t_total, gamma, epsilon, epsilon_self, obs_factor,
           obs_factor_self, jointPadding;

    //n: the initial size of the trajectory,
    //n_max: the final size,
    //max_global_iter: the max # of global chomp interations,
    //max_local_iter: the max # of local smoothing iterations
    size_t n, n_max, min_global_iter, max_global_iter,
                     min_local_iter, max_local_iter;

    //whether or not global and/or local chomp should
    //  be done.
    bool doGlobal, doLocal, doObserve, noFactory, noCollider, 
         noSelfCollision, noEnvironmentalCollision;

    //a basic constructor to initialize values
    ChompInfo() :
        alpha(0.1), obstol(0.01), t_total(1.0), gamma(0.1),
        epsilon( 0.1 ), epsilon_self( 0.01 ), obs_factor( 200 ),
        obs_factor_self( 5 ), jointPadding( 0.05 ),
        n(0), n_max(0),
        min_global_iter( 0 ), max_global_iter( size_t(-1) ), 
        min_local_iter( 0 ), max_local_iter( size_t(-1)), doGlobal( true ),
        doLocal( false), doObserve( false ), noFactory (false),
        noCollider( false ), noSelfCollision( false ),
        noEnvironmentalCollision( false )
        {}
};


class SphereCollisionHelper : public chomp::ChompCollisionHelper{
public:
    
    // a pointer to the module for acces to stuff like the collision
    //  geometry
    mod * module;

    //the first <active_spheres> amount of spheres in the 
    // above vector are active.
    int n_active;

    /* obstacle parameters */
    //environmental collisions
    double epsilon;
    double obs_factor;
   
    //self-collisions
    double epsilon_self;
    double obs_factor_self;
    
    //________________________Public Member Functions____________________//
    
    //the constuctor needs a pointer to the robot in addition to the spaces.
    SphereCollisionHelper( size_t ncspace, size_t nwkspace, size_t nbodies, 
                          mod * module) :
            ChompCollisionHelper( ncspace, nwkspace, nbodies ), module(module),
            epsilon( 0.1 ), obs_factor(200.0), epsilon_self( 0.01 ), 
            obs_factor_self( 5.0 )
    {
    }

    OpenRAVE::KinBodyPtr createCube( double cost,
                                    double size,
                                     OpenRAVE::EnvironmentBasePtr & env,
                                    const OpenRAVE::Vector & pos
                                         );
    OpenRAVE::KinBodyPtr createCube( 
                                 const OpenRAVE::Vector & color,
                                 double size,
                                 OpenRAVE::EnvironmentBasePtr & env,
                                 const OpenRAVE::Vector & pos
                                 );
    void visualizeSDFSlice( size_t sdf_index, size_t axis,
                            size_t slice_index, double time);

    // q - the current configuration
    // body_index - the index of the body to get the gradient, cost and jacobains
    //              for.
    // dx_dq - jacobain of workspace position
    //         dims: nwkspace-by-ncspace
    // cgrad - gradient (jacobian transpose of cost wrt workspace position
    //         dims : ncspace-by-1
    virtual double getCost(const chomp::MatX& q, size_t body_index,
                           chomp::MatX& dx_dq,  chomp::MatX& cgrad); 

private:

    OpenRAVE::dReal getSDFCollisions( const Sphere & sphere,
                                      const OpenRAVE::Vector & position, 
                                      vec3 & gradient );
    OpenRAVE::dReal getSelfCollisions( size_t body_index,
                                       const Sphere & current_sphere,
                                       const OpenRAVE::Vector & position, 
                                       vec3 & gradient,
                                       std::vector<OpenRAVE::dReal> & other);


};



/* the module itself */
class mod : public OpenRAVE::ModuleBase
{
public:

    //____________________PUBLIC MEMBER VARIABLES____________________//
    OpenRAVE::EnvironmentBasePtr environment; /* filled on module creation */
 
    // NOTE : Most of the below variables were in the run structure.

    //the trajectory, start, and endpoint.
    chomp::MatX trajectory, q0, q1;
    ORConstraintFactory * factory;
   
    //This holds basic info relating to an individual 
    //   run of chomp
    ChompInfo info;
  
    // This observes chomp for general debugging purposes
    chomp::DebugChompObserver observer;

    //these are useful for interfacing with the robot.
    OpenRAVE::RobotBasePtr robot;
    std::string robot_name;

    OpenRAVE::RobotBase::ManipulatorPtr active_manip;

    OpenRAVE::KinBodyPtr kinbody;
    std::vector< int > active_indices;
    size_t n_dof;

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
    
    std::vector< ORTSRConstraint * > tsrs;

    //this is a pointer to the chomp class that will pull most of the
    //   weight.
    chomp::Chomp * chomper;
    
    OpenRAVE::TrajectoryBasePtr trajectory_ptr;

    //_______________________PUBLIC MEMBER FUNCTIONS___________________//
    //visualize the collision geometry 
    int viewspheres(std::ostream & sout, std::istream& sinput);

    //compute the distance field for use in collision detection, and
    //   descending the gradient out of collision
    int computedistancefield(std::ostream & sout, std::istream& sinput);
    
    int visualizeslice(std::ostream& sout, std::istream& sinput);

    // NOTE : Find out what this is supposed to do
    int addfield_fromobsarray(std::ostream & sout, std::istream& sinput);

    //
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
    
    int viewtsr(std::ostream & sout, std::istream& sinput);
    
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

    // A small helper function for creating a straight line trajectory between
    //  two endpoints:
    inline void createInitialTrajectory();
    //get the ith state in the trajectory matrix and turn it into an openrave
    //  state vector.
    inline void getIthStateAsVector( size_t i,
                std::vector< OpenRAVE::dReal > & state  );
    //get a state matrix, and turn it into an openrave state vector
    inline void getStateAsVector( const chomp::MatX & state,
                std::vector< OpenRAVE::dReal > & vec  );
    //get a random state that is within the robot's joint limits
    inline void getRandomState( chomp::MatX & vec );
    //take a state matrix, and if any of the values exceed the joint limits,
    //  then clamp that value to the nearest joint limit
    inline void clampToLimits( chomp::MatX & state);
    //returns true if the given state matrix is within the joint limits
    bool isWithinLimits( const chomp::MatX & mat ) const;
    bool isWithinPaddedLimits( const chomp::MatX & mat ) const;
    //Take a state matrix, and set the robot's active DOF values
    inline void setActiveDOFValues( const chomp::MatX & qt );
    
    //Get the collision geometry from the XML files
    void getSpheres();

    //print out the trajectory 
    void coutTrajectory() const;
    //Checks to see if all of the points in the trajectory are within the
    //  joint limits. Print out the status of each point.
    void isTrajectoryWithinLimits() const;
    
    //Returns true if two joints are adjacent
    bool areAdjacent( int first, int second ) const ;


};



inline void vectorToMat(const std::vector< OpenRAVE::dReal > & vec,
                             chomp::MatX & mat )
{
    assert( vec.size() > 0 );
    mat.resize(1, vec.size() );

    for ( size_t i = 0; i < vec.size() ; i ++ ){ mat(i) = vec[i]; }
}


inline void mod::getStateAsVector( const chomp::MatX & state,
                                   std::vector< OpenRAVE::dReal > & vec ){

    vec.resize( n_dof );
    assert( state.size() == int( n_dof ) );
    
    for ( size_t i = 0; i < n_dof; i ++ ){
        vec[i] = state(i);
    }
}


inline void mod::getIthStateAsVector( size_t i, 
                      std::vector< OpenRAVE::dReal > & state )
{
    
    const int width = trajectory.cols();
    state.resize( width );

    for ( int j = 0; j < width; j ++ ){
        state[j] = trajectory(i, j );
    }

}

inline void mod::setActiveDOFValues( const chomp::MatX & qt ){

    std::vector< OpenRAVE::dReal > vec;
    getStateAsVector( qt, vec );

    robot->SetActiveDOFValues( vec );
}



OpenRAVE::KinBodyPtr createBox( const OpenRAVE::Vector & pos,
                                const OpenRAVE::Vector & extents,
                                const OpenRAVE::Vector & color,
                                OpenRAVE::EnvironmentBasePtr & env,
                                float transparency = 0)
{

    //create a cube to be used for collision detection in the world.
    //  create an object and name that object 'cube'
    OpenRAVE::KinBodyPtr cube = RaveCreateKinBody( env);
    
    std::stringstream ss;
    ss   << pos[0] << "_"
         << pos[1] << "_"
         << pos[2] ; 

    const std::string name = ss.str();

    if( env->GetKinBody( name.c_str() ).get() ){return cube; }
    cube->SetName( name.c_str() );

    //set the dimensions of the cube object
    std::vector< OpenRAVE::AABB > vaabbs(1);

    /* extents = half side lengths */
    vaabbs[0].extents = extents;
    vaabbs[0].pos = pos;
    cube->InitFromBoxes(vaabbs, true);
    
    //add the cube to the environment
    env->Add( cube );

    cube->GetLinks()[0]->GetGeometries()[0]->SetAmbientColor( color );
    cube->GetLinks()[0]->GetGeometries()[0]->SetDiffuseColor( color );
    cube->GetLinks()[0]->GetGeometries()[0]->SetTransparency(transparency);
    
    return cube;

}


} /* namespace orchomp */
#endif 
