/** \file orchomp_mod.cpp
 * \brief Implementation of the orchomp module, an implementation of CHOMP
 *        using libcd.
 * \author Christopher Dellin
 * \date 2012
 */

/* (C) Copyright 2012-2013 Carnegie Mellon University */

/* This module (orchomp) is part of libcd.
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

#include <time.h>
#include "orchomp_mod.h"


#include "orchomp_kdata.h"

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


//constructor that registers all of the commands to the openRave
//   command line interface.
mod::mod(OpenRAVE::EnvironmentBasePtr penv) :
    OpenRAVE::ModuleBase(penv), environment( penv ), 
    factory( NULL ), robot( NULL ), sphere_collider( NULL ),
    collisionHelper( NULL ), chomper( NULL )
{
      __description = "orchomp: implementation multigrid chomp";
      RegisterCommand("viewspheres",
              boost::bind(&mod::viewspheres,this,_1,_2),
              "view spheres");
      RegisterCommand("computedistancefield",
               boost::bind(&mod::computedistancefield,this,_1,_2),
               "compute distance field");
      RegisterCommand("addfield_fromobsarray",
            boost::bind( &mod::addfield_fromobsarray,this,_1,_2),
            "compute distance field");
      RegisterCommand("create", 
            boost::bind(&mod::create,this,_1,_2),
            "create a chomp run");
      RegisterCommand("iterate",
            boost::bind(&mod::iterate,this,_1,_2),
            "create a chomp run");
      RegisterCommand("gettraj",
            boost::bind(&mod::gettraj,this,_1,_2),
            "create a chomp run");
      RegisterCommand("destroy",
            boost::bind(&mod::destroy,this,_1,_2),
            "create a chomp run");
      
}

/* ======================================================================== *
 * module commands
 */

// NOTES : viewspheres looks to be fine without editing.
//    - There may be hidden dependency issues
int mod::viewspheres(std::ostream& sout, std::istream& sinput)
{
    parseViewSpheres( sout,  sinput);

    return 0;
}



/* computedistancefield robot Herb2
 * computes a distance field in the vicinity of the passed kinbody
 * 
 * note: does aabb include disabled bodies? probably, but we might hope not ...
 * */

 // NOTES : this is likely to not work for several reasons:
 //         1. It is lacking the necessary libraries for computing
 //         2. Even if the libraries were correct, it is unlikely to work
 //             with the current chomp style of gradients
int mod::computedistancefield(std::ostream& sout, std::istream& sinput)
{
    parseComputeDistanceField( sout,  sinput);
    return 0;
}


int mod::addfield_fromobsarray(std::ostream& sout, std::istream& sinput)
{
    parseAddFieldFromObsArray( sout,  sinput);
   
    return 0;

}

double SphereCollisionHelper::getCost(const chomp::MatX& q,
                                              size_t body_index,
                                              chomp::MatX& dx_dq,
                                              chomp::MatX& cgrad)
{

    //resize the matrices:
    dx_dq.conservativeResize( nwkspace, ncspace );
    cgrad.conservativeResize( ncspace, 1 );

    std::vector< OpenRAVE::dReal > jacobian;

    double cost = 0.0;
    // this should be correct:q.data();
    
    //NOTE : remove the following test code from the final.
    //TODO REMOVE THIS !!!
    double * pose_ptr = (double*)q.data();
    for ( int i = 0; i < q.size(); i ++ ){ assert( q( i ) == pose_ptr[i] ); }

    //TODO make sure that the bounds on this are valid
    //enable the eigen array format to interface with the current position, q
    std::vector<OpenRAVE::dReal> vec( pose_ptr, pose_ptr + q.size()  );

    robot->SetActiveDOFValues(vec);
    
    //loop through all of the active spheres, and check their collision status.
    for (int sphere_index = 0; sphere_index < n_active; sphere_index ++ )
    {
        //extract the current sphere
        const Sphere & current_sphere = spheres[ sphere_index ];
        //get the transformation of the body that the sphere is on.
        const OpenRAVE::Transform t = current_sphere.robot_link->GetTransform();
        //get the transformation from the body to the sphere.
        const OpenRAVE::Vector v = t * OpenRAVE::Vector( current_sphere.pose );
        
        //calculate the jacobians - this is the jacobian of the workspace ... i think
        robot->CalculateJacobian( current_sphere.robot_linkindex, v, jacobian);
        
        //NOTE : I have the jacobian for each of the spheres, how do I transform this into
        // a cost, a jacobian, and a gradient?
        
        // TODO : Make this average the spheres and get stuff done.


    } 

    

    return cost;
}



/* runchomp robot Herb2
 * run chomp from the current config to the passed goal config
 * uses the active dofs of the passed robot
 * initialized with a straight-line trajectory
 * */
int mod::create(std::ostream& sout, std::istream& sinput)
{
    std::cout << "Creating" << std::endl;

    parseCreate( sout,  sinput);
    
    //after the arguments have been collected, pass them to chomp
    createInitialTrajectory();
    
    //now that we have a trajectory, make a chomp object
    chomper = new chomp::Chomp( factory, trajectory,
                                q0, q1, info.n_max, 
                                info.alpha, info.obstol,
                                info.max_global_iter,
                                info.max_local_iter,
                                info.t_total);

     
    //TODO Compute signed distance field

    //TODO setup collision geometry

    //TODO setup momentum stuff ?? maybe ?? 

    //TODO give chomp access to the distance field, and to the
    //  collision helper

    //give chomp a collision helper to 
    //deal with collision and gradients.
    if (sphere_collider){
        collisionHelper = new chomp::ChompCollGradHelper(
                                sphere_collider,info.gamma);
        chomper->ghelper = collisionHelper;
    }

    if ( info.doObserve ){
        chomper->observer = &observer;
    }
    
    std::cout << "Done Creating" << std::endl;
    
    return 0;
}

int mod::iterate(std::ostream& sout, std::istream& sinput)
{
    std::cout << "Iterating" << std::endl;
    
    parseIterate( sout,  sinput);
    chomper->solve( info.doGlobal, info.doLocal );

    std::cout << "Done Iterating" << std::endl;
    return 0;
}

int mod::gettraj(std::ostream& sout, std::istream& sinput)
{
    
    std::cout << "Getting Trajectory" << std::endl;
    OpenRAVE::EnvironmentMutex::scoped_lock lockenv;
    OpenRAVE::TrajectoryBasePtr traj_ptr;

    parseGetTraj( sout,  sinput);
    
    //get the lock for the environment
    lockenv = OpenRAVE::EnvironmentMutex::scoped_lock(
              environment->GetMutex() );
   
    //setup the openrave trajectory pointer to receive the
    //  found trajectory.
    traj_ptr = OpenRAVE::RaveCreateTrajectory(environment);
    traj_ptr->Init(
        robot->GetActiveConfigurationSpecification());
    debugStream << "Done locking environment,"
                << " begun extracting trajectory" << std::endl;
    
    
    //For every state (each row is a state), extract the data,
    //  and turn it into a vector, then give it to 
    // TODO : make sure that this pointer arithmetic
    //  actually works.

    for ( int i = 0; i < trajectory.rows(); i ++ ){
        traj_ptr->Insert( i, getIthStateAsVector( i ) );
    }

    //this is about changing the timing or something
    /*
       OpenRAVE::planningutils::RetimeActiveDOFTrajectory
                (traj_ptr, boostrobot, false, 0.2, 0.2,
                 "LinearTrajectoryRetimer","");
    */
    
    //TODO : check for collisions

    debugStream << "Serializing trajectory output" << std::endl; 
    //serialize the trajectory and send it over the 
    //  output stream.
    traj_ptr->serialize( sout );


    std::cout << "Done Getting Trajectory" << std::endl;
    return 0;
}




int mod::destroy(std::ostream& sout, std::istream& sinput){
    
    if (chomper){ delete chomper; }
    if (sphere_collider){ delete sphere_collider; }
    if (factory){ delete factory;}
    if (collisionHelper){ delete collisionHelper; }

    return 0;
}



//takes the two endpoints and fills the trajectory matrix by
//  linearly interpolating betweeen the two.
inline void mod::createInitialTrajectory()
{
    assert( info.n != 0 );
    assert( q0.size() == q1.size() );

    trajectory.resize(info.n, q0.size() );

    //fill the trajectory matrix 
    for (size_t i=0; i<info.n; ++i) {
        trajectory.row(i) = (i+1)*(q1-q0)/(info.n+1) + q0;
    }
}

inline std::vector< OpenRAVE::dReal >
mod::getIthStateAsVector( size_t i )
{

    const int width = trajectory.cols();
    
    return  std::vector< OpenRAVE::dReal > ( 
                trajectory.row( i ).data(), 
                trajectory.row( i ).data() + width );

}


} /* namespace orchomp */


