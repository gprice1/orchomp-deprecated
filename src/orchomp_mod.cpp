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

#include "orchomp_mod.h"
#include "orchomp_constraint.h"
#include "orchomp_kdata.h"
#include "orchomp_collision.h"
#include "orchomp_observer.h"

#include "chomp-multigrid/chomp/HMC.h"

namespace orchomp
{


//constructor that registers all of the commands to the openRave
//   command line interface.
mod::mod(OpenRAVE::EnvironmentBasePtr penv) :
    OpenRAVE::ModuleBase(penv), environment( penv ),
    chomper( NULL ),
    factory( NULL ), sphere_collider( NULL ),
    collisionHelper( NULL ), observer( NULL ), hmc( NULL )
{
    RAVELOG_INFO( "Constructing\n");
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
       RegisterCommand("execute",
            boost::bind(&mod::execute,this,_1,_2),
            "play a trajectory on a robot");
       RegisterCommand("playback",
            boost::bind(&mod::playback,this,_1,_2),
            "playback a trajectory on a robot");
       RegisterCommand("visualizeslice",
            boost::bind(&mod::visualizeslice,this,_1,_2),
            "playback a trajectory on a robot");
       RegisterCommand("addtsr",
            boost::bind(&mod::addtsr,this,_1,_2),
            "playback a trajectory on a robot");
       RegisterCommand("createtsr",
            boost::bind(&mod::createtsr,this,_1,_2),
            "playback a trajectory on a robot");
       RegisterCommand("removeconstraint",
            boost::bind(&mod::removeconstraint,this,_1,_2),
            "remove a tsr constraint");
       RegisterCommand("viewtsr",
            boost::bind(&mod::viewtsr,this,_1,_2),
            "playback a trajectory on a robot");
}

/* ======================================================================== *
 * module commands
 */


bool mod::playback(std::ostream& sout, std::istream& sinput)
{
    
    double time = -1;
    if ( !sinput.eof() ){
        sinput >> time;
    }
    
    if (time > 10 || time < 0.00001 ){
        time = 0.07;
    }
    
    if ( !chomper ){
        RAVELOG_ERROR( "Iterate must be called before a playback" );
        return false;
    }

    for ( int i = 0; i < chomper->xi.rows(); i ++ ){

        std::vector< OpenRAVE::dReal > vec;
        getStateAsVector( chomper->xi.row(i), vec );
        
        viewspheresVec( chomper->xi.row(i), vec, time );

    }

    return true;
}

void getColor( OpenRAVE::Vector & color, double cost, double eps ){
    if ( cost <= 0.5*eps ){
        float val = cost/(0.5*eps);
        color = OpenRAVE::Vector( 1,val,0 );
    }
    else if ( cost <= eps){
        float val = cost/(0.5*eps) - 1.0;
        color = OpenRAVE::Vector( 1,1,val);
    }else {
        color = OpenRAVE::Vector( 1,1,1 );
    }
}

//view the collision geometry. .
bool mod::viewspheresVec(const chomp::MatX & q,
                        const std::vector< OpenRAVE::dReal > & vec, 
                        double time)
{
    
    robot->SetActiveDOFValues( vec, true );

    if ( !sphere_collider ) { 
        std::cout << "There is no sphere collider, so viewing the" 
                  << " spheres is impossible" << std::endl;
        return true;
    }
    
    char text_buf[1024];
    
    const size_t n_active = active_spheres.size();

    std::vector< OpenRAVE::KinBodyPtr > bodies;
    

    for ( size_t i=0; i < n_active; i ++ )
    {

        double cost( 0 );
        chomp::MatX dxdq, cgrad;

        if( sphere_collider ){ 
            cost = sphere_collider->getCost( q, i, dxdq, cgrad );
            if ( cost <= 0 ){ continue; }
        }
        
        //extract the current sphere
        //  if i is less than n_active, get a sphere from active_spheres,
        //  else, get a sphere from inactive_spheres
        const Sphere & current_sphere = active_spheres[i] ;

        //make a kinbody sphere object to correspond to this sphere.
        OpenRAVE::KinBodyPtr sbody = 
                OpenRAVE::RaveCreateKinBody( environment );
        sprintf( text_buf, "orcdchomp_sphere_%d", int(i));
        sbody->SetName(text_buf);
        
        //set the dimensions and transform of the sphere.
        std::vector< OpenRAVE::Vector > svec;

        //get the position of the sphere in the world 
        OpenRAVE::Transform t =  current_sphere.body->GetLink(
                                 current_sphere.linkname)->GetTransform();
        OpenRAVE::Vector position = t * current_sphere.position;
        
        //set the radius of the sphere
        position.w = current_sphere.radius *
                     ( 1 + std::min(cost/(sphere_collider->epsilon ), 1.0)); 

        //give the kinbody the sphere parameters.
        svec.push_back( position );
        sbody->InitFromSpheres(svec, true);
        
        bodies.push_back( sbody );
        environment->Add( sbody );

        OpenRAVE::Vector color;  
        getColor( color, cost, sphere_collider->epsilon );
        sbody->GetLinks()[0]->GetGeometries()[0]->SetAmbientColor( color );
        sbody->GetLinks()[0]->GetGeometries()[0]->SetDiffuseColor( color );
    }
    
    timer.wait( time );

    for ( size_t i = 0; i < bodies.size() ; i ++ ){
        environment->Remove( bodies[i] );
    }
    
    return true;
}

//view the collision geometry.
bool mod::viewspheres(std::ostream& sout, std::istream& sinput)
{

    
    OpenRAVE::EnvironmentMutex::scoped_lock lockenv(environment->GetMutex());
    parseViewSpheres( sout,  sinput);
    
    //if there are no spheres, get them
    if ( active_spheres.size() + inactive_spheres.size() == 0) {
        getSpheres() ; 
    }
     
    char text_buf[1024];
    
    const size_t n_active = active_spheres.size();
    const size_t n_inactive = inactive_spheres.size();

    for ( size_t i=0; i < n_active + n_inactive ; i ++ )
    {
        //extract the current sphere
        //  if i is less than n_active, get a sphere from active_spheres,
        //  else, get a sphere from inactive_spheres
        const Sphere & sphere = (i < n_active ?
                                 active_spheres[i] :
                                 inactive_spheres[i-n_active]);

        //make a kinbody sphere object to correspond to this sphere.
        OpenRAVE::KinBodyPtr sbody = 
                OpenRAVE::RaveCreateKinBody( environment );
        ;
        sbody->SetName(text_buf);
        

        //set the dimensions and transform of the sphere.
        std::vector< OpenRAVE::Vector > svec;

        //get the position of the sphere in the world 
        OpenRAVE::Transform t = 
                sphere.body->GetLink(sphere.linkname)->GetTransform();
        OpenRAVE::Vector v = t * sphere.position;
        
        //set the radius of the sphere
        v.w = sphere.radius; 

        //give the kinbody the sphere parameters.
        svec.push_back(v);
        sbody->InitFromSpheres(svec, true);

        environment->Add( sbody );
    }
    return true;
}


bool mod::removeconstraint(std::ostream& sout, std::istream& sinput){

    size_t index = 0;

    if ( !sinput.eof() ){
        sinput >> index;
    }

    //if a factory exists, remove the constraint
    if ( factory ){
        factory->removeConstraint( index );
    }

    else {
        std::cout << "There is no factory, so removal cannot happen" <<
                  std::endl;
    }

    return true;
}

bool mod::createtsr( std::ostream& sout, std::istream& sinput){
    
    ORTSRConstraint * c = parseTSR( sinput );
    
    //default to a trajectory wide constraint.
    double starttime(0), endtime(1);

    while ( !sinput.eof() ){
        std::string cmd;
        sinput >> cmd;

        if (!sinput.eof()){ break; }
        if ( cmd == "starttime" ){
            sinput >> starttime;
        }else if ( cmd == "endtime" ){
            sinput >> endtime;
        }
    }
    
    //add the constraint to the factory, and to the list of constraints.
    factory->addConstraint( c, starttime, endtime );
    tsrs.push_back( c );

    return true;
}

bool mod::addtsr(std::ostream& sout, std::istream& sinput){
    

    OpenRAVE::EnvironmentMutex::scoped_lock
                                    lockenv(environment->GetMutex());

    chomp::Transform::quat pose_0_w_rot, pose_w_e_rot;
    chomp::Transform::vec3 pose_0_w_trans, pose_w_e_trans;
    chomp::MatX bounds(6,2);
    OpenRAVE::dReal starttime, endtime;

    while ( !sinput.eof() ){
        std::string cmd;
        sinput >> cmd;
        RAVELOG_INFO( "\tExecuting SubCommand: %s\n", cmd.c_str() );
        if ( cmd == "pose_0_w" ){
            sinput >> pose_0_w_trans[0];
            sinput >> pose_0_w_trans[1];
            sinput >> pose_0_w_trans[2];
            sinput >> pose_0_w_rot[0];
            sinput >> pose_0_w_rot[1];
            sinput >> pose_0_w_rot[2];
        }
        else if ( cmd == "pose_w_e" ){
            sinput >> pose_w_e_trans[0];
            sinput >> pose_w_e_trans[1];
            sinput >> pose_w_e_trans[2];

            sinput >> pose_w_e_rot[0];
            sinput >> pose_w_e_rot[1];
            sinput >> pose_w_e_rot[2];
        }
        else if ( cmd == "bounds" ){
            sinput >> bounds( 0, 0); sinput >> bounds( 0, 1);
            sinput >> bounds( 1, 0); sinput >> bounds( 1, 1);
            sinput >> bounds( 2, 0); sinput >> bounds( 2, 1);
            sinput >> bounds( 3, 0); sinput >> bounds( 3, 1);
            sinput >> bounds( 4, 0); sinput >> bounds( 4, 1);
            sinput >> bounds( 5, 0); sinput >> bounds( 5, 1);
        }
        else if (cmd == "starttime"){ 
            sinput >> starttime;
            if (starttime < 0 || starttime > 1 ){
                RAVELOG_ERROR("Starttime is out of bounds" );
            }
        }
        else if (cmd == "endtime"  ){
            sinput >> endtime;
            if (endtime < 0 || endtime > 1 ){
                RAVELOG_ERROR("Endtime is out of bounds" );
            }
        }
        else if (cmd == "trajwide" || 
                 cmd == "trajectorywide" ){
            starttime = 0; 
            endtime = 1;
        }
        else {
            while ( !sinput.eof() ){
                sinput >> cmd;
                RAVELOG_ERROR("argument %s not known!\n",
                              cmd.c_str() );
            }
            throw OpenRAVE::openrave_exception("Bad arguments!");
        }

    }



    chomp::Transform pose_0_w(pose_0_w_rot, pose_0_w_trans);
    chomp::Transform pose_w_e(pose_w_e_rot, pose_w_e_trans);

    if ( !active_manip.get() ){
        RAVELOG_ERROR( "There is no active manip, needed for addtsr()" );
        throw OpenRAVE::openrave_exception("Bad arguments!");
    }
    
    std::string link_name = active_manip->GetEndEffector()->GetName();
    std::string body_name = robot->GetName();
    std::cout << "Creating TSR" << std::endl;
    ORTSRConstraint * c = new ORTSRConstraint( this, 
                                               -1, 
                                               pose_0_w,
                                               bounds, pose_w_e,
                                               body_name, link_name);
    std::cout << "Done creating TSR " << c << std::endl;
    factory->addConstraint( c, starttime, endtime );
    
    tsrs.push_back( c );
    
    std::cout << "Done with tsr parsing" << std::endl;
    return true;
}


bool mod::viewtsr(std::ostream & sout, std::istream& sinput){
     OpenRAVE::EnvironmentMutex::scoped_lock
                                    lockenv(environment->GetMutex());

    size_t index = 0;
    if ( !sinput.eof() ){
        sinput >> index;
    }

    assert( tsrs.size() > index );

    //TSRs will be green because that is the color that they are.
    OpenRAVE::Vector color( 0,1,0);
    OpenRAVE::Vector size( 
                tsrs[index]->_Bw(0,1) - tsrs[index]->_Bw(0,0),
                tsrs[index]->_Bw(1,1) - tsrs[index]->_Bw(1,0),
                tsrs[index]->_Bw(2,1) - tsrs[index]->_Bw(2,0) );
    OpenRAVE::Vector position( 0,0,0 );

    debugStream << "Making cube" << std::endl;
    OpenRAVE::KinBodyPtr cube = createBox( position, size, color, 0.5 );
    
    chomp::Transform::quat rot = tsrs[index]->_pose_0_w.rotation();
    chomp::Transform::vec3 pos = tsrs[index]->_pose_0_w.translation();

    debugStream << "Setting OR positions" <<
                "\trot: " << rot << "\tpos: " << pos << std::endl;
    OpenRAVE::Vector or_rot( rot[0], rot[1], rot[2], rot[3]);
    OpenRAVE::Vector or_pos( pos[0], pos[1], pos[2] );
    
    debugStream << "Making Transform" << std::endl;
    debugStream << "OR displacements: "  <<
                "\trot: " << or_rot << "\tpos: " << or_pos << std::endl;

    OpenRAVE::Transform xform( or_rot, or_pos );

    debugStream << "Setting xform" << std::endl;
    cube->SetTransform( xform );

    debugStream << "Done" << std::endl;

    return true;
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
bool mod::computedistancefield(std::ostream& sout, std::istream& sinput)
{
    
    //lock the environment
    OpenRAVE::EnvironmentMutex::scoped_lock lock(environment->GetMutex());
    
    //Parse the arguments
    parseComputeDistanceField( sout,  sinput);
    //currently, parse compute does all of the actual work for this function.
    //  That seems peculiar.
    

    return true;
}


bool mod::visualizeslice(std::ostream& sout, std::istream& sinput)
{
    
    //lock the environment
    //OpenRAVE::EnvironmentMutex::scoped_lock lock(environment->GetMutex());
    
    size_t sdf_index(0), axis(2), slice_index(20);
    double time( 3 );
    bool getwhole = false;

    while ( !sinput.eof() ){
        std::string cmd;
        sinput >> cmd;

        if( cmd == "sdf" ){ sinput >> sdf_index; }
        if( cmd == "axis" ){ sinput >> axis; }
        if( cmd == "slice" ){ sinput >> slice_index; }
        if( cmd == "time" ){ sinput >> time; }
        if( cmd == "getwhole" ){ getwhole = true; }
    }
    
    if( sphere_collider ){
        if (getwhole ){
            for ( size_t i = 0; i < sdfs[sdf_index].grid.dims()[axis]; i ++ ){
                sphere_collider->visualizeSDFSlice( sdf_index, axis,
                                                    i, time );
            }
        }else{
            sphere_collider->visualizeSDFSlice( sdf_index, axis,
                                            slice_index, time );
        }
    }

    return true;
}


bool mod::addfield_fromobsarray(std::ostream& sout, std::istream& sinput)
{
    parseAddFieldFromObsArray( sout,  sinput);
   
    return true;

}




/* runchomp robot Herb2
 * run chomp from the current config to the passed goal config
 * uses the active dofs of the passed robot
 * initialized with a straight-line trajectory
 * */
bool mod::create(std::ostream& sout, std::istream& sinput)
{
    //get the lock for the environment
    OpenRAVE::EnvironmentMutex::scoped_lock lockenv(
              environment->GetMutex() );
    
    RAVELOG_INFO( "Creating" );

    parseCreate( sout,  sinput);
    

    //create a padded upper and lower joint limits vectors. These will be
    //  Used for constraining the trajectory to the joint limits.
    //  Since sometimes, chomp allows constraints to be minutely off,
    //  this will prevent that from happening.
    paddedUpperJointLimits.resize( upperJointLimits.size() );
    paddedLowerJointLimits.resize( lowerJointLimits.size() );

    for ( size_t i = 0; i < upperJointLimits.size(); i ++ ){
        OpenRAVE::dReal interval = upperJointLimits[i] - lowerJointLimits[i];
        interval *= info.jointPadding;

        //fill the padded joint limits
        paddedUpperJointLimits[i] = upperJointLimits[i] - interval;
        paddedLowerJointLimits[i] = lowerJointLimits[i] + interval;
    }
    
    clampToLimits( q0 );
    clampToLimits( q1 );

    assert( isWithinLimits( q0 ) );
    assert( isWithinLimits( q1 ) );
    
    //create the sphere collider to actually 
    //  use the sphere information, and pass in its parameters.
    if ( !info.noCollider && !sphere_collider ){
        
        //Setup the collision geometry
        getSpheres();

        //create a collider object to use the geometry
        sphere_collider = new SphereCollisionHelper(
                          n_dof, 3, active_spheres.size(), this);
        sphere_collider->epsilon = info.epsilon;
        sphere_collider->epsilon_self = info.epsilon_self;
        sphere_collider->obs_factor = info.obs_factor;
        sphere_collider->obs_factor_self = info.obs_factor_self;
    }

    if ( !info.noFactory ){
        if (factory ){ delete factory; }
        factory = new ORConstraintFactory( this );
    }

    RAVELOG_INFO( "Done Creating" );
    
    return true;
}


void mod::printChompInfo(){


    RAVELOG_INFO( "Chomp.n = %d\n", info.n );
    RAVELOG_INFO( "Chomp.n_max = %d\n", info.n_max );
    RAVELOG_INFO( "Chomp.alpha = %f\n", info.alpha );
    RAVELOG_INFO( "Chomp.obstol = %f\n", info.obstol );
    RAVELOG_INFO( "Chomp.max_global_iter = %d\n", info.max_global_iter );
    RAVELOG_INFO( "Chomp.min_global_iter = %d\n", info.min_global_iter );
    RAVELOG_INFO( "Chomp.min_local_iter = %d\n", info.min_local_iter );
    RAVELOG_INFO( "Chomp.max_local_iter = %d\n", info.max_local_iter );
    RAVELOG_INFO( "Chomp.t_total = %f\n", info.t_total );
    RAVELOG_INFO( "Chomp.max_time = %f\n", info.timeout_seconds );

    std::stringstream ss;
    std::string configuration;

    ss << q0;
    configuration = ss.str();
    RAVELOG_INFO( "Chomp q0 = %s\n", configuration.c_str() );

    //clear the string stream.
    ss.str( std::string() ) ;
    ss << q1;
    configuration = ss.str();

    RAVELOG_INFO( "Chomp.q1 = %s\n", configuration.c_str() );
}

bool mod::iterate(std::ostream& sout, std::istream& sinput)
{
    RAVELOG_INFO( "Iterating\n" );

    //get the arguments
    parseIterate( sout,  sinput);

    chomp::MatX initialTrajectory;
    //after the arguments have been collected, pass them to chomp
    createInitialTrajectory( initialTrajectory );


    //now that we have a trajectory, make a chomp object
    // if there is an old chomp object, delete it.
    if (chomper){ delete chomper; } 
    chomper = new chomp::Chomp( factory, initialTrajectory,
                                q0, q1, info.n_max, 
                                info.alpha, info.obstol,
                                info.max_global_iter,
                                info.max_local_iter,
                                info.t_total, info.timeout_seconds,
                                info.use_momentum);
    
    if ( info.use_hmc ){
        if (hmc){ delete hmc; }
        hmc = new chomp::HMC( info.hmc_lambda, info.do_not_reject );
        
        if ( info.seed != 0 ){
            hmc->setSeed( info.seed );
        }

        chomper->hmc = hmc;
    }
    
    printChompInfo();

    //give chomp a collision helper to 
    //deal with collision and gradients.
    if (sphere_collider){
        //if a collision helper already exists, delete it, and make a new
        //  one.
        if ( collisionHelper ){ delete collisionHelper; }
        collisionHelper = new chomp::ChompCollGradHelper(
                                sphere_collider, info.gamma);
        chomper->ghelper = collisionHelper;
    }

    
    if ( info.doObserve ){
        observer = new chomp::DebugChompObserver();
        chomper->observer = observer;
    }

    //setup the mins
    chomper->min_global_iter = info.min_global_iter;
    chomper->min_local_iter = info.min_local_iter;


    //get the lock for the environment
    OpenRAVE::EnvironmentMutex::scoped_lock lock(environment->GetMutex() );

    if (!robot.get() ){
        robot = environment->GetRobot( robot_name.c_str() );
    }
    
    timer.start( "CHOMP run" );
    //solve chomp
    chomper->solve( info.doGlobal, info.doLocal );
    
    double elapsedTime = timer.stop( "CHOMP run" );
    double wallTime = timer.getWallElapsed("CHOMP run");

    RAVELOG_INFO( "Chomp process time %fs\n", elapsedTime );
    RAVELOG_INFO( "Chomp wall time    %fs\n", wallTime );
    
    RAVELOG_INFO( "Collision process time:      %fs\n",
                   sphere_collider->timer.getTotal( "collision"));
    RAVELOG_INFO( "Collision FK process time:   %fs\n",
                   sphere_collider->timer.getTotal( "FK"));
    RAVELOG_INFO( "Collision Jacobian process time: %fs\n",
                   sphere_collider->timer.getTotal( "jacobian"));
    RAVELOG_INFO( "xform Collision process time:  %fs\n",
                   sphere_collider->timer.getTotal( "xform"));

    RAVELOG_INFO( "Collision wall time:      %fs\n",
                   sphere_collider->timer.getWallTotal( "collision"));
    RAVELOG_INFO( "Collision FK wall time:   %fs\n",
                   sphere_collider->timer.getWallTotal( "FK"));
    RAVELOG_INFO( "Collision Jacobian wall time: %fs\n",
                   sphere_collider->timer.getWallTotal( "jacobian"));
    RAVELOG_INFO( "SDF Collision wall time:  %fs\n",
                   sphere_collider->timer.getWallTotal( "sdf collision"));
    RAVELOG_INFO( "Self Collision wall time: %fs\n",
                   sphere_collider->timer.getWallTotal( "self collision"));
    RAVELOG_INFO( "Xform Collision wall time: %fs\n",
                   sphere_collider->timer.getWallTotal( "xform"));
   
    RAVELOG_INFO( "Done Iterating" ); 
    return true;
}


void mod::checkTrajectoryForCollision(){

    RAVELOG_INFO("checking trajectory for collision ...\n");
    bool collides = false;

    timer.start( "collision check" );

    OpenRAVE::CollisionReportPtr report(new OpenRAVE::CollisionReport());

    double total_dist = 0.0;

    //get the length of the trajectory
    for ( int i = 0; i < chomper->xi.rows()-1; i ++ ) 
    {
        
        const chomp::MatX & point1 = chomper->xi.row( i );
        const chomp::MatX & point2 = chomper->xi.row( i+1 );
        
        const chomp::MatX diff = point1 - point2;

        const double val = (diff.array() * diff.array()).sum();
            
        total_dist += sqrt( val );
    }
    

    const double step_dist = 0.04;
    const double step_time = trajectory_ptr->GetDuration()
                             *step_dist/total_dist;

    for ( double time = 0.0;
          time < trajectory_ptr->GetDuration();
          time+=step_time)
    {
        std::vector< OpenRAVE::dReal > point;
        trajectory_ptr->Sample(point, time);
        robot->SetActiveDOFValues(point);
        
        if (environment->CheckCollision( robot, report) ||
            robot->CheckSelfCollision(report))
        {
            collides = true;

            if (!info.no_collision_details){
                RAVELOG_ERROR("Collision: %s\n",
                              report->__str__().c_str());
            }
            if (!info.no_collision_exception){
                throw OpenRAVE::openrave_exception(
                    "Resulting trajectory is in collision!");
            }
        }
    }

    if (collides){ RAVELOG_ERROR("   trajectory collides!\n"); }
}


bool mod::gettraj(std::ostream& sout, std::istream& sinput)
{
    
    RAVELOG_INFO( "Getting Trajectory\n" );
    OpenRAVE::EnvironmentMutex::scoped_lock lockenv;

    parseGetTraj( sout,  sinput);
    
    if ( !chomper ){
        RAVELOG_ERROR( "There is no trajectory to get. There must be a"
                       " call to iterate\n" );
    }
                       
    //get the lock for the environment
    lockenv = OpenRAVE::EnvironmentMutex::scoped_lock(
              environment->GetMutex() );
   
    
    //construct the trajectory pointer object, and 
    trajectory_ptr = OpenRAVE::RaveCreateTrajectory(environment);

    if (!robot.get() ){
        robot = environment->GetRobot( robot_name.c_str() );
    }
    trajectory_ptr->Init(
        robot->GetActiveConfigurationSpecification());

    

    //get the start state as an openrave vector
    std::vector< OpenRAVE::dReal > startState;
    getStateAsVector( chomper->q0, startState );

    //insert the start state into the trajectory
    trajectory_ptr->Insert( 0, startState );

    //get the rest of the trajectory
    for ( int i = 0; i < chomper->xi.rows(); i ++ ){
        
        std::vector< OpenRAVE::dReal > state;
        getIthStateAsVector( i, state );
        trajectory_ptr->Insert( i + 1, state );
        
    }
    
    //get the start state as an openrave vector
    std::vector< OpenRAVE::dReal > endState;
    getStateAsVector( chomper->q1, endState );

    //insert the end state into the trajectory
    trajectory_ptr->Insert( chomper->xi.rows(), endState );


    RAVELOG_INFO( "Retiming Trajectory\n" );
    //this times the trajectory so that it can be
    //  sent to a planner
    OpenRAVE::planningutils::RetimeActiveDOFTrajectory(
                             trajectory_ptr, robot, false, 0.2, 0.2,
                             "LinearTrajectoryRetimer","");
    
    if( !info.no_collision_check ){ checkTrajectoryForCollision(); }

    //TODO : check for collisions
    trajectory_ptr->serialize( sout );

    return true;
}

bool mod::execute(std::ostream& sout, std::istream& sinput){

    std::cout << "Executing" << std::endl;
    //get the lock for the environment
    OpenRAVE::EnvironmentMutex::scoped_lock lockenv(
              environment->GetMutex() );


    if ( trajectory_ptr.get() ){
        robot->GetController()->SetPath(trajectory_ptr);
        //robot->WaitForController(0);
    }
    else {
        RAVELOG_ERROR("There is no trajectory to run.\n");
    }
    
    std::cout << "Done Executing" << std::endl;
    return true;
        
}

void mod::delete_items(){
    if (chomper){
        delete chomper;
        chomper = NULL;
    }
    if (sphere_collider){ 
        delete sphere_collider;
        sphere_collider = NULL;
    }
    if (factory){
        delete factory;
        factory = NULL;
    }
    if (collisionHelper){ 
        delete collisionHelper;
        collisionHelper = NULL;
    }
    if (observer){
        delete observer;
        observer = NULL;
    }
    if (hmc){
        delete hmc;
        hmc = NULL;
    }
}

bool mod::destroy(std::ostream& sout, std::istream& sinput){
    
    RAVELOG_INFO( "Deleting orchomp module\n" );
    delete_items();

    return true;
}



//takes the two endpoints and fills the trajectory matrix by
//  linearly interpolating between the two.
inline void mod::createInitialTrajectory( chomp::MatX & trajectory )
{

    //make sure that the number of points is not zero
    assert( info.n != 0 );
    //make sure that the start and endpoint are the same size
    assert( q0.size() == q1.size() );

    //resize the trajectory to hold the current endpoint
    trajectory.resize(info.n, q0.size() );

    //fill the trajectory matrix 
    for (size_t i=0; i<info.n; ++i) {
        trajectory.row(i) = (i+1)*(q1-q0)/(info.n+1) + q0;

        chomp::MatX test = trajectory.row(i);

        //make sure that all of the points are within the joint limits
        //  this is unnecessary, but nice for now.
        //  TODO remove this.
        assert( isWithinLimits( test ));
    }
}


//get the spheres for collision detection
void mod::getSpheres()
{
    
    active_spheres.clear();
    inactive_spheres.clear();

    //a vector holding all of the pertinent bodies in the scene
    std::vector<OpenRAVE::KinBodyPtr> bodies;

    /* consider the robot kinbody, as well as all grabbed bodies */
    robot->GetGrabbed(bodies);
    bodies.push_back(environment->GetRobot( robot->GetName() ));
    
    //iterate over all of the bodies.
    for (size_t i=0; i < bodies.size(); i++)
    {
        OpenRAVE::KinBodyPtr body = bodies[i];

        //get the spheres of the body by creating an xml reader to
        //  extract the spheres from the xml files of the objects
        boost::shared_ptr<orchomp::kdata> data_reader = 
            boost::dynamic_pointer_cast<orchomp::kdata>
                (body->GetReadableInterface("orchomp"));
         
        //bail if there is no orcdchomp data.
        if (data_reader.get() == NULL ) {
            std::string error = "Kinbody called "
                              + body->GetName() 
                              + " does not have a <orchomp> tag defined.";
            
            throw OpenRAVE::openrave_exception(error);
        }
        
        
        for (size_t j = 0; j < data_reader->spheres.size(); j++ )
        {
            
            Sphere & sphere = data_reader->spheres[j];
            
            sphere.body = body.get();
            /* what robot link is this sphere attached to? */
            if (body.get() == robot.get()){
                
                //TODO THIS IS AN OUTRAGEOUS HACK
                //  PLEASE make it not necessary
                if ( sphere.linkname == "/right/wam0" ){
                    sphere.linkname = "/right/wam2";
                }else if ( sphere.linkname == "/left/wam0" ){
                    sphere.linkname = "/left/wam2";
                }

                sphere.link = robot->GetLink(sphere.linkname).get();
            }
            //the sphere is attached to a grabbed kinbody
            else{
                sphere.link = robot->IsGrabbing(body).get();
            }

            //if the link does not exist, throw an exception
            if(!sphere.link){ 
               throw OPENRAVE_EXCEPTION_FORMAT(
                     "link %s in <orcdchomp> does not exist.",
                     sphere.linkname, OpenRAVE::ORE_Failed);
            }
            
            sphere.linkindex = sphere.link->GetIndex();
            
            //if the body is not the robot, then get the transform
            //TODO find out if this is necessary or useful??
            if ( body.get() != robot.get() )
            {
                OpenRAVE::Transform T_w_klink = 
                    body->GetLink(sphere.linkname)->GetTransform();
                OpenRAVE::Transform T_w_rlink = sphere.link->GetTransform();
                OpenRAVE::Vector v = T_w_rlink.inverse()
                                     * T_w_klink 
                                     * sphere.position;
                sphere.position[0] = v.x;
                sphere.position[1] = v.y;
                sphere.position[2] = v.z;
            }
             
            /* is this link affected by the robot's active dofs? */
            for (size_t k = 0; k < n_dof; k++){
                if ( robot->DoesAffect( active_indices[k],
                                        sphere.linkindex    )){
                    active_spheres.push_back( sphere );
                    break;
                }
                //if the sphere is not active, add it to the inactive
                //  vector
                else if ( k == n_dof - 1 ) {
                    inactive_spheres.push_back( sphere);
                }
            }
        }
    }    
}


} /* namespace orchomp */


