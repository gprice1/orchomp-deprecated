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


bool mod::isWithinLimits( const chomp::MatX & mat ) const{
    assert( upperJointLimits.size() > 0 );
    assert( lowerJointLimits.size() > 0 );
    assert( mat.cols() > 0 );

    for ( int i = 0; i < mat.cols(); i ++ ){
        if ( mat(i) > upperJointLimits[i] ||
             mat(i) < lowerJointLimits[i] ){
            return false;
        }
    }
    return true;
}

void mod::coutTrajectory() const
{
    for ( int i = 0; i < trajectory.rows(); i ++ ){


        for ( int j = 0; j < trajectory.cols() ; j ++ ){

            std::cout << trajectory( i, j ) << "\t";
        }
        std::cout << "\n";
    }
}
void mod::isTrajectoryWithinLimits() const {
    for( int i = 0; i < trajectory.rows(); i ++ ){
        chomp::MatX test = trajectory.row(i);
        assert( isWithinLimits( test ) );
        //debugStream << "Point " << i << " is within limits" << std::endl;
    }
}

void ORConstraintFactory::evaluate(
                const std::vector<chomp::Constraint*>& constraints, 
                const chomp::MatX& xi, chomp::MatX& h_tot,
                chomp::MatX& H_tot, int step)
{

    debugStream << "Evaluating Constraints" <<std::endl; 

    size_t DoF = xi.cols();

    assert(size_t(xi.rows()) == constraints.size());

    //the number of rows in the complete matrices.
    size_t numCons = 0;
    
    //the number of total constraints,
    // and the number of timesteps we are actually looking at.
    const size_t size = constraints.size();
    size_t count = 0;
    
    std::vector< chomp::MatX > H_vec, h_vec;

    //annoyingly, with the use of the step, this is the
    //  correct size of the vectors.
    H_vec.resize( (size - 1)/step + 1  );
    h_vec.resize( (size - 1)/step + 1  );
    
    //keeps track of the timestep, while i keeps track of the
    //  vector index.
    std::vector <size_t> constrained_timesteps;
    
    //get all of the jacobians and cost vectors
    for (size_t t=0, i=0; t< size; t+=step, i++) {
        chomp::Constraint* c = constraints[t];
        c->evaluateConstraints( xi.row(t), h_vec[i], H_vec[i] );
        numCons += h_vec[i].rows();
        
        if ( c->numOutputs() != 0 ){
            assert( h_vec[i].rows() == h_vec[i].rows() );
            assert( h_vec[i].cols() == 1);

            constrained_timesteps.push_back( i );
        }

        count ++;
    }
    

    //bail out if there are no constraints.
    if ( constrained_timesteps.size() == 0 ){
        h_tot.resize( 0,0 );
        H_tot.resize( 0,0 );

        debugStream << "Done Evaluating Constraints" <<std::endl; 
        return;
    }


    assert( count == (size - 1)/step + 1  );
    
    // make h
    h_tot.resize( numCons, 1 );
    H_tot.resize( numCons, DoF*constrained_timesteps.size() );
    H_tot.setZero(); 
   
    
    
    //since we don't which of the steps is 0, 
    int row_start = 0;
    for (size_t i=0; i < constrained_timesteps.size(); i ++) {
        
        const size_t index = constrained_timesteps[i];
        const int height = h_vec[index].rows();
        
        h_tot.block(row_start, 0, height, 1 ) = h_vec[index];
        H_tot.block(row_start, i*DoF, height, DoF) = H_vec[index];

        row_start += height;

    }
    debugStream << "Done Evaluating Constraints" <<std::endl; 
    
}

void ORConstraint::evaluateConstraints(const chomp::MatX& qt, 
                                             chomp::MatX& h, 
                                             chomp::MatX& H){
    //jacobain columns must be equal to n_dof
    // the rows must be equal to the number of constraints
    int k = module->n_dof;

    if ( h.cols() != 1 || h.rows() != k ){
        h.resize( k, 1 );
    }
    if ( H.cols() != k || H.rows() != k ){
        H.resize( k, k );
    }
    H.setZero();
    
    int dims = 0;
    for ( int i = 0; i < qt.cols(); i ++ ){
        if ( qt(i) > module->upperJointLimits[i] ){
            h(dims) = qt(i) - module->upperJointLimits[i];
            H( i, dims ) = 1;
            dims++;
        }else if ( qt(i) < module->lowerJointLimits[i] ){
            h(dims) = qt(i) - module->lowerJointLimits[i];
            H( i, dims ) = 1;
            dims++;
        }
    }
    
    n_outputs = dims;
    if (dims != k){ 
        H.conservativeResize( k , dims );
        h.conservativeResize( dims, 1 );
    }
}

//constructor that registers all of the commands to the openRave
//   command line interface.
mod::mod(OpenRAVE::EnvironmentBasePtr penv) :
    OpenRAVE::ModuleBase(penv), environment( penv ), 
    factory( NULL ), sphere_collider( NULL ),
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
      RegisterCommand("execute",
            boost::bind(&mod::execute,this,_1,_2),
            "play a trajectory on a robot");

}

/* ======================================================================== *
 * module commands
 */

// NOTES : viewspheres looks to be fine without editing.
//    - There may be hidden dependency issues
int mod::viewspheres(std::ostream& sout, std::istream& sinput)
{

    
    OpenRAVE::EnvironmentMutex::scoped_lock lockenv(environment->GetMutex());
    parseViewSpheres( sout,  sinput);

    if ( !sphere_collider ) { getSpheres() ; }
    
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
        sprintf( text_buf, "orcdchomp_sphere_%d", int(i));
        sbody->SetName(text_buf);
        

        //set the dimensions and transform of the sphere.
        std::vector< OpenRAVE::Vector > svec;

        //get the position of the sphere in the world 
        OpenRAVE::Transform t = 
                sphere.body->GetLink(sphere.linkname)->GetTransform();
        OpenRAVE::Vector v = t * OpenRAVE::Vector(sphere.pose);
        
        //set the radius of the sphere
        v.w = sphere.radius; 

        //give the kinbody the sphere parameters.
        svec.push_back(v);
        sbody->InitFromSpheres(svec, true);

        environment->Add( sbody );
    }
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
    // Add A distance field to the vector of distance fields.
    sdfs.resize( sdfs.size() + 1 );

    //lock the environment
    OpenRAVE::EnvironmentMutex::scoped_lock(environment->GetMutex());
    
    //Parse the arguments
    parseComputeDistanceField( sout,  sinput);
    
    //extract the distance field to be created.
    DistanceField & sdf = sdfs.back();
    
    sdf.createField( environment );

    return 0;
}


int mod::addfield_fromobsarray(std::ostream& sout, std::istream& sinput)
{
    parseAddFieldFromObsArray( sout,  sinput);
   
    return 0;

}




/* runchomp robot Herb2
 * run chomp from the current config to the passed goal config
 * uses the active dofs of the passed robot
 * initialized with a straight-line trajectory
 * */
int mod::create(std::ostream& sout, std::istream& sinput)
{
    //get the lock for the environment
    OpenRAVE::EnvironmentMutex::scoped_lock lockenv(
              environment->GetMutex() );

    std::cout << "Creating" << std::endl;

    parseCreate( sout,  sinput);
    
    clampToLimits(q0);
    clampToLimits(q1);
    //after the arguments have been collected, pass them to chomp
    createInitialTrajectory();
    
    factory = new ORConstraintFactory( this );
    //now that we have a trajectory, make a chomp object
    chomper = new chomp::Chomp( factory, trajectory,
                                q0, q1, info.n_max, 
                                info.alpha, info.obstol,
                                info.max_global_iter,
                                info.max_local_iter,
                                info.t_total);

     
    //TODO Compute signed distance field
    

    //Setup collision geometry
    getSpheres();

    //create the sphere collider to actually 
    //  use the sphere information
    //TODO make sure that nbodies is correct
    sphere_collider = new SphereCollisionHelper(
                          n_dof, 3, active_spheres.size(), this);


    //TODO setup momentum stuff ?? maybe ?? 

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
    
    //get the arguments
    parseIterate( sout,  sinput);

    //get the lock for the environment
    OpenRAVE::EnvironmentMutex::scoped_lock lock(environment->GetMutex() );

    if (!robot.get() ){
        robot = environment->GetRobot( robot_name.c_str() );
    }

    //solve chomp
    chomper->solve( info.doGlobal, info.doLocal );

    std::cout << "Done Iterating" << std::endl;
    return 0;
}

int mod::gettraj(std::ostream& sout, std::istream& sinput)
{
    
    //get the trajectory from the chomp.
    trajectory = chomper->xi;

    std::cout << "Getting Trajectory" << std::endl;
    OpenRAVE::EnvironmentMutex::scoped_lock lockenv;

    parseGetTraj( sout,  sinput);
    
    //get the lock for the environment
    lockenv = OpenRAVE::EnvironmentMutex::scoped_lock(
              environment->GetMutex() );
   

    debugStream << "Checking Trajectory" << std::endl;
    //check and or print the trajectory
    //isTrajectoryWithinLimits();
    //coutTrajectory();

    //setup the openrave trajectory pointer to receive the
    //  found trajectory.
    trajectory_ptr = OpenRAVE::RaveCreateTrajectory(environment);

    if (!robot.get() ){
        robot = environment->GetRobot( robot_name.c_str() );
    }
    trajectory_ptr->Init(
        robot->GetActiveConfigurationSpecification());

    debugStream << "Done locking environment,"
                << " begun extracting trajectory" << std::endl;
    
    
    //For every state (each row is a state), extract the data,
    //  and turn it into a vector, then give it to


    //get the start state as an openrave vector
    std::vector< OpenRAVE::dReal > startState;
    getStateAsVector( q0, startState );

    //insert the start state into the trajectory
    trajectory_ptr->Insert( 0, startState );

    //get the rest of the trajectory
    for ( int i = 0; i < trajectory.rows(); i ++ ){
        
        std::vector< OpenRAVE::dReal > state;
        getIthStateAsVector( i, state );
        trajectory_ptr->Insert( i + 1, state );
        
    }
    
    //get the start state as an openrave vector
    std::vector< OpenRAVE::dReal > endState;
    getStateAsVector( q1, endState );

    //insert the start state into the trajectory
    trajectory_ptr->Insert( trajectory.rows(), endState );


    //this is about changing the timing or something
    /*
       OpenRAVE::planningutils::RetimeActiveDOFTrajectory
                (traj_ptr, boostrobot, false, 0.2, 0.2,
                 "LinearTrajectoryRetimer","");
    */
    
    //TODO : check for collisions
    
    debugStream << "Retiming Trajectory" << std::endl;
    //this times the trajectory so that it can be
    //  sent to a planner
    OpenRAVE::planningutils::RetimeActiveDOFTrajectory(
                             trajectory_ptr, robot);

    debugStream << "Serializing trajectory output" << std::endl; 
    //serialize the trajectory and send it over the 
    //  output stream.
    trajectory_ptr->serialize( sout );

    std::cout << "Done Getting Trajectory" << std::endl;


    return 0;
}
int mod::execute(std::ostream& sout, std::istream& sinput){

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

        chomp::MatX test = trajectory.row(i);
        assert( isWithinLimits( test ));
    }
}

inline void mod::clampToLimits( chomp::MatX & state ){
    for ( int i = 0; i < state.cols() ; i ++ ){
        if ( state(i) > upperJointLimits[i] ){
            state(i) = upperJointLimits[i];
        }
        else if ( state(i) < lowerJointLimits[i] ){
            state(i) = lowerJointLimits[i];
        }
    }
}



void mod::getSpheres()
{
    
    
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
        boost::shared_ptr<kdata> data_reader = 
            boost::dynamic_pointer_cast<kdata>
                (body->GetReadableInterface("orcdchomp"));
         
        //bail if there is no orcdchomp data.
        if (data_reader.get() == NULL ) {
            debugStream << "Failed to get: " << body->GetName() << std::endl;
            continue;

            throw OpenRAVE::openrave_exception(
                "kinbody does not have a <orcdchomp> tag defined!");
        }
        
        
        for (size_t j = 0; j < data_reader->spheres.size(); j++ )
        {
            
            Sphere & sphere = data_reader->spheres[j];
            
            sphere.body = body.get();
            /* what robot link is this sphere attached to? */
            if (body.get() == robot.get()){
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
            if ( body.get() != robot.get() )
            {
                OpenRAVE::Transform T_w_klink = 
                    body->GetLink(sphere.linkname)->GetTransform();
                OpenRAVE::Transform T_w_rlink = sphere.link->GetTransform();
                OpenRAVE::Vector v = T_w_rlink.inverse()
                                     * T_w_klink 
                                     * OpenRAVE::Vector(sphere.pose);
                sphere.pose[0] = v.x;
                sphere.pose[1] = v.y;
                sphere.pose[2] = v.z;
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


