#include "orchomp_collision.h"
#include "orchomp_mod.h"
#include "orchomp_collision_pruner.h"


namespace orchomp {

typedef Eigen::Map< Eigen::Matrix< double, 
                                   Eigen::Dynamic, 
                                   Eigen::Dynamic, 
                                   Eigen::RowMajor > > MatrixMap;
typedef Eigen::Map< const Eigen::Matrix< double, 
                                         Eigen::Dynamic, 
                                         Eigen::Dynamic, 
                                         Eigen::RowMajor > > ConstMatrixMap;

template <class Derived>
void assertJacobianIsEquivalent(const Eigen::MatrixBase<Derived> & mat,
                                const std::vector<OpenRAVE::dReal> & vec ){
    for ( int i = 0; i < mat.rows(); i ++ ){
        for ( int j = 0; j < mat.cols(); j ++ ){
            assert( mat(i,j) == vec[i*mat.cols() + j]);
        }
    }
}


//the constuctor needs a pointer to the robot in addition to the spaces.
SphereCollisionHelper::SphereCollisionHelper( 
                       size_t ncspace,
                       mod * module,
                       double gamma, 
                       double epsilon, 
                       double obs_factor,
                       double epsilon_self,
                       double obs_factor_self) :
        ncspace(ncspace), nwkspace(3),
        pruner( NULL ),
        module(module),
        inactive_spheres_have_been_set( false ),
        gamma( gamma ),
        epsilon( epsilon ), obs_factor(obs_factor),
        epsilon_self( epsilon_self ), obs_factor_self( obs_factor_self)
{
    //fill the ignorables set with the adjacent links
    ignorables.insert( module->robot->GetAdjacentLinks().begin(),
                       module->robot->GetAdjacentLinks().end() );
    getSpheres();
    initPruner();

}

SphereCollisionHelper::~SphereCollisionHelper(){
    if ( pruner ){ delete pruner; }
}

OpenRAVE::dReal computeCostFromDist( OpenRAVE::dReal dist,
                                     double epsilon,
                                     Eigen::Vector3d & gradient ){

    //compute the cost and gradient function based on the returned
    //  of get dist, and the epsilon value of the planner
    //if the cost is negative, invert the gradient, and calculate the
    //  cost modified by epsilon




    if (dist < 0) {
        gradient *= -1;
        return -dist + 0.5*epsilon;

    } 
    
    if (dist <= epsilon) {

        const double f = dist - epsilon;
        gradient *= f*0.5/epsilon;
        return f*f * 0.5/epsilon;
    }
    //if the gradient is far away enough from the object,
    //  then set the costs and gradient to zero
    
    return 0;

}


double SphereCollisionHelper::addToGradient(const chomp::MatX& xi,
                                            const chomp::MatX& pinit,
                                            const chomp::MatX& pgoal,
                                            double dt,
                                            chomp::MatX& g) {
    timer.start( "collision" );

    q1 = chomp::getTickBorderRepeat(-1, xi, pinit, pgoal, dt).transpose();
    q2 = chomp::getTickBorderRepeat(0, xi, pinit, pgoal, dt).transpose();
    
    chomp::MatX g_sdf( 1, ncspace);
    chomp::MatX g_self( 1, ncspace);

    inv_dt = 1/dt;
    const double inv_dt_squared = inv_dt * inv_dt;
    double total_cost = 0.0;
    
    for (int current_time=0; current_time<xi.rows(); ++current_time) {
        g_sdf.setZero();
        g_self.setZero();
        
        double sdf_cost(0), self_cost(0);

        q0 = q1;
        q1 = q2;
        q2 = chomp::getTickBorderRepeat(current_time+1,
                                        xi, pinit, pgoal, dt).transpose();

        cspace_vel = 0.5 * (q2 - q0) * inv_dt;        
        cspace_accel = (q0 - 2.0*q1 + q2) * inv_dt_squared;
       
        
        timer.start( "FK" );
        //Set the positions of all of the spheres,
        //  for the current configuration.
        setSpherePositions( q1, !inactive_spheres_have_been_set );
        timer.stop( "FK" );

        //clear the jacobian matrices
        jacobians.clear();
    
        for (size_t body_index=0; body_index < nbodies; body_index++) {
            
            timer.start( "self collision" );
            self_cost += getActiveCost( body_index, g_self );
                        
            //get the contribution from the inactive_spheres
            Eigen::Vector3d inactive_gradient( 0,0,0);
            double inactive_cost = getInactiveCost( body_index,
                                                    inactive_gradient);
            inactive_gradient /= inactive_cost;

            
            timer.stop( "self collision" );

            timer.start( "sdf collision" );
            Eigen::Vector3d sdf_gradient(0,0,0);
            double current_sdf_cost = getSDFCollisions( body_index,
                                                        sdf_gradient );
            
            timer.stop( "sdf collision" );
            if ( current_sdf_cost > 0 || inactive_cost > 0 ){
                //turn the jacobian into a matrix
                const std::vector< double > & jacobian = 
                                              getJacobian( body_index ); 

                ConstMatrixMap dx_dq ( jacobian.data(), nwkspace, ncspace);
                

                if ( inactive_cost > 0 ) {

                    timer.start( "self collision" );
                    self_cost += addInWorkspaceGradient( inactive_cost,
                                                         inactive_gradient,
                                                         dx_dq, g_self );
                    timer.stop( "self collision" );
                }
                if ( current_sdf_cost > 0 ){
                    timer.start( "sdf collision" );
                    sdf_cost += addInWorkspaceGradient( current_sdf_cost,
                                                        sdf_gradient,
                                                        dx_dq, g_sdf);
                    timer.stop( "sdf collision" );
                }
            }
        }

        g.row( current_time ) += g_sdf * obs_factor + g_self * obs_factor_self;
        total_cost += sdf_cost * obs_factor + self_cost * obs_factor_self;

    }


    timer.stop( "collision" );
    return total_cost;

}

double SphereCollisionHelper::getActiveCost( size_t body_index, 
                                             chomp::MatX & g_self ){
    
    double total_cost = 0;

    //get the self collisions with the active spheres.
    //We only need to go up, because we only need to check all collisions once.
    //  Furthermore, collisions from sphere_i to sphere_j and from sphere_j
    //  to sphere_i result in the same c-space gradient update, so,
    //  we only need to compute this update for one of the spheres. 
    for( size_t coll_index=body_index+1; coll_index < nbodies; coll_index++ ){
        
        //get the cost and gradient of the collision.
        Eigen::Vector3d gradient;
        double cost = sphereOnSphereCollision( body_index,
                                               coll_index,
                                               gradient);
        if ( cost <= 0.0 ){ continue; }
        
        //get the jacobians of the spheres.
        std::vector< double > jacobian;
        jacobian.resize( nwkspace * ncspace );
        const std::vector< double > & jac1 = getJacobian( body_index ); 
        const std::vector< double > & jac2 = getJacobian( coll_index ); 
        
        for ( size_t i = 0; i < jac1.size() ; i ++ ){
            jacobian[i] = jac1[i] - jac2[i];
        }
        
        //map the jacobian vector into an eigen matrix.
        MatrixMap dx_dq ( jacobian.data(), nwkspace, ncspace );

        //add in the cost of the collision.
        total_cost += addInWorkspaceGradient( cost, gradient, dx_dq, g_self);
    }

    return total_cost;
}


double SphereCollisionHelper::getInactiveCost( size_t body_index,
                                               Eigen::Vector3d & gradient_total ){

    double total_cost = 0.0;

    for ( size_t coll_index=nbodies; coll_index < spheres.size(); coll_index++ )
    {

        Eigen::Vector3d gradient;
        double cost =sphereOnSphereCollision( body_index,
                                               coll_index,
                                               gradient);
        if ( cost <= 0.0 ){ continue; }
        
        gradient_total += gradient * cost;
        total_cost += cost;
    }

    return total_cost;
}


template <class Derived>
double SphereCollisionHelper::addInWorkspaceGradient(
                                   double cost, const Eigen::Vector3d & grad,
                                   const Eigen::MatrixBase<Derived> & dx_dq,
                                   chomp::MatX & g_current)

{
    
    timer.start( "projection" );
    wkspace_vel = dx_dq * cspace_vel;

    //this prevents nans from propagating. Several lines below, 
    //    wkspace_vel /= wv_norm if wv_norm is zero, nans propogate.
    if (wkspace_vel.isZero()){ return 0; }

    wkspace_accel = dx_dq * cspace_accel;
    
    float wv_norm = wkspace_vel.norm();
    wkspace_vel /= wv_norm;

    // add to total
    double scl = wv_norm * gamma / inv_dt;

    P = chomp::MatX::Identity(nwkspace, nwkspace)
        - (wkspace_vel * wkspace_vel.transpose());

    K = (P * wkspace_accel) / (wv_norm * wv_norm);
   
    //          scalar * M-by-W        * (WxW * Wx1   - scalar * Wx1)
    //  The .noalias() call is simply an optimization for eigen.
    g_current.noalias() += (scl * (dx_dq.transpose() *
                 (P * grad - cost * K) ).transpose());

    timer.stop( "projection" );
    return cost * scl;
}




OpenRAVE::dReal SphereCollisionHelper::getSDFCollisions(
                                      size_t body_index,
                                      Eigen::Vector3d & gradient,
                                      bool viewDists)
{
    
    OpenRAVE::dReal dist = HUGE_VAL;
    size_t sdf_index = -1;

    
    //check all of the sdfs, and get the one with the least dist 
    for ( size_t i = 0; i < module->sdfs.size(); i ++ ){
        

        OpenRAVE::dReal current_dist;
        vec3 current_gradient; 

        //get the distance and gradient.
        current_dist = module->sdfs[i].getDist( sphere_positions[body_index], 
                                                current_gradient );
        
        if (current_dist < dist) {
            dist = current_dist;

            gradient[0] = current_gradient[0];
            gradient[1] = current_gradient[1];
            gradient[2] = current_gradient[2];

            sdf_index = i;
        }
    }

    if (dist != HUGE_VAL){

        //adjust the value of the distance by the radius of the
        //  sphere to account for the size of the sphere
        dist -= spheres[body_index].radius;
        
        if ( viewDists && dist > epsilon ){
            createCube( dist, module->sdfs[sdf_index].cube_extent,
                        sphere_positions[body_index], sdf_index);
        }

        //debugStream << "Finished SDF Collisions" << std::endl;
        return computeCostFromDist( dist, epsilon, gradient );
        

    }

    return 0.0;

}

bool SphereCollisionHelper::checkCollision( size_t body1, size_t body2 )
{

    if ( body1 > body2 ){ std::swap( body1, body2 ); }

    //the first element is not an active sphere, so continue
    if ( body1 >= nbodies ){ return false; }
      
    //If Body2 is an SDF
    if ( body2 >= spheres.size() )
    {
        if ( getSDFCollision( body1, body2 - spheres.size() ) ){
            return true;
        }
    }

    // If Body2 is a sphere, and
    //  the pair is not in the ignore set, then check for
    //  a collision.
    else if (sphereOnSphereCollision(body1, body2) ){ return true; }
    return false;
}

bool SphereCollisionHelper::isCollided()
{
    timer.start( "pruner");
    bool collision = pruner->checkPotentialCollisions( this );
    timer.stop( "pruner");

    return collision;

    /*
    assert( pruner );
    std::vector< std::pair<size_t,size_t> > collision_list;

    timer.start( "get coll" );
    pruner->getPotentialCollisions( collision_list );
    timer.stop( "get coll" );

    timer.start( "test coll");
    for ( size_t i=0; i < collision_list.size(); i ++ )
    {
        
        if ( checkCollision( collision_list[i].first,
                             collision_list[i].second ) )
        {

            timer.stop( "test coll");
            return true;
        }
    }
    timer.stop( "test coll");
    return false;
    */
}


bool SphereCollisionHelper::getSDFCollisions(size_t body_index)
{
    
    //check all of the sdfs, and get the one with the least dist 
    const double radius = spheres[body_index].radius;

    for ( size_t i = 0; i < module->sdfs.size(); i ++ ){
        OpenRAVE::dReal dist =
                module->sdfs[i].getDist( sphere_positions[body_index]);
        
        // get the distance from the sdf to the surface of the sphere.
        //  if the distance is less than 0, the surface of the sphere is
        //  inside of an object
        if (dist - radius < 0 ){ return true; }
    }

    return false;
}

bool SphereCollisionHelper::getSDFCollision(size_t body_index, 
                                            size_t sdf_index)
{
    
    const OpenRAVE::dReal dist = module->sdfs[sdf_index].getDist(
                                         sphere_positions[body_index]);
        
    return (dist - spheres[body_index].radius < 0 );
}

OpenRAVE::dReal SphereCollisionHelper::sphereOnSphereCollision(
                                size_t index1, size_t index2,
                                Eigen::Vector3d & gradient,
                                bool ignore){

    const Sphere & sphere1 = spheres[index1];
    const Sphere & sphere2 = spheres[index2];

    if ( ignore && ignoreSphereCollision( sphere1, sphere2 ) ){ 
        return 0.0;
    }

    //calculate the distance between the two centers of the spheres,
    //  also get the vector from the collision sphere to the
    //  current sphere.
    const OpenRAVE::Vector diff = sphere_positions[index1] 
                                - sphere_positions[index2];

    const OpenRAVE::dReal dist_sqrd = diff[0]*diff[0] + 
                                      diff[1]*diff[1] + 
                                      diff[2]*diff[2] ; 

    // get the distance between the spheres, with padding corresponding to
    //  epsilon.
    
    const double radius1 = sphere1.radius;
    const double radius2 = sphere2.radius;
    
    const double total_radius = radius1 + radius2;
    const double total_range  = total_radius + epsilon_self;

    //get the actual distance between the spheres.
    
    //if the spheres are too far away, do not compute any cost.
    if ( dist_sqrd > total_range*total_range ){ return 0.0; }
    
    const OpenRAVE::dReal dist_between_centers = sqrt( dist_sqrd );
    const OpenRAVE::dReal dist = dist_between_centers - total_radius;

    //get a normalized gradient vector, by dividing diff by its length
    //  and add it into the previously computed gradient
    gradient[0] = diff[0] / dist_between_centers;
    gradient[1] = diff[1] / dist_between_centers; 
    gradient[2] = diff[2] / dist_between_centers;

    //compute the cost from the distance.
    return computeCostFromDist( dist, epsilon_self, gradient );

}


bool SphereCollisionHelper::sphereOnSphereCollision( size_t index1,
                                                     size_t index2,
                                                     bool ignore){
    
    const Sphere & sphere1 = spheres[index1];
    const Sphere & sphere2 = spheres[index2];

    if ( ignore && ignoreSphereCollision( sphere1, sphere2 ) ){ 
        return false;
    }

    //calculate the distance between the two centers of the spheres,
    //  also get the vector from the collision sphere to the
    //  current sphere.
    const OpenRAVE::Vector diff = sphere_positions[index1] 
                                - sphere_positions[index2];

    const OpenRAVE::dReal dist_sqrd = diff[0]*diff[0] + 
                                      diff[1]*diff[1] + 
                                      diff[2]*diff[2] ; 

    
    const double total_radius = sphere1.radius + sphere2.radius;

    //if the distance is less than 0, it is in collision;
    if (dist_sqrd <= total_radius * total_radius){  return true; }
    
    return false;

}

void SphereCollisionHelper::setSpherePositions( const chomp::MatX & q,
                                                bool setInactive){
    if ( sphere_positions.size() != spheres.size()){
        sphere_positions.resize( spheres.size() );
    }

    std::vector< OpenRAVE::dReal > vec;
    module->getStateAsVector( q, vec );
    
    setSpherePositions( vec, setInactive );
}
 
void SphereCollisionHelper::setSpherePositions(
                            const std::vector<OpenRAVE::dReal> & state,
                            bool setInactive)
{   

    module->robot->SetActiveDOFValues(state, false);
    
    OpenRAVE::Transform t;
    int current_link_index = -1;
    OpenRAVE::KinBody * current_body = NULL;
    
    //if setInactive is true, then set all of the spheres,
    //  otherwise, only set the active spheres.
    const size_t size = ( setInactive ? spheres.size() : nbodies );

    //get the positions of all of the spheres
    for ( size_t i=0; i < size; i ++ )
    {
        //extract the current sphere
        //  if i is less than n_active, get a sphere from active_spheres,
        //  else, get a sphere from inactive_spheres
        const Sphere & sphere = spheres[i];
        
        //only get a new xform, if the new xform is from a different link.
        //  if they are the same link, just use the only transform.
        if ( sphere.linkindex != current_link_index ||
             sphere.body != current_body ){
            //get the transformation of the body that the sphere is on.
            t = sphere.link->GetTransform();

            //store the current body and link.
            current_link_index = sphere.linkindex;
            current_body = sphere.body;
        }

        //get the transformation from the body to the sphere.
        sphere_positions[i] = t * sphere.position;
    }

    if (pruner){
        timer.start( "sort");
        pruner->sort( sphere_positions, size);
        timer.stop( "sort");
    }

    if ( setInactive ){ inactive_spheres_have_been_set = true;}

}

std::vector< OpenRAVE::dReal > const &
SphereCollisionHelper::getJacobian( size_t sphere_index )
{
    
    timer.start( "jacobian" );
    map::iterator it = jacobians.find( sphere_index );

    //if the object exists, return the thing.
    if ( it != jacobians.end() ){ return it->second; }
    
    const Sphere & sphere = spheres[ sphere_index ];

    //if the jacobian does not exist,create and return it.
    key_value_pair new_jacobian;
    new_jacobian.first = sphere_index;
    std::pair< map::iterator, bool> inserted_element = 
                   jacobians.insert( new_jacobian );

    //actually get the jacobian
    module->robot->CalculateActiveJacobian(
                   sphere.linkindex, 
                   sphere_positions[sphere_index],
                   inserted_element.first->second);


    timer.stop( "jacobian" );
    return inserted_element.first->second;
}


OpenRAVE::KinBodyPtr SphereCollisionHelper::createCube( 
                                            double dist,
                                            double size,
                                            const OpenRAVE::Vector & pos,
                                            size_t sdf_index)
{
    
    OpenRAVE::Vector color, extents( size, size, size );
    colorFromDist( dist, sdf_index, color );
    OpenRAVE::KinBodyPtr cube = module->createBox( pos, extents, color );

    return cube;

}

void SphereCollisionHelper::colorFromDist( double dist,
                                           size_t sdf_index,
                                           OpenRAVE::Vector & color ){

    const DistanceField & df = module->sdfs[ sdf_index ];

    const double min = df.grid.minDist();
    const double max = df.grid.maxDist();
    
    const double cutoff1 = (max - min) / 3;
    const double cutoff2 = cutoff1 * 2;
    
    dist -= min;
    
    //As distance increases, the color goes from red to yellow to green
    // to blue, to black. 
    if ( dist < cutoff1 ){
        double val = dist/cutoff1;
        color = OpenRAVE::Vector( 1-val, val , 0);

    } else if ( dist < cutoff2 ){
        double val = (dist - cutoff1)/cutoff1;
        color = OpenRAVE::Vector( 0, 1-val, val  );

    }else {
        double val = (dist - cutoff2)/cutoff1;
        color = OpenRAVE::Vector( val, val , 1 );
    }

}


void SphereCollisionHelper::visualizeSDFSlice( size_t sdf_index,
                                               size_t axis,
                                               size_t slice_index,
                                               double time)
{

    assert( axis < 3 && axis >=0 );
    assert( module->sdfs.size() > sdf_index && sdf_index >= 0 );
    
    const DistanceField & df = module->sdfs[ sdf_index ];
    
    bool was_visible = false;
    if ( df.kinbody->IsVisible() ){
        df.kinbody->SetVisible( false );
        was_visible = true;
    }

    
    size_t bounds[6] = { 0,0,0, df.grid.nx(), df.grid.ny(), df.grid.nz() };
    bounds[ axis ] = slice_index;
    bounds[axis + 3] = slice_index + 1;

    std::vector< OpenRAVE::KinBodyPtr > cubes;

    for( size_t i = bounds[0]; i < bounds[3]; i ++ ){
    for( size_t j = bounds[1]; j < bounds[4]; j ++ ){
    for( size_t k = bounds[2]; k < bounds[5]; k ++ ){

        double dist = df.grid( i, j, k );
        OpenRAVE::Transform center;
        df.getCenterFromIndex( i,j,k, center );

        cubes.resize( cubes.size() + 1 );
        
        cubes.back() = createCube( dist, df.cube_extent,
                                   center.trans, sdf_index );

    }
    }
    }

    //wait for given amount of time
    struct timespec ticks_tic;
    struct timespec ticks_toc;

    /* start timing voxel grid computation */
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ticks_tic);


    while ( true ){
      /* stop timing voxel grid computation */
      clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ticks_toc);
      CD_OS_TIMESPEC_SUB(&ticks_toc, &ticks_tic);
      if ( time < CD_OS_TIMESPEC_DOUBLE(&ticks_toc) ){ break ; }
    }

    for ( size_t i = 0; i < cubes.size() ; i ++ ){
        module->environment->Remove( cubes[i] );
    }

    if ( was_visible ) { df.kinbody->SetVisible(true); }

}

void SphereCollisionHelper::benchmark( int num_trials,
                                       bool check_all,
                                       double pause_time,
                                       bool print){
    

    int collisions_or(0), collisions_sphere(0),
        collisions_or_self(0), collisions_sphere_self(0), 
        collisions_or_env(0), collisions_sphere_sdf(0);

    int missed_collisions(0), fake_collisions( 0 );
    int missed_collisions_self( 0), missed_collisions_env(0);
    int fake_collisions_self(0), fake_collisions_sdf(0), 
        pruning_missed(0);
    
    for ( int i = 0; i < num_trials; i ++ ){
        if ( i % 1000 == 0 ){
            RAVELOG_INFO( "Completed %d of %d trials.\n", i, num_trials );
        }

        if ( print ) { std::cout << "\nGetting new Config\n";}

        chomp::MatX mat;
        std::vector< OpenRAVE::dReal > state_vector;
        module->getRandomState(mat);
        module->getStateAsVector( mat , state_vector );

        timer.start( "or fk" );
        module->robot->SetActiveDOFValues( state_vector, false );
        timer.stop( "or fk");

        bool is_collided_or( false ), is_collided_sphere(false),
             is_collided_or_self(false), is_collided_sphere_self(false),
             is_collided_or_env(false), is_collided_sphere_sdf(false),
             is_collided_pruner(false);

        OpenRAVE::CollisionReportPtr env_report(new OpenRAVE::CollisionReport());
        OpenRAVE::CollisionReportPtr self_report(new OpenRAVE::CollisionReport());
        
        timer.start( "openrave env" );
        is_collided_or_env = module->environment->CheckCollision(
                                                  module->robot,
                                                  env_report );
        timer.stop( "openrave env" );
        
        timer.start("openrave self");
        is_collided_or_self = module->robot->CheckSelfCollision(
                                                            self_report );
        timer.stop( "openrave self");


        timer.start( "sphere fk" );
        //this will set the inactive spheres, if they have not
        //  already been set.
        setSpherePositions( mat, !inactive_spheres_have_been_set );
        timer.stop( "sphere fk" );

        timer.start("sphere sdf");
        is_collided_sphere_sdf = isCollidedSDF(check_all);
        timer.stop( "sphere sdf");

        timer.start("sphere self");
        is_collided_sphere_self = isCollidedSelf(check_all);
        timer.stop( "sphere self");

        is_collided_pruner = isCollided();

        bool fake( false ), missed( false );
        
        if( is_collided_or_self || is_collided_or_env ){
            collisions_or ++;
            is_collided_or = true;
        }

        if( is_collided_or_self){ 
            collisions_or_self ++;
        }
        if( is_collided_or_env ){
            collisions_or_env ++;
        }
        
        if( is_collided_sphere_sdf || is_collided_sphere_self ){
            collisions_sphere ++;
            is_collided_sphere = true;
        }

        if( is_collided_sphere_self ){ collisions_sphere_self ++;}
        if( is_collided_sphere_sdf ){ collisions_sphere_sdf ++;}

        if( !is_collided_pruner && is_collided_sphere ){
            
            std::string type; 
            if ( is_collided_sphere_self){ type = "self"; }
            if ( is_collided_sphere_sdf){ type = "sdf"; }
            RAVELOG_ERROR( "ERROR in pruning, missed %s collision\n",
                            type.c_str() );
            pruning_missed ++;
        }
        
        if( is_collided_pruner && !is_collided_sphere ){
            RAVELOG_ERROR( "ERROR in pruning, fake collision\n" );
        }
        
        if( (is_collided_or_self || is_collided_or_env ) && 
           !(is_collided_sphere_sdf || is_collided_sphere_self) ){
            missed = true;
            missed_collisions ++;
            if ( is_collided_or_self ){
                missed_collisions_self ++ ;
            }
            if ( is_collided_or_env ){
                missed_collisions_env ++;
            }
        }
        if( !(is_collided_or_self || is_collided_or_env ) && 
             (is_collided_sphere_sdf || is_collided_sphere_self) ){
            fake = true;
            fake_collisions ++;
            if ( is_collided_sphere_sdf ){
                fake_collisions_sdf ++ ;
            }
            if ( is_collided_sphere_self ){
                fake_collisions_self ++;
            }
        }
        
        //print out stuff if the collision was missed:
        if ( missed ){
            
            if (missed) { std::clog << "#Missed collision: "; }
            if( is_collided_or_self){ 
                std::clog << "\tself: "
                          << self_report->plink1->GetName()
                          << ", " 
                          << self_report->plink2->GetName()
                          << "\n";
            }
            if( is_collided_or_env ){
                std::clog << "\tenvironmental: "
                          << env_report->plink1->GetName()
                          << ", " 
                          << env_report->plink2->GetName()
                          << "\n";
            }
            
            if ( pause_time > 0.0 ){
                timer.wait( pause_time );
            }

        }
    }

    double total_sphere, total_or, total_pruner,
           sphere_self, sphere_sdf, sphere_fk,
           or_self, or_env, or_fk, pruner_getting,
           pruner_testing, pruner_sorting;

    sphere_self = timer.getTotal( "sphere self" );
    sphere_sdf =  timer.getTotal( "sphere sdf");
    sphere_fk =   timer.getTotal( "sphere fk" );
    or_self =   timer.getTotal( "openrave self" );
    or_env =    timer.getTotal( "openrave env" );
    or_fk =     timer.getTotal( "or fk" );

    total_sphere = timer.getTotal( "sphere self" ) +
                   timer.getTotal( "sphere sdf") +
                   timer.getTotal( "sphere fk" ) -
                   timer.getTotal( "sort2" ) -
                   timer.getTotal( "sort" );
    total_or = timer.getTotal( "openrave self" ) +
               timer.getTotal( "openrave env" ) +
               timer.getTotal( "or fk" );

    pruner_sorting = timer.getTotal( "sort" );
    pruner_testing = timer.getTotal( "pruner" ); 
    total_pruner = pruner_testing + sphere_fk;

    if (print){
        RAVELOG_INFO("Total collisions:  %d\n" , collisions_or ); 
        RAVELOG_INFO("Sphere collisions: %d\n" , collisions_sphere );

        RAVELOG_INFO("Missed collisions:  %d\n" , missed_collisions ); 
        RAVELOG_INFO("Missed self collisions:  %d\n", missed_collisions_self );        RAVELOG_INFO("Missed env collisions:  %d\n" , missed_collisions_env ); 
        
        RAVELOG_INFO("Fake collisions: %d\n" , fake_collisions );
        RAVELOG_INFO("Fake self collisions: %d\n" , fake_collisions_self );
        RAVELOG_INFO("Fake sdf collisions: %d\n" , fake_collisions_sdf );

        RAVELOG_INFO("Total self collisions:  %d\n" , collisions_or_self ); 
        RAVELOG_INFO("Sphere self collisions: %d\n" , collisions_sphere_self );

        RAVELOG_INFO("Total env collisions:  %d\n" , collisions_or_env ); 
        RAVELOG_INFO("Sphere sdf collisions: %d\n" , collisions_sphere_sdf );

        RAVELOG_INFO("Total sphere collision time:  %f\n" , total_sphere ); 
        RAVELOG_INFO("Total OR collision time:  %f\n" , total_or); 

        RAVELOG_INFO("Total or fk time:  %f\n" , or_fk); 
        RAVELOG_INFO("Total sphere fk time:  %f\n" , sphere_fk ); 

        RAVELOG_INFO("Total OR self collision time:  %f\n" , or_self); 
        RAVELOG_INFO("Total OR env collision time:  %f\n" , or_env); 

        RAVELOG_INFO("Total sphere self collision time:  %f\n" , sphere_self); 
        RAVELOG_INFO("Total sphere sdf collision time:  %f\n" , sphere_sdf); 
    }
    else {
        std::clog << missed_collisions
                  << "\t" << missed_collisions_self
                  << "\t" << missed_collisions_env
                  << "\t" << fake_collisions
                  << "\t" << fake_collisions_self
                  << "\t" << fake_collisions_sdf
                  << "\t" << collisions_or 
                  << "\t" << collisions_or_self
                  << "\t" << collisions_or_env
                  << "\t" << total_or
                  << "\t" << or_fk
                  << "\t" << or_self
                  << "\t" << or_env
                  << "\t" << collisions_sphere
                  << "\t" << collisions_sphere_self
                  << "\t" << collisions_sphere_sdf
                  << "\t" << total_sphere
                  << "\t" << sphere_fk
                  << "\t" << sphere_self
                  << "\t" << sphere_sdf
                  << "\t" << pruning_missed
                  << "\n" << total_pruner 
                  << "\t" << pruner_sorting
                  << "\t" << pruner_testing 
                  << "\n";

    }

}


bool SphereCollisionHelper::isCollidedSDF( bool checkAll ){
    bool isInCollision = false; 

    for ( size_t i = 0; i < nbodies; i ++ ){
        if ( getSDFCollisions( i ) ){
            if ( !checkAll ) { return true; }
            isInCollision = true;
        }
    }

    return isInCollision;
}

bool SphereCollisionHelper::isCollidedSelf( bool checkAll ){

        
    bool isInCollision = false; 

    for ( size_t i = 0; i < nbodies; i ++ ){
        for ( size_t j = i+1; j < spheres.size(); j ++ ){
            
            if ( sphereOnSphereCollision( i, j ) ){
                if ( !checkAll ) { return true; }
                isInCollision = true;
            }
        }
    }

    return isInCollision;
}

void SphereCollisionHelper::getSpheres(){
    
    std::vector< Sphere > inactive_spheres;
    
    //a vector holding all of the pertinent bodies in the scene
    std::vector<OpenRAVE::KinBodyPtr> bodies;

    /* consider the robot kinbody, as well as all grabbed bodies */
    module->robot->GetGrabbed(bodies);
    bodies.push_back(module->environment->GetRobot( module->robot->GetName() ));
    
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
        
        //only get ignorables if the robot is the body
        if ( body.get() == module->robot.get() ){


            for (size_t j = 0; j < data_reader->ignorables.size(); j ++ ){
                int index1 = module->robot->GetLink(
                             data_reader->ignorables[j].first
                             )->GetIndex();
                int index2 = module->robot->GetLink(
                             data_reader->ignorables[j].second
                             )->GetIndex();

                ignorables.insert( getKey( index1, index2 ));
                }
        }
        

        for (size_t j = 0; j < data_reader->spheres.size(); j++ )
        {
            
            Sphere & sphere = data_reader->spheres[j];
            
            sphere.body = body.get();
            /* what robot link is this sphere attached to? */
            if (body.get() == module->robot.get()){
                sphere.link = module->robot->GetLink(sphere.linkname).get();
            }
            //the sphere is attached to a grabbed kinbody
            else{
                sphere.link = module->robot->IsGrabbing(body).get();
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
            if ( body.get() != module->robot.get() )
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
             
            bool is_active = false;
            
            //is the sphere designated as inactive in the config file?
            if ( !sphere.inactive ){
                /* is this link affected by the robot's active dofs? */
                for (size_t k = 0; k < module->n_dof; k++){
                    if ( module->robot->DoesAffect( module->active_indices[k],
                                                    sphere.linkindex )){
                        spheres.push_back( sphere );
                        is_active = true;
                        break;
                    }
                }
            }

            //if the sphere is not active, add it to the inactive
            //  vector
            if( !is_active ) { inactive_spheres.push_back( sphere); }
        }
    }
    
    //set the number of bodies equivalent to the number of active spheres.
    nbodies = spheres.size();
    //insert the inactive spheres into the spheres vector.
    spheres.insert( spheres.end(), inactive_spheres.begin(),
                                   inactive_spheres.end() );

}

inline int SphereCollisionHelper::getKey( int linkindex1,
                                          int linkindex2 ) const
{
    if (linkindex1 > linkindex2){ std::swap( linkindex1, linkindex2 );}
    return ( linkindex2 << 16 ) | (linkindex1 & 0x0000FFFF);
}

inline bool SphereCollisionHelper::ignoreSphereCollision(
                                           const Sphere & sphere1,
                                           const Sphere & sphere2) const
{
    
    //only check if the bodies are the same.
    if (sphere1.body == sphere2.body )
    {
        //extract the link indices. 
        const int first = sphere1.linkindex;
        const int second = sphere2.linkindex;

        //if the links are the same, ignore the collision.
        if( first == second ){ return true; }

        //if the body is the robot, find out if some
        //  links can be ignored.
        if ( sphere1.body == module->robot.get() ){

            const int key = getKey( first, second );
            boost::unordered_set<int>::const_iterator it =
                                               ignorables.find(key);

            //if the iterator points to a valid item,
            //  return true.
            if ( it != ignorables.end() ){ return true; }
        }
    }

    return false;
}

void SphereCollisionHelper::initPruner(){
    if ( !pruner ){
        pruner = new ArrayCollisionPruner( 0, spheres, module->sdfs );
    }  
}

inline bool SphereCollisionHelper::ignoreSphereCollision(
                                           size_t sphere_index1,
                                           size_t sphere_index2 ) const
{

    
    const Sphere & sphere1 = spheres[sphere_index1]; 
    const Sphere & sphere2 = spheres[sphere_index2];
    
    return ignoreSphereCollision( sphere1, sphere2 );
}


} //namespace
