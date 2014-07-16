#include "orchomp_collision.h"
#include "orchomp_mod.h"

namespace orchomp {


OpenRAVE::dReal computeCostFromDist( OpenRAVE::dReal dist,
                                     double epsilon,
                                     vec3 & gradient ){

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


double SphereCollisionHelper::addToGradient(const chomp::Chomp& c,
                                            chomp::MatX& g) {

    q1 = c.getTickBorderRepeat(-1).transpose();
    q2 = c.getTickBorderRepeat(0).transpose();
    
    chomp::Matx g_sdf( 1, ncspace);
    chomp::Matx g_self( 1, ncspace);    
    total_cost = 0.0;
    
    for (current_time=0; current_time<c.N; ++current_time) {
        g_sdf.setZero();
        g_self.setZero();
        
        double sdf_cost(0), self_cost(0);

        q0 = q1;
        q1 = q2;
        q2 = c.getTickBorderRepeat(current_time+1).transpose();

        cspace_vel = 0.5 * (q2 - q0) * c.inv_dt;        
        cspace_accel = (q0 - 2.0*q1 + q2) * (c.inv_dt * c.inv_dt);
        
        //Set the positions of all of the spheres, for the current configuration.
        setSpherePositions( q1 );

        //clear the jacobian matrices
        jacobians.clear();
    
        for (size_t body_index=0; body_index < nbodies; body_index++) {

            
            Sphere & current_sphere = module->active_spheres[body_index];

            //get the self collisions with the active spheres.
            for( size_t coll_index = body_index+1;
                 coll_index < nbodies;
                 coll_index++ )
            {

                Sphere & collision_sphere = module->active_spheres[coll_index];
                
                //if the spheres are on adjacent bodies, do not test them for collisions
                if( current_sphere.body == collision_sphere.body &&
                    module->areAdjacent( current_sphere.linkindex, 
                                         collision_sphere.linkindex))
                {
                    continue;
                }

                vec3 & gradient;
                double cost = sphereOnSphereCollision( body_index, coll_index, gradient);
                if ( cost <= 0.0 ){ continue; }
                
                std::vector< double > jacobian;
                jacobian.resize( nwkspace * ncspace );
                const std::vector< double > & jac1 = getJacobian( body_index ); 
                const std::vector< double > & jac2 = getJacobian( coll_index ); 

                Eigen::Map<const Eigen::MatrixXd> dx_dq
                            ( jacobian.data(), nwkspace, ncspace );

                sdf_cost += addInWorkspaceGradient( cost, gradient, dx_dq, g_self);
            }
            
            //get the contribution from the inactive_spheres
            vec3 inactive_gradient( 0,0,0);
            double inactive_cost = getInactiveCost( body_index, inactive_gradient);
            inactive_gradient /= inactive_cost;


            vec3 sdf_gradient;
            double current_sdf_cost = getSDFCollisions( body_index, sdf_gradient );
            
            if ( sdf_cost > 0 || inactive_cost > 0 ){
                //turn the jacobian into a matrix
                const std::vector< double > & jacobian = getJacobian( body_index ); 
                Eigen::Map<const Eigen::MatrixXd> dx_dq ( jacobian.data(), 
                                                          nwkspace, 
                                                          ncspace );
                if ( self_cost > 0 ) {
                    self_cost += addInWorkspaceGradient( inactive_cost,
                                                         inactive_gradient,
                                                         dx_dq, g_self );
                }
                if ( sdf_cost > 0 ){
                    sdf_cost += addInWorkspaceGradient( 
                                    current_sdf_cost, sdf_gradient, dx_dq, g_sdf);
                }
            }

        }

        g.row( current_time ) += g_sdf * obs_factor + g_self * obs_factor_self;
        total_cost += sdf_cost * obs_factor + self_cost * obs_factor_self;

    }
    return total_cost;

}

double SphereCollisionHelper::getInactiveCost( size_t body_index, vec3 & gradient_total ){

    double cost_total = 0.0;
    const Sphere & current_sphere = module->active_spheres[body_index];

    for ( size_t coll_index = nbodies;
          coll_index < nbodies + module->inactive_spheres.size();
          coll_index ++ )
    {
        const int sphere_index = coll_index - nbodies;
        const Sphere & collision_sphere = module->inactive_spheres[sphere_index];

        //if the spheres are on adjacent bodies, do not test them for collisions
        if( current_sphere.body == collision_sphere.body &&
            module->areAdjacent( current_sphere.linkindex, 
                                 collision_sphere.linkindex))
        {
            continue;
        }

        vec3 & gradient;
        double cost = sphereOnSphereCollision( body_index, coll_index, gradient);
        if ( cost <= 0.0 ){ continue; }
        
        
        gradient_total += gradient * cost;
        cost_total += cost;
    }

    return cost_total;
}

void SphereCollisionHelper::addInWorkspaceGradient(
                                    double cost, const vec3 & grad,
                                    Eigen::Map<const Eigen::MatrixXd> & dx_dq,
                                    chomp::MatX & g_current)
{

    wkspace_vel = dx_dq * cspace_vel;

    //this prevents nans from propagating. Several lines below, 
    //    wkspace_vel /= wv_norm if wv_norm is zero, nans propogate.
    if (wkspace_vel.isZero()){ return; }

    wkspace_accel = dx_dq * cspace_accel;
    
    float wv_norm = wkspace_vel.norm();
    wkspace_vel /= wv_norm;

    // add to total
    scl = wv_norm * gamma / c.inv_dt;

    P = chomp::MatX::Identity(nwkspace, nwkspace)
        - (wkspace_vel * wkspace_vel.transpose());

    K = (P * wkspace_accel) / (wv_norm * wv_norm);
   
    //          scalar * M-by-W        * (WxW * Wx1   - scalar * Wx1)
    g_current += (scl * (dx_dq.transpose() *
                 (P * cgrad - cost * K)).transpose());

    return cost * scl;
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


OpenRAVE::dReal SphereCollisionHelper::getSDFCollisions(
                                      size_t body_index,
                                      vec3 & gradient,
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
            gradient = current_gradient;
            sdf_index = i;
        }
    }

    if (dist != HUGE_VAL){

        //adjust the value of the distance by the radius of the
        //  sphere to account for the size of the sphere
        dist -= module->active_spheres[body_index].radius;
        
        if ( viewDists && dist > epsilon ){
            createCube( dist, module->sdfs[sdf_index].cube_extent,
                        sphere_positions[body_index], sdf_index);
        }

        //debugStream << "Finished SDF Collisions" << std::endl;
        return computeCostFromDist( dist, epsilon, gradient );
        

    }

    return 0.0;

}

OpenRAVE::dReal SphereCollisionHelper::sphereOnSphereCollision(
                                size_t index1, size_t index2,
                                vec3 & gradient ){

    const OpenRAVE::Vector & pos1 = sphere_positions[index1];
    const OpenRAVE::Vector & pos2 = sphere_positions[index2];

    //calculate the distance between the two centers of the spheres,
    //  also get the vector from the collision sphere to the
    //  current sphere.
    const OpenRAVE::Vector diff = sphere_positions[index1] - sphere_positions[index2];
    const OpenRAVE::dReal dist_sqrd = diff[0]*diff[0] + 
                                      diff[1]*diff[1] + 
                                      diff[2]*diff[2] ; 

    // get the distance between the spheres, with padding corresponding to
    //  epsilon.
    const OpenRAVE::dReal dist_between_centers = sqrt( dist_sqrd );
    
    const double radius1 = active_spheres[ index1 ].radius;
    const double radius2 = (index2 < active_spheres.size() ? 
                            active_spheres[index2] :
                            inactive_spheres[ index2 - active_spheres.size() ] )

    //get the actual distance between the spheres.
    const OpenRAVE::dReal dist = dist_between_centers - radius1 - radius2;
    
    //if the spheres are too far away, do not compute any cost.
    if ( dist_self > epsilon_self ){ return 0.0; }

    //get a normalized gradient vector, by dividing diff by its length
    //  and add it into the previously computed gradient
    gradient[0] = diff[0] / dist_between_centers;
    gradient[1] = diff[1] / dist_between_centers; 
    gradient[2] = diff[2] / dist_between_centers;

    //compute the cost from the distance.
    return computeCostFromDist( dist_self, epsilon_self, gradient );

}


OpenRAVE::dReal SphereCollisionHelper::getSelfCollisions(
                        size_t body_index, vec3 & gradient,
                        std::vector< OpenRAVE::dReal > & jacobian)
{

    gradient = vec3(0,0,0);
    OpenRAVE::dReal cost = 0;

    Sphere & current_sphere = module->active_spheres[ body_index ];
    
    //if the collision cache is not equal to null, get the data, and then
    //  delete the data.
    if ( current_sphere.cache != NULL ){
        gradient = current_sphere.cache->gradient;
        cost = current_sphere.cache->cost;
        jacobian = current_sphere.cache->jacobian;
        
        delete current_sphere.cache;
        current_sphere.cache = NULL;
    }

    //index over all the other spheres, and check for collisions.
    const size_t n_active = module->active_spheres.size();
    const size_t n_inactive = module->inactive_spheres.size();
    
    //do not repeat work, accomplish this by not going backwards.
    for ( size_t i = body_index+1; i < n_active; i ++ )
    {
        
        //extract the current sphere
        //  if i is less than n_active, get a sphere from active_spheres,
        //  else, get a sphere from inactive_spheres
        Sphere & collision_sphere = module->active_spheres[i];
        
        //if the spheres are attached to the same link, do not compute
        //  collisions
        if (current_sphere.linkindex == collision_sphere.linkindex){
            continue;
        }
        
        //if the spheres are on adjacent bodies, do not test them for collisions
        if( current_sphere.body == collision_sphere.body &&
            module->areAdjacent( current_sphere.linkindex, 
                                 collision_sphere.linkindex))
        {
            continue;
        }

        //get the cost and gradient of a collision;
        vec3 gradient_collision;
        double current_cost = sphereOnSphereCollision( body_index, i, gradient_collision);

        if ( current_cost <= 0 ){ continue; }
        
        //get the collision jacobian.
        const std::vector<OpenRAVE::dReal> & collision_jacobian = 
                getJacobian( collision_sphere, i);

       //add the current jacobian to the total
        for ( size_t i = 0; i < collision_jacobian.size(); i ++ ){
            jacobian[i] += current_cost * collision_jacobian[i];
        }

        cost += current_cost;

        gradient_collision *= current_cost;
        gradient += gradient_collision;
        
        //if the current collision sphere is an active sphere, cache the collision 
        //  information.
        if ( ! collision_sphere.cache ){ 
            collision_sphere.cache =
                    new CollisionDataCache( current_cost, -gradient_collision,
                                            getJacobian( current_sphere, body_index ));
        } else {
            collision_sphere.cache->addData( 
                    current_cost, -gradient_collision,
                    getJacobian( current_sphere, body_index ) );
        }
    }
    
    //debugStream << "Finished Self Collisions" << std::endl;
    return cost;
}





double SphereCollisionHelper::getCost(const chomp::MatX& q,
                                              size_t body_index,
                                              chomp::MatX& dx_dq,
                                              chomp::MatX& cgrad)
{
    
    timer.start("collision");
    
    //debugStream << "Getting Gradient Costs" << std::endl;

    //resize the matrices:
    dx_dq.resize( nwkspace, ncspace );
    cgrad.resize( nwkspace, 1 );
    
    //Set the positions of all of the spheres.
    if( body_index == 0 )
    {
        jacobians.clear();
        timer.start( "FK" );
        setSpherePositions( q );
        timer.stop( "FK" );
    }

    //COLLISION DETECTION: distance field collisions, and self collisions
    vec3 gradient_sdf(0,0,0), gradient_self(0,0,0);
    OpenRAVE::dReal cost_sdf(0), cost_self(0);
    

    timer.start("sdf collision");
    if (!module->info.noEnvironmentalCollision){
        cost_sdf = getSDFCollisions( body_index, gradient_sdf, false);
    }
    timer.stop("sdf collision");

    timer.start("self collision");
    //fill the jacobian with zeros.
    std::vector<OpenRAVE::dReal> otherJacobian( nwkspace * ncspace, 0.0);
    if (!module->info.noSelfCollision){
        cost_self = getSelfCollisions( body_index, gradient_self,
                                       otherJacobian);
    }
    timer.stop("self collision");

    if ( cost_sdf + cost_self == 0 ){ return 0.0; }
    

    timer.start( "jacobian" );

    //create the structure for the jacobian computation
    const std::vector< OpenRAVE::dReal > & jacobian =
             getJacobian( module->active_spheres[body_index], body_index );
    
    //make sure that the size of the jacobian is the correct size
    assert( jacobian.size() == ncspace * nwkspace );
    
    const double total_factor = obs_factor + obs_factor_self;

    //copy over data
    for (  size_t i = 0; i < nwkspace; i ++ ){

        //copy the gradient information
        if ( cost_self > 0.0000001 ){

            cgrad(i) = obs_factor*gradient_sdf[i]
                     + obs_factor_self*gradient_self[i] / cost_self;
        }else {
            cgrad(i) = obs_factor*gradient_sdf[i];
        }
        
        for ( size_t j = 0; j < ncspace; j ++ ){

            //copy the jacobian information

            if (cost_self > 0.0000001 ){
                dx_dq( i, j ) = (jacobian[ i * ncspace + j ]
                                  * total_factor )
                                - ( otherJacobian[ i * ncspace + j ]
                                    / cost_self * obs_factor_self );
                            
            }else {
                dx_dq( i, j ) = jacobian[ i * ncspace + j ] * obs_factor;
            }
        }
    }
    
    timer.stop( "jacobian" );

    //debugStream << "Done Gradient Cost: " <<  cost_sdf + cost_self << std::endl;
    timer.stop("collision");
    
    //setup the gradient and cost for return
    return obs_factor*cost_sdf + obs_factor_self*cost_self;

}


void SphereCollisionHelper::setSpherePositions( const chomp::MatX & q ){

    const size_t n_active = module->active_spheres.size();
    const size_t n_inactive = module->inactive_spheres.size();
    const size_t total = n_active + n_inactive;

    sphere_positions.resize( total );

    std::vector< OpenRAVE::dReal > vec;
    module->getStateAsVector( q, vec );
    module->robot->SetActiveDOFValues(vec, false);
    
    OpenRAVE::Transform t;
    int current_link_index = -1;
    OpenRAVE::KinBody * current_body = NULL;

    //get the positions of all of the spheres
    for ( size_t i=0; i < total ; i ++ )
    {
        //extract the current sphere
        //  if i is less than n_active, get a sphere from active_spheres,
        //  else, get a sphere from inactive_spheres
        const Sphere & sphere = (i < n_active ?
                                 module->active_spheres[i] :
                                 module->inactive_spheres[i-n_active]);
        
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
}


std::vector< OpenRAVE::dReal > const &
SphereCollisionHelper::getJacobian( size_t sphere_index )
{
    
    map::iterator it = jacobians.find( sphere_index );

    //if the object exists, return the thing.
    if ( it != jacobians.end() ){ return it->second; }
    
    const Sphere & sphere = ( i < module->inactive_spheres.size() ?
                         module->active_spheres[i] :
                         module->inactive_spheres[ i - module->active_spheres.size()]);

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
    return inserted_element.first->second;
}

} //namespace
