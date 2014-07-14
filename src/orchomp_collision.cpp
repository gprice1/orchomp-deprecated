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
                                      const Sphere & sphere,
                                      const OpenRAVE::Vector & position, 
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
        current_dist = module->sdfs[i].getDist( position,current_gradient );
        
        if (current_dist < dist) {
            dist = current_dist;
            gradient = current_gradient;
            sdf_index = i;
        }
    }

    if (dist != HUGE_VAL){

        //adjust the value of the distance by the radius of the
        //  sphere to account for the size of the sphere
        dist -= sphere.radius;
        
        if ( viewDists && dist > epsilon ){
            createCube( dist, module->sdfs[sdf_index].cube_extent,
                        position, sdf_index);
        }

        //debugStream << "Finished SDF Collisions" << std::endl;
        return computeCostFromDist( dist, epsilon, gradient );
        

    }

    return 0.0;

}


OpenRAVE::dReal SphereCollisionHelper::getSelfCollisions(
                                       size_t body_index,
                                       const Sphere & current_sphere,
                                       const OpenRAVE::Vector & position, 
                                       vec3 & gradient,
                        std::vector< OpenRAVE::dReal > & otherJacobian)
{
    
   
    
    gradient = vec3(0,0,0);
    OpenRAVE::dReal cost = 0;
    //index over all the other spheres, and check for collisions.
    const size_t n_active = module->active_spheres.size();
    const size_t n_inactive = module->inactive_spheres.size();

    for ( size_t i=0; i < n_active + n_inactive ; i ++ )
    {
        
        //if the current index would give the current sphere, skip it.
        if ( i == body_index ) { continue;}

        //extract the current sphere
        //  if i is less than n_active, get a sphere from active_spheres,
        //  else, get a sphere from inactive_spheres
        const Sphere & collision_sphere = (i < n_active ?
                                           module->active_spheres[i] :
                                           module->inactive_spheres[i-n_active]);
        
        //if the spheres are attached to the same link, do not compute
        //  collisions
        if (current_sphere.linkindex == collision_sphere.linkindex){
            continue;
        }
        
        //TODO : is this something that I want?
        //if the spheres are on adjacent bodies, do not test them for collisions
        if( current_sphere.body == collision_sphere.body &&
            module->areAdjacent( current_sphere.linkindex, 
                                 collision_sphere.linkindex))
        {
            continue;
        }

        //get the transformation of the body that the sphere is on.
        const OpenRAVE::Transform link_xform =
                collision_sphere.link->GetTransform();

        //get the position of the collision sphere in the world.
        const OpenRAVE::Vector collision_pos =
                  link_xform * collision_sphere.position;
        
        //calculate the distance between the two centers of the spheres,
        //  also get the vector from the collision sphere to the
        //  current sphere.
        const OpenRAVE::Vector diff = (position - collision_pos);
        const OpenRAVE::dReal dist_sqrd = diff[0]*diff[0] + 
                                          diff[1]*diff[1] + 
                                          diff[2]*diff[2] ; 

        // get the distance between the spheres, with padding corresponding to
        //  epsilon.
        const OpenRAVE::dReal dist_between_centers = sqrt( dist_sqrd );

        //get a normalized gradient vector, by dividing diff by its length
        //  and add it into the previously computed gradient
        vec3 gradient_collision ( diff[0] / dist_between_centers,
                                  diff[1] / dist_between_centers, 
                                  diff[2] / dist_between_centers );

        //get the actual distance between the spheres.
        OpenRAVE::dReal dist_self = dist_between_centers
                                    - collision_sphere.radius
                                    - current_sphere.radius;
        
        //compute the cost from the distance.
        double current_cost = computeCostFromDist( dist_self, 
                                     epsilon_self,
                                     gradient_collision );

        std::vector<OpenRAVE::dReal> collision_jacobian;
        if ( current_cost > 0 ){
            module->robot->CalculateActiveJacobian(
                                    collision_sphere.linkindex,
                                    collision_pos,
                                    collision_jacobian);
        }
        for ( size_t i = 0; i < collision_jacobian.size(); i ++ ){
            otherJacobian[i] += collision_jacobian[i] * current_cost;
        }

        
        cost += current_cost;
        gradient += gradient_collision * current_cost;

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
    

    if( body_index == 0 )
    {
        timer.start( "FK" );
        std::vector< OpenRAVE::dReal > vec;
        module->getStateAsVector( q, vec );
        
        /*
        //debugStream << q << std::endl;
        //assert( module->isWithinLimits( q ));
        
        //check the state vector for nan's.
        // TODO remove this eventually
        for ( int i = 0; i < q.size(); i ++ ){
            if ( q(i) != q(i) ){
                debugStream << q <<std::endl;
                break;
            }
        }
        */
        module->robot->SetActiveDOFValues(vec, false);
        timer.stop( "FK" );
    }
    
    timer.start( "xform" );
    //loop through all of the active spheres,
    //and check their collision status.
    
    //extract the current sphere
    const Sphere & current_sphere = module->active_spheres[ body_index ];


    //get the transformation of the body that the sphere is on.
    const OpenRAVE::Transform t =
            current_sphere.link->GetTransform();
    
    //get the transformation from the body to the sphere.
    const OpenRAVE::Vector current_pos = t * current_sphere.position;
    
    //COLLISION DETECTION: distance field collisions, and self collisions
    vec3 gradient_sdf(0,0,0), gradient_self(0,0,0);
    OpenRAVE::dReal cost_sdf(0), cost_self(0);
    
    timer.stop( "xform" );

    timer.start("sdf collision");
    if (!module->info.noEnvironmentalCollision){
        cost_sdf = getSDFCollisions( current_sphere, current_pos,
                                     gradient_sdf, false);
    }
    timer.stop("sdf collision");

    timer.start("self collision");
    //fill the jacobian with zeros.
    std::vector<OpenRAVE::dReal> otherJacobian( nwkspace * ncspace, 0.0);
    if (!module->info.noSelfCollision){
        cost_self = getSelfCollisions( body_index, current_sphere,
                                   current_pos, gradient_self,
                                   otherJacobian);
    }
    timer.stop("self collision");
    
    timer.start( "jacobian" );
    //create the structure for the jacobian computation
    std::vector< OpenRAVE::dReal > jacobian;

    //calculate the jacobian 
    module->robot->CalculateActiveJacobian( current_sphere.linkindex,
                                      current_pos,
                                      jacobian);
    
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


} //namespace
