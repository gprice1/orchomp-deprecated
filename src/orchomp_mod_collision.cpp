#include "orchomp_mod.h"
#include "orchomp_kdata.h"

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
    gradient = vec3(0,0,0);
    return 0;

}

//TODO : Make this a real function
bool mod::areAdjacent( int first, int second ) const {

    //second must be larger than first, so if that is not the case, swap
    //  them.
    if ( first > second ){ std::swap( first, second ); }
    
    //for some reason, openRAVE uses this idiotic structure
    //  to store adjacency values.
    int value = ( second << 16 ) | (first & 0x0000FFFF);
    std::set<int>::iterator it = robot->GetAdjacentLinks().find(value);
    
    //if the iterator does not hold an item, 
    //  then the links are not adjacent
    if ( robot->GetAdjacentLinks().end() == it ){
        return false;
    }
    return true;
}


OpenRAVE::dReal SphereCollisionHelper::getSDFCollisions(
                                      const Sphere & sphere,
                                      const OpenRAVE::Vector & position, 
                                      vec3 & gradient )
{
    
    //if there are no sdfs for collision detection, do not calculate this
    if ( module->sdfs.size() == 0 ){
        return 0.0;
    }


    //debugStream << "Computing SDF Collisions" << std::endl;
    
    vec3 trans( position[0], position[1], position[2] );
    OpenRAVE::dReal dist = HUGE_VAL;

    //check all of the sdfs, and get the one with the least dist 
    for ( size_t i = 0; i < module->sdfs.size(); i ++ ){
        
        //if the point is not within the field, do not get the dist
        if ( !module->sdfs[i].grid.isInside( trans )){
            continue;
        }

        OpenRAVE::dReal current_dist;
        vec3 current_gradient; 

        //get the distance and gradient.
        current_dist = module->sdfs[i].getDist( trans, current_gradient );
        
        if (current_dist < dist) {
            dist = current_dist;
            gradient = current_gradient;
        }
    }
    
    //adjust the value of the distance by the radius of the
    //  sphere to account for the size of the sphere
    dist -= sphere.radius;

    //debugStream << "Finished SDF Collisions" << std::endl;
    
    return computeCostFromDist( dist, epsilon, gradient );

}


OpenRAVE::dReal SphereCollisionHelper::getSelfCollisions(
                                       size_t body_index,
                                       const Sphere & current_sphere,
                                       const OpenRAVE::Vector & position, 
                                       vec3 & gradient )
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
                  link_xform * OpenRAVE::Vector( collision_sphere.pose );
        
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
        cost += computeCostFromDist( dist_self, 
                                     epsilon_self,
                                     gradient_collision );

        gradient += gradient_collision;

    }
    
    //debugStream << "Finished Self Collisions" << std::endl;
    return cost;
}

double SphereCollisionHelper::getCost(const chomp::MatX& q,
                                              size_t body_index,
                                              chomp::MatX& dx_dq,
                                              chomp::MatX& cgrad)
{
    
    //debugStream << "Getting Gradient Costs" << std::endl;

    //resize the matrices:
    //TODO make sure that dx_dq dimensions are correct
    dx_dq.conservativeResize( nwkspace, ncspace );
    cgrad.conservativeResize( nwkspace, 1 );
    
    if( body_index == 0 ){
        std::vector< OpenRAVE::dReal > vec;
        module->getStateAsVector( q, vec );

        //debugStream << q << std::endl;
        //assert( module->isWithinLimits( q ));
        
        //check the state vector for nan's.
        for ( int i = 0; i < q.size(); i ++ ){
            if ( q(i) != q(i) ){
                debugStream << q <<std::endl;
                break;
            }
        }
        module->robot->SetActiveDOFValues(vec, 0);
    }

    //loop through all of the active spheres,
    //and check their collision status.
    
    //extract the current sphere
    const Sphere & current_sphere = module->active_spheres[ body_index ];


    //get the transformation of the body that the sphere is on.
    const OpenRAVE::Transform t =
            current_sphere.link->GetTransform();
    
    //get the transformation from the body to the sphere.
    const OpenRAVE::Vector current_pos =
                           t * OpenRAVE::Vector( current_sphere.pose );
    
    //COLLISION DETECTION: distance field collisions, and self collisions
    vec3 gradient_sdf(0,0,0), gradient_self(0,0,0);
    OpenRAVE::dReal cost_sdf(0), cost_self(0);
    
    if (!module->info.noEnvironmentalCollision){
        cost_sdf = getSDFCollisions( current_sphere, current_pos, gradient_sdf );
    }
    if (!module->info.noSelfCollision){
        cost_self = getSelfCollisions( body_index, current_sphere,
                                   current_pos, gradient_self);
    }

    //create the structure for the jacobian computation
    std::vector< OpenRAVE::dReal > jacobian;

    //calculate the jacobian 
    module->robot->CalculateActiveJacobian( current_sphere.linkindex,
                                      current_pos,
                                      jacobian);
    
    //make sure that the size of the jacobian is the correct size
    assert( jacobian.size() == ncspace * nwkspace );
     
    //copy over data
    for (  size_t i = 0; i < nwkspace; i ++ ){

        //copy the gradient information
        cgrad(i) = gradient_sdf[i] + gradient_self[i];

        for ( size_t j = 0; j < ncspace; j ++ ){

            //copy the jacobian information
            dx_dq( i, j ) = jacobian[ i * ncspace + j ];
        }
    }
    
    //debugStream << "Done Gradient Cost: " <<  cost_sdf + cost_self << std::endl;
    
    //setup the gradient and cost for return
    return obs_factor*cost_sdf + obs_factor_self*cost_self;
}


} //namespace
