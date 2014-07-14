

#include "orchomp_observer.h"

namespace orchomp{

//this should get all of the positions at an iteration. 
int ORObserver::notify(const chomp::Chomp& chomper, 
                       chomp::ChompEventType event,
                       size_t iter,
                       double curObjective,
                       double lastObjective,
                       double constraintViolation)

{
    jacobians.clear();

    int n_times = chomper.xi.rows();
    
    const size_t n_active = module->active_spheres.size();
    
    //make sure that the positions and transforms are the correct size.
    if (positions.size() != n_times * n_spheres){
        positions.resize( n_times * n_spheres );
    }

    for ( int i = 0; i < n_times; i ++ ){
        
        module->setActiveDOFValues( chomper.xi.row( i ));

        //get all of the positions of the spheres.
        for ( size_t j=0; j < n_spheres ; j ++ )
        {
            //extract the current sphere
            //  if i is less than n_active,
            //  get a sphere from active_spheres,
            //  else, get a sphere from inactive_spheres
            const Sphere & sphere = (j < n_active ?
                                    module->active_spheres[j] :
                                    module->inactive_spheres[j-n_active]);
            const OpenRAVE::Transform link_xform = 
                                    sphere.link->GetTransform();

            positions[ i * n_spheres + j ] = link_xform * sphere.position;
        }
    }
    return 0;
}


std::vector< OpenRAVE::dReal > const &
ORObserver::getJacobian( Sphere & sphere, size_t time,
                                        size_t sphere_index )
{
    size_t key = getKey(time, sphere_index); 
    map::iterator it = jacobians.find( key );

    //if the object exists, return the thing.
    if ( it != jacobians.end() ){ return it->second; }
    
    //if the jacobian does not exist,create and return it.
    key_value_pair new_jacobian;
    new_jacobian.first = key;
    std::pair< map::iterator, bool> inserted_element = 
                   jacobians.insert( new_jacobian );

    //actually get the jacobian
    module->robot->CalculateActiveJacobian(
                   sphere.linkindex, 
                   getPosition(time, sphere_index),
                   inserted_element.first->second);
    return inserted_element.first->second;
}


OpenRAVE::Vector const & ORObserver::getPosition(
                                   size_t time, size_t sphere)
{
    return positions[ time * n_spheres + sphere ];
}

  
size_t ORObserver::getKey( size_t time, size_t object )
{
    return (object << sizeof(size_t)/2) + time;
}

}//namespace
