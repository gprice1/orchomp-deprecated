#ifndef __ORCHOMP_HELPER_H_
#define __ORCHOMP_HELPER_H_

#include "orchomp_mod.h"
#include <boost/unordered_map.hpp>

namespace orchomp{

class ORObserver : public chomp::ChompObserver{

    typedef std::pair< unsigned long int, std::vector<OpenRAVE::dReal> > 
                key_value_pair;
    typedef boost::unordered_map< unsigned long int,
                                  std::vector< OpenRAVE::dReal > > map;

  public: 
    ORObserver( mod * module, size_t n_spheres ) : 
             module( module ), n_spheres( n_spheres){}
    

    virtual int notify(const chomp::Chomp& chomper, 
                       chomp::ChompEventType event,
                       size_t iter,
                       double curObjective,
                       double lastObjective,
                       double constraintViolation);

    //gets the jacobian of the sphere.
    std::vector< OpenRAVE::dReal > const& 
            getJacobian( Sphere & sphere, size_t time,
                                          size_t sphere_index);    
    OpenRAVE::Vector const & getPosition(size_t time,
                                         size_t sphere);

    //get the position of the sphere
    void getPosition( size_t time, size_t object,
                      const OpenRAVE::Vector & position );

  private:
    
    //a pointer to the module
    mod * module;

    //the number of spheres, and the number of timesteps
    size_t n_spheres, n_times;
    //the map holding jacobians
    map jacobians;
    //the vector of positions
    std::vector< OpenRAVE::Vector > positions;

    size_t getKey( size_t time, size_t object );
};

}//namespace
#endif
