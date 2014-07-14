
#ifndef _ORCHOMP_COLLISION_H_
#define _ORCHOMP_COLLISION_H_

#include "orchomp_distancefield.h"
#include "chomp-multigrid/chomp/Chomp.h"
#include <openrave/openrave.h>
#include "orchomp_kdata.h"

namespace orchomp{

class mod;
//typedef vec3_t< OpenRAVE::dReal > vec3;

OpenRAVE::dReal computeCostFromDist( OpenRAVE::dReal dist,
                                     double epsilon,
                                     vec3 & gradient );

class SphereCollisionHelper : public chomp::ChompCollisionHelper{
public:
    
    // a pointer to the module for acces to stuff like the collision
    //  geometry
    mod * module;
    
    Timer timer;

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
            ChompCollisionHelper( ncspace, nwkspace, nbodies ),
            module(module),
            epsilon( 0.1 ), obs_factor(0.7), epsilon_self( 0.01 ), 
            obs_factor_self( 0.3)
    {
    }

    OpenRAVE::KinBodyPtr createCube( double dist,
                                    double size,
                                    const OpenRAVE::Vector & pos,
                                    size_t sdf_index);

    void colorFromDist( double dist, size_t sdf_index,
                        OpenRAVE::Vector & color );

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
    
    //get collisions with the environment from a list of signed distance
    //  fields.
    OpenRAVE::dReal getSDFCollisions( const Sphere & sphere,
                                      const OpenRAVE::Vector & position, 
                                      vec3 & gradient,
                                      bool viewDists=false);

    //get collisions with the current robot using its collision geometry.
    OpenRAVE::dReal getSelfCollisions( size_t body_index,
                                       const Sphere & current_sphere,
                                       const OpenRAVE::Vector & position, 
                                       vec3 & gradient,
                                       std::vector<OpenRAVE::dReal> & other);


};

} // namespace end
#endif
