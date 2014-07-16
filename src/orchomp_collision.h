
#ifndef _ORCHOMP_COLLISION_H_
#define _ORCHOMP_COLLISION_H_

#include "orchomp_distancefield.h"
#include "chomp-multigrid/chomp/Chomp.h"
#include <openrave/openrave.h>
#include "orchomp_kdata.h"
#include <boost/unordered_map.hpp>

namespace orchomp{

class mod;
//typedef vec3_t< OpenRAVE::dReal > vec3;

OpenRAVE::dReal computeCostFromDist( OpenRAVE::dReal dist,
                                     double epsilon,
                                     vec3 & gradient );

class SphereCollisionHelper : public chomp::ChompGradientHelper{
    typedef std::pair< unsigned long int, std::vector<OpenRAVE::dReal> > 
                key_value_pair;
    typedef boost::unordered_map< unsigned long int,
                                  std::vector< OpenRAVE::dReal > > map;
public:
    
    
    size_t ncspace, nwkspace, nbodies;
    
    //all of this is from the chomp collision gradient helper.
    chomp::MatX q0, q1, q2;
    chomp::MatX cspace_vel, cspace_accel;
    chomp::MatX wkspace_vel, wkspace_accel;
    chomp::MatX P;
    chomp::MatX K;

    double gamma;
    double scl;
    // a pointer to the module for acces to stuff like the collision
    //  geometry
    mod * module;
     
    Timer timer;

    size_t current_time;
    double total_cost;
    /* obstacle parameters */
    //environmental collisions
    double epsilon;
    double obs_factor;
   
    //self-collisions
    double epsilon_self;
    double obs_factor_self;
    
    std::vector< OpenRAVE::Vector > sphere_positions;
    map jacobians; //an unordered map of the jacobians.
    
    //________________________Public Member Functions____________________//
    
    //the constuctor needs a pointer to the robot in addition to the spaces.
    SphereCollisionHelper( size_t ncspace, size_t nwkspace, size_t nbodies, 
                          mod * module) :
            ncspace(ncspace), nwkspace(nwkspace), nbodies(nbodies),
            gamma( 0.1 ),
            module(module),
            epsilon( 0.1 ), obs_factor(0.7), epsilon_self( 0.01 ), 
            obs_factor_self( 0.3)
    {
        //timer.print = true;
    }
    
    virtual double addToGradient(const chomp::Chomp& c, chomp::MatX& g);

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
    void addInWorkspaceGradient( double cost, const vec3 & grad,
                                 Eigen::Map<const Eigen::MatrixXd> & dx_dq,
                                 chomp::MatX & g_current)

    OpenRAVE::dReal sphereOnSphereCollision( size_t index1, size_t index2,
                                             OpenRAVE::Vector & direction );

    //gets the jacobian of the sphere.
    std::vector< OpenRAVE::dReal > const& 
            getJacobian( const Sphere & sphere, size_t sphere_index); 

    void setSpherePositions( const chomp::MatX & q );

    //get collisions with the environment from a list of signed distance
    //  fields.
    OpenRAVE::dReal getSDFCollisions( size_t body_index,
                                      vec3 & gradient,
                                      bool viewDists=false);

    //get collisions with the current robot using its collision geometry.
    OpenRAVE::dReal getSelfCollisions( size_t body_index,
                                       vec3 & gradient,
                                       std::vector<OpenRAVE::dReal> & jacobian);

    void workspaceToCspace(const vec3 & workspace_gradient );
    
};


class CollisionDataCache{

  public:
    
    double cost;
    vec3 gradient;
    std::vector< OpenRAVE::dReal > jacobian;

    CollisionDataCache() : cost( 0 ), gradient( 0,0,0) {}
    CollisionDataCache( double new_cost, const vec3 & new_gradient,
                        const std::vector< OpenRAVE::dReal > & new_jacobian) :
            cost( new_cost ), gradient( new_gradient )
    {
        this->jacobian.resize( new_jacobian.size() );

        for( size_t i = 0; i < new_jacobian.size(); i ++ ){
            this->jacobian[i] = new_jacobian[i] * new_cost;
        }
    }

    void addData( double new_cost, const vec3 & new_gradient,
                  const std::vector< OpenRAVE::dReal > & new_jacobian ){
        
        this->cost += new_cost;
        this->gradient += new_gradient;
        
        if ( this->jacobian.size() != new_jacobian.size() ){
            this->jacobian = new_jacobian;
        }

        for( size_t i = 0; i < new_jacobian.size(); i ++ ){
            this->jacobian[i] += new_jacobian[i] * new_cost;
        }
    }
};

} // namespace end
#endif
