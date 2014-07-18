
#ifndef _ORCHOMP_COLLISION_H_
#define _ORCHOMP_COLLISION_H_

#include "orchomp_distancefield.h"
#include "chomp-multigrid/chomp/Chomp.h"
#include <openrave/openrave.h>
#include "orchomp_kdata.h"
#include <boost/unordered_map.hpp>

namespace orchomp{

class mod;
//typedef Eigen::Vector3d_t< OpenRAVE::dReal > Eigen::Vector3d;

OpenRAVE::dReal computeCostFromDist( OpenRAVE::dReal dist,
                                     double epsilon,
                                     Eigen::Vector3d & gradient );

class SphereCollisionHelper : public chomp::ChompGradientHelper{
    typedef std::pair< unsigned long int, std::vector<OpenRAVE::dReal> > 
                key_value_pair;
    typedef boost::unordered_map< unsigned long int,
                                  std::vector< OpenRAVE::dReal > > map;
public:
    
    //the dimensions of c-space, 
    //the dimensions of workspace
    //and the number of bodies.
    size_t ncspace, nwkspace, nbodies;
    
    //all of this is from the chomp collision gradient helper.
    chomp::MatX q0, q1, q2;
    chomp::MatX cspace_vel, cspace_accel;
    chomp::MatX wkspace_vel, wkspace_accel;
    chomp::MatX P;
    chomp::MatX K;

    double gamma;
    double inv_dt;

    // a pointer to the module for acces to stuff like the collision
    //  geometry
    mod * module;
     
    Timer timer;

    /* obstacle parameters */
    //environmental collisions
    double epsilon;
    double obs_factor;
   
    //self-collisions
    double epsilon_self;
    double obs_factor_self;
    
    //the positions of the spheres for a given configuration.
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
    }

    
    //The main call for this class. Find the workspace collision gradient for the 
    //  current trajectory.
    virtual double addToGradient(const chomp::Chomp& c, chomp::MatX& g);


    //these are mostly helper functions for addToGradient. 

    //get the cost of the active sphere on active sphere collisions. 
    //  Add it to the C-space gradient.
    double getActiveCost( size_t body_index, chomp::MatX & g_self );
    
    //get the cost of collisions from active to inactive spheres.
    //  get a gradient in workspace.
    double getInactiveCost( size_t body_index, Eigen::Vector3d & gradient_total );
    
    //Multiply the workspace gradient through the jacobian, and add it into
    //   the c-space gradient.
    template <class Derived>
    double addInWorkspaceGradient( double cost, const Eigen::Vector3d & grad,
                                   const Eigen::MatrixBase<Derived> & dx_dq,
                                   chomp::MatX & cgrad);
    
    //calculate the cost and direction for a collision between two spheres.
    OpenRAVE::dReal sphereOnSphereCollision( size_t index1, size_t index2,
                                             Eigen::Vector3d & direction );
    bool sphereOnSphereCollision( size_t index1, size_t index2);

    //get collisions with the environment from a list of signed distance
    //  fields.
    OpenRAVE::dReal getSDFCollisions( size_t body_index,
                                      Eigen::Vector3d & gradient,
                                      bool viewDists=false);
    bool getSDFCollisions( size_t body_index );

    //gets the jacobian of the sphere.
    std::vector< OpenRAVE::dReal > const& 
            getJacobian( size_t sphere_index); 
  
    //for a given configuration q, set the sphere_positions vector, to the
    //  positions of the spheres for the configuration.
    void setSpherePositions( const chomp::MatX & q );


  public:
  //Public methods for visualization and testing purposes:

    OpenRAVE::KinBodyPtr createCube( double dist,
                                    double size,
                                    const OpenRAVE::Vector & pos,
                                    size_t sdf_index);

    void colorFromDist( double dist, size_t sdf_index,
                        OpenRAVE::Vector & color );

    void visualizeSDFSlice( size_t sdf_index, size_t axis,
                            size_t slice_index, double time);
  
  public: 
    bool isCollidedSDF( bool checkAll=true);
    bool isCollidedSelf( bool checkAll=true);

    void benchmark();

};


} // namespace end
#endif
