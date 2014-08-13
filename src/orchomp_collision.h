
#ifndef _ORCHOMP_COLLISION_H_
#define _ORCHOMP_COLLISION_H_

#include "orchomp_distancefield.h"
#include "chomp-multigrid/chomp/Chomp.h"
#include <openrave/openrave.h>
#include "orchomp_kdata.h"
#include <boost/unordered_map.hpp>
#include "boost/unordered_set.hpp"

namespace orchomp{

class CollisionPruner;
class ArrayCollisionPruner;
class mod;

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
    
    ArrayCollisionPruner * pruner;

    // a pointer to the module for acces to stuff like the collision
    //  geometry
    mod * module;
    
    bool inactive_spheres_have_been_set;

    //the magnitude of the gradient update
    double gamma_sdf, gamma_self;

    /* obstacle parameters */
    //environmental collisions
    double epsilon;
    //self-collisions
    double epsilon_self;
    
    //all of this is from the chomp collision gradient helper.
    chomp::MatX q0, q1, q2;
    chomp::MatX cspace_vel, cspace_accel;
    chomp::MatX wkspace_vel, wkspace_accel;
    chomp::MatX P;
    chomp::MatX K;
    
    double inv_dt;

    //the positions of the spheres for a given configuration.
    std::vector< OpenRAVE::Vector > sphere_positions;
    map jacobians; //an unordered map of the jacobians.

    std::vector< Sphere > spheres; // the container holding spheres.
 
    //This is a set, used to hold joint pairs that can be ignored
    //  during collision checking.
    boost::unordered_set<int> ignorables;

    //used to time stuff.
    Timer timer;

    
    //________________________Public Member Functions____________________//
    
    //the constuctor needs a pointer to the robot in addition to the spaces.
    SphereCollisionHelper( size_t ncspace, mod * module,
                           double gamma=0.1, 
                           double epsilon=0.1, 
                           double obs_factor=0.7,
                           double epsilon_self=0.01,
                           double obs_factor_self=0.3);
    ~SphereCollisionHelper();

    //The main call for this class.
    //Find the workspace collision gradient for the 
    //  current trajectory.
    virtual double addToGradient(const chomp::MatX& xi,
                                 const chomp::MatX& pinit,
                                 const chomp::MatX& pgoal,
                                 double dt,
                                 chomp::MatX& g);

    //The main call for this class.
    //Find the workspace collision gradient for the 
    //  current trajectory.
    virtual double addToGradient(chomp::ConstMatMap& xi,
                                 const chomp::MatX& pinit,
                                 const chomp::MatX& pgoal,
                                 double dt,
                                 chomp::MatMap& g);


    template <class Derived>
    double getCollisionCostAndGradient( int index1, int index2,
                          const Eigen::MatrixBase<Derived> & g);

    //get the cost of the active sphere on active sphere collisions. 
    //  Add it to the C-space gradient.
    template <class Derived>
    double getSphereCost( int index1, int index2,
                          const Eigen::MatrixBase<Derived> & g);

    template <class Derived>
    double getSDFCost( int index1, int index2,
                       const Eigen::MatrixBase<Derived> & g);
    
    //Multiply the workspace gradient through the jacobian, and add it into
    //   the c-space gradient.
    template <class Derived1, class Derived2>
    double addInWorkspaceGradient(double cost, 
                                  const Eigen::Vector3d & grad,
                                  const Eigen::MatrixBase<Derived1> & dx_dq,
                                  const Eigen::MatrixBase<Derived2> & g,
                                  bool is_self_collision=false);
    
    
    //get collisions with the environment from a list of signed distance
    //  fields.
    double getSDFCollision( int sphere_index, int sdf_index,
                            Eigen::Vector3d & gradient  );
    //return true if the sphere corresponding to body_index,
    //  and the sdf corresponding to sdf_index are in collision
    bool getSDFCollision(size_t body_index, size_t sdf_index);
    bool getSDFCollisions( size_t body_index );

    //calculate the cost and direction for a collision between two spheres.
    double sphereOnSphereCollision( size_t index1, size_t index2,
                                    Eigen::Vector3d & direction,
                                    bool ignore=true);
    bool sphereOnSphereCollision( size_t index1, size_t index2,
                                  bool ignore = true);
    
    //returns true if the sphere is in collision with either a sphere
    //  or an sdf.
    bool isCollided();

    //gets the jacobian of the sphere.
    std::vector< OpenRAVE::dReal > const& getJacobian( size_t sphere_index); 
    //for a given configuration q, set the sphere_positions vector, to the
    //  positions of the spheres for the configuration.
    virtual void setSpherePositions( const chomp::MatX & q,
                                     bool setInactive=false);
    virtual void setSpherePositions(
                            const std::vector<OpenRAVE::dReal> & state,
                            bool setInactive=false);

    bool checkCollision( size_t body1, size_t body2 );

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
    void benchmark( int num_trials = 100,
                    bool check_all=false,
                    double pause_time = 0.0,
                    bool print = false);


  private:
    void getSpheres();
    void initPruner();

    //inline methods for ignoring sphere collisions.
    int getKey( int linkindex1, int linkindex2 ) const;

    virtual bool ignoreSphereCollision( const Sphere & sphere1,
                                        const Sphere & sphere2) const; 
    virtual bool ignoreSphereCollision( size_t sphere_index1,
                                        size_t sphere_index2 ) const;

};


} // namespace end
#endif
