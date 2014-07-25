
#ifndef _ORCHOMP_COLLISION_H_
#define _ORCHOMP_COLLISION_H_

#include "orchomp_distancefield.h"
#include "chomp-multigrid/chomp/Chomp.h"
#include <openrave/openrave.h>
#include "orchomp_kdata.h"
#include <boost/unordered_map.hpp>
#include "boost/unordered_set.hpp"

namespace orchomp{

class mod;
//typedef Eigen::Vector3d_t< OpenRAVE::dReal > Eigen::Vector3d;

OpenRAVE::dReal computeCostFromDist( OpenRAVE::dReal dist,
                                     double epsilon,
                                     Eigen::Vector3d & gradient );

class CollisionPruner{
  public:

    class Node{
      public:
        Node *next, *prev;
        double value;
        size_t index;
    };

    class PruningIterator{
      public:
        size_t index;
        Node * next_node;

        //go forward until, the index has been reached.
        bool hasNext(){ 
            if ( next_node->index == index ){ return false; }
            return true;
        }
        size_t getNext(){
            size_t sphere_index = next_node->index;
            next_node = next_node->next;
            return sphere_index;
        }
    };
    
    Node * head;
    std::vector< std::pair<Node, Node> > nodes;
    const int axis;
    
    
    CollisionPruner( int axis, const std::vector< Sphere > & spheres ) :
        axis( axis )
    {

        nodes.resize( spheres.size() );

        head = &(nodes[0].first);

        for ( size_t i = 0; i < spheres.size(); i ++ ){
            nodes[i].first.next = &( nodes[i].second );
            if ( i == 0 ){
                nodes[i].first.prev = NULL;
            }else { 
                nodes[i].first.prev = &( nodes[i-1].second );
            }

            nodes[i].first.index = i;
            nodes[i].second.index = i;
            
            nodes[i].second.prev = &(nodes[i].first);
            if ( i + 1 == spheres.size() ){
                nodes[i].second.next = NULL;
            }else { 
                nodes[i].second.next = &(nodes[i+1].first);
            }
        }
    }
    

    void assertSorted(){
        Node * current = head;
        size_t count = 0;
        while ( current->next != NULL ){
            assert( current->value < current->next->value );
            current = current->next;
            count ++;
        }
        assert( count == nodes.size() - 1 );
    }

    void sort( const std::vector< OpenRAVE::Vector > & positions, 
               const std::vector< Sphere > & spheres,
               int n_set_positions)
    {
        for ( size_t i = 0; i < n_set_positions; i ++ ){
            const double pos = positions[i][axis];
            const double rad = spheres[i].radius;

            nodes[i].first.value  = pos - rad;
            nodes[i].second.value = pos + rad;
        }

        insertion_sort();
        
        RAVELOG_ERROR( "Remove this assertion" );
        assertSorted();

    }

    void insertion_sort(){
        Node * sorted_list_tail = head;
        Node * current_element = head->next;

        while ( current_element != NULL ){

            //if the element is out of place then sort it, if it is in place,
            //  setup the next iteration
            if ( sorted_list_tail->value > current_element->value ){

                //save the next element to be processed.
                Node * next_element = current_element->next;
                
                //perform insertion
                insert( sorted_list_tail, current_element );

                //set the current element for the next run.
                current_element = next_element;

            }else {
                //setup the next step of insertion sort;
                sorted_list_tail->next = current_element;
                sorted_list_tail = current_element;
                current_element = current_element->next;
            }
        }
    }
    

    void insert( Node * list_tail, Node * current_element ){
        
        //if the prev is null, then we have hit the end of the list, so 
        //  we stop.
        const double value = current_element->value;

        while ( list_tail->prev != NULL ){
            list_tail = list_tail->prev;
            
            //we have found the correct spot for the current element,
            //  it should be placed, after list_tail, because 
            //  its value is larger.
            if ( list_tail->value <= value ){

                //set the pointers of the current element.
                current_element->next = list_tail->next;
                current_element->prev = list_tail;
                
                //set the pointers of the elements surrounding 
                //  current_element
                list_tail->next->prev = current_element;
                list_tail->next = current_element;
            }
        }
        head = current_element;
        current_element->prev = NULL;
        current_element->next = list_tail;
        list_tail->prev = current_element;
    }

    PruningIterator getIterator( size_t sphere_index ){
        PruningIterator iter;
        iter.index = sphere_index;
        //get the start node corresponding to the current sphere.
        iter.next_node = nodes[sphere_index].first.next;

        return iter;
    }

};


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
    
    CollisionPruner * pruner;

    // a pointer to the module for acces to stuff like the collision
    //  geometry
    mod * module;
    
    bool inactive_spheres_have_been_set;

    //the magnitude of the gradient update
    double gamma;

    /* obstacle parameters */
    //environmental collisions
    double epsilon, obs_factor;
    //self-collisions
    double epsilon_self, obs_factor_self;
    
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

    //The main call for this class.
    //Find the workspace collision gradient for the 
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
    virtual void setSpherePositions( const chomp::MatX & q,
                                     bool setInactive=false);
    virtual void setSpherePositions(
                            const std::vector<OpenRAVE::dReal> & state,
                            bool setInactive=false);

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

    //inline methods for ignoring sphere collisions.
    int getKey( int linkindex1, int linkindex2 ) const;

    virtual bool ignoreSphereCollision( const Sphere & sphere1,
                                        const Sphere & sphere2) const; 
    virtual bool ignoreSphereCollision( size_t sphere_index1,
                                        size_t sphere_index2 ) const;

};


} // namespace end
#endif
