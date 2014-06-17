/** \file orcdchomp_mod.h
 * \brief Interface to the orcdchomp module, an implementation of CHOMP
 *        using libcd.
 * \author Christopher Dellin
 * \edited by Temple Price for use with multigrid chomp
 * \date 2012
 */

/* (C) Copyright 2012-2013 Carnegie Mellon University */

/* This module (orcdchomp) is part of libcd.
 *
 * This module of libcd is free software: you can redistribute it
 * and/or modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This module of libcd is distributed in the hope that it will be
 * useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * A copy of the GNU General Public License is provided with libcd
 * (license-gpl.txt) and is also available at <http://www.gnu.org/licenses/>.
 */

/* requires:
 *  - openrave/openrave.h
 * */

#ifndef _ORCHOMP_MOD_H_
#define _ORCHOMP_MOD_H_

#include "chomp-multigrid/chomp/Chomp.h"
#include "chomp-multigrid/chomp/Constraint.h"
#include "chomp-multigrid/chomp/ConstraintFactory.h"
#include "orchomp_distancefield.h"

#include <openrave/openrave.h>
#include <openrave/planningutils.h>

extern "C" {
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "utils/os.h"
}

#define DEBUG_COUT 0
#define debugStream \
    if (DEBUG_COUT) {} \
    else std::cout


namespace orchomp
{
class mod;
class Sphere;

//this is a structure used to hold and initialize values
//  for an eventual call to chomp.
class ChompInfo {
public:
    // al: alpha
    // obstol: Relative error of the objective function - 
    //          once the error gets to a low enough point,
    //          reltive to the objective function, quit chomping.
    // t_total: the total amount of timesteps that
    //          N+1 timesteps will take            
    // gamma: how much the gradient term counts for.
    double alpha, obstol, t_total, gamma;

    //n: the initial size of the trajectory,
    //n_max: the final size,
    //max_global_iter: the max # of global chomp interations,
    //max_local_iter: the max # of local smoothing iterations
    size_t n, n_max, max_global_iter, max_local_iter;

    //whether or not global and/or local chomp should
    //  be done.
    bool doGlobal, doLocal, doObserve;

    //a basic constructor to initialize values
    ChompInfo() :
        alpha(0.1), obstol(0.01), t_total(1.0), gamma(0.1),
        n(0), n_max(0), max_global_iter( size_t(-1) ), 
        max_local_iter( size_t(-1)), doGlobal( true ),
        doLocal( false), doObserve( false )
        {}
};


class SphereCollisionHelper : public chomp::ChompCollisionHelper{
public:
    
    //a pointer to the robot model. 
    OpenRAVE::RobotBase * robot;
   

    //a vector containing the collision geometry
    std::vector< Sphere > spheres;
    std::vector< Sphere > inactive_spheres;

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
   
    
    //the constuctor needs a pointer to the robot in addition to the spaces.
    SphereCollisionHelper( size_t ncspace, size_t nwkspace, size_t nbodies, 
                          OpenRAVE::RobotBase * robot) :
            ChompCollisionHelper( ncspace, nwkspace, nbodies ), robot( robot )
    {
    }

   
    //________________________Public Member Functions____________________//
    virtual double getCost(const chomp::MatX& q,         // configuration
                          size_t body_index,     // which body
                          chomp::MatX& dx_dq, // Jacobian of workspace pos (nwkspace-by-ncspace)
                          chomp::MatX& cgrad); // gradient (Jacobian transpose)
                                        // of cost wrt workspace pos (ncspace-by-1)

    int viewspheres( int argc, char * argv[], std::ostream& sout);

};



/* the module itself */
class mod : public OpenRAVE::ModuleBase
{
public:

   //____________________PUBLIC MEMBER VARIABLES____________________//
   OpenRAVE::EnvironmentBasePtr environment; /* filled on module creation */

   // NOTE : Most of the below variables were in the run structure.

   //the trajectory, start, and endpoint.
   chomp::MatX trajectory, q0, q1;
   chomp::ConstraintFactory * factory;
   
   //This holds basic info relating to an individual 
   //   run of chomp
   ChompInfo info;
  
   // This observes chomp for general debugging purposes
   chomp::DebugChompObserver observer;

   //these are useful for interfacing with the robot.
   OpenRAVE::RobotBasePtr robot;
   std::string robot_name;
   OpenRAVE::KinBodyPtr kinbody;
   std::vector< int > active_indices;
   size_t n_dof;

   //This holds information pertinent to the collision geometry, and
   //   assists chomp in computing collision gradients.
   SphereCollisionHelper * sphere_collider;
   chomp::ChompCollGradHelper * collisionHelper;

    //the upper and lower limits for the joints.
   std::vector< OpenRAVE::dReal > upperJointLimits,
                                  lowerJointLimits;

   //This vector holds all of the sdf's.
   std::vector< DistanceField > sdfs;

   //this is a pointer to the chomp class that will pull most of the
   //   weight.
   chomp::Chomp * chomper;
    
   OpenRAVE::TrajectoryBasePtr trajectory_ptr;

   //_______________________PUBLIC MEMBER FUNCTIONS___________________//
   //visualize the collision geometry 
   int viewspheres(std::ostream & sout, std::istream& sinput);

   //compute the distance field for use in collision detection, and
   //   descending the gradient out of collision
   int computedistancefield(std::ostream & sout, std::istream& sinput);

   // NOTE : Find out what this is supposed to do
   int addfield_fromobsarray(std::ostream & sout, std::istream& sinput);

   //
   int create(std::ostream & sout, std::istream& sinput);

   //GO through one iteration of chomp
   int iterate(std::ostream & sout, std::istream& sinput);

   //Get the current trajectory
   int gettraj(std::ostream & sout, std::istream& sinput);

   //destroy the current chomp iteration.
   int destroy(std::ostream & sout, std::istream& sinput);
   
   //execute a construsted trajectory
   int execute(std::ostream & sout, std::istream& sinput);
 
   //constructor that registers all of the commands to the openRave
   //   command line interface.
   mod(OpenRAVE::EnvironmentBasePtr penv);

   //Destructor
   virtual ~mod() {}
   void Destroy() { RAVELOG_INFO("module unloaded from environment\n"); }

   /* This is called on e.LoadProblem(m, 'command') */
   int main(const std::string& cmd)
       { RAVELOG_INFO("module init cmd: %s\n", cmd.c_str()); return 0; }


    // this is all helper code for the parsing.
    // The source code for these functions is in orchomp_mod_parse.cpp, 
    //      not the same file that contains many of the other functions.
    void parseCreate(std::ostream & sout, std::istream& sinput);
    void parseViewSpheres(std::ostream & sout, std::istream& sinput);
    void parseIterate(std::ostream & sout, std::istream& sinput);
    void parseGetTraj(std::ostream & sout, std::istream& sinput);
    void parseDestroy(std::ostream & sout, std::istream& sinput);
    void parseComputeDistanceField(std::ostream & sout,
                                   std::istream& sinput);
    void parseAddFieldFromObsArray(std::ostream & sout,
                                   std::istream& sinput);
    void parsePoint( std::istream & sinput, chomp::MatX & point);
    void parseExecute( std::ostream & sout , std::istream & sinput );
    // A small helper function for creating a straight line trajectory between
    //  two endpoints:
    inline void createInitialTrajectory();
    
    inline void getIthStateAsVector( size_t i,
                std::vector< OpenRAVE::dReal > & state  );
    inline void getStateAsVector( const chomp::MatX & state,
                std::vector< OpenRAVE::dReal > & vec  );
    inline void getRandomState( chomp::MatX & vec );
    inline void clampToLimits( chomp::MatX & state);
    
    void getSpheres();
    bool isWithinLimits( const chomp::MatX & mat ) const;

    void coutTrajectory() const;
    void isTrajectoryWithinLimits() const;

};


class ORConstraint : public chomp::Constraint {
  public:
    mod * module;
    int n_outputs;

    ORConstraint( mod * module) : module( module ), n_outputs(1) {}
    virtual void evaluateConstraints(const chomp::MatX& qt, 
                                     chomp::MatX& h, 
                                     chomp::MatX& H);
    virtual size_t numOutputs(){
        return n_outputs;
    }
};


//this does nothing right now.
class ORConstraintFactory : public chomp::ConstraintFactory {
  public: 

    mod * module;

    virtual chomp::Constraint* getConstraint(size_t t, size_t total){
        return new ORConstraint( module );
    }

    ORConstraintFactory( mod * module ) : module( module ){}

    virtual void evaluate( const std::vector<chomp::Constraint*>& constraints, 
                   const chomp::MatX& xi, chomp::MatX& h_tot,
                   chomp::MatX& H_tot, int step);

};



void run_destroy(struct run * r);


inline void vectorToMat(const std::vector< OpenRAVE::dReal > & vec,
                             chomp::MatX & mat )
{
    assert( vec.size() > 0 );
    mat.resize(1, vec.size() );

    for ( size_t i = 0; i < vec.size() ; i ++ ){ mat(i) = vec[i]; }
}


} /* namespace orchomp */


#endif 
