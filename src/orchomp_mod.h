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

#include "orwrap.h"
#include "chomp-multigrid/chomp/Chomp.h"
#include "chomp-multigrid/chomp/Constraint.h"
#include "chomp-multigrid/chomp/ConstraintFactory.h"
#include <openrave/openrave.h>

extern "C" {
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "utils/os.h"
#include "utils/util_shparse.h"
}


namespace orchomp
{

class sdf{
    int garbage;
};

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
    bool doGlobal, doLocal;

    //a basic constructor to initialize values
    ChompInfo() :
        alpha(0.1), obstol(0.01), t_total(1.0), gamma(0.1),
        n(0), n_max(0), max_global_iter( size_t(-1) ), 
        max_local_iter( size_t(-1)), doGlobal( true ),
        doLocal( false)
        {}
};

//this does nothing right now.
class ORConstraintFactory : public chomp::ConstraintFactory {
    
    virtual chomp::Constraint* getConstraint(size_t t, size_t total) {
        return new chomp::NullConstraint();
    }

};


class Sphere{
    public: 
        // The radius of the sphere
        double radius;

        //a pointer to the kinBody that
        //the sphere comes off of.
        OpenRAVE::KinBody::Link * robot_link;

        //The index of the link
        int robot_linkindex;   

        //The transform between the coordinates of the
        //  robot link and the sphere. Since it is a sphere,
        //  rotation is meaningless, so it should just be an
        //  xyz vector.
        double pose[3];

};


class SphereCollisionHelper : public chomp::ChompCollisionHelper{
public:
    
   //a pointer to the robot model. 
   OpenRAVE::RobotBase * robot;
   
   //a vector containing the collision geometry
   std::vector<Sphere> spheres;

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
   ORConstraintFactory * factory;
   
   //This holds basic info relating to an individual 
   //   run of chomp
   ChompInfo info;
    
   //these are useful for interfacing with the robot.
   OpenRAVE::RobotBase * robot;
   OpenRAVE::KinBodyPtr kinbody;
   std::vector< int > active_indices;
   
   //This holds information pertinent to the collision geometry, and
   //   assists chomp in computing collision gradients.
   SphereCollisionHelper * sphere_collider;
   chomp::ChompCollGradHelper * collisionHelper;


   //This vector holds all of the sdf's.
   std::vector< sdf > sdfs;

   //this is a pointer to the chomp class that will pull most of the
   //   weight.
   chomp::Chomp * chomper;
    


   //_______________________PUBLIC MEMBER FUNCTIONS___________________//
   //visualize the collision geometry 
   int viewspheres(int argc, char * argv[], std::ostream& sout);

   //compute the distance field for use in collision detection, and
   //   descending the gradient out of collision
   int computedistancefield(int argc, char * argv[], std::ostream& sout);

   // NOTE : Find out what this is supposed to do
   int addfield_fromobsarray(int argc, char * argv[], std::ostream& sout);

   //
   int create(int argc, char * argv[], std::ostream& sout);

   //GO through one iteration of chomp
   int iterate(int argc, char * argv[], std::ostream& sout);

   //Get the current trajectory
   int gettraj(int argc, char * argv[], std::ostream& sout);

   //destroy the current chomp iteration.
   int destroy(int argc, char * argv[], std::ostream& sout);
    
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
    void parseCreate(int argc, char * argv[], std::ostream& sout);
    void parseViewSpheres(int argc, char * argv[], std::ostream& sout);
    void parseIterate(int argc, char * argv[], std::ostream& sout);
    void parseGetTraj(int argc, char * argv[], std::ostream& sout);
    void parseDestroy(int argc, char * argv[], std::ostream& sout);
    void parseComputeDistanceField(int argc, char * argv[], std::ostream& sout);
    void parseAddFieldFromObsArray(int argc, char * argv[], std::ostream& sout);

    // A small helper function for creating a straight line trajectory between
    //  two endpoints:
    inline void createInitialTrajectory();
    
    inline std::vector< OpenRAVE::dReal > getIthStateAsVector( size_t i );
};

void run_destroy(struct run * r);


} /* namespace orchomp */


#endif 
