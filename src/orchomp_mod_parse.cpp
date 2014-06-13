/** \file orchomp_mod_parse.cpp
 * \brief Implementation of the orchomp module, an implementation of CHOMP
 *        using libcd.
 * \author Christopher Dellin
 * \date 2012
 */

/* (C) Copyright 2012-2013 Carnegie Mellon University */

/* This module (orchomp) is part of libcd.
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
 *
 * This handles all of the parsing for the mod class from ordchomp_mod.h
 */
 
#include <time.h>

#include "orchomp_mod.h"
#include "orchomp_kdata.h"

#define DOPARSE 0
#define CREATEPARSE 1
namespace orchomp{


inline void mod::getRandomState( chomp::MatX & state ){
    assert( n_dof > 0 );
    if ( size_t( state.cols()) != n_dof ){
        state.resize( 1, n_dof );
    }

    for ( size_t i = 0; i < n_dof; i ++ ){
         const double rand_val = double( rand() ) / double(RAND_MAX);
         state(i) = lowerJointLimits[i] 
                    + (upperJointLimits[i]-lowerJointLimits[i])
                    * rand_val;
         debugStream << "Got random joint value - index: "
                     << i << "\tLower: " << lowerJointLimits[i]
                     << "\tFound: " << state(i) 
                     << "\tUpper: " << upperJointLimits[i] << std::endl;
    }
}
         
        
void mod::parsePoint( std::istream& sinput, chomp::MatX & point ){
    if ( active_indices.size() <= 0 ){
        RAVELOG_ERROR("n_dof must be set before states can be added" );
        throw OpenRAVE::openrave_exception("Bad arguments!");
    }

    point.resize( 1, active_indices.size() );
    for ( size_t i = 0; i <active_indices.size() ; i ++ ){
        sinput >> point(i);
    }   

    debugStream << "\t\t-Added Point: " << point << std::endl;
}

void mod::parseCreate(std::ostream & sout, std::istream& sinput)
{
#if DOPARSE || CREATEPARSE
    
   std::string cmd;
   /* parse command line arguments */
   while (!sinput.eof () ){
      sinput >> cmd;
      debugStream << "\t-ExecutingCommand: " << cmd << std::endl;
      
      if (!sinput){ break; }
      if (cmd == "loadrobot"){
            std::string robot_location;
            sinput >> robot_location;
            debugStream << "Location: " << robot_location << std::endl;
            if (!environment.get() ){
                debugStream << "There is no Environment" << std::endl;
            }
            environment->Load( robot_location.c_str() );
            debugStream << "Done loading robot" << std::endl;
      }

      else if (cmd == "robot")
      {
         sinput >> robot_name;

         if (robot.get()) { 
             throw OpenRAVE::openrave_exception(
                        "Only one robot can be passed!");
         }
         robot = environment->GetRobot( robot_name.c_str() );
         active_indices = robot->GetActiveDOFIndices();
         robot->GetDOFLimits( lowerJointLimits, upperJointLimits,
                              active_indices);
         n_dof = active_indices.size();

         if (!robot.get()) {
                throw OpenRAVE::openrave_exception(
                        "Robot name not valid");
         }
      }else if (cmd =="n") {
            sinput >> info.n;
      }else if (cmd =="n_max"){
            sinput >> info.n_max;
      }else if (cmd =="alpha") {
            sinput >> info.alpha;
      }else if (cmd =="max_global_iter") {
            sinput >> info.max_global_iter;
      }else if (cmd =="max_local_iter") {
            sinput >> info.max_local_iter;
      }else if (cmd == "q0" ){
            parsePoint( sinput, q0 );
      }else if ( cmd == "q1" ){
            parsePoint( sinput, q1);
      }else if ( cmd == "randomstart" ){
            getRandomState( q0 );
            debugStream << "\t\tRandom Start: " << q0 << std::endl;
      }else if ( cmd == "randomend" ){
            getRandomState( q1 );
            debugStream << "\t\tRandom end: " << q1 << std::endl;
      }
    
      else if ( cmd == "dolocal"  ){ info.doLocal   = true;  }
      else if ( cmd == "nolocal"  ){ info.doLocal   = false; }
      else if ( cmd == "doglobal" ){ info.doGlobal  = true;  } 
      else if ( cmd == "doobserve"){ info.doObserve = true;  }
      //else if ( cmd =="noglobal") { info.doGlobal == false; }

      //error case
      else{ 
          while ( !sinput.eof() ){
              sinput >> cmd;
              RAVELOG_ERROR("argument %s not known!\n", cmd.c_str() );
          }
          throw OpenRAVE::openrave_exception("Bad arguments!");
      }
   }

   if (size_t( q0.cols() )!= active_indices.size() ){
       std::vector< OpenRAVE::dReal > values;
       robot->GetDOFValues( values, active_indices );
       vectorToMat( values, q0 );
   }


   
#endif
}

void mod::parseViewSpheres(std::ostream & sout, std::istream& sinput)
{
}
void mod::parseIterate(std::ostream & sout, std::istream& sinput)
{
}
void mod::parseGetTraj(std::ostream & sout, std::istream& sinput)
{
}
void mod::parseDestroy(std::ostream & sout, std::istream& sinput)
{
}
void mod::parseComputeDistanceField(std::ostream & sout, std::istream& sinput)
{

}
void mod::parseAddFieldFromObsArray(std::ostream & sout, std::istream& sinput)
{
}
void mod::parseExecute(std::ostream & sout, std::istream& sinput)
{
}


} /* orchomp namespace */
