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
#include <openrave/planningutils.h>

#include "orchomp_kdata.h"
#include "orchomp_mod.h"

#define DOPARSE 0
#define CREATEPARSE 1
namespace orchomp{


void mod::parseCreate(int argc, char * argv[], std::ostream& sout)
{
#if DOPARSE || CREATEPARSE
    
   /* parse command line arguments */
   for (int i=1; i<argc; i++)
   {
      if (strcmp(argv[i],"robot")==0 && i+1<argc)
      {
         
         if (robot) { 
             throw OpenRAVE::openrave_exception(
                        "Only one robot can be passed!");
         }
         OpenRAVE::RobotBasePtr rob = environment->GetRobot(argv[++i]);
         robot = rob.get();

         if (!robot) {
                throw OpenRAVE::openrave_exception(
                        "Only one robot can be passed!");
         }
      }else if (strcmp(argv[i],"n")==0 && i+1<argc) {
            info.n = atoi( argv[++i]);
      }else if (strcmp(argv[i],"n_max")==0 && i+1<argc) {
            info.n_max = atoi( argv[++i]);
      }else if (strcmp(argv[i],"alpha")==0 && i+1<argc) {
            info.alpha = atof( argv[ ++i]);
      }else if (strcmp(argv[i],"max_global_iter")==0 && i+1<argc) {
            info.max_global_iter = atoi( argv[ ++i]);
      }else if (strcmp(argv[i],"max_local_iter")==0 && i+1<argc) {
            info.max_local_iter = atoi( argv[ ++i]);
      }else if (strcmp(argv[i],"dolocal")==0 ) {
            info.doLocal = true;
      }else if (strcmp(argv[i],"noglobal")==0 ) {
            info.doGlobal = false;
      }
      else if (i<argc){
          for (; i<argc; i++) RAVELOG_ERROR("argument %s not known!\n", argv[i]);
          throw OpenRAVE::openrave_exception("Bad arguments!");
          break;
      }
   }
   
#endif
}

void mod::parseViewSpheres(int argc, char * argv[], std::ostream& sout)
{

#if DOPARSE
   /* parse command line arguments */
   for (int i = 1; i<argc; i++)
   {
      if (strcmp(argv[i],"robot")==0 && i+1<argc)
      {
         if (robot.get()){
             throw OpenRAVE::openrave_exception("Only one robot can be passed!");
         }

         RAVELOG_INFO("Getting robot named |%s|.\n", argv[i+1]);
         robot = environment->GetRobot(argv[++i]);

         if (!robot.get()) {
             throw OpenRAVE::openrave_exception("Could not find robot with that name!");
         }
         RAVELOG_INFO("Using robot %s.\n", robot->GetName().c_str());
      }

      //an error in the parsing
      else if ( i < argc ){
          for (; i<argc; i++) RAVELOG_ERROR("argument %s not known!\n", argv[i]);
          throw OpenRAVE::openrave_exception("Bad arguments!");
          break;
      }
   }

   /* check that we have everything */
   if (!robot.get()) { 
       throw OpenRAVE::openrave_exception("Did not pass all required args!");
   }

}
void mod::parseIterate(int argc, char * argv[], std::ostream& sout)
{
}
void mod::parseGetTraj(int argc, char * argv[], std::ostream& sout)
{
}
void mod::parseDestroy(int argc, char * argv[], std::ostream& sout)
{
}
void mod::parseComputeDistanceField(int argc, char * argv[], std::ostream& sout)
{

       /* parse command line arguments */
   for (int i = 1; i<argc; i++)
   {
      if (strcmp(argv[i],"kinbody")==0 && i+1<argc)
      {
         if (kinbody.get()){
             throw OpenRAVE::openrave_exception("Only one kinbody can be passed!");
         }
         kinbody = environment->GetKinBody(argv[++i]);
         if (!kinbody.get()){
             throw OpenRAVE::openrave_exception("Could not find kinbody with that name!");
         }
      } else if (strcmp(argv[i],"aabb_padding")==0 && i+1<argc){
         aabb_padding = atof(argv[++i]);
      }else if (strcmp(argv[i],"cube_extent")==0 && i+1<argc){
         cube_extent = atof(argv[++i]);
      }else if (strcmp(argv[i],"cache_filename")==0 && i+1<argc){
         cache_filename = argv[++i];

      //if there is an error in the parsing such as an unknown arg or a bad one.
      } else if (i<argc){
         for (; i<argc; i++){
             RAVELOG_ERROR("argument %s not known!\n", argv[i]);
         }
         throw OpenRAVE::openrave_exception("Bad arguments!");
         break;
      }
   }
   
   RAVELOG_INFO("Using kinbody %s.\n", kinbody->GetName().c_str());
   RAVELOG_INFO("Using aabb_padding |%f|.\n", aabb_padding);
   RAVELOG_INFO("Using cube_extent |%f|.\n", cube_extent);
   RAVELOG_INFO("Using cache_filename |%s|.\n", cache_filename ? cache_filename : "none passed");

   /* check that we have everything */
   if (!kinbody.get()) throw OpenRAVE::openrave_exception("Did not pass all required args!");
   
   /* make sure we don't already have an sdf loaded for this kinbody */
   for (int i=0; i < sdfs.size(); i++){
      if (strcmp( sdfs[i].kinbody_name, kinbody->GetName().c_str()) == 0){
          throw OpenRAVE::openrave_exception("We already have an sdf for this kinbody!");
          break;
      }
   }
#endif 
}
void mod::parseAddFieldFromObsArray(int argc, char * argv[], std::ostream& sout)
{

#if DOPARSE
    /* parse command line arguments */
   for (i=1; i<argc; i++)
   {
      if (strcmp(argv[i],"kinbody")==0 && i+1<argc)
      {
         if (kinbody.get()) throw OpenRAVE::openrave_exception("Only one kinbody can be passed!");
         kinbody = this->e->GetKinBody(argv[++i]);
         if (!kinbody.get()) throw OpenRAVE::openrave_exception("Could not find kinbody with that name!");
      }
      else if (strcmp(argv[i],"obsarray")==0 && i+1<argc)
         sscanf(argv[++i], "%p", &obsarray);
      else if (strcmp(argv[i],"sizes")==0 && i+1<argc)
      {
         char ** sizes_argv;
         int sizes_argc;
         cd_util_shparse(argv[++i], &sizes_argc, &sizes_argv);
         if (sizes_argc != 3) { free(sizes_argv); throw OpenRAVE::openrave_exception("sizes must be length 3!"); }
         for (j=0; j<3; j++)
            sizes[j] = atoi(sizes_argv[j]);
         free(sizes_argv);
      }
      else if (strcmp(argv[i],"lengths")==0 && i+1<argc)
      {
         char ** lengths_argv;
         int lengths_argc;
         cd_util_shparse(argv[++i], &lengths_argc, &lengths_argv);
         if (lengths_argc != 3) { free(lengths_argv); throw OpenRAVE::openrave_exception("lengths must be length 3!"); }
         for (j=0; j<3; j++)
            lengths[j] = atof(lengths_argv[j]);
         free(lengths_argv);
      }
      else if (strcmp(argv[i],"pose")==0 && i+1<argc)
      {
         char ** pose_argv;
         int pose_argc;
         cd_util_shparse(argv[++i], &pose_argc, &pose_argv);
         if (pose_argc != 7) { free(pose_argv); throw OpenRAVE::openrave_exception("pose must be length 7!"); }
         for (j=0; j<7; j++)
            pose[j] = atof(pose_argv[j]);
         free(pose_argv);
      }
      else break;
   }
   if (i<argc)
   {
      for (; i<argc; i++) RAVELOG_ERROR("argument %s not known!\n", argv[i]);
      throw OpenRAVE::openrave_exception("Bad arguments!");
   }
   
   RAVELOG_INFO("Using kinbody %s.\n", kinbody->GetName().c_str());
   RAVELOG_INFO("Using obsarray %p.\n", obsarray);
   RAVELOG_INFO("Using sizes %d %d %d.\n", sizes[0], sizes[1], sizes[2]);
   RAVELOG_INFO("Using lengths %f %f %f.\n", lengths[0], lengths[1], lengths[2]);
   RAVELOG_INFO("Using pose %f %f %f %f %f %f %f.\n",
      pose[0], pose[1], pose[2], pose[3], pose[4], pose[5], pose[6]);
   
   /* check that we have everything */
   if (!kinbody.get()) throw OpenRAVE::openrave_exception("Did not pass a kinbody!");
   if (!obsarray) throw OpenRAVE::openrave_exception("Did not pass an obsarray!");
   for (i=0; i<3; i++) if (sizes[i] <= 0) break;
   if (i<3) throw OpenRAVE::openrave_exception("Didn't pass non-zero sizes!");
   for (i=0; i<3; i++) if (lengths[i] <= 0.0) break;
   if (i<3) throw OpenRAVE::openrave_exception("Didn't pass non-zero lengths!");
   cd_kin_pose_normalize(pose);
   
   /* make sure we don't already have an sdf loaded for this kinbody */
   for (i=0; i<this->n_sdfs; i++)
      if (strcmp(this->sdfs[i].kinbody_name, kinbody->GetName().c_str()) == 0)
         break;
   if (i<this->n_sdfs)
      throw OpenRAVE::openrave_exception("We already have an sdf for this kinbody!");
   
   /* copy in name */
   strcpy(sdf_new.kinbody_name, kinbody->GetName().c_str());

   /* set pose of grid wrt kinbody frame */
   cd_mat_memcpy(sdf_new.pose, pose, 7, 1);

#endif 
}
void mod::parseDestroy(int argc, char * argv[], std::ostream& sout)
{

#if DOPARSE
       /* parse arguments */
   for (i=1; i<argc; i++)
   {
      if (strcmp(argv[i],"run")==0 && i+1<argc)
      {
         if (r) throw OpenRAVE::openrave_exception("Only one run can be passed!");
         nscan = sscanf(argv[++i], "%p", &r);
         if (nscan != 1) throw OpenRAVE::openrave_exception("Could not parse r!");
      }
      else break;
   }
   if (i<argc)
   {
      for (; i<argc; i++) RAVELOG_ERROR("argument %s not known!\n", argv[i]);
      throw OpenRAVE::openrave_exception("Bad arguments!");
   }
   if (!r) throw OpenRAVE::openrave_exception("you must pass a created run!");

#endif 
}


} /* orchomp namespace */
