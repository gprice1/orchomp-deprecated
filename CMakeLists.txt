cmake_minimum_required(VERSION 2.6)
include($ENV{ROS_ROOT}/core/rosbuild/rosbuild.cmake)
rosbuild_init()
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")

add_subdirectory( src/chomp-multigrid )


rosbuild_add_library(orchomp
   src/orchomp.cpp
   src/orchomp_mod.cpp
   src/orchomp_kdata.cpp
   src/orwrap.cpp

   src/chomp-multigrid/chomp/Chomp.cpp   
   src/chomp-multigrid/chomp/ConstraintFactory.cpp   
   src/chomp-multigrid/chomp/Constraint.cpp   

   src/utils/os.c
   src/utils/util_shparse.c
   
   )


#target_link_libraries(orcdchomp blas lapacke lapack gsl)
