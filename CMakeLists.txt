cmake_minimum_required(VERSION 2.6)
project( orchomp )

set( CMAKE_BUILD_TYPE Debug)

#$ cmake -DCMAKE_BUILD_TYPE="Debug" .
set(CMAKE_C_FLAGS_DEBUG "-O -g")
set(CMAKE_CXX_FLAGS_DEBUG "-O -g")

#$ cmake -DCMAKE_BUILD_TYPE="Release" 
set(CMAKE_C_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

include(src/chomp-multigrid/mzsettings.cmake)

find_package(OpenRAVE REQUIRED)
include_directories(${OpenRAVE_INCLUDE_DIRS})
link_directories(${OpenRAVE_LIBRARY_DIRS})

add_subdirectory( src/chomp-multigrid )

#This is all of the ros stuff.
set( ROS TRUE )
if (ROS )
    include($ENV{ROS_ROOT}/core/rosbuild/rosbuild.cmake)
    rosbuild_init()
    set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
    set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
endif( ROS )

# end of the ros stuff

set( ORCHOMPLIBS   
    src/orchomp.cpp

    src/orchomp_mod_utils.cpp
    src/orchomp_mod_parse.cpp
    src/orchomp_mod.cpp
    
    src/orchomp_kdata.cpp
    src/orchomp_distancefield.cpp
    src/orchomp_collision.cpp
    src/orchomp_constraint.cpp

    src/utils/os.c
    src/utils/util_shparse.c
    src/utils/timer.h
    )

#uncomment the below line to build normally.
if( ROS )
    rosbuild_add_library( orchomp ${ORCHOMPLIBS} )
else( ROS )
    add_library(orchomp SHARED ${ORCHOMPLIBS})
endif( ROS )

set_target_properties(orchomp PROPERTIES COMPILE_FLAGS
                      "${OpenRAVE_CXX_FLAGS}" LINK_FLAGS 
                      "${OpenRAVE_LINK_FLAGS}")
target_link_libraries(orchomp chomp mzcommon
                      gsl ${OpenRAVE_LIBRARIES})


ADD_CUSTOM_TARGET(debug
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Debug"
  )

ADD_CUSTOM_TARGET(release
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Release"
  )
