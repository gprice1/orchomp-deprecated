# CMakeLists.txt for hubomz
cmake_minimum_required(VERSION 2.6)

find_package( Threads )

set(CMAKE_BUILD_TYPE Debug)
project(CHOMP_MULTIGRID)

include(mzsettings.cmake)

set(BASE_LIBS
  ${OPENGL_LIBRARY}
)

set(GUI_LIBS
  ${BASE_LIBS}
  ${GLUT_LIBRARY}
  ${PNG_LIBRARY}
)

add_subdirectory(mzcommon)
include_directories(mzcommon)

add_subdirectory(chomp)
include_directories( chomp )

add_subdirectory( nlopt )
include_directories( nlopt )

add_subdirectory( demo )

link_directories ( ${CMAKE_CURRENT_SOURCE_DIR}/mzcommon
                   ${CMAKE_CURRENT_SOURCE_DIR}/chomp
                   ${CMAKE_CURRENT_SOURCE_DIR}/nlopt
                   )
