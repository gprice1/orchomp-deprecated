cmake_minimum_required(VERSION 2.8.3)
project(pr_constraint_tests)
find_package(catkin REQUIRED COMPONENTS
  pr_constraint pr_constraint_tsr pr_constraint_openrave pr_constraint_ompl)
find_package(Boost REQUIRED COMPONENTS thread)
find_package(OpenRAVE)
find_package(OMPL REQUIRED)
catkin_package()
include_directories(${catkin_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  ${OpenRAVE_INCLUDE_DIRS}
  ${OMPL_INCLUDE_DIRS})
link_directories(${OpenRAVE_LIBRARY_DIRS})
link_libraries(${catkin_LIBRARIES}
  ${Boost_LIBRARIES}
  ${OpenRAVE_LIBRARIES}
  ${OpenRAVE_CORE_LIBRARIES}
  ${OMPL_LIBRARIES})
add_executable(test-ompl-circle src/test-ompl-circle.cpp)
add_executable(test-ompl-or-tsr src/test-ompl-or-tsr.cpp)
install(TARGETS test-ompl-circle
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)
