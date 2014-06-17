cmake_minimum_required(VERSION 2.8.3)
project(pr_constraint_ompl)
find_package(catkin REQUIRED COMPONENTS pr_constraint)
find_package(OMPL REQUIRED)
catkin_package(
   INCLUDE_DIRS include
   LIBRARIES ${PROJECT_NAME}
)
include_directories(include ${catkin_INCLUDE_DIRS} ${OMPL_INCLUDE_DIRS})
link_libraries(${catkin_LIBRARIES} ${OMPL_LIBRARIES})
add_library(${PROJECT_NAME}
  src/geometric/planners/rrt/CRRT.cpp
  src/geometric/planners/rrt/CRRTConnect.cpp)
install(TARGETS ${PROJECT_NAME}
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)
install(DIRECTORY include/${PROJECT_NAME}/
  DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
)
