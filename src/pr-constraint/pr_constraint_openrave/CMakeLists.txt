cmake_minimum_required(VERSION 2.8.3)
project(pr_constraint_openrave)
find_package(catkin REQUIRED COMPONENTS pr_constraint)
find_package(OpenRAVE)
catkin_package(
   INCLUDE_DIRS include
   LIBRARIES ${PROJECT_NAME}
)
include_directories(include ${catkin_INCLUDE_DIRS} ${OpenRAVE_INCLUDE_DIRS})
link_directories(${OpenRAVE_LIBRARY_DIRS})
link_libraries(${catkin_LIBRARIES} ${OpenRAVE_LIBRARIES})
add_library(${PROJECT_NAME} src/robot-link.cpp)
install(TARGETS ${PROJECT_NAME}
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)
install(DIRECTORY include/${PROJECT_NAME}/
  DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
)
