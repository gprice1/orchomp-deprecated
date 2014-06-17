cmake_minimum_required(VERSION 2.8.3)
project(pr_constraint)
find_package(catkin REQUIRED)
catkin_package(
   INCLUDE_DIRS include
   LIBRARIES ${PROJECT_NAME}
)
include_directories(include)
link_libraries(blas lapacke lapack)
add_library(${PROJECT_NAME} src/holonomic.cpp src/holonomic-circle.cpp)
install(TARGETS ${PROJECT_NAME}
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)
install(DIRECTORY include/${PROJECT_NAME}/
  DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
)
