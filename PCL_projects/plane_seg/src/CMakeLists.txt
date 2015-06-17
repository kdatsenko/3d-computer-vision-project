cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(plane_segmentation)

find_package(PCL 1.2 REQUIRED)
find_package(Boost 1.55 COMPONENTS program_options REQUIRED)

set(INCLUDE_DIR ../../include) # relative to source directory (CMakeLists.txt location)
set(SOURCES planar_segmentation.cpp ${INCLUDE_DIR}/point_types.h ${INCLUDE_DIR}/visualization.cpp)

include_directories(${PCL_INCLUDE_DIRS} ${Boost_INCLUDE_DIR} ${INCLUDE_DIR})
link_directories(${PCL_LIBRARY_DIRS} ${Boost_LIBRARY_DIR})

add_definitions(${PCL_DEFINITIONS})

add_executable (plane_segmentation ${SOURCES})
target_link_libraries (plane_segmentation ${PCL_LIBRARIES} ${Boost_LIBRARIES})
