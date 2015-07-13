cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(region_growing_segmentation)

find_package(PCL 1.2 REQUIRED)
find_package(Boost 1.55 COMPONENTS program_options REQUIRED)

set(INCLUDE_DIR ../../include)  # relative to source directory (CMakeLists.txt location)
set(INCLUDE_HELP ./)  # source directory
set(SOURCES region_growing_segmentation.cpp ${INCLUDE_DIR}/point_types.h ${INCLUDE_DIR}/visualization.cpp helpers.h compute_hausdorff.cpp)

include_directories(${PCL_INCLUDE_DIRS} ${Boost_INCLUDE_DIR} ${INCLUDE_DIR} ${INCLUDE_HELP})
link_directories(${PCL_LIBRARY_DIRS} ${Boost_LIBRARY_DIR})
add_definitions(${PCL_DEFINITIONS})

add_executable (region_growing_segmentation ${SOURCES})
target_link_libraries (region_growing_segmentation ${PCL_LIBRARIES} ${Boost_LIBRARIES})
