cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
set (CMAKE_CXX_STANDARD 11)
set( CMAKE_BUILD_TYPE Release )
project( background_subtract )

# Find OpenCV
find_package( OpenCV REQUIRED )

# Import headers
include_directories( ${PROJECT_SOURCE_DIR} )
include_directories( ${OpenCV_INCLUDE_DIRS} )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

file( GLOB SOURCES "*.cpp" )

add_executable( background_subtract ${SOURCES})
target_link_libraries( background_subtract ${OpenCV_LIBS} )
