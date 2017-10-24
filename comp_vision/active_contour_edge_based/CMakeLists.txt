cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
set (CMAKE_CXX_STANDARD 11)
set( CMAKE_BUILD_TYPE Release )
project( active_contour )

# Find OpenCV
find_package( OpenCV REQUIRED )
# Find Boost
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.58.0 COMPONENTS filesystem regex REQUIRED)

# Import headers
include_directories( ${PROJECT_SOURCE_DIR} )
include_directories( ${OpenCV_INCLUDE_DIRS} )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories(${Boost_INCLUDE_DIRS})

file( GLOB SOURCES "*.cpp" )

add_executable( active_contour ${SOURCES})
target_link_libraries( active_contour ${OpenCV_LIBS} ${Boost_LIBRARIES} )
