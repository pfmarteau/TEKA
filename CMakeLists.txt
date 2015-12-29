#
# This is a CMake makefile.  You can find the cmake utility and
# information about it at http://www.cmake.org
#

cmake_minimum_required(VERSION 2.6)

PROJECT(eKATS)
set(Project_HOME_DIR "./")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${Project_HOME_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${Project_HOME_DIR}/lib")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")

add_subdirectory(src)

include_directories(${Project_HOME_DIR}/src)
LINK_DIRECTORIES(${Project_HOME_DIR}/lib)

MACRO(add_bin name)
   ADD_EXECUTABLE(${name} ${Project_HOME_DIR}/${name}/${name}.cpp)
   TARGET_LINK_LIBRARIES(${name} ekats )
ENDMACRO()

#here we apply our macros 
add_bin(eKats)

