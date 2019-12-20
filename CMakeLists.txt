cmake_minimum_required(VERSION 3.3)
project(openfoam)

add_definitions("-Wall" "-g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# All directories with H files
#MACRO(HEADER_DIRECTORIES return_list)
#    FILE(GLOB_RECURSE new_list *.H)
#    SET(dir_list "")
#    FOREACH(file_path ${new_list})
#        GET_FILENAME_COMPONENT(dir_path ${file_path} PATH)
#        SET(dir_list ${dir_list} ${dir_path})
#    ENDFOREACH()
#    LIST(REMOVE_DUPLICATES dir_list)
#    SET(${return_list} ${dir_list})
#ENDMACRO()
#
#HEADER_DIRECTORIES(all_H_include_dirs)
#
#include_directories(${PROJECT_SOURCE_DIR} ${all_H_include_dirs})
include_directories(${PROJECT_SOURCE_DIR} src/OpenFOAM/lnInclude src/finiteVolume/lnInclude)

# Sources
file(GLOB_RECURSE SRC_FILES *.C)
#file(GLOB_RECURSE HDR_FILES *.H)

#message(${SRC_FILES})

add_executable(openfoam ${SRC_FILES})
