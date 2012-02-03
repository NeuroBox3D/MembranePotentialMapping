################################################################################
# Your description here...
################################################################################

cmake_minimum_required(VERSION 2.6)

project(MEMBRANE_POTENTIAL_MAPPING)

# include the definitions and dependencies for ug-plugins.
include("../../../ug_plugin_cmake_includes.txt")


# sources
set(SOURCES	plugin_main.cpp
			vm2ug.cpp
			bg.cpp
			my_commands.cpp
			solve.cpp
			gating.cpp
      unit_tests/src/basic_tests.cpp)
			
# get the path to annlib
option(ANNLIB "Sets the path in which annlib shall be searched." OFF)

if(NOT ANNLIB)
	message(FATAL_ERROR "ERROR: libann is required. Please set the option -DANNLIB=path_to_libann.")
endif(NOT ANNLIB)

find_library(ANNLIB_LIBS NAMES ANN PATHS ${ANNLIB}/lib)
if(${ANNLIB_LIBS} STREQUAL "ANNLIB_LIBS-NOTFOUND")
	message(FATAL_ERROR "ERROR: Couldn't find ANNLIB in the specified path: " ${ANNLIB}/lib)
else(${ANNLIB_LIBS} STREQUAL "ANNLIB_LIBS-NOTFOUND")
	include_directories(${ANNLIB}/include)		
endif(${ANNLIB_LIBS} STREQUAL "ANNLIB_LIBS-NOTFOUND")


# create a shared library from the sources and link it against ug4.
add_library(MembranePotentialMappingPlugin SHARED ${SOURCES})

if(STATIC)
	target_link_libraries(MembranePotentialMappingPlugin ug4_s ${ANNLIB_LIBS})
else(STATIC)
	target_link_libraries(MembranePotentialMappingPlugin ug4 ${ANNLIB_LIBS})
endif(STATIC)
