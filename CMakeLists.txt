################################################################################
# MEMBRANE_POTENTIAL_MAPPING ugplugin
#
#
# adopted by martin scherer
# created on: 10.02.2012
# author: smgrein
#
################################################################################

cmake_minimum_required(VERSION 2.6)

project(MEMBRANE_POTENTIAL_MAPPING)

option(MPM "Build MPMTestsuite" OFF)
message(STATUS "Info: MPMTestsuite:      "${MPMTestsuite}" (options are: ON, OFF)")

# include the definitions and dependencies for ug-plugins.
include("../../../ug_plugin_cmake_includes.txt")


# sources
set(SOURCES	plugin_main.cpp
			vm2ug.cpp
			bg.cpp
			my_commands.cpp
			solve.cpp
			gating.cpp
            spannung.cpp)
			
set(SOURCES_TEST unit_tests/src/basic_tests.cpp)
			
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

if(MPMTestsuite)
	add_executable(MPMTestsuite ${SOURCES_TEST})
endif(MPMTestsuite)

if(STATIC)
	target_link_libraries(MembranePotentialMappingPlugin ug4_s ${ANNLIB_LIBS})
	if(MPMTestsuite)
		target_link_libraries (MPMTestsuite MembranePotentialMappingPlugin ug4_s)
	endif(MPMTestsuite)
else(STATIC)
	target_link_libraries(MembranePotentialMappingPlugin ug4 ${ANNLIB_LIBS})
	if(MPMTestsuite)
		target_link_libraries (MPMTestsuite MembranePotentialMappingPlugin ug4)
	endif(MPMTestsuite)
endif(STATIC)
