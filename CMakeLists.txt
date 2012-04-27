################################################################################
# MEMBRANE_POTENTIAL_MAPPING ugplugin
#
#
# adopted by sebastian reiter (and martin scherer)
# created on: 10.02.2012
# author: smgrein
#
################################################################################

cmake_minimum_required(VERSION 2.6)

project(MEMBRANE_POTENTIAL_MAPPING)

option(MPMTestsuite "Build MPMTestsuite" ON)
message(STATUS "Info: MPMTestsuite:      "${MPMTestsuite}" (options are: ON, OFF)")

option(FLAVOR "Build FLAVOR" ON)
message(STATUS "Info: FLAVOR:            "${FLAVOR}" (options are: ON, OFF)")
 
# include the definitions and dependencies for ug-plugins.
include("../../../ug_plugin_cmake_includes.txt")

# sources
set(SOURCES_simple	       bg_simple/plugin_main.cpp
			               bg_simple/vm2ug.cpp
			               bg_simple/bg.cpp
			               bg_simple/my_commands.cpp
			               bg_simple/solve.cpp
			               bg_simple/gating.cpp
                           bg_simple/spannung.cpp)

set(SOURCES_default  bg_default/plugin_main.cpp
                     bg_default/spannung.cpp
                     bg_default/my_commands.cpp
                     #bg_default/vm2ug.cpp
                     bg_default/bg.cpp
                     bg_default/solve.cpp
                     bg_default/gating.cpp
                     #experimental
                     vm2ug.cpp)
                    # mvec.cpp)
                     
					    
set(SOURCES_TEST unit_tests/src/tests.cpp)
			
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

# decide which flavor of the plugin you want to build
if(FLAVOR)
    add_library(MembranePotentialMappingPlugin SHARED ${SOURCES_default})
else(FLAVOR)
    add_library(MembranePotentialMappingPlugin SHARED ${SOURCES_simple})
endif(FLAVOR)

#if(${FLAVOR} STREQUAL "vm2ug")
#	add_library(MembranePotentialMappingPlugin SHARED ${SOURCES_vm2ug})
#else(${FLAVOR} STREQUAL "vm2ug")
#   if (${FLAVOR} STREQUAL "default")
#	   add_library(MembranePotentialMappingPlugin SHARED ${SOURCES_default})
#   else(${FLAVOR} STREQUAL "default")
#     add_library(MembranePotentialMappingPlugin SHARED ${SOURCES_simple})
#  endif(${FLAVOR} STREQUAL "default")
#endif(${FLAVOR} STREQUAL "vm2ug")


# decide if you want to build the boost testsuite
if(MPMTestsuite)
	add_executable(MPMTestsuite ${SOURCES_TEST})
endif(MPMTestsuite)

# create a shared library from the sources and link it against ug4.
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

# ENDOF CMAKE
