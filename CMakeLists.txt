################################################################################
# MEMBRANE_POTENTIAL_MAPPING ugplugin
#
#
# adopted by martin scherer & sebastian reiter
# created on: 10.02.2012
# author: smgrein
#
################################################################################

cmake_minimum_required(VERSION 2.6)

project(MEMBRANE_POTENTIAL_MAPPING)

option(MPM "Build MPMTestsuite" OFF)
message(STATUS "Info: MPMTestsuite:      "${MPMTestsuite}" (options are: ON, OFF)")

option(FLAVOR "Build FLAVOR" simple)
message(STATUS "Info: FLAVOR:            "${FLAVOR}" (options are: simple, vm2ug, userdef)")
 
# include the definitions and dependencies for ug-plugins.
include("../../../ug_plugin_cmake_includes.txt")

# sources
set(SOURCES_simple	 plugin_main.cpp
			               vm2ug.cpp
			               bg.cpp
			               my_commands.cpp
			               solve.cpp
			               gating.cpp
                     spannung.cpp)

set(SOURCES_vm2ug  bg_vm2ug/plugin_main.cpp
						             bg_vm2ug/vm2ug.cpp
                         bg_vm2ug/bg.cpp
						             bg_vm2ug/my_commands.cpp
					               bg_vm2ug/solve.cpp
					               bg_vm2ug/gating.cpp
					               bg_vm2ug/spannung.cpp)

set(SOURCES_userdef  bg_userdef/plugin_main.cpp
                     bg_userdef/spannung.cpp
                     vm2ug.cpp
                     bg.cpp
                     solve.cpp
                     gating.cpp)
					    
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

# decide which flavor of the plugin you want to build
if(${FLAVOR} STREQUAL "vm2ug")
	add_library(MembranePotentialMappingPlugin SHARED ${SOURCES_vm2ug})
else(${FLAVOR} STREQUAL "vm2ug")
   if (${FLAVOR} STREQUAL "userdef")
	   add_library(MembranePotentialMappingPlugin SHARED ${SOURCES_userdef})
   else(${FLAVOR} STREQUAL "userdef")
     add_library(MembranePotentialMappingPlugin SHARED ${SOURCES_simple})
  endif(${FLAVOR} STREQUAL "userdef")
endif(${FLAVOR} STREQUAL "vm2ug")

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
