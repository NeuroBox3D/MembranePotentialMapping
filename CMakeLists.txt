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

set(FLAVOR "default")

set(MPMTestsuite OFF)

set(STD_C++0x OFF) #currently compiles only with gcc (when using c++0x language features!)

IF(STD_C++0x)                                                          
  SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-std=c++0x")
ENDIF(STD_C++0x)

option(MPM "Build MPMTestsuite" OFF)
message(STATUS "Info: MPMTestsuite:      "${MPMTestsuite}" (options are: ON, OFF)")

option(FLAVOR "Build FLAVOR" default)
message(STATUS "Info: FLAVOR:            "${FLAVOR}" (options are: default, simple, vm2ug)")
 
option(STD_C++0x "Build with C++0x language features" defaultSTD_C++0x)
message(STATUS "Info: C++0x              "${STD_C++0x}" (options are: ON, OFF)")

# include the definitions and dependencies for ug-plugins.
include("../../../ug_plugin_cmake_includes.txt")

# sources (same as in bg_simple) -> need to change that
#set(SOURCES_simple	 plugin_main.cpp
#			               vm2ug.cpp
#			               bg.cpp
#			               my_commands.cpp
#			               solve.cpp
#			               gating.cpp
#                     spannung.cpp)

set(SOURCES_simple	 bg_simple/plugin_main.cpp
			               bg_simple/vm2ug.cpp
			               bg_simple/bg.cpp
			               bg_simple/my_commands.cpp
			               bg_simple/solve.cpp
			               bg_simple/gating.cpp
                     bg_simple/spannung.cpp)

set(SOURCES_vm2ug    bg_vm2ug/plugin_main.cpp
						         bg_vm2ug/vm2ug.cpp
                     bg_vm2ug/bg.cpp
						         bg_vm2ug/my_commands.cpp
					           bg_vm2ug/solve.cpp
					           bg_vm2ug/gating.cpp
					           bg_vm2ug/spannung.cpp)

set(SOURCES_default  bg_default/plugin_main.cpp
                     bg_default/spannung.cpp
                     bg_default/my_commands.cpp
                     #bg_default/vm2ug.cpp
                     bg_default/bg.cpp
                     bg_default/solve.cpp
                     bg_default/gating.cpp
                     vm2ug.cpp #experimental
                     mvec.cpp) #experimental
	
if(${FLAVOR} STREQUAL simple)
   set(SOURCES_TEST unit_tests/src/basic_tests_simple.cpp)
else(${FLAVOR} STREQUAL simple)
   set(SOURCES_TEST unit_tests/src/basic_tests_default_vm2ug.cpp)
endif(${FLAVOR} STREQUAL simple)
			
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
   if (${FLAVOR} STREQUAL "default")
	   add_library(MembranePotentialMappingPlugin SHARED ${SOURCES_default})
   else(${FLAVOR} STREQUAL "default")
     add_library(MembranePotentialMappingPlugin SHARED ${SOURCES_simple})
  endif(${FLAVOR} STREQUAL "default")
endif(${FLAVOR} STREQUAL "vm2ug")

# decide if you want to build the boost testsuite
if(MPMTestsuite)
	add_executable(MPMTestsuite ${SOURCES_TEST})
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${UG_ROOT_PATH}/bin/)
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
