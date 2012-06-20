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

message(STATUS "Info: Additional options:")

option(MPMTestsuite "Build MPMTestsuite" ON)
message(STATUS "Info: MPMTestsuite:      "${MPMTestsuite}" (options are: ON, OFF)")

option(MPMDefault "Build MPMDefault MPM" ON)
message(STATUS "Info: MPMDefault:        "${MPMDefault}" (options are: ON, OFF)")
 
# include the definitions and dependencies for ug-plugins.
include("../../../ug_plugin_cmake_includes.txt")

# sources
set(SOURCES_simple	 plugin_main.cpp
			         vm2ug.cpp 
			         bg_simple/bg.cpp
			         bg_simple/my_commands.cpp
			         bg_simple/solve.cpp
			         bg_simple/gating.cpp
                     bg_simple/spannung.cpp)

set(SOURCES_default  plugin_main.cpp
                     vm2ug.cpp 
                     bg_default/spannung.cpp
                     bg_default/my_commands.cpp
                     bg_default/bg.cpp
                     bg_default/solve.cpp
                     bg_default/gating.cpp)
                     
# set sources for unit tests
set(SOURCES_TEST unit_tests/src/tests.cpp)
		
# get the path to annlib and set annlib include path
option(ANNLIB "Sets the path in which annlib shall be searched." OFF)

if(NOT ANNLIB)
	message(FATAL_ERROR "ERROR: libann is required. Please set the option -DANNLIB=path_to_libann.")
endif(NOT ANNLIB)

find_library(ANNLIB_LIBS NAMES ANN PATHS ${ANNLIB}/lib)

if(${ANNLIB_LIBS} STREQUAL "ANNLIB_LIBS-NOTFOUND")
	message(FATAL_ERROR "ERROR: Couldn't find ANNLIB in the sporecified path: " ${ANNLIB}/lib)
else(${ANNLIB_LIBS} STREQUAL "ANNLIB_LIBS-NOTFOUND")
	include_directories(${ANNLIB}/include)		
endif(${ANNLIB_LIBS} STREQUAL "ANNLIB_LIBS-NOTFOUND")

# decide if you want to build the boost testsuite
if(MPMTestsuite)
	add_executable(MPMTestsuite ${SOURCES_TEST})
	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${UG_ROOT_PATH}/bin/)
endif(MPMTestsuite)

# create a shared library from the sources and link it against ug4.
if(buildEmbeddedPlugins)
# decide which build of the plugin you want to build
    if(${MPMDefault} STREQUAL "ON")
	    EXPORTSOURCES(plugins/experimental/membrane_potential_mapping ${SOURCES_MPMDefault})
	else(${MPMDefault} STREQUAL "ON")
	    EXPORTSOURCES(plugins/experimental/membrane_potential_mapping ${SOURCES_simple})
	endif(${MPMDefault})
	
	if(MPMTestsuite)
		target_link_libraries (MPMTestsuite MembranePotentialMappingPlugin ug4_s)
	endif(MPMTestsuite)
else(buildEmbeddedPlugins)
# decide which build of the plugin you want to build
    if(${MPMDefault} STREQUAL "ON")
        add_library(MembranePotentialMappingPlugin SHARED ${SOURCES_default})
        add_definitions(-DMPMDEFAULT) 
    else(${MPMDefault} STREQUAL "ON")
        add_library(MembranePotentialMappingPlugin SHARED ${SOURCES_simple})
    endif(${MPMDefault} STREQUAL "ON")
    
	target_link_libraries(MembranePotentialMappingPlugin ug4 ${ANNLIB_LIBS})
	
	if(MPMTestsuite)
		target_link_libraries (MPMTestsuite MembranePotentialMappingPlugin ug4)
	endif(MPMTestsuite)
endif(buildEmbeddedPlugins)

# ENDOF CMAKE