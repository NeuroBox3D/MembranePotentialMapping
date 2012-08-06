################################################################################
# CMake for the membrane_potential_mapping ug plugin 
#
#
# created by: Stephan Grein, adopted by Sebastian Reiter and Martin Scherer
# on 02. Feb, 2012
# see: docs for additional informations
#
################################################################################

# minimum required CMake version for project below
cmake_minimum_required(VERSION 2.6)

# specify project name
project(MEMBRANE_POTENTIAL_MAPPING)

# default values
set(MPMDim 3)
set(MPMTestsuite ON)
set(MPMDefault ON)
set(MPMParallel OFF)

message(STATUS "Info: Additional options for membrane_potential_mapping:")

# options for building membrane_potential_mapping (e. g. -DMPMTestsuite=OFF)
option(MPMTestsuite "Build MPMTestsuite" ${MPMTestsuite})
message(STATUS "Info: MPMTestsuite:      "${MPMTestsuite}" (options are: ON, OFF)")

option(MPMDefault "Build MPMDefault    " ${MPMDefault})
message(STATUS "Info: MPMDefault:        "${MPMDefault}" (options are: ON, OFF)")

option(MPMDim  "Build MPMDim:           " ${MPMDim})
message(STATUS "Info: MPMDim:            "${MPMDim}" (options are: 1, 2, 3)") 

option(MPMParallel "Build MPMParallel:  " ${MPMParallel})
message(STATUS "Info: MPMParallel:       " ${MPMParallel}" (options are: ON, OFF)")
 
# include the definitions and dependencies for ug-plugins.
include("../../../ug_plugin_cmake_includes.txt")

# python includes
include(FindPythonLibs   REQUIRED)
include(FindPythonInterp REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

# boost includes
find_package(Boost COMPONENTS filesystem REQUIRED)
find_package(Boost COMPONENTS system REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# set the sources 
set(SOURCES_simple	 plugin_main.cpp
			         vm2ug.cpp
			         transform.cpp
			         bg_simple/bg.cpp
			         bg_simple/my_commands.cpp
			         bg_simple/solve.cpp
			         bg_simple/gating.cpp
                     bg_simple/spannung.cpp)

set(SOURCES_default  plugin_main.cpp
                     vm2ug.cpp 
                     transform.cpp
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

# decide if you want to build the boost testsuite executable (MPMTestsuite)
if(${MPMTestsuite} STREQUAL "ON")
	add_executable(MPMTestsuite ${SOURCES_TEST})
	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${UG_ROOT_PATH}/bin/)
endif(${MPMTestsuite} STREQUAL "ON")

# decide which dimension you want to build (currently only DIM=3 is supported!)
if (${MPMDim} EQUAL 3)
    add_definitions(-DMPMDIM=${MPMDim})
else(${MPMDim} EQUAL 3)
    message(FATAL_ERROR "ERROR: Can currently only build with dimension 3, not ${MPMDim}.")
endif(${MPMDim} EQUAL 3)

# decide if you want to build membrane_potential_mapping with OMP support
if (${MPMParallel} STREQUAL "ON")
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        message(FATAL_ERROR "ERROR: Clang currently does not support OpenMP.")
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        SET(CMAKE_CXX_FLAGS “${CMAKE_CXX_FLAGS} -fopenmp”)
        SET(CMAKE_C_FLAGS “${CMAKE_C_FLAGS} -fopenmp”)
    endif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
endif(${MPMParallel} STREQUAL "ON")

# create a shared library from the sources and link it against ug4.
if(buildEmbeddedPlugins)
# decide which build of the plugin you want to build
    if(${MPMDefault} STREQUAL "ON")
	    EXPORTSOURCES(plugins/experimental/membrane_potential_mapping ${SOURCES_MPMDefault})
	else(${MPMDefault} STREQUAL "ON")
	    EXPORTSOURCES(plugins/experimental/membrane_potential_mapping ${SOURCES_simple})
	endif(${MPMDefault})
	
	if(${MPMTestsuite} STREQUAL "ON")
		target_link_libraries (MPMTestsuite MembranePotentialMapping ug4_s ${PYTHON_LIBRARIES}  ${Boost_LIBRARIES})
	endif(${MPMTestsuite} STREQUAL "ON")
else(buildEmbeddedPlugins)


# decide which style of the plugin you want to build: note soon only one style will be available, BG channels will be put fully into either plasma_membrane or voltage_gated_channels.
    if(${MPMDefault} STREQUAL "ON")
        add_library(MembranePotentialMapping SHARED ${SOURCES_default})
        add_definitions(-DMPMDEFAULT) 
    else(${MPMDefault} STREQUAL "ON")
        add_library(MembranePotentialMapping SHARED ${SOURCES_simple})
    endif(${MPMDefault} STREQUAL "ON")
    
	target_link_libraries(MembranePotentialMapping ug4 ${ANNLIB_LIBS} ${PYTHON_LIBRARIES}  ${Boost_LIBRARIES})
	
	if(${MPMTestsuite} STREQUAL "ON")
		target_link_libraries (MPMTestsuite MembranePotentialMapping ug4 ${PYTHON_LIBRARIES}  ${Boost_LIBRARIES})
	endif(${MPMTestsuite} STREQUAL "ON")
	
endif(buildEmbeddedPlugins)

# ENDOF CMAKE