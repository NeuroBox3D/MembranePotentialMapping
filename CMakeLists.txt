################################################################################
# CMake for the membrane_potential_mapping ug plugin 
# Please note: to activate use -DMembranePotentialMapping=ON during CMake config
#
# created by: Stephan Grein, adopted by Sebastian Reiter and Martin Scherer
# on 02. Feb, 2012
# see: docs for additional informations
# prerequisites for usage with VRL/ugshell: 
#      a working NEURON, NeuGen and NeuTria installation is needed on runtime,
#      additional make sure that ugshell is linked against libpython.so 
#      (see ldd ugshell for further)
#
################################################################################

# minimum required CMake version for project below
cmake_minimum_required(VERSION 2.6)

# specify project name
project(MEMBRANE_POTENTIAL_MAPPING)

# default values
set(MPMTestsuite ON)
set(MPMDefault ON)
set(MPMParallel OFF)

message(STATUS "")
message(STATUS "Info: Additional options for the MembranePotentialMapping (MPM) plugin:")

# options for building membrane_potential_mapping (e. g. -DMPMTestsuite=OFF)
option(MPMTestsuite "Build Testsuite" ${MPMTestsuite})
message(STATUS "Info: Testsuite:   "${MPMTestsuite}" (options are: ON, OFF)")

option(MPMDefault "Build VGCC:    " ${MPMDefault})
message(STATUS "Info: VGCC:        "${MPMDefault}" (options are: ON, OFF)")

option(MPMParallel "Build OpenMP:  " ${MPMParallel})
message(STATUS "Info: OpenMP:      " ${MPMParallel}" (options are: ON, OFF)")

message(STATUS "")
 
# include the definitions and dependencies for ug-plugins.
include("../../../cmake/ug_plugin_includes.cmake")

include(FindPythonLibs  REQUIRED)
include(FindPythonInterp  REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
# fixes include issue (Python.h is NOT found) on some MacOS X systems (>= 10.7.4 (Lion)) 
include_directories(${PYTHON_INCLUDE_PATH})

# boost includes (for now required - FIXME: remove that non-header-only components!!!)
find_package(Boost COMPONENTS filesystem REQUIRED)
find_package(Boost COMPONENTS system REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
message(STATUS "")

# set the sources 
set(SOURCES_simple	 transform.cpp
			         vm2ug.cpp
			         plugin_main.cpp
			         bg_simple/bg.cpp
			         bg_simple/my_commands.cpp
			         bg_simple/solve.cpp
			         bg_simple/gating.cpp
                     bg_simple/spannung.cpp
                     aux/aux_bridge.cpp
                     )

set(SOURCES_default  transform.cpp
                     vm2ug.cpp 
                     plugin_main.cpp
                     bg_default/spannung.cpp
                     bg_default/my_commands.cpp
                     bg_default/bg.cpp
                     bg_default/solve.cpp
                     bg_default/gating.cpp
                     aux/aux_bridge.cpp
                     )
                     
# set sources for unit tests
set(SOURCES_TEST unit_tests/src/tests.cpp)
		
# get the path to annlib and set annlib include path
option(ANNLIB "Sets the path in which annlib shall be searched." OFF)

if(NOT ANNLIB)
	message(FATAL_ERROR "ERROR: libann is required. Please set the option -DANNLIB=path_to_libann.")
endif(NOT ANNLIB)

find_library(ANNLIB_LIBS NAMES ANN PATHS ${ANNLIB}/lib)
# fixes include issue (ANN.h is NOT found) on some MacOS X systems (>= 10.7.4 (Lion)) 
include_directories(${ANNLIB}/include)

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

# decide if you want to build membrane_potential_mapping with OMP support
if (${MPMParallel} STREQUAL "ON")
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        message(FATAL_ERROR "ERROR: Clang currently does not support OpenMP.")
    else ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
    endif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
endif(${MPMParallel} STREQUAL "ON")

# create a shared library from the sources and link it against ug4.
if(buildEmbeddedPlugins)
# decide which build of the plugin you want to build
    if(${MPMDefault} STREQUAL "ON")
	    EXPORTSOURCES(plugins/experimental/membrane_potential_mapping ${SOURCES_default})
	else(${MPMDefault} STREQUAL "ON")
	    EXPORTSOURCES(plugins/experimental/membrane_potential_mapping ${SOURCES_simple})
	endif(${MPMDefault} STREQUAL "ON")
	
	if(${MPMTestsuite} STREQUAL "ON")
		target_link_libraries (MPMTestsuite ug4 ${PYTHON_LIBRARIES} ${Boost_LIBRARIES}) # was: MPMTestSuite MembranePotentialMapping ug4_s ...
	endif(${MPMTestsuite} STREQUAL "ON")
else(buildEmbeddedPlugins)


# decide which style of the plugin you want to build: note soon only one style will be available, all other BG channels will be maintained in voltage_gated_channels.
    if(${MPMDefault} STREQUAL "ON")
       add_library(MembranePotentialMapping SHARED ${SOURCES_default})
       add_definitions(-DMPMDEFAULT) 
    else(${MPMDefault} STREQUAL "ON")
        add_library(MembranePotentialMapping SHARED ${SOURCES_simple})
    endif(${MPMDefault} STREQUAL "ON")
    
	target_link_libraries(MembranePotentialMapping ug4 ${ANNLIB_LIBS} ${PYTHON_LIBRARIES} ${Boost_LIBRARIES})
	
	if(${MPMTestsuite} STREQUAL "ON")
		target_link_libraries (MPMTestsuite MembranePotentialMapping ug4 ${PYTHON_LIBRARIES} ${Boost_LIBRARIES})
	endif(${MPMTestsuite} STREQUAL "ON")
	
endif(buildEmbeddedPlugins)

# decide which dimension you want to build (currently only DIM=3 is supported!)
if (NOT ${DIM} STREQUAL 3 AND NOT ${DIM} STREQUAL 2 AND NOT ${DIM} STREQUAL 1)
    message(WARNING "The MPM plugin requires DIM to be set to one of {1, 2, 3}.")
endif(NOT ${DIM} STREQUAL 3 AND NOT ${DIM} STREQUAL 2 AND NOT ${DIM} STREQUAL 1)


# ENDOF CMAKE
