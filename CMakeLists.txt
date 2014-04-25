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
set(MPMVGCC ON)
set(MPMOpenMP OFF)

if(buildEmbeddedPlugins)
   set(MPMTestsuite OFF)
endif(buildEmbeddedPlugins)


message(STATUS "")
message(STATUS "Info: Additional options for the MembranePotentialMapping (-DMPM) plugin:")

# options for building membrane_potential_mapping (e. g. -DMPMTestsuite=OFF)
option(MPMTestsuite "Build Testsuite" ${MPMTestsuite})
message(STATUS "Info: Testsuite:   " ${MPMTestsuite} " (options are: ON, OFF)")

option(MPMVGCC "Build VGCC:    " ${MPMVGCC})
message(STATUS "Info: VGCC:        " ${MPMVGCC} " (options are: ON, OFF)")

option(MPMOpenMP "Build OpenMP:  " ${MPMOpenMP})
message(STATUS "Info: OpenMP:      " ${MPMOpenMP} " (options are: ON, OFF)")

message(STATUS "")
 
# include the definitions and dependencies for ug-plugins.
include("../../../cmake/ug_plugin_includes.cmake")


## ANN (requirement for vm2ug) ###############################################

# get the path to annlib and set annlib include path and linker options
option(ANN "Sets the path in which ann shall be searched." OFF)

if(ANN)
    # the user can suply a ANN for his own compiled ANN   
    set(INTERNAL_ANN_PATH "${ANN}")
    message(STATUS "Info: ANN is ON, using user provided path ${ANN}")
else(ANN)    
    UGInstall("ANN" INTERNAL_ANN_PATH)
	message(STATUS "Info: ANN is ON, using uginstall provided path ${INTERNAL_ANN_PATH}")    
endif(ANN)

# find includes
FIND_PATH(ANN_INCLUDE_DIR ANN/ANN.h PATHS ${INTERNAL_ANN_PATH}/include)
IF (ANN_INCLUDE_DIR STREQUAL "ANN_INCLUDE_DIR-NOTFOUND")
	MESSAGE(FATAL_ERROR "ERROR: Could not find ANN header file ANN.h in ${INTERNAL_ANN_PATH}/include.")
ELSE (ANN_INCLUDE_DIR STREQUAL "ANN_INCLUDE_DIR-NOTFOUND")
	MESSAGE(STATUS "Info: ANN include:   ${ANN_INCLUDE_DIR}/ANN/ANN.h")
ENDIF (ANN_INCLUDE_DIR STREQUAL "ANN_INCLUDE_DIR-NOTFOUND")

# find libraries
SET(ANN_NAMES ${ANN_NAMES} ANN)
FIND_LIBRARY(ANN_LIBRARY NAMES ${ANN_NAMES} PATHS ${INTERNAL_ANN_PATH}/lib)
IF (ANN_LIBRARY STREQUAL "ANN_LIBRARY-NOTFOUND")
	MESSAGE(FATAL_ERROR "ERROR: Could not find ANN libraries.")
ELSE (ANN_LIBRARY STREQUAL "ANN_LIBRARY-NOTFOUND")
	MESSAGE(STATUS "Info: ANN library:   ${ANN_LIBRARY}")
ENDIF (ANN_LIBRARY STREQUAL "ANN_LIBRARY-NOTFOUND")

IF (ANN_LIBRARY AND ANN_INCLUDE_DIR)
	SET(ANN_LIBRARIES ${ANN_LIBRARY})
	SET(ANN_FOUND "YES")
ELSE (ANN_LIBRARY AND ANN_INCLUDE_DIR)
	SET(ANN_FOUND "NO")
ENDIF (ANN_LIBRARY AND ANN_INCLUDE_DIR)

IF (ANN_FOUND)
	include_directories(${ANN_INCLUDE_DIR})
	set(linkLibraries ${linkLibraries} ${ANN_LIBRARY})
endif(ANN_FOUND)

##############################################################################


# set the sources 
set(SOURCES_simple	 
			         #vm2ug.cpp
			         plugin_main.cpp
			         bg_simple/bg.cpp
			         bg_simple/my_commands.cpp
			         bg_simple/solve.cpp
			         bg_simple/gating.cpp
                     bg_simple/spannung.cpp
                     auxiliary/auxiliary_bridge.cpp
                     )

set(SOURCES_default  
                     #vm2ug.cpp 
                     plugin_main.cpp
                     bg_default/spannung.cpp
                     bg_default/my_commands.cpp
                     bg_default/bg.cpp
                     bg_default/solve.cpp
                     bg_default/gating.cpp
                     auxiliary/auxiliary_bridge.cpp
                     )
                     
# set sources for unit tests
set(SOURCES_TEST unit_tests/src/tests.cpp)
		

# decide if you want to build the boost testsuite executable (MPMTestsuite)
if(${MPMTestsuite} STREQUAL "ON")
	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${UG_ROOT_PATH}/bin/)
	add_executable(MPMTestsuite ${SOURCES_TEST})
endif(${MPMTestsuite} STREQUAL "ON")

# decide if you want to build membrane_potential_mapping with OMP support
if (${MPMOpenMP} STREQUAL "ON")
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        message(FATAL_ERROR "ERROR: Clang currently does not support OpenMP.")
    else ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
    endif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
endif(${MPMOpenMP} STREQUAL "ON")

# create a shared library from the sources and link it against ug4.
if(buildEmbeddedPlugins)
    exportDependencies(${ANN_LIBRARY})
        exportIncludes(${ANN_LIBRARY}/include)
    if(${MPMVGCC} STREQUAL "ON")
      exportDefinitions("-DMPMVGCC")
	    EXPORTSOURCES(plugins/experimental/membrane_potential_mapping ${SOURCES_default})
	else(${MPMVGCC} STREQUAL "ON")
	    EXPORTSOURCES(plugins/experimental/membrane_potential_mapping ${SOURCES_simple})
	endif(${MPMVGCC} STREQUAL "ON")
      EXPORTDEPENDENCIES(${ANNLIB_LIBS})
#target_link_libraries(ug4 ${ANNLIB_LIBS})
		#target_link_libraries (MPMTestsuite ug4 ${PYTHON_LIBRARIES} ${Boost_LIBRARIES}) # was: MPMTestSuite MembranePotentialMapping ug4_s ...
	#if(${MPMTestsuite} STREQUAL "ON")
#		target_link_libraries (MPMTestsuite ug4 ${PYTHON_LIBRARIES} ${Boost_LIBRARIES}) # was: MPMTestSuite MembranePotentialMapping ug4_s ...
#	endif(${MPMTestsuite} STREQUAL "ON")
else(buildEmbeddedPlugins)


# decide which style of the plugin you want to build: note soon only one style will be available, all other BG channels will be maintained in voltage_gated_channels.
    if(${MPMVGCC} STREQUAL "ON")
       add_library(MembranePotentialMapping SHARED ${SOURCES_default})
       add_definitions(-DMPMVGCC) 
    else(${MPMVGCC} STREQUAL "ON")
        add_library(MembranePotentialMapping SHARED ${SOURCES_simple})
    endif(${MPMVGCC} STREQUAL "ON")
    
	target_link_libraries(MembranePotentialMapping ug4 ${ANN_LIBRARY} ${PYTHON_LIBRARIES} ${Boost_LIBRARIES})
	
	if(${MPMTestsuite} STREQUAL "ON")
		target_link_libraries (MPMTestsuite MembranePotentialMapping ug4 ${PYTHON_LIBRARIES} ${Boost_LIBRARIES})
	endif(${MPMTestsuite} STREQUAL "ON")
	
endif(buildEmbeddedPlugins)

# warn of future obsolete VGCC components in this plugin
message(WARNING "The MembranePotentialMapping's VGCC implementation and unit tests will be unified in voltage_gated_channels in future.")

# ENDOF CMAKE
