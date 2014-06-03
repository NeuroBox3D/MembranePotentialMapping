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

# check for NEURON api
if (NOT MPMNEURON)
   message(STATUS "NEURON API is off, you can turn it on by -DMPMNEURON=/path/to/NEURON/install")
else (NOT MPMNEURON)
   # check for architecture (64 bit or 32 bit on x86 systems). hack: works only if target == build host ...
      set(MPMNEURON_ARCH "i686")
      

   if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
      set(NEURON_LIBRARY_PATH ${MPMNEURON}/src/oc/modlreg.o
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/nvkludge.o 
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libnrnoc.dll
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libnrniv.dll
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libivoc.dll
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/liboc.dll
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libneuron_gnu.dll
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libscopmath.dll
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libsparse13.dll
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libsundials.dll
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libnrnmpi.dll
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libmemacs.dll
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libmeschach.dll
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libivos.dll
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libocxt.dll)
   elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
      set(NEURON_LIBRARY_PATH ${MPMNEURON}/src/oc/modlreg.o
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/nvkludge.o 
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libnrnoc.so
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libnrniv.so
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libivoc.so
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/liboc.so
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libneuron_gnu.so
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libscopmath.so
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libsparse13.so
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libsundials.so
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libnrnmpi.so
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libmemacs.so
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libmeschach.so
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libivos.so
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libocxt.so)
   elseif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
      set(NEURON_LIBRARY_PATH ${MPMNEURON}/src/oc/modlreg.o
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/nvkludge.o 
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libnrnoc.dylib
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libnrniv.dylib
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libivoc.dylib
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/liboc.dylib
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libneuron_gnu.dylib
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libscopmath.dylib
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libsparse13.dylib
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libsundials.dylib
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libnrnmpi.dylib
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libmemacs.dylib
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libmeschach.dylib
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libivos.dylib
                              ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libocxt.dylib)
   endif()
   message(STATUS "Detected architecture ${MPMNEURON_ARCH} on ${CMAKE_SYSTEM_NAME} for the NEURON API. Don't forget that SHELLTYPE readline is required.")
   # additional include directories
   include_directories(${MPMNEURON} ${MPMNEURON}/include ${MPMNEURON}/src)
   # definition for usage in source code, if we have no NEURON API then we just use the membrane potential mapping plugin without NEURON API dependent features for convenience
   add_definitions(-DMPMNEURON)
   # output information on success
   message(STATUS "NEURON API is on (by -DMPMNEURON), you can use the full NEURON hoc language in Lua now.")
   message(STATUS "")
endif(NOT MPMNEURON)

# include the definitions and dependencies for ug-plugins.
include("../../../cmake/ug_plugin_includes.cmake")

if(NOT MPMNEURON)
# set the sources 
set(SOURCES_simple
			         vm2ug.cpp
			         plugin_main.cpp
			         bg_simple/bg.cpp
			         bg_simple/my_commands.cpp
			         bg_simple/solve.cpp
			         bg_simple/gating.cpp
                     bg_simple/spannung.cpp
                     a_u_x/a_u_x_bridge.cpp
                     )

set(SOURCES_default  
                     vm2ug.cpp 
                     plugin_main.cpp
                     bg_default/spannung.cpp
                     bg_default/my_commands.cpp
                     bg_default/bg.cpp
                     bg_default/solve.cpp
                     bg_default/gating.cpp
                     a_u_x/a_u_x_bridge.cpp
                     )
                     
else(NOT MPMNEURON)
set(SOURCES_simple
			         vm2ug.cpp
			         plugin_main.cpp
			         bg_simple/bg.cpp
			         bg_simple/my_commands.cpp
			         bg_simple/solve.cpp
			         bg_simple/gating.cpp
                     bg_simple/spannung.cpp
                     a_u_x/a_u_x_bridge.cpp
                     transformator.cpp
                     hoc_command.cpp
                     )

set(SOURCES_default  
                     vm2ug.cpp 
                     plugin_main.cpp
                     bg_default/spannung.cpp
                     bg_default/my_commands.cpp
                     bg_default/bg.cpp
                     bg_default/solve.cpp
                     bg_default/gating.cpp
                     a_u_x/a_u_x_bridge.cpp
                     transformator.cpp
                     hoc_command.cpp
                     )

endif(NOT MPMNEURON)
# set sources for unit tests
set(SOURCES_TEST unit_tests/src/tests.cpp)
		
# get the path to annlib and set annlib include path
option(ANNLIB "Sets the path in which annlib shall be searched." OFF)

if(NOT ANNLIB)
	message(FATAL_ERROR "ERROR: libann is required. Please set the option -DANNLIB=path_to_libann.")
endif(NOT ANNLIB)

find_library(ANNLIB_LIBS NAMES ANN PATHS ${ANNLIB}/lib)
include_directories(${ANNLIB}/include)

if(${ANNLIB_LIBS} STREQUAL "ANNLIB_LIBS-NOTFOUND")
	message(FATAL_ERROR "ERROR: Couldn't find ANNLIB in the sporecified path: " ${ANNLIB}/lib)
else(${ANNLIB_LIBS} STREQUAL "ANNLIB_LIBS-NOTFOUND")
	include_directories(${ANNLIB}/include)		
endif(${ANNLIB_LIBS} STREQUAL "ANNLIB_LIBS-NOTFOUND")

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
# decide which build of the plugin you want to build
    if(${MPMVGCC} STREQUAL "ON")
	    EXPORTSOURCES(plugins/experimental/membrane_potential_mapping ${SOURCES_default})
	   ExportDependencies(${ANNLIB}/lib/libANN.a)
   ExportDefinitions(-DMPMVGCC)
   ExportIncludes(${ANNLIB}/include)
	else(${MPMVGCC} STREQUAL "ON")
	    EXPORTSOURCES(plugins/experimental/membrane_potential_mapping ${SOURCES_simple})
	  ExportDependencies(${ANNLIB}/lib/libANN.a)
   ExportDefinitions(-DMPMVGCC)
   ExportIncludes(${ANNLIB}/include)
	endif(${MPMVGCC} STREQUAL "ON")
	
	if(${MPMTestsuite} STREQUAL "ON")
		target_link_libraries (MPMTestsuite ug4 ${NEURON_LIBRARY_PATH}) 
	endif(${MPMTestsuite} STREQUAL "ON")
else(buildEmbeddedPlugins)


# decide which style of the plugin you want to build: note soon only one style will be available, all other BG channels will be maintained in voltage_gated_channels.
    if(${MPMVGCC} STREQUAL "ON")
       add_library(MembranePotentialMapping SHARED ${SOURCES_default})
       add_definitions(-DMPMVGCC) 
    else(${MPMVGCC} STREQUAL "ON")
        add_library(MembranePotentialMapping SHARED ${SOURCES_simple})
    endif(${MPMVGCC} STREQUAL "ON")
    
	target_link_libraries(MembranePotentialMapping ug4 ${ANNLIB_LIBS} ${NEURON_LIBRARY_PATH})
	
	if(${MPMTestsuite} STREQUAL "ON")
		target_link_libraries (MPMTestsuite MembranePotentialMapping ug4)
	endif(${MPMTestsuite} STREQUAL "ON")
	
endif(buildEmbeddedPlugins)

# warn of future obsolete VGCC components in this plugin
message(WARNING "The MembranePotentialMapping's VGCC implementation and unit tests will be unified in voltage_gated_channels in future.")

# ENDOF CMAKE membrane_potential_mapping
