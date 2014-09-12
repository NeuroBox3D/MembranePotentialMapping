################################################################################
# Purpose:
# CMakeLists.txt for the experimental ug plugin membrane_potential_mapping, 
# here fter we refer to the plugin by the abbreviation mpm in general.
#
# Remark: 
# Adapted by a CMakeLists file from Sebastian Reiter and Andreas Vogel.

# Notes:
#   I) Activate plugin by -DMembranePotentialMapping=ON 
#  II) Activate the boost testsuite by -DMPMTestsuite=ON
# III) Activate NEURON API by -DMPMNEURON=/path/to/your/local/neuron/install
#  IV) Decide which NEURON version to use (either HEAD or pre 7.3 branch) by
#      -DMPMNEURON_REVISION=ON or -DMPMNEURON_REVISION=OFF
#   V) For now one needs to activate at a minimum the ANN C++ external library 
#      by -DANNLIB=/path/to/your/local/annlib/install

# Known bugs: 
#  I) static build won't allow mod files for now, figure out why this is the case

# Author: Stephan Grein
# Created: 2nd of February, 2012
# See: docs directory for additional informations on the mpm plugin
#
################################################################################

# minimum required CMake version for our project below
cmake_minimum_required(VERSION 2.6)

# specify project name
project(MEMBRANE_POTENTIAL_MAPPING)

# default values for cmake options
set(MPMTestsuite OFF)
set(MPMVGCC ON)
set(MPMOpenMP OFF)
set(MPMNEURON_REVISION ON)
set(MPMNEURON_STATIC ON)

# show optional arguments for the plugin
message(STATUS "")
message(STATUS "Info: Additional options for the MembranePotentialMapping (-DMPM) plugin:")

# enable the boost testsuite
option(MPMTestsuite "Build Testsuite" ${MPMTestsuite})
message(STATUS "Info: Testsuite:       " ${MPMTestsuite} " (options are: ON, OFF)")

# which channel type
option(MPMVGCC "Build VGCC:    " ${MPMVGCC})
message(STATUS "Info: VGCC:            " ${MPMVGCC} " (options are: ON, OFF)")

# build with openmp 
option(MPMOpenMP "Build OpenMP:  " ${MPMOpenMP})
message(STATUS "Info: OpenMP:          " ${MPMOpenMP} " (options are: ON, OFF)")

# use head NEURON revision or pre 7.3
option(MPMNEURON_REVISION "Use head revision of NEURON or pre 7.3 build? " ${MPMNEURON_REVISION})
message(STATUS "Info: NEURON head rev: " ${MPMNEURON_REVISION} " (options are: ON, OFF)")

option(MPMNEURON_STATIC "Use static NEURON build?" ${MPMNEURON_STATIC})
message(STATUS "Info: NEURON static:   " ${MPMNEURON_STATIC} " (options are: ON, OFF)")
message(STATUS "")

# check for the presence of the NEURON API
if (NOT MPMNEURON)
   message(STATUS "NEURON API is off, you can turn it on by -DMPMNEURON=/path/to/NEURON/install")
else (NOT MPMNEURON)
   # check for 64 or 32 bit - works only if target and build host are identically 
   if(CMAKE_SIZEOF_VOID_P EQUAL 8)
      set(MPMNEURON_ARCH "x86_64")
   else(CMAKE_SIZEOF_VOID_P EQUAL 8)
      set(MPMNEURON_ARCH "i686")
   endif(CMAKE_SIZEOF_VOID_P EQUAL 8)
   
   # check for NEURON API revision (head or pre 7.3)
   if (MPMNEURON_REVISION)
      add_definitions(-DMPMNEURON_REVISION)
       set(MPMNEURON_REVISION_NAME "post 7.3 (head)")
   else(MPMNEURON_REVISION)
       set(MPMNEURON_REVISION_NAME "pre 7.3")
   endif(MPMNEURON_REVISION)
   
   # determine system, i. e. windows, linux or mac and set appropriate linking order 
   if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
      set(MPMNEURON_SYSTEM_EXT "dll")
   elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
      set(MPMNEURON_SYSTEM_EXT "so")
   elseif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
      set(MPMNEURON_SYSTEM_EXT "dylib")
   endif()
   
   # static libraries to use? (remember to compile NEURON static!) -> override in case we have static build...
   if (MPMNEURON_STATIC)
      set(MPMNEURON_SYSTEM_EXT "a")
      # for now set the library manually in the static case (this is just for testing purposes, in general you dont need this line) remove this later
      # set(MPMNEURON "/Users/stephan/Code/hg/neuron/nrn/") -> my static build with head revision (without memacs please otherwise we need ncurses at compile/linking time
      #                                                     -> this means build NEURON without memacs support at configure (this is just some terminal cursor support, which we dont need anyways)
      #                                                     -> TODO: if build static, NEURON won't load compiled mechanisms (i. e. mod files), we need a way to figure out why
      # set(MPMNEURON "/Users/stephan/Code/git/neuron/nrn-7.3/") -> my dynamic build with pre head revision (with memacs and ncurses at run time must be present), note this okay for non-VRL builds
      # note we could check for static build, then remove from NEURON_LIBRARY_PATH the memacs, but this is not necessary since a non-available lib/file will be ignored by cmake
   endif(MPMNEURON_STATIC)

   set(NEURON_LIBRARY_PATH 
       ${MPMNEURON}/src/oc/modlreg.o
       ${MPMNEURON}/${MPMNEURON_ARCH}/lib/nvkludge.o 
       ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libnrnoc.${MPMNEURON_SYSTEM_EXT}
       ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libnrniv.${MPMNEURON_SYSTEM_EXT}
       ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libivoc.${MPMNEURON_SYSTEM_EXT}
       ${MPMNEURON}/${MPMNEURON_ARCH}/lib/liboc.${MPMNEURON_SYSTEM_EXT}
       ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libneuron_gnu.${MPMNEURON_SYSTEM_EXT}
       ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libscopmath.${MPMNEURON_SYSTEM_EXT}
       ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libsparse13.${MPMNEURON_SYSTEM_EXT}
       ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libsundials.${MPMNEURON_SYSTEM_EXT}
       ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libnrnmpi.${MPMNEURON_SYSTEM_EXT}
       ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libmemacs.${MPMNEURON_SYSTEM_EXT}
       ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libmeschach.${MPMNEURON_SYSTEM_EXT}
       ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libivos.${MPMNEURON_SYSTEM_EXT}
       ${MPMNEURON}/${MPMNEURON_ARCH}/lib/libocxt.${MPMNEURON_SYSTEM_EXT})
   
   # set the includes for NEURON
   set(NEURON_INCLUDE_PATH
         ${MPMNEURON} 
         ${MPMNEURON}/include
         ${MPMNEURON}/src
         ${MPMNEURON}/src/ivoc
         ${MPMNEURON}/include/nrn 
         ${MPMNEURON}/src/ivos
         ${MPMNEURON}/src/oc)
         
   # NEURON required linking to readline (this fixes the issue, that SHELLTYPE=readline was required previously!)
   include(FindReadline.cmake)

   # use additional includes for NEURON
   include_directories(${NEURON_INCLUDE_PATH})

   # definition for usage in source code, if we have no NEURON API then we just use the membrane potential mapping plugin without NEURON API dependent features for convenience
   add_definitions(-DMPMNEURON)
   
   # status report - NEURON arch and OS
   message(STATUS "Detected architecture ${MPMNEURON_ARCH} on ${CMAKE_SYSTEM_NAME} for the NEURON API.")

   # output information on success or failure
   message(STATUS "NEURON API is on (by -DMPMNEURON:${MPMNEURON}), you can use the full hoc language of NEURON (as provided by ${MPMNEURON_REVISION_NAME}).")
   
   # check for readline support
   if (READLINE_FOUND)
       message(STATUS "Readline library was found (${READLINE_LIBRARY})")
   else(READLINE_FOUND)
       message(STATUS "Readline library not found - you need SHELLTYPE readline now.")
   endif(READLINE_FOUND)
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
                     a_u_x/a_u_x_bridge.cpp)

set(SOURCES_default  
                     vm2ug.cpp 
                     plugin_main.cpp
                     bg_default/spannung.cpp
                     bg_default/my_commands.cpp
                     bg_default/bg.cpp
                     bg_default/solve.cpp
                     bg_default/gating.cpp
                     a_u_x/a_u_x_bridge.cpp)
                     
else(NOT MPMNEURON)
   # simple sources for VGCCs
   set(SOURCES_simple
			         vm2ug.cpp
			         plugin_main.cpp
			         bg_simple/bg.cpp
			         bg_simple/my_commands.cpp
			         bg_simple/solve.cpp
			         bg_simple/gating.cpp
                     bg_simple/spannung.cpp
                     a_u_x/a_u_x_bridge.cpp
                     transformator.cpp)

   # default sources for VGCCs 
   set(SOURCES_default  
               vm2ug.cpp 
               plugin_main.cpp
               bg_default/spannung.cpp
               bg_default/my_commands.cpp
               bg_default/bg.cpp
               bg_default/solve.cpp
               bg_default/gating.cpp
               a_u_x/a_u_x_bridge.cpp
               transformator.cpp)
endif(NOT MPMNEURON)

# set sources for unit tests
set(SOURCES_TEST unit_tests/src/tests.cpp)
		
# get the path to annlib and set annlib include path
option(ANNLIB "Sets the path in which annlib shall be searched." OFF)

if(NOT ANNLIB)
	message(FATAL_ERROR "ERROR: libann is required. Please set the option -DANNLIB=path_to_libann.")
endif(NOT ANNLIB)

# add include directoris for ANNLIB
find_library(ANNLIB_LIBS NAMES ANN PATHS ${ANNLIB}/lib)
include_directories(${ANNLIB}/include)

if(${ANNLIB_LIBS} STREQUAL "ANNLIB_LIBS-NOTFOUND")
	message(FATAL_ERROR "ERROR: Couldn't find ANNLIB in the sporecified path: " ${ANNLIB}/lib)
else(${ANNLIB_LIBS} STREQUAL "ANNLIB_LIBS-NOTFOUND")
	include_directories(${ANNLIB}/include)		
	message(STATUS "ANN Library API is on (by -ANNLIB=${ANNLIB}).")
	message(STATUS "")
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
  # detect if NEURON API should be used
   if (MPMNEURON)
      ExportIncludes(${MPMNEURON} ${MPMNEURON}/include ${MPMNEURON}/src ${MPMNEURON}/src/ivoc ${MPMNEURON}/include/nrn ${MPMNEURON}/src/ivos ${MPMNEURON}/src/oc)

      ExportDependencies(${NEURON_LIBRARY_PATH})
      ExportDependencies(${READLINE_LIBRARY})

      ExportDefinitions(-DMPMNEURON)

      # if we need head revision set it (this needs to be fixed in mpm plugin code!)
      if(${MPMNEURON_REVISION} STREQUAL "ON")
          ExportDefinitions(-DMPMNEURON_REVISION)
      endif(${MPMNEURON_REVISION} STREQUAL "ON")

   endif(MPMNEURON)
 
   # exports for vm2ug (libann) and borg graham style
   ExportDependencies(${ANNLIB}/lib/libANN.a)
   ExportIncludes(${ANNLIB}/include)
   
# decide which build of the plugin you want to build
    if(${MPMVGCC} STREQUAL "ON")
	    ExportSources(plugins/experimental/membrane_potential_mapping ${SOURCES_default})
      ExportDefinitions(-DMPMVGCC)
	else(${MPMVGCC} STREQUAL "ON")
	    ExportSources(plugins/experimental/membrane_potential_mapping ${SOURCES_simple})
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

   if (${SHELLTYPE} STREQUAL "readline")
	   target_link_libraries(MembranePotentialMapping ug4 ${ANNLIB_LIBS} ${NEURON_LIBRARY_PATH})
   else(${SHELLTYPE} STREQUAL "readline")
	   target_link_libraries(MembranePotentialMapping ug4 ${ANNLIB_LIBS} ${NEURON_LIBRARY_PATH} ${READLINE_LIBRARY})
   endif(${SHELLTYPE} STREQUAL "readline")
	
	if(${MPMTestsuite} STREQUAL "ON")
		target_link_libraries (MPMTestsuite MembranePotentialMapping ug4)
	endif(${MPMTestsuite} STREQUAL "ON")
	
endif(buildEmbeddedPlugins)
