##
 # Copyright (c) 2009-2020: G-CSC, Goethe University Frankfurt
 #
 # Author: Stephan Grein
 # Creation date: 2012-02-02
 #
 # This file is part of NeuroBox, which is based on UG4.
 #
 # NeuroBox and UG4 are free software: You can redistribute it a #/or modify it
 # under the terms of the GNU Lesser General Public License version 3
 # (as published by the Free Software Foundation) with the following additional
 # attribution requirements (according to LGPL/GPL v3 §7):
 #
 # (1) The following notice must be displayed in the appropriate legal notices
 # of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 #
 # (2) The following notice must be displayed at a prominent place in the
 # terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 #
 # (3) The following bibliography is recommended for citation and must be
 # preserved in all covered files:
 # "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 #   parallel geometric multigrid solver on hierarchically distributed grids.
 #   Computing and visualization in science 16, 4 (2013), 151-164"
 # "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 #   flexible software system for simulating PDE based models on high performance
 #   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 # "Stepniewski, M., Breit, M., Hoffer, M. and Queisser, G.
 #   NeuroBox: computational mathematics in multiscale neuroscience.
 #   Computing and visualization in science (2019).
 # "Breit, M. et al. Anatomically detailed and large-scale simulations studying
 #   synapse loss and synchrony using NeuroBox. Front. Neuroanat. 10 (2016), 8"
 #
 # This program is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 # GNU Lesser General Public License for more details.
 ##
 
################################################################################
# Import notes:
# Minimum required NEURON version is v7.3. The package nrn is required from the
# following website: http://www.neuron.yale.edu/ftp/neuron/versions/v7.3/
# This will make MembranePotentialMapping work with OSX, Linux and Windows
################################################################################

# minimum required CMake version for our project below
cmake_minimum_required(VERSION 2.6)

# specify project name
project(MEMBRANE_POTENTIAL_MAPPING)


# show optional arguments for the plugin
message(STATUS "")
message(STATUS "Additional options for the MembranePotentialMapping (-DMPM) plugin:")

# enable the boost testsuite
option(MPMTestsuite "Build Testsuite" OFF)
message(STATUS "      Testsuite:       " ${MPMTestsuite} " (options are: ON, OFF)")

# which channel type
option(MPMVGCC "Build VGCC:    " 0)
message(STATUS "      VGCC:            " ${MPMVGCC} " (options are: 0 (none), 1 (simple), 2 (default))")

# build with openmp 
option(MPMOpenMP "Build OpenMP:  " OFF)
message(STATUS "      OpenMP:          " ${MPMOpenMP} " (options are: ON, OFF)")

# use head NEURON revision or pre 7.3
option(MPMNEURON_REVISION "Use head revision of NEURON or pre 7.3 build? " ON)
message(STATUS "      NEURON_REVISION: " ${MPMNEURON_REVISION} " (options are: ON, OFF)")

option(MPMNEURON_STATIC "Use static NEURON build?" ON)
message(STATUS "      NEURON_STATIC:   " ${MPMNEURON_STATIC} " (options are: ON, OFF)")

option(SHELLTYPE "The type of shell that should be used: (options are none, readline)" "none")
message(STATUS "      SHELLTYPE:       " ${SHELLTYPE} " (options are none, readline)")


# check for the presence of the NEURON API
if (NOT MPMNEURON)
  message(STATUS "      NEURON:          OFF (turn Neuron API on using -DMPMNEURON=/path/to/NEURON/installation")
else (NOT MPMNEURON)
  message(STATUS "      NEURON:      " $(MPMNEURON))

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
   
  # For now static libraries for NEURON are used (Remember to compile NEURON fully static!)
  # If ug4 build is STATIC as well, NEURON .modfiles will not be loaded, thus we need DYNAMIC build for ug4 instead.
  if (MPMNEURON_STATIC)
    set(MPMNEURON_SYSTEM_EXT "a")
      # My static NEURON build with head revision (Without memacs (terminal cursosr) 
      # support (otherwise ncurses will ne required at compile/link time and this 
      # may not be available on all platforms we want to support)
      # set(MPMNEURON "/Users/stephan/Code/hg/neuron/nrn/")

      # My dynamic build of NEURON with pre v7.3 version (memacs and ncurses must be installed and present at runtime).
      # Note, that this is okay for non-VRL builds in general but not for VRL builds: Issues could possibly happen.
      # set(MPMNEURON "/Users/stephan/Code/git/neuron/nrn-7.3/")
  
      # Note: One could check for static build, then remove from NEURON_LIBRARY_PATH memacs, 
      #       but this is not necessary since a non-avaialble lib or file will be ignored
      #       by CMake and throw at most a warning: Anyway this should be fixed in the future.
  endif(MPMNEURON_STATIC)


  # Libraries are set manually for now, should provide a CMakeFind command for NEURON installations
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

  message(STATUS "Note: NEURON v7.3 is recommended!")
  message(STATUS "Note: SHELLTYPE=readline is recommended!")
         
  # NEURON required linking to readline (this fixes the issue, that SHELLTYPE=readline was required previously!)
  include(FindReadline.cmake)

  # use additional includes for NEURON
  include_directories(${NEURON_INCLUDE_PATH})
   
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
include(${UG_ROOT_CMAKE_PATH}/ug_plugin_includes.cmake)

# set the sources
set(SOURCES
  kdtree/kd_tree.cpp
  plugin_main.cpp
  #a_u_x/a_u_x_bridge.cpp
)

if("${MPMVGCC}" STREQUAL "1")
  set (SOURCES ${SOURCES}
    bg_simple/bg.cpp
    bg_simple/my_commands.cpp
    bg_simple/solve.cpp
    bg_simple/gating.cpp
    bg_simple/voltage.cpp
  )
elseif("${MPMVGCC}" STREQUAL "2")
  set (SOURCES ${SOURCES}
    bg_default/voltage.cpp
    bg_default/my_commands.cpp
    bg_default/bg.cpp
    bg_default/solve.cpp
    bg_default/gating.cpp
  )
endif("${MPMVGCC}" STREQUAL "1")

# add NEURON-specific sources
if (MPMNEURON)
  set(SOURCES ${SOURCES} transformator.cpp)
endif (MPMNEURON)

# set sources for unit tests
set(SOURCES_TEST unit_tests/src/tests.cpp)

################################################################################

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



# This will copy the configurable file config_cmake.h to
# the final configuration header neti_config.h,
# setting the #defines as given by the current CMake variables.
configure_file(config_cmake.h ${CMAKE_BINARY_DIR}/mpm_config.h)


# create a shared library from the sources and link it against ug4.
if(buildEmbeddedPlugins)
  # detect if NEURON API should be used
  if (MPMNEURON)
    ExportIncludes(${MPMNEURON} ${MPMNEURON}/include ${MPMNEURON}/src ${MPMNEURON}/src/ivoc ${MPMNEURON}/include/nrn ${MPMNEURON}/src/ivos ${MPMNEURON}/src/oc)

    ExportDependencies(${NEURON_LIBRARY_PATH})
    ExportDependencies(${READLINE_LIBRARY})

    # if we need head revision set it (this needs to be fixed in mpm plugin code!)
  endif(MPMNEURON)
    
  ExportSources(${CMAKE_CURRENT_SOURCE_DIR} ${SOURCES})
  
  if(${MPMTestsuite} STREQUAL "ON")
    target_link_libraries (MPMTestsuite ug4 ${NEURON_LIBRARY_PATH}) 
  endif(${MPMTestsuite} STREQUAL "ON")

else(buildEmbeddedPlugins)
  add_library(MembranePotentialMapping SHARED ${SOURCES})
 
  if(MPMNEURON)
   if(${SHELLTYPE} STREQUAL "readline")
     target_link_libraries(MembranePotentialMapping ug4 ${linkLibraries} ${NEURON_LIBRARY_PATH} ${READLINE_LIBRARY})
   else(${SHELLTYPE} STREQUAL "readline")
     target_link_libraries(MembranePotentialMapping ug4 ${linkLibraries} ${NEURON_LIBRARY_PATH})
   endif(${SHELLTYPE} STREQUAL "readline")
  else(MPMNEURON)
     target_link_libraries(MembranePotentialMapping ug4 ${linkLibraries})
  endif(MPMNEURON)
  
  if(${MPMTestsuite} STREQUAL "ON")
    target_link_libraries (MPMTestsuite MembranePotentialMapping ug4 ${linkLibraries})
  endif(${MPMTestsuite} STREQUAL "ON")
  
endif(buildEmbeddedPlugins)
