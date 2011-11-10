################################################################################
# Your description here...
################################################################################

cmake_minimum_required(VERSION 2.6)

project(MEMBRANE_POTENTIAL_MAPPING)

# include the definitions and dependencies for ug-plugins.
include("../../../ug_plugin_cmake_includes.txt")
include_directories("annlib/include")

# sources
set(SOURCES	plugin_main.cpp
			#vm2ug.cc
			bg.cc
			my_commands.cc
			solve.cc
			gating.cc)

set(ANN_SOURCES	annlib/src/ANN.cpp
				annlib/src/brute.cpp
				annlib/src/kd_tree.cpp
				annlib/src/kd_util.cpp
				annlib/src/kd_split.cpp
				annlib/src/kd_dump.cpp
				annlib/src/kd_search.cpp
				annlib/src/kd_pr_search.cpp
				annlib/src/kd_fix_rad_search.cpp
				annlib/src/bd_tree.cpp
				annlib/src/bd_search.cpp
				annlib/src/bd_pr_search.cpp
				annlib/src/bd_fix_rad_search.cpp
				annlib/src/perf.cpp)

# create a shared library from the sources and link it against ug4.
add_library(MembranePotentialMappingPlugin SHARED ${SOURCES} ${ANN_SOURCES})

if(STATIC)
	target_link_libraries (MembranePotentialMappingPlugin ug4_s)
else(STATIC)
	target_link_libraries (MembranePotentialMappingPlugin ug4)
endif(STATIC)
