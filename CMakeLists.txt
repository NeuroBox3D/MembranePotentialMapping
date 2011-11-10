################################################################################
# Your description here...
################################################################################

cmake_minimum_required(VERSION 2.6)

project(MEMBRANE_POTENTIAL_MAPPING)

# include the definitions and dependencies for ug-plugins.
include("../../../ug_plugin_cmake_includes.txt")


# sources
set(SOURCES	plugin_main.cpp)


# create a shared library from the sources and link it against ug4.
add_library(MembranePotentialMappingPlugin SHARED ${SOURCES})

if(STATIC)
	target_link_libraries (MembranePotentialMappingPlugin ug4_s)
else(STATIC)
	target_link_libraries (MembranePotentialMappingPlugin ug4)
endif(STATIC)
