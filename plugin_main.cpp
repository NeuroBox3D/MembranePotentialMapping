/**
 * \file plugin_main.cpp
 * \brief registry of the ug plugin membrane_potential_mapping
 *
 *  \date created on Apr 30, 2012
 *  \author stephan grein
 */

// includes
#include <string>

#include "vm2ug.h"
#ifdef MPMDEFAULT
#include "bg_default/bg.h"
#else
#include "bg_simple/bg.h"
#endif
#include "transform.h"

#include <common/log.h>
#include "registry/registry.h"
#include "registry/error.h"
#include "common/ug_config.h"
#include "common/error.h"
#include "lib_disc/domain.h"
#include "lib_grid/lib_grid.h"


// using directives
using namespace ug;


extern "C" UG_API void
InitUGPlugin_MembranePotentialMapping(ug::bridge::Registry* reg, std::string parentGroup) {
	try {
	   // group membership of the membrane_potential_mapping plugin
	   std::string grp(parentGroup); grp.append("Neuro/");

	   // typedefs
	   typedef membrane_potential_mapping::Vm2uG<std::string> TVm2uG;
	   typedef membrane_potential_mapping::bg::BG TBG;
	   typedef membrane_potential_mapping::Transform TTransform;

	   /** registry Vm2uG (\see vm2ug) */
	   reg->add_class_<TVm2uG>("MembranePotentialMapper", grp)
		   .add_constructor<void (*)(std::string, std::string, bool)>("Initial timestep|load-dialog|endings=[\"csv\", \"txt\"];description=\"Timestep files\"#Suffix|selection|value=[\"csv\", \"txt\"]#Static Nodes|selection|value=[True, False]")
		   .add_method("get_potential", &TVm2uG::get_potential, "Potential|default", "x|default#y|default#z|default#Timestep|default")
		   .add_method("build_tree", &TVm2uG::buildTree, grp)
		   .add_method("get_potential_lin", &TVm2uG::get_potential_lin, "Potential|default", "x|default#y|default#z|default#Timestep|default#k|default")
		   .add_method("get_potential_bililn", &TVm2uG::get_potential_bilin, "Potential|default", "x|default#y|default#z|default#Timestep|default#k|default");

	   /** registry BG (\see bg.h), but soon obsolete. TODO: remove BG (see meeting.pdf) */
	   	   reg->add_class_<TBG>("BorgGraham", grp)
		   .add_constructor()
		   .add_method("install_can_gates", &TBG::install_can_gates, grp)
			#ifdef MPMDEFAULT
		   .add_method("get_current", (double (TBG::*)(const double, const double, const double)) (&TBG::timestepping_of_gates_and_calc_current), "t [s] |default#delta t [s]|default#custom membrane potential [mV]|default",  grp)
			#else
		   .add_method("get_current", (double (TBG::*)(const double, const double)) (&TBG::timestepping_of_gates_and_calc_current), "t [s] |default#delta t [s]|default",  grp)
			#endif
		   .add_method("calc_current_at_start", &TBG::calc_current_at_start, grp)
		   .add_method("get_Neumann_Flux", &TBG::get_Neumann_Flux, grp);


	   /** registry Transform (\see Transform.h) */
	   	reg->add_class_<TTransform>("TransformHocToObj", grp)
	   		.add_constructor<void (*)(std::string, std::string)>("hocfile|load-dialog|endings=[\"hoc\"];description=\"Hoc file for transformation\"#timestep directory|load-dialog;description=\"Location to store extracted timesteps\"#Delta t|default#steps|default#initial membrane potential|default")
	   		.add_method("get_hocfile", &TTransform::get_hocfile, "hocfile|default", "", grp)
	   		.add_method("get_objfile", &TTransform::get_objfile, "objfile|default", "", grp)
	   	    .add_method("get_xmlfile", &TTransform::get_xmlfile, "xmlfile|default", "", grp)
	   	    .add_method("get_timestepdirectory", &TTransform::get_timestepdirectory, "timestep directory|default", "", grp)
	   	    .add_method("get_dt", &TTransform::get_dt, "Delta t", "", grp)
	   	    .add_method("extract_timesteps_and_obj", (void (TTransform::*)(bool, std::string, std::string)) &TTransform::extract_timesteps_and_obj, "", "generate object|default#neugen exe|load-dialog|endings=[\"jar\"];description=\"Neugen executable\"#neutria exe|load-dialog;description=\"Neutria executable\"", grp);

		} catch (const ug::bridge::UGRegistryError& error) {
			UG_LOG("### ERROR in UGRegistry: InitUGPlugin_MembranePotentialMappingPlugin failed registering with message: " << error.name << "" << std::endl);
		}
}
