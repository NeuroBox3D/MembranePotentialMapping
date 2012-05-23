/*
 * plugin_main.cpp
 *
 *  Created on: Apr 30, 2012
 *      Author: stephan grein
 */

#include "registry/registry.h"
#include "common/ug_config.h"
#include "common/error.h"
#include "lib_disc/domain.h"
#include "lib_grid/lib_grid.h"

#include <string>

#include "vm2ug.h"
#ifdef MPMDEFAULT
#include "bg_default/bg.h"
#else
#include "bg_simple/bg.h"
#endif

using namespace ug;

extern "C" UG_API void
InitUGPlugin_MembranePotentialMappingPlugin(ug::bridge::Registry* reg, string parentGroup) {
   // group to which the plugin shall belong
   std::string grp(parentGroup); grp.append("Neuro/");

   // register functionalities of the Vm2uG tool
   typedef vug::Vm2uG<std::string> TVm2uG;

   reg->add_class_<TVm2uG>("MembranePotentialMapper", grp)
       .add_constructor<void (*)(std::string, std::string, bool)>("Timestep-File#File-Extension#NodesCanChange")
       .add_method("get_potential", &TVm2uG::get_potential, grp)
       .add_method("build_tree", &TVm2uG::buildTree, grp)
       .add_method("get_potential_lin", &TVm2uG::get_potential_lin, grp)
       .add_method("get_potential_bilin", &TVm2uG::get_potential_bilin, grp);

   // register functionalities of the Borg Graham model
   typedef bg::BG TBG;

   reg->add_class_<TBG>("BorgGraham", grp)
       .add_constructor()
       .add_method("install_can_gates", &TBG::install_can_gates, grp)
    //   .add_method("get_current", &TBG::timestepping_of_gates_and_calc_current, grp)
       .add_method("calc_current_at_start", &TBG::calc_current_at_start, grp)
       .add_method("get_Neumann_Flux", &TBG::get_Neumann_Flux, grp);
}

