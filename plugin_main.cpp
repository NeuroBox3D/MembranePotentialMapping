/*!
 * \file plugins/experimental/membrane_potential_mapping/plugin_main.cpp
 * \brief registry of the ug plugin membrane_potential_mapping
 *
 * \date created on Apr 30, 2012
 * \author stephan grein
 */

#include <string>

#include <bridge/util.h>
#include <bridge/util_domain_dependent.h>
#include <common/error.h>
#include <common/ug_config.h>
#include <common/log.h>
#include <registry/registry.h>
#include <registry/error.h>
#include <lib_disc/domain.h>
#include <lib_grid/lib_grid.h>

#include "vm2ug.h"
#ifdef MPMVGCC
#include "bg_default/bg.h"
#else
#include "bg_simple/bg.h"
#endif
#include "transform.h"

#include "aux/edge_utilities.h"
#include "aux/aux_bridge.cpp"

// begin namespace ug
namespace ug {
	// begin namespace mpm
	namespace membrane_potential_mapping {
		/*!
		 * \defgroup mpm_plugin Membrane Potential Mapping plugin
		 * \ingroup plugins_experimental
		 * \{
		 *
		 */
		// the functionality which is to be registered
		struct Functionality {
			template <typename TDomain> static void
			Domain(ug::bridge::Registry& reg, const std::string& parentGroup)
			{
			   // group membership of the membrane_potential_mapping plugin
			   std::string grp(parentGroup);
			   grp.append("Neuro/");

			   // typedefs
			   typedef membrane_potential_mapping::Vm2uG<std::string> TVm2uG;
			   typedef membrane_potential_mapping::bg::BG TBG;
			   typedef membrane_potential_mapping::Transform TTransform;

			   // registry of auxiliarities (\see aux/edge_utilities.h)
			   reg.add_function("get_edge_sum", static_cast<number (*)(TDomain&,  ISubsetHandler&, int, int)> (&membrane_potential_mapping::aux::EdgeSum<TDomain>), "", grp);
			   reg.add_function("get_edge_sum_sq", static_cast<number (*)(TDomain&, ISubsetHandler&, int, int)> (&membrane_potential_mapping::aux::EdgeSumSq<TDomain>), "", grp);

			   // registry of Vm2uG (\see vm2ug.h)
			   reg.add_class_<TVm2uG>("MembranePotentialMapper", grp)
				  .add_constructor<void (*)(const std::string&, const std::string&, bool)>("Initial timestep|load-dialog|endings=[\"csv\", \"txt\"];description=\"Timestep files\"#Suffix|selection|value=[\"csv\", \"txt\"]#Static Nodes|selection|value=[True, False]")
				  .add_method("get_potential", &TVm2uG::get_potential, "Potential|default", "x|default#y|default#z|default#Timestep|default")
				  .add_method("build_tree", &TVm2uG::buildTree, grp)
				  .add_method("get_potential_lin", &TVm2uG::get_potential_lin, "Potential|default", "x|default#y|default#z|default#Timestep|default#k|default")
				  .add_method("get_potential_bilin", &TVm2uG::get_potential_bilin, "Potential|default", "x|default#y|default#z|default#Timestep|default#k|default");

			   // registry of BG (\see bg.h)
			   reg.add_class_<TBG>("BorgGraham", grp)
				  .add_constructor()
				  .add_method("install_can_gates", &TBG::install_can_gates, grp)
				  #ifdef MPMVGCC
				  .add_method("install_can_gates_cfp", &TBG::install_can_gates_cfp, grp)
				  .add_method("get_current", (double (TBG::*)(const double, const double, const double)) (&TBG::timestepping_of_gates_and_calc_current), "t [s] |default#delta t [s]|default#custom membrane potential [mV]|default",  grp)
				  .add_method("get_current", (double (TBG::*)(const double, const double, const double, const double, const double)) (&TBG::timestepping_of_gates_and_calc_current), "t [s] |default#delta t [s]|default#custom membrane potential [mV]|default#IC Calcium [Mol]|default#EC Calcium [Mol]|default",  grp)
				  #else
				  .add_method("get_current", (double (TBG::*)(const double, const double)) (&TBG::timestepping_of_gates_and_calc_current), "t [s] |default#delta t [s]|default",  grp)
				  #endif
				  .add_method("calc_current_at_start", (double (TBG::*)(const double)) (&TBG::calc_current_at_start), grp)
				  .add_method("calc_current_at_start", (double (TBG::*)(const double, const double, const double, const double, const double)) (&TBG::calc_current_at_start), grp)
				  .add_method("get_Neumann_Flux", &TBG::get_Neumann_Flux, grp)
				  .add_method("get_Neumann_Flux_As_Concentration", &TBG::get_Neumann_Flux_as_Concentration, grp)
				  #ifdef MPMVGCC
				  .add_method("dCadCa_i", &TBG::dCa_dCa_i, "derivative w.r.t. internal calcium concentration|default", "", grp)
				  .add_method("dCadCa_o", &TBG::dCa_dCa_i, "derivative w.r.t. external calcium concentration|default", "", grp)
				  #else
				  .add_method("dCa", &TBG::dCa, "derivative if no calcium concentration is considered|default", "", grp)
				  #endif
				  ;

			   // registry of Transform (\see transform.h)
				reg.add_class_<TTransform>("TransformHocToObj", grp)
				   .add_constructor<void (*)(const std::string&, const std::string&)>("hocfile|load-dialog|endings=[\"hoc\"];description=\"Hoc file for transformation\"#timestep directory|load-dialog;description=\"Location to store extracted timesteps\"#Delta t|default#steps|default#initial membrane potential|default")
				   .add_method("get_hocfile", &TTransform::get_hocfile, "hocfile|default", "", grp)
				   .add_method("get_objfile", &TTransform::get_objfile, "objfile|default", "", grp)
				   .add_method("get_xmlfile", &TTransform::get_xmlfile, "xmlfile|default", "", grp)
				   .add_method("get_timestepdirectory", &TTransform::get_timestepdirectory, "timestep directory|default", "", grp)
				   .add_method("get_dt", &TTransform::get_dt, "Delta t", "", grp)
				   .add_method("extract_timesteps_and_obj", (void (TTransform::*)(bool, const std::string&, const std::string&)) &TTransform::extract_timesteps_and_obj, "", "generate object|default#neugen exe|load-dialog|endings=[\"jar\"];description=\"Neugen executable\"#neutria exe|load-dialog;description=\"Neutria executable\"", grp);
			}
		// end of functionality which is to be exported
		};
	/// \}
	// end namespace mpm
	}

	/// \addtogroup mpm_plugin
	extern "C" void
	InitUGPlugin_MembranePotentialMapping(bridge::Registry& reg, std::string& parentGroup)
	{
		parentGroup.append("MembranePotentialMapping");
		typedef membrane_potential_mapping::Functionality Functionality;

		#if defined(UG_DIM_3)
			try {
				bridge::RegisterDomain3dDependent<Functionality>(reg, parentGroup);
			} UG_REGISTRY_CATCH_THROW(parentGroup);
		#else
			UG_WARNING("Reasonable usage of MPM plugin assured with DIM=3.")
		#endif
	}
// end namespace ug
}
