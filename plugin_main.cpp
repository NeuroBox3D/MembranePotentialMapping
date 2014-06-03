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
#ifdef MPMNEURON
#include "transformator.h"
#include "hoc_command.h"
#endif

#include "a_u_x/edge_utilities.h"
#include "a_u_x/a_u_x_bridge.cpp"

#include <lib_disc/function_spaces/grid_function.h>
#include <bridge/util_domain_algebra_dependent.h>
#include "lib_disc/domain.h"
#include "lib_disc/common/function_group.h"

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
			template <typename TDomain, typename TAlgebra> static void
			DomainAlgebra(ug::bridge::Registry& reg, const std::string& parentGroup)
			{
			   // group membership of the membrane_potential_mapping plugin
			   std::string grp(parentGroup);
			   grp.append("Neuro/");

				// define algebra
				typedef typename TAlgebra::vector_type vector_type;
				typedef typename TAlgebra::matrix_type matrix_type;

				// define the grid function depending on domain and algebra
				typedef GridFunction<TDomain, TAlgebra> TGridFunction;

			   // typedefs
			   typedef membrane_potential_mapping::Vm2uG<std::string> TVm2uG;
			   typedef membrane_potential_mapping::bg::BG TBG;
			   typedef membrane_potential_mapping::Transformator TTransformator;
			   typedef membrane_potential_mapping::HocCommand THC;

				// do we want to use smartpointers?
				const bool bSmartPointer = true;

			   // registry of auxiliarities (\see aux/edge_utilities.h)
			   reg.add_function("get_edge_sum", static_cast<number (*)(TDomain&,  ISubsetHandler&, int, int)> (&membrane_potential_mapping::aux::EdgeSum<TDomain>), "", grp);
			   reg.add_function("get_edge_sum_sq", static_cast<number (*)(TDomain&, ISubsetHandler&, int, int)> (&membrane_potential_mapping::aux::EdgeSumSq<TDomain>), "", grp);

			   // registry of Vm2uG (\see vm2ug.h)
			   reg.add_class_<TVm2uG>("MembranePotentialMapper", grp)
#ifdef MPMNEURON
				   .add_constructor<void (*)(SmartPtr<Transformator>)>("Transformation setup")
#endif
				  .add_constructor<void (*)(const std::string&, const std::string&, bool)>("Initial timestep|load-dialog|endings=[\"csv\", \"txt\"];description=\"Timestep files\"#Suffix|selection|value=[\"csv\", \"txt\"]#Static Nodes|selection|value=[True, False]")
			#ifdef MPMNEURON
				  .add_method("build_tree", (void (TVm2uG::*)()) (&TVm2uG::buildTree), grp)
			#endif
#ifndef MPMNEURON
				  .add_method("get_potential", (number (TVm2uG::*)(number, number, number, const std::string&)) (&TVm2uG::get_potential), "Potential|default", "x|default#y|default#z|default#Timestep|default", grp)
#else
				  .add_method("get_potential", (number (TVm2uG::*)(number, number, number)) (&TVm2uG::get_potential), grp)
#endif
#ifndef MPMNEURON
				  .add_method("build_tree", &TVm2uG::buildTree, grp)
				  .add_method("get_potential_lin", &TVm2uG::get_potential_lin, "Potential|default", "x|default#y|default#z|default#Timestep|default#k|default")
				  .add_method("get_potential_bilin", &TVm2uG::get_potential_bilin, "Potential|default", "x|default#y|default#z|default#Timestep|default#k|default")
#endif
				  .set_construct_as_smart_pointer(bSmartPointer);


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
				  .set_construct_as_smart_pointer(bSmartPointer);

#ifdef MPMNEURON
					/// registry of Transformator
				reg.add_class_<TTransformator>("Transformator", grp)
					.add_constructor<void (*)()>("")
					.add_method("load_geom", (void (TTransformator::*)(const std::string&))(&TTransformator::load_geom), "", "hoc geometry", grp)
					.add_method("load_stim", (void (TTransformator::*)(const std::string&))(&TTransformator::load_stim), "", "hoc stimulation protocol", grp)
					.add_method("get_sections", (size_t (TTransformator::*)())(&TTransformator::get_sections), "number of sections", "", grp)
					.add_method("get_total_points", (size_t (TTransformator::*)())(&TTransformator::get_total_points), "number of total points", "", grp)
					.add_method("get_t", (number (TTransformator::*)())(&TTransformator::get_t), "t", "", grp)
					.add_method("get_dt", (number (TTransformator::*)())(&TTransformator::get_dt), "dt", "", grp)
					.add_method("get_tstop", (number (TTransformator::*)())(&TTransformator::get_tstop), "tstop", "", grp)
					.add_method("get_tstart", (number (TTransformator::*)())(&TTransformator::get_tstart), "tstart", "", grp)
					.add_method("get_finitialize", (number (TTransformator::*)())(&TTransformator::get_finitialize), "finitialize", "", grp)
					.add_method("setup_hoc", (void (TTransformator::*)(number, number, number, number))(&TTransformator::setup_hoc), "", "tstart|tstop|dt|finitialize", grp)
					.add_method("extract_vms", (void (TTransformator::*)(size_t, size_t))(&TTransformator::extract_vms), "", "limit", grp)
					.add_method("purge", (void (TTransformator::*)())(&TTransformator::purge), "", "reinit", grp)
					.add_method("print_setup", (void (TTransformator::*)(bool))(&TTransformator::print_setup), "", "prints setup", grp)
	    			.add_method("adjust_resistivity_obstacle", (void (TTransformator::*)(const char* , const char* , double ))(&TTransformator::adjust_resistivity_obstacle), "", "", grp)
	    			.add_method("feedback", (void (TTransformator::*)(const char* uCmp, const char* , const char* ,double ))(&TTransformator::feedback), "", "", grp)
	    			.add_method("feedbacks", (void (TTransformator::*)(const char* uCmp, const char* , const char* ,double ))(&TTransformator::feedbacks), "", "", grp)
	    			.add_method("set_feedback", (void (TTransformator::*)())(&TTransformator::set_feedback), "", "", grp)
	    			.add_method("set_feedbacks", (void (TTransformator::*)())(&TTransformator::set_feedback), "", "", grp)
	  			    .add_method("set_grid", (void (TTransformator::*)(TGridFunction&, const char*))(&TTransformator::set_grid<TGridFunction, TDomain>), "", "", grp)
	    			.add_method("init_feedback", (void (TTransformator::*)(const char* uCmp, const char* , const char* ,double, const char*))(&TTransformator::init_feedback), "", "", grp)
	    			.add_method("init_feedbacks", (void (TTransformator::*)(const char* uCmp, const char* , const char* ,double, const char* ))(&TTransformator::init_feedbacks), "", "", grp)
					.set_construct_as_smart_pointer(bSmartPointer);

				reg.add_class_<THC>("HocCmd", grp)
					.add_constructor<void (*)()>("");
#endif
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
				bridge::RegisterDomain3dAlgebraDependent<Functionality>(reg, parentGroup);
			} UG_REGISTRY_CATCH_THROW(parentGroup);
		#else
			UG_WARNING("Reasonable usage of MPM plugin assured with DIM=3.")
		#endif
	}
// end namespace ug
}
