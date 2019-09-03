/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Stephan Grein
 * Creation date: 2012-04-30
 *
 * This file is part of NeuroBox, which is based on UG4.
 *
 * NeuroBox and UG4 are free software: You can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3
 * (as published by the Free Software Foundation) with the following additional
 * attribution requirements (according to LGPL/GPL v3 §7):
 *
 * (1) The following notice must be displayed in the appropriate legal notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 *
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 *
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating PDE based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * "Stepniewski, M., Breit, M., Hoffer, M. and Queisser, G.
 *   NeuroBox: computational mathematics in multiscale neuroscience.
 *   Computing and visualization in science (2019).
 * "Breit, M. et al. Anatomically detailed and large-scale simulations studying
 *   synapse loss and synchrony using NeuroBox. Front. Neuroanat. 10 (2016), 8"
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
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

#include "mpm_config.h"  // for project-specific defines

//#include "vm2ug.h"
#include "vm2ug.h"
#if (defined(MPMVGCC) && MPVGCC == 1)
#include "bg_simple/bg.h"
#endif
#if (defined(MPMVGCC) && MPVGCC == 2)
#include "bg_default/bg.h"
#endif
#ifdef MPMNEURON
#include "transformator.h"
#include "neuron_mpm.h"
#endif

#include "vm2ug_mpm.h"

#include "kdtree/kd_tree.h"

//#include "a_u_x/edge_utilities.h"
//#include "a_u_x/a_u_x_bridge.cpp"

#include <lib_disc/function_spaces/grid_function.h>
#include <bridge/util_domain_algebra_dependent.h>
#include "lib_disc/domain.h"
#include "lib_disc/common/function_group.h"
#include <bridge/util.h>
#include <bridge/util_domain_dependent.h>
#include <registry/registry.h>

#include <bridge/util.h>
#include <bridge/util_domain_dependent.h>
#include <common/error.h>
#include <common/ug_config.h>
#include <common/log.h>
#include <registry/registry.h>
#include <registry/error.h>


using namespace ug::bridge;

namespace ug {
namespace membrane_potential_mapping {

/*!
 * \defgroup mpm_plugin Membrane Potential Mapping plugin
 * \ingroup plugins_experimental
 * \{
 *
 */
// the functionality which is to be registered
struct Functionality {
	/**
	 * Function called for the registration of Domain and Algebra independent parts.
	 * All Functions and Classes not depending on Domain and Algebra
	 * are to be placed here when registering.
	 *
	 * @param reg   registry
	 * @param grp   group for sorting of functionality
	 */
	static void Common(Registry& reg, std::string grp)
	{
		// Mapper
		{
			typedef membrane_potential_mapping::Mapper<3, number> TMapper;
			std::string name("Mapper");
			reg.add_class_<TMapper>(name, grp)
				.add_constructor<void (*)()>("", "", "")
				.add_method("build_tree", static_cast<void (TMapper::*)()> (&TMapper::build_tree), grp)
				.add_method("build_tree_from_file", static_cast<void (TMapper::*)(const std::string&)> (&TMapper::build_tree), grp)
				//.add_method("build_tree_from_memory", static_cast<void (TMapper::*)(const std::vector<std::pair<std::vector<number>, number> >&)> (&TMapper::build_tree), grp)
				//.add_method("add_node_with_meta", static_cast<void (TMapper::*)(const std::pair<std::vector<number>, number>&)> (&TMapper::add_node), grp)
				.add_method("add_node_with_meta", static_cast<void (TMapper::*)(const std::vector<number>&, const number&)> (&TMapper::add_node), grp)
				.add_method("get_data_from_nn", static_cast<number (TMapper::*)(const std::vector<number>&) const> (&TMapper::get_data_from_nearest_neighbor), grp)
#ifndef UG_FOR_VRL
				//.add_method("add_node_with_meta_mv", static_cast<void (TMapper::*)(const std::pair<MathVector<3, number>, number>&)> (&TMapper::add_node), grp)
				.add_method("get_data_from_nn_mv", static_cast<number (TMapper::*)(const MathVector<3, number>&) const> (&TMapper::get_data_from_nearest_neighbor), grp)
#endif
			;
		}

		// MPM Mapper
#ifdef MPMNEURON
		{
			typedef membrane_potential_mapping::NeuronMPM TMPMMapper;
			std::string name("NeuronMPM");

			reg.add_class_<TMPMMapper>(name, grp)
			   .add_constructor()
			   .add_method("set_transformator", static_cast<void (TMPMMapper::*)(SmartPtr<Transformator>)> (&TMPMMapper::set_transformator), "", "", "")
			   .add_method("get_vm", static_cast<number (TMPMMapper::*)(number, number, number)> (&TMPMMapper::get_vm), "", "", "")
			   .add_method("build_tree", static_cast<void (TMPMMapper::*)()> (&TMPMMapper::build_tree), "", "", "")
			   .add_method("get_mapper", static_cast<SmartPtr<Mapper<3, number> > (TMPMMapper::*)() const> (&TMPMMapper::get_mapper), "", "", "")
			   .set_construct_as_smart_pointer(true);
		}
#endif
		/*
		// Vm2uG (\see vm2ug.h)
		{
		   typedef membrane_potential_mapping::Vm2uG<std::string> TVm2uG;
		   std::string name("MembranePotentialMapper");
		   reg.add_class_<TVm2uG>(name, grp)
				.add_constructor<void (*)(const std::string&, const std::string&, bool)>
				   ("Initial timestep|load-dialog|endings=[\"csv\", \"txt\"];description=\"Timestep files\"#Suffix|selection|value=[\"csv\", \"txt\"]#Static Nodes|selection|value=[True, False]")
#ifdef MPMNEURON
				.add_constructor<void (*)(SmartPtr<Transformator>)>("Transformation setup")
#endif
				.add_method("build_tree", static_cast<void (TVm2uG::*)(const std::string&)>(&TVm2uG::buildTree), grp)
				.add_method("get_potential", static_cast<number (TVm2uG::*)(number, number, number, const std::string&)> (&TVm2uG::get_potential),
					"Potential|default", "x|default#y|default#z|default#Timestep|default", grp)
				.add_method("get_potential_lin", static_cast<number (TVm2uG::*)(number, number, number, const std::string&, number, size_t)> (&TVm2uG::get_potential_lin),
					"Potential|default", "x|default#y|default#z|default#Timestep|default#cutoff|default#k|default")
				.add_method("get_potential_bilin", static_cast<number (TVm2uG::*)(number, number, number, const std::string&, number, size_t)> (&TVm2uG::get_potential_bilin),
					"Potential|default", "x|default#y|default#z|default#Timestep|default#cutoff|default#k|default")
#ifdef MPMNEURON
				.add_method("build_tree", (void (TVm2uG::*)()) (&TVm2uG::buildTree), grp)
				.add_method("get_potential", (number (TVm2uG::*)(number, number, number)) (&TVm2uG::get_potential), grp)
#endif
				.set_construct_as_smart_pointer(true);
		}*/

#if (defined(MPMVGCC) && (MPVGCC == 1 || MPVGCC == 2))
		// registry of BG (\see bg.h)
		{
			typedef membrane_potential_mapping::bg::BG TBG;
			std::string name("BorgGraham");
			reg.add_class_<TBG>(name, grp)
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
				.set_construct_as_smart_pointer(true);
		}
#endif
	}

	/**
	 * Function called for the registration of Dimension dependent parts.
	 * All Functions and Classes depending on the Dimension
	 * are to be placed here when registering. The method is called for all
	 * available Dimension types, based on the current build options.
	 *
	 * @param reg   registry
	 * @param grp   group for sorting of functionality
	 */
	template <size_t dim>
	static void Dimension(Registry& reg, std::string grp)
	{
		std::string suffix = GetDimensionSuffix<dim>();
		std::string tag = GetDimensionTag<dim>();

		// KD node
		{
			typedef kd_node<dim, number> TKDNode;
			std::string name = std::string("KDNode").append(suffix);
			reg.add_class_<TKDNode>(name, grp)
				.template add_constructor<void (*)()>("", "", "");
			reg.add_class_to_group(name, std::string("KDNode"), tag);
		}

		// KD tree
		{
			typedef kd_tree<dim, number> TKDTree;
			std::string name = std::string("KDTree").append(suffix);
			reg.add_class_<TKDTree>(name, grp)
				.template add_constructor<void (*)()>("", "", "")
				/// UG_FOR_VRL activated iff TARGET=vrl
#ifndef UG_FOR_VRL
				.add_method("add_node_meta", (void (TKDTree::*)(const MathVector<dim, number>&, number))(&TKDTree::add_node_meta), grp)
				.add_method("add_node", (void (TKDTree::*)(const MathVector<dim, number>&))(&TKDTree::add_node), grp)
				.add_method("build_tree", (bool (TKDTree::*)()) &TKDTree::build_tree)
				.add_method("query", (number (TKDTree::*)(const MathVector<dim>&)) &TKDTree::query)
#endif
			;
			reg.add_class_to_group(name, std::string("KDTree"), tag);
		}
	}

	/**
	 * Function called for the registration of Domain dependent parts.
	 * All Functions and Classes depending on the Domain
	 * are to be placed here when registering. The method is called for all
	 * available Domain types, based on the current build options.
	 *
	 * @param reg   registry
	 * @param grp   group for sorting of functionality
	 */
	template <typename TDomain>
	static void Domain(Registry& reg, std::string grp)
	{
		std::string suffix = GetDomainSuffix<TDomain>();
		std::string tag = GetDomainTag<TDomain>();
	}

	/**
	 * Function called for the registration of Domain and Algebra dependent parts.
	 * All Functions and Classes depending on both Domain and Algebra
	 * are to be placed here when registering. The method is called for all
	 * available Domain and Algebra types, based on the current build options.
	 *
	 * @param reg   registry
	 * @param grp   group for sorting of functionality
	 */
	template <typename TDomain, typename TAlgebra> static void
	DomainAlgebra(ug::bridge::Registry& reg, const std::string& grp)
	{
		std::string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
		std::string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();

		typedef GridFunction<TDomain, TAlgebra> TGridFunction;

#ifdef MPMNEURON
		// Transformator
		{
			typedef membrane_potential_mapping::Transformator TTransformator;
			std::string name = std::string("Transformator").append(suffix);
			reg.add_class_<TTransformator>(name, grp)
				.add_constructor<void (*)()>("")
				.add_method("load_geom", (void (TTransformator::*)(const std::string&))(&TTransformator::load_geom),
					"", "hoc geometry|load-dialog", grp)
				.add_method("load_stim", (void (TTransformator::*)(const std::string&))(&TTransformator::load_stim),
					"", "hoc stimulation protocol|load-dialog", grp)
				.add_method("get_sections", (size_t (TTransformator::*)())(&TTransformator::get_sections),
					"number of sections", "", grp)
				.add_method("get_total_points", (size_t (TTransformator::*)())(&TTransformator::get_total_points),
					"number of total points", "", grp)
				.add_method("get_t", (number (TTransformator::*)())(&TTransformator::get_t), "t", "", grp)
				.add_method("get_dt", (number (TTransformator::*)())(&TTransformator::get_dt), "dt", "", grp)
				.add_method("get_tstop", (number (TTransformator::*)())(&TTransformator::get_tstop), "tstop", "", grp)
				.add_method("get_tstart", (number (TTransformator::*)())(&TTransformator::get_tstart), "tstart", "", grp)
				.add_method("get_finitialize", (number (TTransformator::*)())(&TTransformator::get_finitialize),
					"finitialize", "", grp)
				.add_method("setup_hoc", (void (TTransformator::*)(number, number, number, number))(&TTransformator::setup_hoc),
					"", "tstart#tstop#dt#finitialize", grp)
				.add_method("extract_vms", (void (TTransformator::*)(size_t, size_t))(&TTransformator::extract_vms),
					"", "limit", grp)
				.add_method("purge", (void (TTransformator::*)())(&TTransformator::purge), "", "reinit", grp)
				.add_method("print_setup", (void (TTransformator::*)(bool))(&TTransformator::print_setup),
					"", "prints setup", grp)
				.add_method("adjust_resistivity_obstacle", (void (TTransformator::*)(const char* , const char* , double ))
					(&TTransformator::adjust_resistivity_obstacle), "", "", grp)
				.add_method("feedback", (void (TTransformator::*)(const char* uCmp, const char* , const char* ,double ))
					(&TTransformator::feedback), "", "", grp)
				.add_method("feedbacks", (void (TTransformator::*)(const char* uCmp, const char* , const char* ,double ))
					(&TTransformator::feedbacks), "", "", grp)
				.add_method("set_feedback", (void (TTransformator::*)())(&TTransformator::set_feedback), "", "", grp)
				.add_method("set_feedbacks", (void (TTransformator::*)())(&TTransformator::set_feedback), "", "", grp)
				.add_method("set_grid", (void (TTransformator::*)(TGridFunction&, const char*))
					(&TTransformator::set_grid<TGridFunction, TDomain>), "", "", grp)
				.add_method("init_feedback", (void (TTransformator::*)(const char* uCmp, const char* , const char* ,double, const char*))
					(&TTransformator::init_feedback), "", "", grp)
				.add_method("init_feedbacks", (void (TTransformator::*)(const char* uCmp, const char* , const char* ,double, const char* ))
					(&TTransformator::init_feedbacks), "", "", grp)
				.add_method("execute_hoc_stmt", (number (TTransformator::*)(const std::string& stmt))
					(&TTransformator::execute_hoc_stmt), "success or failure", "any valid hoc statement", grp)
				.add_method("set_hoc_variable", (bool (TTransformator::*)(const std::string& var, number value))
					(&TTransformator::set_hoc_variable), "sucess or failure", "variable|value", grp)
				.add_method("set_hoc_variable_section", (bool (TTransformator::*)(const std::string& var, number value, const std::string& section))
					(&TTransformator::set_hoc_variable_sec), "sucess or failure", "variable|value|section", grp)
				.add_method("get_hoc_variable", (number (TTransformator::*)(const std::string& var))(&TTransformator::get_hoc_variable),
					"sucess or failure", "variable", grp)
				.add_method("get_hoc_variable_section", (number (TTransformator::*)(const std::string& var, const std::string& section))
					(&TTransformator::get_hoc_variable_sec), "sucess or failure", "variable|section", grp)
				.add_method("fadvance", (bool (TTransformator::*)())(&TTransformator::fadvance), "success or failure", "fadvance", grp)
				.add_method("get_transformator", (SmartPtr<TTransformator> (TTransformator::*)())(&TTransformator::get_transformator),
					"", "HOC Interpreter", grp)
#ifdef MPMNEURON_REVISION
				.add_method("get_section_names_all", (std::vector<std::string> (TTransformator::*)())
					(&TTransformator::get_all_sections), "", "", grp)
				.add_method("get_section_names_as_string", (std::string (TTransformator::*)())
					(&TTransformator::get_all_sections_as_string), "", "", grp)
#endif
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "Transformator", tag);
		}
#endif
	}
}; // end of functionality which is to be exported
/// \}

} // namespace membrane_potential_mapping



/// \addtogroup mpm_plugin
extern "C" void
InitUGPlugin_MembranePotentialMapping(bridge::Registry& reg, std::string& grp)
{
	grp.append("/MembranePotentialMapping");
	typedef membrane_potential_mapping::Functionality Functionality;

	/// Note:
	/// functionality only implemented for 3D: in case of necessity for 2D/1D implementation
	/// one needs to adapt slightly the get_potential and build_tree functions, since they
	/// rely on 3D points and 1D time data in CSV files, or in case of NEURON interpreter
	/// embedded within Vm2uG - if you have other data you may adapt all functions of Vm2uG
	try
	{
		// register domain/algebra first as transformator is needed in vm2ug
		RegisterDomain3dAlgebraDependent<Functionality>(reg, grp); // only 3d supported
		RegisterDimensionDependent<Functionality>(reg, grp);
		RegisterCommon<Functionality>(reg, grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

} // namespace ug

