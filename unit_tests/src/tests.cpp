/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Stephan Grein
 * Creation date: 2012-04-27
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

#define BOOST_TEST_MODULE __CPP__UNIT_TESTS__UG__MEMBRANE_POTENTIAL_MAPPING__

#include <boost/test/included/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <cmath>

#include "../inc/unit_test_helper.h"
#include "../inc/fixtures.cpp"

#include "mpm_config.h"  // for project-specific defines

#if (defined(MPMVGCC) && MPMVGCC == 1)
#include "../../bg_simple/bg.h"
#elif (defined(MPMVGCC) && MPMVGCC == 2)
#include "../../bg_default/bg.h"
#endif


#include "../../vm2ug.h"
#include "../../kdtree/kd_tree.h"


using namespace boost::unit_test;
using namespace ug::membrane_potential_mapping;
using ug::membrane_potential_mapping::SMALL;
using ug::membrane_potential_mapping::VERY_SMALL;



#if (defined(MPMVGCC) && (MPMVGCC == 1 || MPMVGCC == 2))

using namespace ug::membrane_potential_mapping::bg;

// BOOST Test Suite for testing of bg class
BOOST_FIXTURE_TEST_SUITE(bg, FixtureBG);

// test construct BG channel
BOOST_AUTO_TEST_CASE(constructInstanceBG) {

	BOOST_REQUIRE_MESSAGE(bg, "BG instance cannot be constructed");
}

// test install different gate types
BOOST_AUTO_TEST_CASE(install_gates) {

	bg->install_can_gates();
	BOOST_CHECK_MESSAGE(bg->installed_can_gates(), "can Gates could not be installed");
	bg->install_cal_gates();
	BOOST_CHECK_MESSAGE(bg->installed_cal_gates(), "cal Gates could not be installed");
	bg->install_cat_gates();
	BOOST_CHECK_MESSAGE(bg->installed_cat_gates(), "cat Gates could not be installed");
}

// check fluxes for the ohmic BG model
#if MPMVGCC == 1
BOOST_AUTO_TEST_CASE(check_fluxes_ohmic) {

	bg->install_can_gates(1000.0);
	bg->calc_current_at_start(0);
	std::vector<number> results;
	results.push_back(1.6345251904980005);
	results.push_back(1.6345251904980005);
	results.push_back(7.80493);
	results.push_back(41460.8);
	results.push_back(14315.6);
	results.push_back(1946.1);
	results.push_back(240.059);
	results.push_back(36.044762202708);
	number delta_t = 0.001;
	for (size_t i = 0; i < results.size(); i++) {
		number d = bg->timestepping_of_gates_and_calc_current(delta_t * i, delta_t);
		BOOST_CHECK_CLOSE(d, results[i], SMALL);
	}
}

//  check molar fluxes for the ohmic BG model
BOOST_AUTO_TEST_CASE(check_molar_fluxes_ohmic) {

	bg->install_can_gates(1000.0);
	bg->calc_current_at_start(0);
	std::vector<number> results;
	results.push_back(8.4684798976316203e-06);
	results.push_back(8.4684798976316203e-06);
	results.push_back(4.0437351365245749e-05);
	number delta_t = 0.001;
	for (size_t i = 0; i < results.size(); i++) {
		bg->timestepping_of_gates_and_calc_current(delta_t * i, delta_t);
		number d = bg->get_Neumann_Flux_as_Concentration(delta_t);
	    BOOST_CHECK_CLOSE(d, results[i], SMALL);
	}
}

#elif MPMVGCC == 2
// check fluxes for the cfp BG model
BOOST_AUTO_TEST_CASE(check_fluxes_cfp) {

	bg->install_can_gates_cfp();
	bg->calc_current_at_start(0.0);
	std::vector<number> results;
	results.push_back(3.70861628e-09);
	results.push_back(4.48751466e-09);
	results.push_back(4.73544e-09);
	results.push_back(4.80714e-09);
	results.push_back(4.820455125e-09);
	results.push_back(4.81453e-09);
	results.push_back(4.802373071129e-09);
	results.push_back(4.788292607170e-09);
	results.push_back(4.773723314275e-09);
	number delta_t = 0.001;
	for (size_t i = 0; i < results.size(); i++) {
		number d = bg->timestepping_of_gates_and_calc_current(delta_t * i, delta_t, -55.0);
		BOOST_CHECK_CLOSE(d, results[i], SMALL);
	}

}
// check molar fluxes for the cfp BG model
BOOST_AUTO_TEST_CASE(check_molar_fluxes_cfp) {

	bg->install_can_gates_cfp();
	bg->calc_current_at_start(0.0);
	std::vector<number> results;
	results.push_back(1.9214352053711674e-29);
	results.push_back(2.3249826893955825e-29);
	results.push_back(2.4534319739280018e-29);
	number delta_t = 0.001;
	for (size_t i = 0; i < results.size(); i++) {
		bg->timestepping_of_gates_and_calc_current(delta_t * i, delta_t, -55.0);
		number d = bg->get_Neumann_Flux_as_Concentration(delta_t);
		BOOST_CHECK_CLOSE(d, results[i], SMALL);
	}
}
#endif

BOOST_AUTO_TEST_SUITE_END();

#endif


/// BOOST Test Suite for avoid warning
BOOST_AUTO_TEST_SUITE(dummy);
BOOST_AUTO_TEST_CASE(dummy) {
	BOOST_CHECK_SMALL(0.0, SMALL);
}
BOOST_AUTO_TEST_SUITE_END();




/// vm2ug rework tests
BOOST_AUTO_TEST_SUITE(MAPPER_REWORK);

/// tests neuron implementation of mpm
BOOST_AUTO_TEST_CASE(NEURON_MPM) {

}

/// tests vm2ug implementation of mpm
BOOST_AUTO_TEST_CASE(VM2UG_MPM) {

}

/// tests kdtree
BOOST_AUTO_TEST_CASE(KDTREE) {
	using namespace ug;
	typedef kd_tree<3, number> KDTree3;
	typedef std::vector<std::pair<MathVector<3>, number> >::const_iterator ITER;

	// rand setup
	number LO = 0;
	number HI = 10000;
	srand (static_cast <unsigned> (time(0)));

	// tree setup
	SmartPtr<KDTree3> tree(make_sp(new KDTree3()));
	std::vector<std::pair<MathVector<3>, number> >saved;
	size_t treeSize = 1e6;
	size_t numberQuerys = 1e2;
	size_t portion = treeSize - numberQuerys;

	// popoulate tree
	for (size_t i = 0; i < treeSize; i++) {
		MathVector<3> vec(
				LO + static_cast <number> (rand()) /( static_cast <number> (RAND_MAX/(HI-LO))),
				LO + static_cast <number> (rand()) /( static_cast <number> (RAND_MAX/(HI-LO))),
				LO + static_cast <number> (rand()) /( static_cast <number> (RAND_MAX/(HI-LO))));
		number vm = LO + static_cast <number> (rand()) /( static_cast <number> (RAND_MAX/(HI-LO)));
		tree->add_node_meta(vec, vm);
		if ((i % portion) == 0) saved.push_back(std::make_pair(vec, vm));
	}

	// build tree
	tree->build_tree();

	// query tree and assert vms are same
	for (ITER it = saved.begin(); it != saved.end(); ++it)
		BOOST_CHECK_CLOSE(it->second, tree->query(it->first), 1e-6);
}

BOOST_AUTO_TEST_SUITE_END();
