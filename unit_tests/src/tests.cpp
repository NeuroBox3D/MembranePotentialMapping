/*!
 * \file plugins/experimental/membrane_potential_mapping/unit_tests/src/tests.cpp
 * \brief boost unit test for the membrane_potenial_mapping ug plugin
 *
 * \date created on Apr 27, 2012
 * \author Stephan Grein
 */

#define BOOST_TEST_MODULE __CPP__UNIT_TESTS__UG__MEMBRANE_POTENTIAL_MAPPING__

#include <boost/test/included/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <cmath>

#include "../inc/unit_test_helper.h"
#include "../inc/fixtures.cpp"

#ifdef MPMVGCC
#include "../../bg_default/bg.h"
#else
#include "../../bg_simple/bg.h"
#endif

//#include "../../vm2ug.h"
#include "../../mvec.h"
#include "../../common_typedefs.h"

#include "../../vm2ug_rework.h"

using namespace boost::unit_test;
using namespace ug::membrane_potential_mapping;
using namespace ug::membrane_potential_mapping::bg;
using std::string;
using ug::membrane_potential_mapping::SMALL;
using ug::membrane_potential_mapping::VERY_SMALL;

// BOOST Test Suite for testing of mvec class
BOOST_AUTO_TEST_SUITE(vec);

// test basic constructor
BOOST_AUTO_TEST_CASE(test_constructor) {

	mvec<number, 3> m1;
	for (size_t i = 0; i < 3; i++)
			m1.push_back(2.0);

	for (std::vector<number>::const_iterator it = m1.begin(); it < m1.end(); it++)
		BOOST_CHECK_CLOSE(*it, 2.0, SMALL);


}

// test default constructor
BOOST_AUTO_TEST_CASE(test_default_constructor) {

	std::vector<number> a;
	for (size_t i = 0; i < 3; i++)
		a.push_back(2.0);

	mvec<number, 3> m1(a);
	for (std::vector<number>::const_iterator it = m1.begin(); it < m1.end(); it++)
		BOOST_CHECK_CLOSE(*it, 2.0, SMALL);


}

// test determinante
BOOST_AUTO_TEST_CASE(test_determinante) {


	std::vector<number> a;
	std::vector<number> b;

		a.push_back(-2.0);
		a.push_back(5.0);

		b.push_back(3.0);
		b.push_back(0.5);


	mvecd3 mvec1 = a;
	mvecd3 mvec2 = b;


	std::vector<mvecd3> mvecs;
	mvecs.push_back(mvec1);
	mvecs.push_back(mvec2);

	number ret = mvec<number, 3>::det(mvecs);
	BOOST_REQUIRE_CLOSE(-16.0, ret, SMALL);

	a.clear();
	b.clear();
	mvecs.clear();

	std::vector<number> c;

	a.push_back(2.0);
	a.push_back(-1.0);
	a.push_back(3.0);

	b.push_back(-2.0);
	b.push_back(5.0);
	b.push_back(4.0);

	c.push_back(0.0);
	c.push_back(1.0);
	c.push_back(5.0);

	mvec1 = a;
	mvec2 = b;
	mvecd3 mvec3 = c;

	mvecs.push_back(mvec1);
	mvecs.push_back(mvec2);
	mvecs.push_back(mvec3);

	number ret2 = mvec<number, 3>::det(mvecs);
	BOOST_REQUIRE_CLOSE(26.0, ret2, SMALL);


}

// test norm
BOOST_AUTO_TEST_CASE(test_norm) {

	std::vector<number> a;
	for (size_t i = 0; i < 3; i++) {
		a.push_back(2.0);
	}

	mvec<number, 3> m1(a);

	number res = m1.norm(EUCLIDEAN);
	number check = std::sqrt(12.0);

	BOOST_REQUIRE_CLOSE(res, check, SMALL);


}



// test add
BOOST_AUTO_TEST_CASE(test_add) {

	std::vector<number> a;
	std::vector<number> b;
	for (size_t i = 0; i < 3; i++) {
		a.push_back(2.0);
		b.push_back(1.0);
	}

	mvec<number, 3> m1(a);
	mvec<number, 3> m2(b);
	mvec<number, 3> c = m1 + m2;

	for (std::vector<number>::const_iterator it = c.begin(); it < c.end(); it++) {
		BOOST_REQUIRE_CLOSE(*it, 3.0, SMALL);
	}

}

// test sub
BOOST_AUTO_TEST_CASE(test_sub) {

	std::vector<number> a;
	std::vector<number> b;

	for (size_t i = 0; i < 3; i++) {
		a.push_back(2.0);
		b.push_back(1.0);
	}

	mvec<number, 3> m1(a);
	mvec<number, 3> m2(b);
	mvec<number, 3> c = m1 - m2;

	for (std::vector<number>::const_iterator it = c.begin(); it < c.end(); it++) {
		BOOST_REQUIRE_CLOSE(*it, 1.0, SMALL);
	}

}

// test vec
BOOST_AUTO_TEST_CASE(test_vec) {

	std::vector<number> a;
	std::vector<number> b;
	for (size_t i = 0; i < 3; i++) {
		a.push_back(2.0);
		b.push_back(1.0);
	}

	mvec<number, 3> m1(a);
	mvec<number, 3> m2(b);
	mvec<number, 3> c = m1 % m2;

	for (std::vector<number>::const_iterator it = c.begin(); it < c.end(); it++) {
		BOOST_REQUIRE_CLOSE(*it, 0.0, SMALL);
	}

}

// test dot
BOOST_AUTO_TEST_CASE(test_dot) {

	std::vector<number> a;
	for (size_t i = 0; i < 3; i++) {
		a.push_back(2.0);
	}

	mvec<number, 3> m1(a);

	number c = m1 * m1;

	BOOST_REQUIRE_CLOSE(c, 12.0, SMALL);


}

// test neg
BOOST_AUTO_TEST_CASE(test_neg) {

	std::vector<number> a;
	for (size_t i = 0; i < 3; i++) {
		a.push_back(2.0);
	}

	mvec<number, 3> m1(a);
	mvec<number, 3> m2 = -m1;

	for (size_t i = 0; i < 3; i++)
	BOOST_REQUIRE_CLOSE(m2[i], -m1[i], SMALL);

}

// test id
BOOST_AUTO_TEST_CASE(test_id) {

	std::vector<number> a;
	for (size_t i = 0; i < 3; i++) {
		a.push_back(2.0);
	}
	mvec<number, 3> m1(a);
	mvec<number, 3> m2 = +m1;

	for (size_t i = 0; i < 3; i++) {
		BOOST_REQUIRE_CLOSE(m2[i], m1[i], SMALL);
	}


}

// test assignment
BOOST_AUTO_TEST_CASE(test_assignment) {

	std::vector<number> a;
	for (size_t i = 0; i < 3; i++)
	a.push_back(2.0);
	mvec<number, 3> m1 = a;

	for (std::vector<number>::const_iterator cit = m1.begin(); cit < m1.end(); cit++)
	BOOST_REQUIRE_CLOSE(*cit, 2.0, SMALL);


}

// test assignment plus
BOOST_AUTO_TEST_CASE(test_assignment_plus) {

	std::vector<number> a;
	for (size_t i = 0; i < 3; i++)
		a.push_back(2.0);

	std::vector<number> b;
	for (size_t i = 0; i < 3; i++)
		b.push_back(1.0);

	mvec<number, 3> m1 = a;
	mvec<number, 3> m2 = b;

	m2 += m1;

	for (std::vector<number>::const_iterator cit = m2.begin(); cit < m2.end(); cit++)
		BOOST_REQUIRE_CLOSE(*cit, 3.0, SMALL);

	for (std::vector<number>::const_iterator cit = m1.begin(); cit < m1.end(); cit++)
		BOOST_REQUIRE_CLOSE(*cit, 2.0, SMALL);


}

// test assignment minus
BOOST_AUTO_TEST_CASE(test_assignment_minus) {

	std::vector<number> a;
	for (size_t i = 0; i < 3; i++)
		a.push_back(2.0);

	std::vector<number> b;
		for (size_t i = 0; i < 3; i++)
	b.push_back(1.0);

		mvec<number, 3> m1 = a;
	mvec<number, 3> m2 = b;

	m1 -= m2;

	for (std::vector<number>::const_iterator cit = m2.begin(); cit < m2.end(); cit++)
		BOOST_REQUIRE_CLOSE(*cit, 1.0, SMALL);

	for (std::vector<number>::const_iterator cit = m1.begin(); cit < m1.end(); cit++)
		BOOST_REQUIRE_CLOSE(*cit, 1.0, SMALL);


}

//test assignment vec
BOOST_AUTO_TEST_CASE(test_assignment_vec) {

	std::vector<number> a;
	for (size_t i = 0; i < 3; i++)
		a.push_back(2.0);

	std::vector<number> b;
	for (size_t i = 0; i < 3; i++)
		b.push_back(1.0);

	mvec<number, 3> m1 = a;
	mvec<number, 3> m2 = b;

	m1 %= m2;

	for (std::vector<number>::const_iterator cit = m2.begin(); cit < m2.end(); cit++)
		BOOST_REQUIRE_CLOSE(*cit, 1.0, SMALL);

	for (std::vector<number>::const_iterator cit = m1.begin(); cit < m1.end(); cit++)
		BOOST_CHECK_CLOSE(*cit, 0.0, SMALL);


}

BOOST_AUTO_TEST_SUITE_END();

/*
// BOOST Test Suite for testing of vm2ug class
BOOST_FIXTURE_TEST_SUITE(vm2ug, FixtureVUG<string>);

// test constructor
BOOST_AUTO_TEST_CASE(construct_Vm2uG) {

	BOOST_REQUIRE_MESSAGE(vm2ug, "Vm2uG<string> instance cannot be constructed");

}

// test build_tree
#ifndef MPMNEURON
BOOST_AUTO_TEST_CASE(build_tree) {

	BOOST_REQUIRE_MESSAGE(vm2ug, "Vm2uG<string> instance cannot be constructed");
	vm2ug->buildTree("timestep0.000000.csv");
	BOOST_CHECK_MESSAGE(vm2ug->treeBuilt(), "tree could not be rebuild");

}
#endif

#ifndef MPMNEURON
// test get_potential
BOOST_AUTO_TEST_CASE(get_potential) {

	BOOST_REQUIRE_MESSAGE(vm2ug, "Vm2uG<string> instance cannot be constructed");
	vm2ug->buildTree("timestep0.000000.csv");
	BOOST_CHECK_MESSAGE(vm2ug->treeBuilt(), "tree could not be rebuild");
	BOOST_CHECK_MESSAGE(vm2ug->get_potential(0,0,0, "timestep0.000000.csv") == -75.0, "initial potential should be -75.0");

}
#endif

BOOST_AUTO_TEST_SUITE_END();
*/

// BOOST Test Suite for testing of bg class
BOOST_FIXTURE_TEST_SUITE(bg, FixtureBG);

// test construct BG channel
BOOST_AUTO_TEST_CASE(constructInstanceBG) {

	BOOST_REQUIRE_MESSAGE(bg, "BG instance cannot be constructed");

}

// test install different gate types
BOOST_AUTO_TEST_CASE(install_gates) {

	BOOST_REQUIRE_MESSAGE(bg, "BG instance cannot be constructed");
	bg->install_can_gates();
	BOOST_CHECK_MESSAGE(bg->installed_can_gates(), "can Gates could not be installed");
	bg->install_cal_gates();
	BOOST_CHECK_MESSAGE(bg->installed_cal_gates(), "cal Gates could not be installed");
	bg->install_cat_gates();
	BOOST_CHECK_MESSAGE(bg->installed_cat_gates(), "cal Gates could not be installed");

}

// check fluxes for the ohmic BG model
#ifndef MPMVGCC
BOOST_AUTO_TEST_CASE(check_fluxes_ohmic) {

	BG* b = new BG();
	BOOST_CHECK_MESSAGE(b, "BG instance cannot be constructed");
	b->install_can_gates(1000.0);
	b->calc_current_at_start(0);
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
		number d = b->timestepping_of_gates_and_calc_current(delta_t * i, delta_t);
		BOOST_REQUIRE_CLOSE(d, results[i], SMALL);
	}

}

//  check molar fluxes for the ohmic BG model
BOOST_AUTO_TEST_CASE(check_molar_fluxes_ohmic) {

	BG* b = new BG();
	BOOST_CHECK_MESSAGE(b, "BG instance cannot be constructed");
	b->install_can_gates(1000.0);
	b->calc_current_at_start(0);
	std::vector<number> results;
	results.push_back(8.4684798976316203e-06);
	results.push_back(8.4684798976316203e-06);
	results.push_back(4.0437351365245749e-05);
	number delta_t = 0.001;
	for (size_t i = 0; i < results.size(); i++) {
		b->timestepping_of_gates_and_calc_current(delta_t * i, delta_t);
		number d = b->get_Neumann_Flux_as_Concentration(delta_t);
	    BOOST_REQUIRE_CLOSE(d, results[i], SMALL);
	}

}
#else
// check fluxes for the cfp BG model
BOOST_AUTO_TEST_CASE(check_fluxes_cfp) {

	BG* b = new BG();
	BOOST_CHECK_MESSAGE(b, "BG instance cannot be constructed");
	b->install_can_gates_cfp();
	b->calc_current_at_start(0.0);
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
		number d = b->timestepping_of_gates_and_calc_current(delta_t * i, delta_t, -55.0);
		BOOST_REQUIRE_CLOSE(d, results[i], SMALL);
	}

}
// check molar fluxes for the cfp BG model
BOOST_AUTO_TEST_CASE(check_molar_fluxes_cfp) {

	BG* b = new BG();
    BOOST_CHECK_MESSAGE(b, "BG instance cannot be constructed");
	BOOST_CHECK_MESSAGE(b, "BG instance cannot be constructed");
	b->install_can_gates_cfp();
	b->calc_current_at_start(0.0);
	std::vector<number> results;
	results.push_back(1.9214352053711674e-29);
	results.push_back(2.3249826893955825e-29);
	results.push_back(2.4534319739280018e-29);
	number delta_t = 0.001;
	for (size_t i = 0; i < results.size(); i++) {
		b->timestepping_of_gates_and_calc_current(delta_t * i, delta_t, -55.0);
		number d = b->get_Neumann_Flux_as_Concentration(delta_t);
		BOOST_REQUIRE_CLOSE(d, results[i], SMALL);
	}

}
#endif

BOOST_AUTO_TEST_SUITE_END();

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
		if ((i % portion) == 0) saved.push_back(make_pair(vec, vm));
	}

	// build tree
	tree->build_tree();

	// query tree and assert vms are same
	for (ITER it = saved.begin(); it != saved.end(); ++it)
		BOOST_REQUIRE_CLOSE(it->second, tree->query(it->first), 1e-6);
}

BOOST_AUTO_TEST_SUITE_END();
