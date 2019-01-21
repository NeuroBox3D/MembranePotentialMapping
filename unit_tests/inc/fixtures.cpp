/**
 * \file plugins/experimental/membrane_potential_mapping/unit_tests/inc/fixtures.cpp
 * \brief fixtures for testing the membrane_potential_mapping ug plugin
 *
 * \date May 27, 2012
 * \author Stephan Grein
 */

// includes
#include <boost/test/included/unit_test.hpp>

// #include "../../vm2ug.h"
#include "../../mvec.h"
#include "../../common_typedefs.h"

#include "mpm_config.h"  // for project-specific defines

#ifdef MPMVGCC
#include "../../bg_default/bg.h"
#else
#include "../../bg_simple/bg.h"
#endif


// using directives
using namespace boost::unit_test;


/**
 * \brief fixture for class Vm2uG
template <class T> struct FixtureVUG {
	ug::membrane_potential_mapping::Vm2uG<T>* vm2ug;
	FixtureVUG() : vm2ug(new ug::membrane_potential_mapping::Vm2uG<T>("","")) { BOOST_TEST_MESSAGE("setup fixture >>vm2ugG<<");}
	~FixtureVUG() { BOOST_TEST_MESSAGE("teardown fixture >>vm2ug<<"); }

};
*/

/**
 * \brief fixture for class BG
 */
struct FixtureBG {
	ug::membrane_potential_mapping::bg::BG* bg;
	static size_t count;
	/**
	 * \brief construct a fixture for class bg
	 */
	FixtureBG() : bg(new ug::membrane_potential_mapping::bg::BG()) {  BOOST_TEST_MESSAGE("setup fixture >>bg<< no. " << ++count);}
	~FixtureBG() { BOOST_TEST_MESSAGE("teardown fixture >>bg<<"); }

};

// initialization of static members
size_t FixtureBG::count = 0;
