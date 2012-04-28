/*
 * fixtures.cpp
 *
 *  Created on: Apr 27, 2012
 *      Author: stephan grein
 */

#include <boost/test/included/unit_test.hpp>
#include "../../vm2ug.h"
#include "../../mvec.h"

//#ifdef FLAVOR
#include "../../bg.h"
//#else
// #include "../bg_simple/bg.h"
//#endif

#include "../../common_typedefs.h"

using namespace boost::unit_test;
using namespace vug;
using namespace bg;

template <class T> struct Fixture {
	Fixture() : vm2ug(new Vm2uG<T>("","")) { BOOST_TEST_MESSAGE("setup fixture >>vm2ug<<"); }
	~Fixture() { BOOST_TEST_MESSAGE("teardown fixture >>vm2ug<<"); }

	Vm2uG<T>* vm2ug;
};

struct Fixture2 {
	Fixture2() : bg(new BG()) { BOOST_TEST_MESSAGE("setup fixture >>vm2ug<<"); }
	~Fixture2() { BOOST_TEST_MESSAGE("teardown fixture >>bg<<"); }
	BG* bg;
};
