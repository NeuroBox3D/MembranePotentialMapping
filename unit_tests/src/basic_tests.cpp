/*
 * basic_tests.cpp
 *
 * author stephanmg
 * date 01m12d12y
 *
 */

#define BOOST_TEST_MODULE basic_tests

#include <boost/test/included/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>

#include "../../bg.h"
#include "../../vm2ug.h"
#include "../inc/unit_test_helper.h"

using namespace boost::unit_test;
using namespace vug;
using namespace bg;
using std::string;

BOOST_AUTO_TEST_SUITE(VM2UGBG1);

BOOST_AUTO_TEST_CASE(constructVm2uG) {
   BOOST_MESSAGE("Starting test >>constructVm2uG<<");
   Vm2uG<string>* vm2ug = new Vm2uG<string>("", "");
   BOOST_CHECK_MESSAGE(vm2ug, "Vm2uG<string> instance cannot be constructed");
}

BOOST_AUTO_TEST_CASE(buildTree) {
   BOOST_MESSAGE("Starting test >>buildTree<<");
   Vm2uG<string>* vm2ug = new Vm2uG<string>("", "");
   vm2ug->buildTree("timestep0.000000.csv");
   BOOST_CHECK_MESSAGE(vm2ug->treeBuild(), "tree could not be rebuild");
}

BOOST_AUTO_TEST_CASE(get_potential) {
   BOOST_MESSAGE("Starting test >>get_potential<<");
   Vm2uG<string>* vm2ug = new Vm2uG<string>("", "");
   vm2ug->buildTree("timestep0.000000.csv");
   BOOST_CHECK_MESSAGE(vm2ug->treeBuild(), "tree could not be rebuild");
   BOOST_CHECK_MESSAGE(!vm2ug->get_potential(0,0,0, "timestep0.000000.csv") == -75.0, "initial potential should be -75.0");
}

BOOST_AUTO_TEST_CASE(constructInstanceBG) {
   BOOST_MESSAGE("Starting test <<constructBG>>");
   BG* b = new BG(); 
   BOOST_CHECK_MESSAGE(b, "BG instance cannot be constructed"); 
}

BOOST_AUTO_TEST_CASE(installGates) {
   BOOST_MESSAGE("Starting test <<installGates>>");
   BG* b = new BG();
   BOOST_CHECK_MESSAGE(b, "BG instance cannot be constructed"); 
   b->install_can_gates();
   BOOST_CHECK_MESSAGE(b->installed_can_gates(), "can Gates could not be installed");
}

BOOST_AUTO_TEST_CASE(checkFluxesAndCurrents) {
   BOOST_MESSAGE("Starting test <<flux>>");
   BG* b = new BG();
   BOOST_CHECK_MESSAGE(b, "BG instance cannot be constructed");
   b->install_can_gates();
   b->calc_current_at_start(0);
   std::vector<double> results;
   results.push_back(1.63453);
   results.push_back(1.63453);
   results.push_back(7.80493);
   results.push_back(41460.8);
   results.push_back(14315.6);
   results.push_back(1946.1);
   results.push_back(240.059);
   results.push_back(36.0448);
   double delta_t = 0.001;
   for (int i = 0; i < results.size(); i++) {
      BOOST_CHECK_MESSAGE(AreSame(results[i], b->timestepping_of_gates_and_calc_current(delta_t * i, delta_t)), "Currents not equal!");
   }
}

BOOST_AUTO_TEST_SUITE_END();
