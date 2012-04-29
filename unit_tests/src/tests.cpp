/*
 * tests.cpp
 *
 *  Created on: Apr 27, 2012
 *      Author: stephan grein
 */
#define BOOST_TEST_MODULE tests

#include <boost/test/included/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <cmath>

#include "../inc/unit_test_helper.h"
#include "../inc/fixtures.cpp"

// #ifdef FLAVOR
#include "../../bg.h"
// #else
// #include "../bg_default/bg.h"
// #endif

#include "../../vm2ug.h"
#include "../../mvec.h"
#include "../../common_typedefs.h"

using namespace boost::unit_test;
using namespace vug;
using namespace bg;
using std::string;

BOOST_AUTO_TEST_SUITE(vec);

BOOST_AUTO_TEST_CASE(test_add)  {
	       std::vector<double> a;
		   std::vector<double> b;
		   for (size_t i = 0; i < 3; i++) {
			   a.push_back(2.0);
			   b.push_back(1.0);
		   }

		   mvec<double, 3> m1(a);
		   mvec<double, 3> m2(b);
		   mvec<double, 3> c = m1 + m2;

		   for (std::vector<double>::const_iterator it = c.begin(); it < c.end(); it++) {
			   BOOST_REQUIRE_CLOSE(*it, 3.0, SMALL);
		   }
}
BOOST_AUTO_TEST_CASE(test_sub)  {
	   std::vector<double> a;
	   std::vector<double> b;

	  		   for (size_t i = 0; i < 3; i++) {
	  			   a.push_back(2.0);
	  			   b.push_back(1.0);
	  		   }



	   mvec<double, 3> m1(a);
	   mvec<double, 3> m2(b);
	   mvec<double, 3> c = m1 - m2;

	   for (std::vector<double>::const_iterator it = c.begin(); it < c.end(); it++) {
	  			  BOOST_REQUIRE_CLOSE(*it, 1.0, SMALL);
	  		   }

}

BOOST_AUTO_TEST_CASE(test_vec)  {
	       std::vector<double> a;
		   std::vector<double> b;
		   for (size_t i = 0; i < 3; i++) {
			  			   a.push_back(2.0);
			  			   b.push_back(1.0);
			  		   }


		   mvec<double, 3> m1(a);
		   mvec<double, 3> m2(b);
		  mvec<double, 3> c = m1 % m2;

	   for (std::vector<double>::const_iterator it = c.begin(); it < c.end(); it++) {
	  			  BOOST_REQUIRE_CLOSE(*it, 0.0, SMALL);
	  		   }

}

BOOST_AUTO_TEST_CASE(test_dot)  {
	 	  std::vector<double> a;
	 	 for (size_t i = 0; i < 3; i++) {
	 				  			   a.push_back(2.0);
	 	 }


		   mvec<double, 3> m1(a);

		   double c = m1 * m1;

		   BOOST_REQUIRE_CLOSE(c, 12.0, SMALL);

}
BOOST_AUTO_TEST_CASE(test_neg)  {
	 std::vector<double> a;
		 	 for (size_t i = 0; i < 3; i++) {
		 				  			   a.push_back(2.0);
		 	 }

		 mvec<double, 3> m1(a);
		 mvec<double, 3> m2 = -m1;

		for (size_t i = 0; i < 3; i++)
			BOOST_REQUIRE_CLOSE(m2[i], -m1[i], SMALL);
}
BOOST_AUTO_TEST_CASE(test_id)   {
	 std::vector<double> a;
		 	 for (size_t i = 0; i < 3; i++) {
		 				  			   a.push_back(2.0);
		 	 }
	 mvec<double, 3> m1(a);
	 mvec<double, 3> m2 = +m1;

	for (size_t i = 0; i < 3; i++) {
		BOOST_REQUIRE_CLOSE(m2[i], m1[i], SMALL);
	}

}
BOOST_AUTO_TEST_CASE(test_norm) {

	           std::vector<double> a;
	          		 	 for (size_t i = 0; i < 3; i++) {
	          		 				  			   a.push_back(2.0);
	          		 	 }

			   mvec<double, 3> m1(a);

			   double res = m1.norm(EUCLIDEAN);
			   double check = std::sqrt(12.0);

			   BOOST_REQUIRE_CLOSE(res, check, SMALL);


}

BOOST_AUTO_TEST_CASE(test_assignment) {
	std::vector<double> a;
	for (size_t i = 0; i < 3; i++)
			a.push_back(2.0);
	mvec<double, 3> m1 = a;

	for (typename std::vector<double>::const_iterator cit = m1.begin(); cit < m1.end(); cit++)
		BOOST_REQUIRE_CLOSE(*cit, 2.0, SMALL);
}

BOOST_AUTO_TEST_CASE(test_assignment_plus) {
	std::vector<double> a;
		for (size_t i = 0; i < 3; i++)
				a.push_back(2.0);

		std::vector<double> b;
			for (size_t i = 0; i < 3; i++)
					b.push_back(1.0);
		mvec<double, 3> m1 = a;
		mvec<double, 3> m2 = b;

		m2 += m1;

		for (typename std::vector<double>::const_iterator cit = m2.begin(); cit < m2.end(); cit++)
				BOOST_REQUIRE_CLOSE(*cit, 3.0, SMALL);

		for (typename std::vector<double>::const_iterator cit = m1.begin(); cit < m1.end(); cit++)
				BOOST_REQUIRE_CLOSE(*cit, 2.0, SMALL);
}

BOOST_AUTO_TEST_CASE(test_assignment_minus) {
	std::vector<double> a;
		for (size_t i = 0; i < 3; i++)
				a.push_back(2.0);

		std::vector<double> b;
			for (size_t i = 0; i < 3; i++)
					b.push_back(1.0);
		mvec<double, 3> m1 = a;
		mvec<double, 3> m2 = b;

		m1 -= m2;

		for (typename std::vector<double>::const_iterator cit = m2.begin(); cit < m2.end(); cit++)
				BOOST_REQUIRE_CLOSE(*cit, 1.0, SMALL);

		for (typename std::vector<double>::const_iterator cit = m1.begin(); cit < m1.end(); cit++)
				BOOST_REQUIRE_CLOSE(*cit, 1.0, SMALL);

}
BOOST_AUTO_TEST_CASE(test_assignment_vector) {
	std::vector<double> a;
		for (size_t i = 0; i < 3; i++)
				a.push_back(2.0);

		std::vector<double> b;
			for (size_t i = 0; i < 3; i++)
					b.push_back(1.0);

		mvec<double, 3> m1 = a;
		mvec<double, 3> m2 = b;

		m1 %= m2;

		for (typename std::vector<double>::const_iterator cit = m2.begin(); cit < m2.end(); cit++)
				BOOST_REQUIRE_CLOSE(*cit, 1.0, SMALL);

		for (typename std::vector<double>::const_iterator cit = m1.begin(); cit < m1.end(); cit++)
				BOOST_CHECK_CLOSE(*cit, 0.0, SMALL);
}

BOOST_AUTO_TEST_SUITE_END();

BOOST_FIXTURE_TEST_SUITE(vm2ug, Fixture<string>);

BOOST_AUTO_TEST_CASE(construct_Vm2uG) {
   BOOST_MESSAGE("Starting test >>construct_Vm2uG<<");
   BOOST_REQUIRE_MESSAGE(vm2ug, "Vm2uG<string> instance cannot be constructed");
   BOOST_MESSAGE("End test >>construct_Vm2uG<<");
}

BOOST_AUTO_TEST_CASE(build_tree) {
   BOOST_MESSAGE("Starting test >>build_tree<<");
   BOOST_REQUIRE_MESSAGE(vm2ug, "Vm2uG<string> instance cannot be constructed");
   vm2ug->buildTree("timestep0.000000.csv");
   BOOST_CHECK_MESSAGE(vm2ug->treeBuild(), "tree could not be rebuild");
   BOOST_MESSAGE("End test >>build_tree<<");
}

BOOST_AUTO_TEST_CASE(get_potential) {
   BOOST_MESSAGE("Starting test >>get_potential<<");
   BOOST_REQUIRE_MESSAGE(vm2ug, "Vm2uG<string> instance cannot be constructed");
   vm2ug->buildTree("timestep0.000000.csv");
   BOOST_CHECK_MESSAGE(vm2ug->treeBuild(), "tree could not be rebuild");
   BOOST_CHECK_MESSAGE(vm2ug->get_potential(0,0,0, "timestep0.000000.csv") == -75.0, "initial potential should be -75.0");
   BOOST_MESSAGE("End test >>get_potential<<");
}

BOOST_AUTO_TEST_SUITE_END();

BOOST_FIXTURE_TEST_SUITE(bg, Fixture2);

BOOST_AUTO_TEST_CASE(constructInstanceBG) {
   BOOST_MESSAGE("Starting test >>constructBG<<");
   BOOST_REQUIRE_MESSAGE(bg, "BG instance cannot be constructed");
   BOOST_MESSAGE("End test >>constructBG<<");
}

BOOST_AUTO_TEST_CASE(install_gates) {
   BOOST_MESSAGE("Starting test >>install_gates<<");
   BOOST_REQUIRE_MESSAGE(bg, "BG instance cannot be constructed");
   bg->install_can_gates();
   BOOST_CHECK_MESSAGE(bg->installed_can_gates(), "can Gates could not be installed");
   bg->install_cal_gates();
   BOOST_CHECK_MESSAGE(bg->installed_cal_gates(), "cal Gates could not be installed");
   bg->install_cat_gates();
   BOOST_CHECK_MESSAGE(bg->installed_cat_gates(), "cal Gates could not be installed");
   BOOST_MESSAGE("End test >>install_gates<<");
}

#ifndef FLAVOR
BOOST_AUTO_TEST_CASE(check_fluxes) {
   BOOST_MESSAGE("Starting test >>flux<<");
   BG* b = new BG();
   BOOST_CHECK_MESSAGE(b, "BG instance cannot be constructed");
   b->install_can_gates(1000.0);
   b->calc_current_at_start(0);
   std::vector<double> results;
   results.push_back(1.6345251904980005);
   results.push_back(1.6345251904980005);
   results.push_back(7.80493);
   results.push_back(41460.8);
   results.push_back(14315.6);
   results.push_back(1946.1);
   results.push_back(240.059);
   results.push_back(36.044762202708);
   double delta_t = 0.001;
   for (size_t i = 0; i < results.size(); i++) {
      double d = b->timestepping_of_gates_and_calc_current(delta_t * i, delta_t);
      BOOST_REQUIRE_CLOSE(d, results[i], SMALL);
      BOOST_MESSAGE("End test >>check_fluxes<<");
   }
}
#endif

BOOST_AUTO_TEST_SUITE_END();



