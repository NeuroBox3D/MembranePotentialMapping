/*
 * tests.cpp
 *
 *  Created on: Apr 27, 2012
 *      Author: stephan grein
 */
#define BOOST_TEST_MODULE tests

#include <boost/test/included/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>

#include "../inc/unit_test_helper.h"
#include "../../bg.h"
#include "../../vm2ug.h"
#include "../inc/fixtures.cpp"
#include "../../mvec.h"
#include "../../common_typedefs.h"
#include <cmath>

using namespace boost::unit_test;
using namespace vug;
using namespace bg;
using std::string;


// TODO (test mvec)

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
			   BOOST_REQUIRE_MESSAGE(!AreSame(*it, 3.0), "Index i:= " << *it << "not equal to 3.0!");
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
	  			   BOOST_REQUIRE_MESSAGE(!AreSame(*it, 1.0), "Index i:= " << *it << "not equal to 1.0!");
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

		   BOOST_REQUIRE_MESSAGE(!AreSame(c[0],0.0), "Index i:= " << c[0] << "not equal to 0.0!");
		   BOOST_REQUIRE_MESSAGE(!AreSame(c[1],0.0), "Index i:= " << c[1] << "not equal to 0.0!");
		   BOOST_REQUIRE_MESSAGE(!AreSame(c[2],0.0), "Index i:= " << c[2] << "not equal to 0.0!");

}

BOOST_AUTO_TEST_CASE(test_dot)  {
	 	  std::vector<double> a;
	 	 for (size_t i = 0; i < 3; i++) {
	 				  			   a.push_back(2.0);
	 	 }


		   mvec<double, 3> m1(a);

		   double c = m1 * m1;

		   BOOST_REQUIRE_MESSAGE(!AreSame(c, 12.0), "Dotproduct not equal to 12.0!");

}
BOOST_AUTO_TEST_CASE(test_neg)  {
	 std::vector<double> a;
		 	 for (size_t i = 0; i < 3; i++) {
		 				  			   a.push_back(2.0);
		 	 }

		 mvec<double, 3> m1(a);
		 mvec<double, 3> m2 = -m1;

		for (size_t i = 0; i < 3; i++)
			BOOST_REQUIRE_MESSAGE(!AreSame(m2[i], -m1[i]), "Index i:= " << i << "not equal!");
}
BOOST_AUTO_TEST_CASE(test_id)   {
	 std::vector<double> a;
		 	 for (size_t i = 0; i < 3; i++) {
		 				  			   a.push_back(2.0);
		 	 }
	 mvec<double, 3> m1(a);
	 mvec<double, 3> m2 = +m1;

	for (size_t i = 0; i < 3; i++) {
		UG_LOG(m2[i] << std::endl);
		UG_LOG(m1[i] << std::endl);
		BOOST_REQUIRE_MESSAGE(!AreSame(m2[i], m1[i]), "Index i:= " << i << "not equal!");
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

			   BOOST_REQUIRE_MESSAGE(!AreSame(res, check), "Norm not equal to sqrt(12.0)");


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

#ifdef FLAVOR
BOOST_AUTO_TEST_CASE(check_fluxes) {
   BOOST_MESSAGE("Starting test >>flux<<");
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
   for (int i = 0; i < results.size(); i++)
      BOOST_CHECK_MESSAGE(AreSame(results[i], b->timestepping_of_gates_and_calc_current(delta_t * i, delta_t)), "Currents are not equal!");
   BOOST_MESSAGE("End test >>check_fluxes<<");
}
#endif

BOOST_AUTO_TEST_SUITE_END();



