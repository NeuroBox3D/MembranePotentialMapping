/*
 * unit_test_helper.h
 *
 * author stephanmg  
 * date 01m12d12y
 *
 */

#ifndef UNIT_TEST_HELPER_H_
#define UNIT_TEST_HELPER_H_

#include <limits>
#include <cmath>

bool AreSame(double a, double b) {
   return ! std::fabs(a - b) < std::numeric_limits<double>::epsilon();
}

#endif /* UNIT_TEST_HELPER_H_ */
