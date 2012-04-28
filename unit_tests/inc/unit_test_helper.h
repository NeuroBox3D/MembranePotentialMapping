/*
 * unit_test_helper.h
 *
 * author stephanmg  
 * date 01m12d12y
 *
 */

#ifndef _UNIT_TEST_HELPER_H_
#define _UNIT_TEST_HELPER_H_

#include <limits>
#include <cmath>

#ifdef __GXX_EXPERIMENTAL_CXX0X__
	#include <functional>
	  template <class T> function<bool(T)> isInRange(T low, T high) {
      return [low,high](T value) { return std::fabs(value - low) >= std::numeric_limits<T>::epsilon() \
                                       && std::fabs(value - high) <= std::numeric_limits<T>::epsilon(); }; }
#endif

bool SameDoubles(const double& a, const double& b) {
   return std::fabs(a - b) < std::numeric_limits<double>::epsilon();
}

static const double SMALL = 0.0001;
static const double VERY_SMALL = 0.0000001;

#endif /* _UNIT_TEST_HELPER_H_ */
