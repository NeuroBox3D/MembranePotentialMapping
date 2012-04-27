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

#ifdef __GXX_EXPERIMENTAL_CXX0X__
	#include <functional>
	  template <class T> function<bool(T)> isInRange(T low, T high) {
      return [low,high](T value) { return std::fabs(value - low) >= std::numeric_limits<T>::epsilon() \
                                       && std::fabs(value - high) <= std::numeric_limits<T>::epsilon(); }; }
#endif


bool AreSame(double a, double b) {
   return ! std::fabs(a - b) < std::numeric_limits<double>::epsilon();
}

#endif /* UNIT_TEST_HELPER_H_ */
