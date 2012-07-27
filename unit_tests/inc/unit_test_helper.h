/**
 * \file unit_test_helper.h
 * \brief helper functions for unit tests
 *
 * \date created on Apr 27, 2012
 * \author Stephan Grein
 */

#ifndef __H__UG__MEMBRANE_POTENTIAL_MAPPING__UNIT_TEST_HELPER__
#define __H__UG__MEMBRANE_POTENTIAL_MAPPING__UNIT_TEST_HELPER__

// includes
#include <limits>
#include <cmath>


// namespace ug (ug)
namespace ug {
	// namespace mpm (membrane_potential_mapping)
	namespace membrane_potential_mapping {
		/**
		 * \brief checks if a value v is in a given range, i. e. v in [low, high].
		 *
		 * \tparam T the type for the values, e. g. int, double, etc.
		 * \param[in] low: the lower bound
		 * \param[in] high: the upper bound
		 *
		 * \return \c bool indicates if the value is in the range
		 */
		#ifdef __GXX_EXPERIMENTAL_CXX0X__
			#include <functional>
				template <class T> function<bool(T)> isInRange(T low, T high) {
					return [low,high](T value) { return std::fabs(value - low) >= std::numeric_limits<T>::epsilon() && std::fabs(value - high) <= std::numeric_limits<T>::epsilon(); };
				}
		#endif

		// initialization of static member
		static const double SMALL = 10e-5;
		static const double VERY_SMALL = 10e-8;
	// end namespace mpm (membrane_potential_mapping)
	}
// end namespace ug (ug)
}
#endif /* __H__UG__MEMBRANE_POTENTIAL_MAPPING__UNIT_TEST_HELPER__ */
