/**
 * \file common_typedefs.h
 * \brief common typedefs for mvec usage (\see mvec.h)
 *
 * \author Stephan Grein
 * \date Created on: 07 April, 2012
 */

#ifndef __H__UG__MEMBRANE_POTENTIAL_MAPPING__COMMON_TYPEDEFS__
#define __H__UG__MEMBRANE_POTENTIAL_MAPPING__COMMON_TYPEDEFS__

// includes
#include <vector>
#include <cstddef>


// start namespace ug (ug)
namespace ug {
	// start namespace membrane_potential_mapping (mpm)
	namespace membrane_potential_mapping {
		// forward declarations
		template <class T, size_t i> class mvec;

		// common typedefs
		typedef mvec<double, 3> mvecd3;
		typedef mvec<double, 2> mvecd2;
		typedef mvec<double, 1> mvecd1;

		typedef std::vector<double>::const_iterator DITC;
		typedef std::vector<double>::iterator DIT;

		enum NORM {
			INF=0, MANHATTAN, EUCLIDEAN
		};
	// end namespace membrane_potential_mapping (mpm)
	}
// end namespace ug (ug)
}
#endif /* __H__UG__MEMBRANE_POTENTIAL_MAPPING__COMMON_TYPEDEFS__ */
