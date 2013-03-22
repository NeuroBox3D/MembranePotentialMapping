/*!
 * \file common_typedefs.h
 * \brief common typedefs for mvec usage (\see mvec.h)
 *
 * \author Stephan Grein
 * \date Created on: 07 April, 2012
 */

#ifndef __H__UG__MEMBRANE_POTENTIAL_MAPPING__COMMON_TYPEDEFS__
#define __H__UG__MEMBRANE_POTENTIAL_MAPPING__COMMON_TYPEDEFS__

/* standard includes */
#include <vector>
#include <cstddef>

/* begin namespace ug */
namespace ug {
	/* begin namespace mpm */
	namespace membrane_potential_mapping {
		/* forward declarations for mvec template */
		template <class T, size_t i> class mvec;

		/* commonly used typedefs for mvec */
		typedef mvec<double, 3> mvecd3;
		typedef mvec<double, 2> mvecd2;
		typedef mvec<double, 1> mvecd1;

		/* commonly used typedefs for iterators */
		typedef std::vector<double>::iterator DIT;
		typedef std::vector<double>::const_iterator DITC;

		/*!
		 * \enum NORM
		 * \brief currently only three norms supported
		 */
		enum NORM {
			INF=0, MANHATTAN, EUCLIDEAN
		};
	/* end namespace mpm */
	}
/* end namespace ug */
}
#endif /* __H__UG__MEMBRANE_POTENTIAL_MAPPING__COMMON_TYPEDEFS__ */
