/**
 * \file common_typedefs.h
 * \brief common typedefs. TODO: move to namespace ug::membrane_potential_mapping
 *
 * \author Stephan Grein
 * \date Created on: 07 April, 2012
 */

#ifndef __H__COMMON_TYPEDEFS__
#define __H__COMMON_TYPEDEFS__

#include <vector>
#include <cstddef>

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

#endif /* __H__COMMON_TYPEDEFS__ */
