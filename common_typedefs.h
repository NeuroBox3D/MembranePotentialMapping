/**
 * author: stephanmg
 * date: 04m07d12y
 * file: common_typedefs.h
 */

#ifndef _COMMON_TYPEDEFS_H_
#define _COMMON_TYPEDEFS_H_

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

#endif /* _COMMON_TYPEDEFS_H_ */
