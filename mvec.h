/*!
 * \file mvec.h
 * \brief vector functionalities for linear and bilinear interpolation of membrane potentials.
 * \addtogroup mpm_plugin
 *
 * \author: Stephan Grein
 * \date Created on: Apr 07, 2012
 *
 * Note: Migrate to UG's vector class
 */
#ifndef __H__UG__MEMBRANE_POTENTIAL_MAPPING__MVEC__
#define __H__UG__MEMBRANE_POTENTIAL_MAPPING__MVEC__

#include <vector>
#include <cstddef>
#include <cmath>
#include <numeric>

#include "common_typedefs.h"

// begin namespace ug
namespace ug {
	// begin namespace mpm
	namespace membrane_potential_mapping {
			/*!
			 * \brief mvec class for vector calculations
			 * \addtogroup mpm_plugin
			 *
			 * \param[in] T the precision
			 * \param[in] i the vector's length
			 */
			template <class T = number, size_t i = 3> class mvec : public std::vector<T> {
					public:
					/*!
					 * \brief default constructor
					 */
					mvec();

					/*!
					 * \brief main constructor
					 */
					mvec(const std::vector<T>& init);

					/*!
					 * \brief default destructor
					 */
					~mvec();

					/*!
					 * \brief calculates the determinant for a list of m n-dimensional vectors
					 *
					 * \param[in] mvecs the list of m n-dimensional vectors
					 *
					 * \return \c number determinant
					 */
					static number det(const std::vector<mvec<T, i> >& mvecs);

					/*!
					 * \brief calculates the norm of an n-dimensional vector
					 *
					 * \return \c number norm of n-dimensional vector
					 */
					number norm(NORM) const;

					// operators follow below
					inline mvec<T, i> operator+(const mvec<T, i>& rhs) const { return this->add(rhs); }
					inline mvec<T, i> operator-(const mvec<T, i>& rhs) const { return this->sub(rhs); }
					inline const number operator*(const mvec<T, i>& rhs) const { return this->dot(rhs); }
					inline mvec<T, i> operator%(const mvec<T, i>& rhs) const { return this->vec(rhs); }
					inline mvec<T, i> operator-() const { return this->neg(); }
					inline mvec<T, i> operator+() const { return this->id(); }
					inline mvec<T, i> operator=(const mvec<T, i>& rhs);
					inline mvec<T, i>& operator+=(const mvec<T, i>& rhs);
					inline mvec<T, i>& operator-=(const mvec<T, i>& rhs);
					inline mvec<T, i>& operator%=(const mvec<T, i>& rhs);

				private:
					/*!
					 * \brief adds to vectors component-wise
					 *
					 * \param[in] rhs the right hand site vector
					 *
					 * \return \c mvec<T, i> the sum vector
					 */
					mvec<T, i> add(const mvec<T, i>& rhs) const;

					/*!
					 * \brief subs to vectors component-wise
					 *
					 * \param[in] rhs the right hand site vector
					 *
					 * \return \c mvec<T, i> the diff vector
					 */
					mvec<T, i> sub(const mvec<T, i>& rhs) const;

					/*!
					 * \brief calculates the dot product of two vectors
					 *
					 * \param[in] rhs the right hand site vector
					 *
					 * \return \c number the dot product
					 */
					inline number dot(const mvec<T, i>& rhs) const;

					/*!
					 * \brief calculates the cross product of two vectors
					 *
					 * \param[in] rhs the right hand site vector
					 *
					 * \return mvec<T, i> the cross product
					 */
					mvec<T, i> vec(const mvec<T, i>& rhs) const;

					/*!
					 * \brief negates a given vector component-wise
					 *
					 * \return mvec<T, i> the negated vector
					 */
					mvec<T, i> neg() const;

					/*!
					 * \brief identity
					 *
					 * \return mvec<T, i> the unmodified vector
					 */
					mvec<T, i> id() const;
		};
	// end namespace mpm
	}
// end namespace ug
}

// include implementation of mvec
#include "mvec_impl.h"

#endif // __H__UG__MEMBRANE_POTENTIAL_MAPPING__MVEC__
