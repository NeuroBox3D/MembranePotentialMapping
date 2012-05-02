/*
 * mvec.h
 *
 *  Created on: Apr 07, 2012
 *      Author: stephan grein
 */

#ifndef _MVEC_H_
#define _MVEC_H_
#include <vector>
#include <cstddef>
#include <cmath>
#include <numeric>
#include "common_typedefs.h"

template <class T = double, size_t i = 3> class mvec : public std::vector<T> {

		public:
			// constructors
			mvec();

			mvec(const std::vector<T>& init);

			// default destructor
			~mvec();

			// determinant
			static const double det(const std::vector<mvec<T, i> >& mvecs);

			// norm
			const double norm(NORM) const;

			// operators
			inline mvec<T, i> operator+(const mvec<T, i>& rhs) const { return this->add(rhs); }

			inline mvec<T, i> operator-(const mvec<T, i>& rhs) const { return this->sub(rhs); }

			inline const double operator*(const mvec<T, i>& rhs) const { return this->dot(rhs); }

			inline mvec<T, i> operator%(const mvec<T, i>& rhs) const { return this->vec(rhs); }

		    inline mvec<T, i> operator-() const { return this->neg(); }

		    inline mvec<T, i> operator+() const { return this->id(); }

		    inline mvec<T, i> operator=(const mvec<T, i>& rhs);

	        inline mvec<T, i>& operator+=(const mvec<T, i>& rhs);

	        inline mvec<T, i>& operator-=(const mvec<T, i>& rhs);

	        inline mvec<T, i>& operator%=(const mvec<T, i>& rhs);

		private:
			// add
			mvec<T, i> add(const mvec<T, i>& rhs) const;

			// sub
			mvec<T, i> sub(const mvec<T, i>& rhs) const;

			// dot product
			inline const double dot(const mvec<T, i>& rhs) const;

			// cross product
			mvec<T, i> vec(const mvec<T, i>& rhs) const;

			// negate
			mvec<T, i> neg() const;

			// identity
			mvec<T, i> id() const;

};

#include "mvec_impl.h"

#endif /* _MVEC_H_ */


