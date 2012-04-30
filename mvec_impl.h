/*
 * mvec_impl.h
 *
 *  Created on: Apr 27, 2012
 *      Author: stephan grein
 */

#include "mvec.h"
#include "common_typedefs.h"
#include <boost/lexical_cast.hpp>
#include <numeric> // needed for std::inner_product
#include <common/log.h>
#include <iostream>

// default constructor
template <class T, size_t i> mvec<T, i>::mvec() {
	this->reserve(i);
}

// main constructor
template<class T, size_t i> mvec<T, i>::mvec(const std::vector<T>& init) {
   this->reserve(i);

   for (typename std::vector<T>::const_iterator cit = init.begin(); cit < init.end(); cit++)
		this->push_back(boost::lexical_cast<T>(*cit));
}

// default destructor
template<class T, size_t i> mvec<T, i>::~mvec() {
}

// norm
template <class T, size_t i> const double mvec<T, i>::norm(NORM norm) const {
     typedef typename std::vector<T>::const_iterator CIT;
     double ret = 0.0;

     switch(norm) {
     case INF:
    	 break;
     case MANHATTAN:
    	 break;
     case EUCLIDEAN:
    	 for (CIT cit = this->begin(); cit < this->end(); cit++)
    		 ret += boost::lexical_cast<double>((*cit) * (*cit));
    	 return std::sqrt(ret);
     default:
    	 for (CIT cit = this->begin(); cit < this->end(); cit++)
    		 ret += boost::lexical_cast<double>((*cit) * (*cit));
    	 return std::sqrt(ret);
     }
     return ret;
}

// operators
template<class T, size_t i> mvec<T, i> mvec<T, i>::add(const mvec<T, i>& rhs) const {
	const size_t lhs_s = this->size();
	const size_t rhs_s = rhs.size();

	mvec<T, i> ret = std::vector<T>(i);

	if (! (lhs_s == rhs_s) ) {
		UG_LOG("mvec::add size of mvec lhs (" << lhs_s << ")" << "and rhs (" << rhs_s << ") are not equal!");
		return ret;
	}
	 ret.clear();
	 typename mvec<T, i>::const_iterator it2 = rhs.begin();

	 for (typename std::vector<T>::const_iterator it = this->begin(); it < this->end() && it2 < rhs.end(); it++, it2++)
			 ret.push_back(boost::lexical_cast<T>((*it) + (*it2)));

	return ret;
}

template<class T, size_t i> mvec<T, i> mvec<T, i>::sub(const mvec<T, i>& rhs) const {
	const size_t lhs_s = this->size();
	const size_t rhs_s = rhs.size();

	mvec<T, i> ret = std::vector<T>(i);

	if (! (lhs_s == rhs_s) ) {
		UG_LOG("mvec::add size of mvec lhs (" << lhs_s << ")" << "and rhs (" << rhs_s << ") are not equal!");
		return ret;
	}
	 ret.clear();
	 typename mvec<T, i>::const_iterator it2 = rhs.begin();

	 for (typename std::vector<T>::const_iterator it = this->begin(); it < this->end() && it2 < rhs.end(); it++, it2++)
			 ret.push_back(boost::lexical_cast<T>((*it) - (*it2)));

	return ret;
}

template <class T, size_t i> const double mvec<T, i>::dot(const mvec<T, i>& rhs) const {
	return std::inner_product(this->begin(), this->end(), rhs.begin(), 0);
}


template <class T, size_t i> mvec<T, i> mvec<T, i>::vec(const mvec<T, i>& rhs) const {
	const size_t lhs_s = this->size();
	const size_t rhs_s = rhs.size();

	mvec<T, i> ret = std::vector<T>(i);

	if (! (lhs_s == rhs_s) ) {
		UG_LOG("mvec::vec size of mvec lhs (" << lhs_s << ")" << "and rhs (" << rhs_s << ") are not equal!");
		if ( (lhs_s != 3) || (rhs_s != 3) )
			UG_LOG("mvec::vec currently only works for dim := 3")
		return ret;
	}

	std::vector<T> temp;
	ret.clear();
	temp.clear();

	for (typename std::vector<T>::const_iterator it = this->begin(); it < this->end(); it++)
		temp.push_back(boost::lexical_cast<T>(*it));

    ret.push_back(boost::lexical_cast<T>(temp[1] * rhs[2] - temp[2] * rhs[1]));
    ret.push_back(boost::lexical_cast<T>(temp[2] * rhs[0] - temp[0] * rhs[2]));
	ret.push_back(boost::lexical_cast<T>(temp[0] * rhs[1] - temp[1] * rhs[0]));

	return ret;
}

template <class T, size_t i> mvec<T, i> mvec<T, i>::neg() const {

	mvec<T, i> ret = std::vector<T>(i);
	ret.clear();

	for (typename std::vector<T>::const_iterator it = this->begin(); it < this->end(); it++)
				 ret.push_back(boost::lexical_cast<T>(-(*it)));
	return ret;
}

template <class T, size_t i> mvec<T, i> mvec<T, i>::id() const {
	return *this;
}

template <class T, size_t i> mvec<T, i> mvec<T, i>::operator=(const mvec<T, i>& rhs) {
	this->clear();

	for (typename std::vector<T>::const_iterator cit = rhs.begin(); cit < rhs.end(); cit++)
			this->push_back(boost::lexical_cast<T>(*cit));

	return *this;
}

template <class T, size_t i > mvec<T, i>& mvec<T, i>::operator+=(const mvec<T, i>& rhs) {
	    const size_t lhs_s = this->size();
		const size_t rhs_s = rhs.size();

		if (! (lhs_s == rhs_s) ) {
			UG_LOG("mvec::add size of mvec lhs (" << lhs_s << ")" << "and rhs (" << rhs_s << ") are not equal!");
			return *this;
		}

		for (size_t t = 0; t < lhs_s; t++)
			this->data()[t] = boost::lexical_cast<T>(this->data()[t] + rhs[t]);

		return *this;

}

template <class T, size_t i > mvec<T, i>& mvec<T, i>::operator-=(const mvec<T, i>& rhs) {
	    const size_t lhs_s = this->size();
		const size_t rhs_s = rhs.size();

		if (! (lhs_s == rhs_s) ) {
			UG_LOG("mvec::add size of mvec lhs (" << lhs_s << ")" << "and rhs (" << rhs_s << ") are not equal!");
			return *this;
		}

		for (size_t t = 0; t < lhs_s; t++)
			this->data()[t] = boost::lexical_cast<T>(this->data()[t] - rhs[t]);

		return *this;

}

template <class T, size_t i> mvec<T, i>& mvec<T, i>::operator%=(const mvec<T, i>& rhs) {
	const size_t lhs_s = this->size();
	const size_t rhs_s = rhs.size();

	if (! (lhs_s == rhs_s) ) {
		UG_LOG("mvec::vec size of mvec lhs (" << lhs_s << ")" << "and rhs (" << rhs_s << ") are not equal!");
		if ( (lhs_s != 3) || (rhs_s != 3) )
			UG_LOG("mvec::vec currently only works for dim := 3")
		return *this;
	}

	const T first = this->data()[0];
	const T second = this->data()[1];

	this->data()[0] = boost::lexical_cast<T>(this->data()[1] * rhs[2] - this->data()[2] * rhs[1]);
    this->data()[1] = boost::lexical_cast<T>(this->data()[2] * rhs[0] - first * rhs[2]);
	this->data()[2] = boost::lexical_cast<T>(first * rhs[1] - second * rhs[0]);

	return *this;
}
