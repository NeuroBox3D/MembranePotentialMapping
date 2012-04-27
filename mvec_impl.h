/*
 * mvec_impl.h
 *
 *  Created on: Apr 27, 2012
 *      Author: stephan
 */

#include "mvec.h"
#include "common_typedefs.h"
#include <boost/lexical_cast.hpp>
#include <numeric> // needed for std::inner_product
#include <common/log.h>
#include <iostream>

/**
 * default constructors
 */
template <class T, size_t i> mvec<T, i>::mvec() {
	this->reserve(i);
	temp.reserve(i);
}

/*template <class T, size_t i> mvec<T, i>::mvec(size_t t) {
	this->reserve(t);
	temp.reserve(i);
}*/

/**
 * main constructor
 */
template<class T, size_t i> mvec<T, i>::mvec(const std::vector<T>& init) {
   this->reserve(i);
   temp.reserve(i);

	for (typename std::vector<T>::const_iterator cit = init.begin(); cit < init.end(); cit++)
		this->push_back(boost::lexical_cast<T>(*cit));
	// call to this->assign(init); is equivlaent
}
/**
 * default destructor
*/
template<class T, size_t i> mvec<T, i>::~mvec() {
}

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

// TODO: implement vec product for higher dimensions: dim >= 4.
template <class T, size_t i> mvec<T, i> mvec<T, i>::vec(const mvec<T, i>& rhs) const {
	const size_t lhs_s = this->size();
	const size_t rhs_s = rhs.size();

	mvec<T, i> ret = std::vector<T>(i);

	if (! (lhs_s == rhs_s) ) {
		UG_LOG("mvec::vec size of mvec lhs (" << lhs_s << ")" << "and rhs (" << rhs_s << ") are not equal!");
		if ( (lhs_s != 3) || (rhs_s != 3) )
			UG_LOG("mvec::vec currently only works for dimension: 1 <= dim <= 3")
		return ret;
	}

	std::vector<T> temp2;
	ret.clear();
	temp2.clear();

	for (typename std::vector<T>::const_iterator it = this->begin(); it < this->end(); it++)
		temp2.push_back(boost::lexical_cast<T>(*it));

    ret.push_back(boost::lexical_cast<T>(temp2[1] * rhs[2] - temp2[2] * rhs[1]));
    ret.push_back(boost::lexical_cast<T>(temp2[2] * rhs[0] - temp2[0] * rhs[2]));
	ret.push_back(boost::lexical_cast<T>(temp2[0] * rhs[1] - temp2[1] * rhs[0]));

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
	//mvec<T, i> ret = std::vector<T>(i);
	return *this;
}

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

