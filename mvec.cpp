/**
 * mvec.cpp
 *
 *  Created on: Apr 9, 2012
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

template <class T, size_t i> mvec<T, i>::mvec(size_t t) {
	this->reserve(t);
	temp.reserve(i);
}

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


/*
 * main method
 */
/*
int main() {
  std::vector<double> foofoo;
   foofoo.push_back(3.0);
   foofoo.push_back(3.0);
   foofoo.push_back(3.0);
   std::vector<double> foofoo2;
   foofoo2.push_back(5.0);
   foofoo2.push_back(5.0);
   foofoo2.push_back(5.0);

  mvec<double, 3> blabla();
	mvec<double, 3> foo(foofoo);
  mvec<double, 3> foo2(foofoo2);
	mvec<double, 3> a = foofoo;
	mvec<double, 3> b = foofoo2;

	mvec<double, 3> c = a + b;
	for (DITC mit3 = c.begin(); mit3 < c.end(); mit3++)
		std::cout << *mit3;
	std::cout << std::endl;

	//mvec<double, 3> c = a.add(b).add(b);

//	mvec<double, 3> d = c.dot(c);

	/*mvec<double, 3> c = foo + foo2;
  mvecd3 d = c + foo;
   typedef mvec<double, 3>::const_iterator MIT;
   std::cout << "START" << std::endl;
   for (MIT mit = d.begin(); mit < d.end(); mit++)
      std::cout << *mit;


   double d2 = c * foo;
   std::cout << "d2" << d2 << std::endl;
}
   */


