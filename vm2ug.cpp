/*!
 * \file vm2ug.cpp
 * \brief forward declarations for vm2ug
 * \see vm2ug.h
 *
 * \date Created on July, 2011
 * \author Stephan Grein
 */

/* standard includes */
#include <cstdlib>
#include <string>
#include <cstdio>
#include <sstream>
#include <fstream>
#include <locale>
#include <cmath>
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

// mpm includes
#include "vm2ug.h"

namespace ug {
namespace membrane_potential_mapping {


sPoint::sPoint(const std::vector<number>& coordinates, number dist, number Vm, int index, long timestep) {
	if (coordinates.size() != 0) {
		try {
			this->coordinates = coordinates;
			this->dist = dist;
			this->Vm = Vm;
			this->index = index;
			this->timestep = timestep; // was: timestep
		} catch (const char* const str) {
			std::cout << "Initialization failed: " << str << std::endl;
		}
	}
	else {
		std::cout << "You created an EMPTY point (no coords). Check?!" << std::endl;
	}
}

sPoint::sPoint() { }

sPoint::~sPoint() { }

number sPoint::getVm() const {
	return this->Vm;
}

number sPoint::getDist() const {
	return this->dist;
}

number sPoint::getIndex() const {
	return this->index;
}

long sPoint::getTimestep() const {
	return this->timestep;
}

std::vector<number> sPoint::getCoordinates() const {
	return this->coordinates;
}

std::ostream& operator<<(std::ostream& output, const sPoint& p) {
	output << "{hoc} (";
	unsigned int i;
	for (i=0; i < p.coordinates.size()-1; i++) output << p.coordinates[i] << ",";
	output << p.coordinates[i] << ")" << std::endl;
	output << "\tVm: " << p.Vm << " (mV)" << std::endl;
	output << "\tDist to query Point {uG}: " << p.dist << std::endl;
	output << "\tIndex {in timestep (Hashcode=" << p.timestep << ") file}: " << p.index << std::endl;
	return output;
}

/* }}} */

/* uGPoint {{{ */
uGPoint::uGPoint(const std::vector<number>& coordinates, const std::vector<sPoint>& nearestNeighbors) {
	if (coordinates.size() != 0) {
		try {
			this->coordinates = coordinates;
			this->nearestNeighbors = nearestNeighbors;
		} catch (const char* const str) {
			std::cout << "Initialization failed: " << str << std::endl;
		}
	}
	else {
		std::cout << "You created an EMPTY point (no coords). Check?!" << std::endl;
	}
}

uGPoint::uGPoint() { }

uGPoint::~uGPoint() { }

number uGPoint::getVm() {
	return nearestNeighbors[0].getVm();
}

number uGPoint::getDist() {
	return nearestNeighbors[0].getDist();
}


std::vector<number> uGPoint::getCoordinates() const {
	return this->coordinates;
}

std::vector<sPoint> uGPoint::getNearestNeighbors() const {
	return this->nearestNeighbors;
}

std::ostream& operator<<(std::ostream& output, const uGPoint& p) {
	output << "----------------------------------------------------------------------------------------" << std::endl;
	output << "query Point {uG} (";
	unsigned int i;
	for (i=0; i < p.coordinates.size()-1; i++) output << p.coordinates[i] << ", ";
	output << p.coordinates[i] << ")" << std::endl;
	output << "----------------------------------------------------------------------------------------" << std::endl;
	for (i=0; i < p.nearestNeighbors.size(); i++) output << "nearest neighbor #" << i << ": " << p.nearestNeighbors[i] << std::endl;
	output << "****************************************************************************************" << std::endl;
	return output;
}



/* not needed
template class Vm2uG<int>;
template class Vm2uG<number>;
template class Vm2uG<std::string>;

template std::ostream& operator<<(std::ostream& output, const Vm2uG<int>& p);
template std::ostream& operator<<(std::ostream& output, const Vm2uG<number>& p);
template std::ostream& operator<<(std::ostream& output, const Vm2uG<std::string>& p);
*/

#if 0
// template declarations for Vm2uG -- not needed
template Vm2uG<int>::~Vm2uG();
template Vm2uG<number>::~Vm2uG();
template Vm2uG<std::string>::~Vm2uG();

#ifdef MPMNEURON
template Vm2uG<int>::Vm2uG(SmartPtr<Transformator>);
template Vm2uG<number>::Vm2uG(SmartPtr<Transformator>);
template Vm2uG<std::string>::Vm2uG(SmartPtr<Transformator>);
#endif

template Vm2uG<int>::Vm2uG (const std::string& dataFileBaseName, const std::string& dataFileExt, bool promise);
template Vm2uG<number>::Vm2uG (const std::string& dataFileBaseName, const std::string& dataFileExt, bool promise);
template Vm2uG<std::string>::Vm2uG (const std::string& dataFileBaseName, const std::string& dataFileExt, bool promise);

#ifdef MPMNEURON
template void Vm2uG<int>::buildTree();
template void Vm2uG<number>::buildTree();
template void Vm2uG<std::string>::buildTree();
#endif

template void Vm2uG<int>::buildTree(const int& timestep);
template void Vm2uG<number>::buildTree(const number& timestep);
template void Vm2uG<std::string>::buildTree(const std::string& timestep);

template uGPoint Vm2uG<int>::vm_t(const int& timestep, number node[]);
template uGPoint Vm2uG<number>::vm_t(const number& timestep, number node[]);
template uGPoint Vm2uG<std::string>::vm_t(const std::string& timestep, number node[]);

template std::vector<uGPoint> Vm2uG<int>::vm_t_many_k(const int& timestep, number nodes[][DIM], int k);
template std::vector<uGPoint> Vm2uG<number>::vm_t_many_k(const number& timestep, number nodes[][DIM], int k);

#ifdef MPMNEURON
template number Vm2uG<int>::vm_t(number node[]);
template number Vm2uG<number>::vm_t(number node[]);
template number Vm2uG<std::string>::vm_t(number node[]);
#endif

template number Vm2uG<int>::interp_bilin_vms(const int& timestep, number node[], number cutoff, int k);
template number Vm2uG<number>::interp_bilin_vms(const number& timestep, number node[], number cutoff, int k);
template number Vm2uG<std::string>::interp_bilin_vms(const std::string& timestep, number node[], number cutoff, int k);

template number Vm2uG<int>::interp_lin_vms(const int& timestep, number node[], number cutoff, int k);
template number Vm2uG<number>::interp_lin_vms(const number& timestep, number node[],  number cutoff, int k);
template number Vm2uG<std::string>::interp_lin_vms(const std::string& timestep, number node[],  number cutoff, int k);
#endif

} // end namespace mpm
} // end namespace ug
