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

/* mpm includes */
#include "vm2ug.h"

// begin namespace ug 
namespace ug {
	// begin namespace mpm 
	namespace membrane_potential_mapping {
		// template declarations for Vm2uG 
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
#else
		template void Vm2uG<int>::buildTree(const int& timestep);
		template void Vm2uG<number>::buildTree(const number& timestep);
		template void Vm2uG<std::string>::buildTree(const std::string& timestep);
#endif

#ifndef MPMNEURON
		template uGPoint<int> Vm2uG<int>::vm_t(const int& timestep, number node[]);
		template uGPoint<number> Vm2uG<number>::vm_t(const number& timestep, number node[]);
		template uGPoint<std::string> Vm2uG<std::string>::vm_t(const std::string& timestep, number node[]);

		template std::vector<uGPoint<int> > Vm2uG<int>::vm_t_many_k(const int& timestep, number nodes[][DIM], int k);
		template std::vector<uGPoint<number> > Vm2uG<number>::vm_t_many_k(const number& timestep, number nodes[][DIM], int k);
#else
		template number Vm2uG<int>::vm_t(number node[]);
		template number Vm2uG<number>::vm_t(number node[]);
		template number Vm2uG<std::string>::vm_t(number node[]);
#endif

		template std::ostream& operator<<(std::ostream& output, const Vm2uG<int>& p);
		template std::ostream& operator<<(std::ostream& output, const Vm2uG<number>& p);
		template std::ostream& operator<<(std::ostream& output, const Vm2uG<std::string>& p);

		template uGPoint<int>::~uGPoint();
		template uGPoint<number>::~uGPoint();
		template uGPoint<std::string>::~uGPoint();

		template number uGPoint<int>::getVm();
		template number uGPoint<number>::getVm();
		template number uGPoint<std::string>::getVm();

		template number uGPoint<int>::getDist();
		template number uGPoint<number>::getDist();
		template number uGPoint<std::string>::getDist();

		template std::ostream& operator<<(std::ostream& output, const uGPoint<int>& p);
		template std::ostream& operator<<(std::ostream& output, const uGPoint<number>& p);
		template std::ostream& operator<<(std::ostream& output, const uGPoint<std::string>& p);

#ifndef MPMNEURON
		template number Vm2uG<int>::interp_bilin_vms(const int& timestep, number node[], number cutoff, int k);
		template number Vm2uG<number>::interp_bilin_vms(const number& timestep, number node[], number cutoff, int k);
		template number Vm2uG<std::string>::interp_bilin_vms(const std::string& timestep, number node[], number cutoff, int k);

		template number Vm2uG<int>::interp_lin_vms(const int& timestep, number node[], number cutoff, int k);
		template number Vm2uG<number>::interp_lin_vms(const number& timestep, number node[],  number cutoff, int k);
		template number Vm2uG<std::string>::interp_lin_vms(const std::string& timestep, number node[],  number cutoff, int k);
#endif
	// end namespace mpm 
	}
// end namespace ug 
}
