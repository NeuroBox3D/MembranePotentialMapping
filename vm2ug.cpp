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
#include <cstring>
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

using namespace std;

/* begin namespace ug */
namespace ug {
	/* begin namespace mpm */
	namespace membrane_potential_mapping {
		/* template declarations for Vm2uG */
		template Vm2uG<int>::~Vm2uG();
		template Vm2uG<double>::~Vm2uG();
		template Vm2uG<string>::~Vm2uG();

		template Vm2uG<int>::Vm2uG (string dataFileBaseName, string dataFileExt, const bool promise);
		template Vm2uG<double>::Vm2uG (string dataFileBaseName, string dataFileExt, const bool promise);
		template Vm2uG<string>::Vm2uG (string dataFileBaseName, string dataFileExt, const bool promise);

		template void Vm2uG<int>::buildTree(const int& timestep);
		template void Vm2uG<double>::buildTree(const double& timestep);
		template void Vm2uG<string>::buildTree(const string& timestep);

		template uGPoint<int> Vm2uG<int>::vm_t(const int& timestep, const double node[]);
		template uGPoint<double> Vm2uG<double>::vm_t(const double& timestep, const double node[]);
		template uGPoint<string> Vm2uG<string>::vm_t(const string& timestep, const double node[]);

		template std::vector<uGPoint<int> > Vm2uG<int>::vm_t_many_k(const int& timestep, const double nodes[][DIM], const int& k);
		template std::vector<uGPoint<double> > Vm2uG<double>::vm_t_many_k(const double& timestep, const double nodes[][DIM], const int& k);

		template ostream& operator<<(ostream& output, const Vm2uG<int>& p);
		template ostream& operator<<(ostream& output, const Vm2uG<double>& p);
		template ostream& operator<<(ostream& output, const Vm2uG<string>& p);

		template uGPoint<int>::~uGPoint();
		template uGPoint<double>::~uGPoint();
		template uGPoint<string>::~uGPoint();

		template double uGPoint<int>::getVm();
		template double uGPoint<double>::getVm();
		template double uGPoint<string>::getVm();

		template double uGPoint<int>::getDist();
		template double uGPoint<double>::getDist();
		template double uGPoint<string>::getDist();

		template ostream& operator<<(ostream& output, const uGPoint<int>& p);
		template ostream& operator<<(ostream& output, const uGPoint<double>& p);
		template ostream& operator<<(ostream& output, const uGPoint<string>& p);

		template const double Vm2uG<int>::interp_bilin_vms(const int& timestep, const double node[], const double cutoff, const int k);
		template const double Vm2uG<double>::interp_bilin_vms(const double& timestep, const double node[], const double cutoff, const int k);
		template const double Vm2uG<string>::interp_bilin_vms(const string& timestep, const double node[], const double cutoff, const int k);

		template const double Vm2uG<int>::interp_lin_vms(const int& timestep, const double node[], const double cutoff, const int k);
		template const double Vm2uG<double>::interp_lin_vms(const double& timestep, const double node[], const double cutoff, const int k);
		template const double Vm2uG<string>::interp_lin_vms(const string& timestep, const double node[], const double cutoff, const int k);
	/* end namespace mpm */
	}
/* end namespace ug */
}
