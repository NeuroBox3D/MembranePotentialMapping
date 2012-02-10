/* Includes {{{ */
#include "vm2ug.h"
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
/* }}} */

/* Namespaces {{{ */
using namespace std;
/* }}} */

/* Template declarations for vug {{{ */
namespace vug {

#include "vm2ug_impl.h"

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

template std::vector<uGPoint<int> > Vm2uG<int>::vm_t_many_k(const int& timestep, const double nodes[][QDIM], const int& k);
template std::vector<uGPoint<double> > Vm2uG<double>::vm_t_many_k(const double& timestep, const double nodes[][QDIM], const int& k);

template ostream& operator<<(ostream& output, const Vm2uG<int>& p);
template ostream& operator<<(ostream& output, const Vm2uG<double>& p);
template ostream& operator<<(ostream& output, const Vm2uG<string>& p);

template uGPoint<int>::~uGPoint();
template uGPoint<double>::~uGPoint();
template uGPoint<string>::~uGPoint();

template DOUBLE uGPoint<int>::getVm();
template DOUBLE uGPoint<double>::getVm();
template DOUBLE uGPoint<string>::getVm();

template DOUBLE uGPoint<int>::getDist();
template DOUBLE uGPoint<double>::getDist();
template DOUBLE uGPoint<string>::getDist();

template ostream& operator<<(ostream& output, const uGPoint<int>& p);
template ostream& operator<<(ostream& output, const uGPoint<double>& p);
template ostream& operator<<(ostream& output, const uGPoint<string>& p);

}
/* }}} */