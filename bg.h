#ifndef __bg_h__
#define __bg_h__

#include <stddef.h>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <string>
#include "solve.h"
#include <fstream>

using namespace std;
namespace bg {

class BG
{
public:

 std::vector<solve_gating> solgat;

 string output_file_current;
 void install_can_gates(double conductivity = 1000);
 void install_cat_gates(double conductivity = 0.3);
 void install_cal_gates(double conductivity = 0.3*0.6);

 double Neumann_flux;

 BG();

double timestepping_of_gates_and_calc_current( double time, double delta_t, double Vm );

double calc_current_at_start( double Vm = -75);

double getCurrent(double delta_t, double Vm);

double getNeumannFlux(double delta_t, double Vm);

double initVm;

};
	 
}
#endif /* __bg_h__ */ 
