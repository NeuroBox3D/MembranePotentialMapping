/****************************************************************************/
/*                                                                          */
/* File:      my_commands.cc                                            	*/
/*                                                                          */
/* Purpose:   several own functions in order to use the the classes    */
/* 
 * Author:	  Markus M. Knodel                                              */
/*                Goethe Center for Scientific Computing             */
/*                University of Frankfurt                   */
/*                Kettenhofweg 139                           */
/*                60325 Frankfurt              */
/*                Germany              */
/*                email: markus.knodel@gcsc.uni-frankfurt.de    */
/*                                                                  */
/* History:   2009 begin            									   
 *                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

#include "my_commands.h"
#include "bg.h"
#include "solve.h"
#include <fstream>
#include <iostream>
#include "spannung.h"

using namespace std;

using namespace ug::membrane_potential_mapping::bg;
/////////////////////////////////////////////////////////////////////////////

double BG::timestepping_of_gates_and_calc_current(double time, double delta_t,
		double myVm, double Ca_i, double Ca_o) {
	double time_in_ms = time * 1000.;
	double delta_t_in_ms = delta_t * 1000.;

	Neumann_flux = solgat.calculate_current_expli_next_timestep(time_in_ms,
			delta_t_in_ms, myVm, Ca_i, Ca_o);
	return solgat.Current_current(time_in_ms, myVm, Ca_i, Ca_o);
}

double BG::timestepping_of_gates_and_calc_current(double time, double delta_t,
		double myVm) {
	double time_in_ms = time * 1000.;
	double delta_t_in_ms = delta_t * 1000.;

	Neumann_flux = solgat.calculate_current_expli_next_timestep(time_in_ms,
			delta_t_in_ms, myVm);
	return solgat.Current_current(time_in_ms, myVm);
}

/////////////////////////////////////////////////////////////////////////////

double BG::calc_current_at_start(double time) {

	Neumann_flux = solgat.initial_value_for_current();

	return Neumann_flux;
}

double BG::calc_current_at_start(double time, double basic_voltage, double myVm,
		double Ca_i, double Ca_o) {
	Neumann_flux = solgat.initial_value_for_current(basic_voltage, myVm, Ca_i,
			Ca_o);
	return Neumann_flux;
}
/////////////////////////////////////////////////////////////////////////////
double BG::dCa() const {
	return solgat.dCa();
}
double BG::dCa_dCa_i() const {
	return solgat.dCadCa_i();
}
double BG::dCa_dCa_o() const {
	return solgat.dCadCa_o();
}

