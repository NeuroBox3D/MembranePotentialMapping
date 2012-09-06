/****************************************************************************/
/*                                                                          */
/* File:      solve.cc                                          	*/
/*                                                                          */
/* Purpose:   solves the odes of the borg-graham model     */
/*                                                                          */
/* Author:	  Markus M. Knodel                                              */
/*                Goethe Center for Scientific Computing             */
/*                University of Frankfurt                   */
/*                Kettenhofweg 139                           */
/*                60325 Frankfurt              */
/*                Germany              */
/*                email: markus.knodel@gcsc.uni-frankfurt.de    */
/*																			*/
/* History:   2009 begin            									    */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

#include "solve.h"
#include <cmath>
#include "spannung.h"
#include <assert.h>
#include "bg.h"

using namespace ug::membrane_potential_mapping::bg;

double solve_gating::I(double x, double y, double t, double myVm) {
	d_Ca = 0.0; // since no calcium is considered!
	return -g * pow(x, xp) * pow(y, yp) * (myVm - V_ret);
}

double solve_gating::I(double x, double y, double t, double myVm, double Ca_i,
		double Ca_o) {
	d_CadCa_o =
			-(pow(BG::F, 2.0) * myVm * pow(x, xp) * pow(y, yp) * permeability
					* pow(valency, 2.0))
					/ (BG::R * BG::T
							* (std::exp(
									(BG::F * myVm * valency) / (BG::R * BG::T))
									- 1));
	d_CadCa_i = -(pow(BG::F, 2.0) * myVm * pow(x, xp) * pow(y, yp)
			* permeability * pow(valency, 2.0))
			/ (BG::R * BG::T
					* (1 / std::exp((BG::F * myVm * valency) / (BG::R * BG::T))
							- 1));

	return pow(x, xp) * pow(y, yp) * permeability * valency * valency * BG::F
			* BG::F * myVm / BG::R / BG::T
			* (Ca_i - Ca_o * std::exp(-valency * BG::F * myVm / (BG::R * BG::T)))
			/ (1 - std::exp(-valency * BG::F * myVm / (BG::R * BG::T))) * 1e3;
}

/////////////////////////////////////////////////////////////////////////////

double solve_gating::initial_value_for_current(double basic_voltage,
		double myVm) // t = 0
		{
	current_x_0 = gate_x.x_infty(basic_voltage);
	current_y_0 = gate_y.x_infty(basic_voltage);

	return I(current_x_0, current_y_0, 0., myVm);
}

double solve_gating::initial_value_for_current(double basic_voltage,
		double myVm, double Ca_i, double Ca_o) {
	current_x_0 = gate_x.x_infty(basic_voltage);
	current_y_0 = gate_y.x_infty(basic_voltage);

	return I(current_x_0, current_y_0, 0., myVm, Ca_i, Ca_o);
}

double solve_gating::calculate_current_expli_next_timestep(double time,
		double delta_t, double myVm) {

	current_x_0 = gate_x.x_1(time, current_x_0, delta_t, myVm);
	current_y_0 = gate_y.x_1(time, current_y_0, delta_t, myVm);

	return I(current_x_0, current_y_0, time, myVm);
}

double solve_gating::calculate_current_expli_next_timestep(double time,
		double delta_t, double myVm, double Ca_i, double Ca_o) {

	current_x_0 = gate_x.x_1(time, current_x_0, delta_t, myVm);
	current_y_0 = gate_y.x_1(time, current_y_0, delta_t, myVm);

	return I(current_x_0, current_y_0, time, myVm, Ca_i, Ca_o);
}

/////////////////////////////////////////////////////////////////////////////
