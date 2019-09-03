/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus M. Knodel
 * Creation date: 2009-08
 *
 * This file is part of NeuroBox, which is based on UG4.
 *
 * NeuroBox and UG4 are free software: You can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3
 * (as published by the Free Software Foundation) with the following additional
 * attribution requirements (according to LGPL/GPL v3 §7):
 *
 * (1) The following notice must be displayed in the appropriate legal notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 *
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 *
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating PDE based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * "Stepniewski, M., Breit, M., Hoffer, M. and Queisser, G.
 *   NeuroBox: computational mathematics in multiscale neuroscience.
 *   Computing and visualization in science (2019).
 * "Breit, M. et al. Anatomically detailed and large-scale simulations studying
 *   synapse loss and synchrony using NeuroBox. Front. Neuroanat. 10 (2016), 8"
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include "solve.h"
#include <cmath>
#include <assert.h>
#include "bg.h"
#include "voltage.h"

using namespace ug::membrane_potential_mapping::bg;

double solve_gating::I(double x, double y, double t, double myVm) {
	d_Ca = 0.0; // since no calcium is considered!
	return -g * pow(x, xp) * pow(y, yp) * (myVm - V_ret);
}

double solve_gating::I(double x, double y, double t, double myVm, double Ca_i, double Ca_o) {
// old derivatives
	/*	d_CadCa_o =
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
							- 1)); */

// new derivatives
	d_CadCa_i = -(247955322265625*std::pow(BG::F, 2.0)*myVm*permeability*std::pow(valency, 2.0))/(23929278055826824*BG::R*BG::T*valency*(1/std::exp((BG::F*myVm*valency)/(BG::R*BG::T)) - 1));
	d_CadCa_o = -(247955322265625*std::pow(BG::F, 2.0)*myVm*permeability*std::pow(valency, 2.0))/(23929278055826824*BG::R*BG::T*valency*(1/std::exp((BG::F*myVm*valency)/(BG::R*BG::T)) - 1));

	return pow(x, xp) * pow(y, yp) * permeability * valency * valency * BG::F
			* BG::F * myVm / BG::R / BG::T
			* (Ca_i - Ca_o * std::exp(-valency * BG::F * myVm / (BG::R * BG::T)))
			/ (1 - std::exp(-valency * BG::F * myVm / (BG::R * BG::T))); // was: *1e3;
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
