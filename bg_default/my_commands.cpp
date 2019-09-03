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

#include "my_commands.h"
#include "bg.h"
#include "solve.h"
#include <fstream>
#include <iostream>

#include "voltage.h"

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
