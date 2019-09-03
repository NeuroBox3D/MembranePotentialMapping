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

#include "bg.h"
using namespace ug::membrane_potential_mapping::bg;


void BG::install_cal_gates_cfp(double permeability) {
	gating_parameter gate_cal_m(3, -36, 1.5);
	gating_parameter gate_cal_h(0, 0, 0);

	solgat = solve_gating(gate_cal_m, gate_cal_h, permeability, 2, 0, 2.0);

	inst_cal_gates = true;
	inst_cat_gates = false;
	inst_can_gates = false;

	this->permeability = permeability;
}

void BG::install_cat_gates_cfp(double permeability) {
	gating_parameter gate_cat_m(3, -36, 1.5);
	gating_parameter gate_cat_h(-5.2, -68, 10);

	solgat = solve_gating(gate_cat_m, gate_cat_h, permeability, 2, 1, 2.0);

	inst_cal_gates = false;
	inst_cat_gates = true;
	inst_can_gates = false;

	this->permeability = permeability;
}

void BG::install_can_gates_cfp(double permeability) {

	gating_parameter gate_can_m(3.4, -21, 1.5);
	gating_parameter gate_can_h(-2, -40, 75);
	solgat = solve_gating(gate_can_m, gate_can_h, permeability, 2, 1, 2.0);

	inst_cal_gates = false;
	inst_cat_gates = false;
	inst_can_gates = true;

	this->permeability = permeability;
}

void BG::install_cal_gates(double conductivity) {
	gating_parameter gate_cal_m(3, -36, 1.5);
	gating_parameter gate_cal_h(0, 0, 0);
	solgat = solve_gating(gate_cal_m, gate_cal_h, conductivity, 135, 2, 0);

	inst_cal_gates = true;
	inst_cat_gates = false;
	inst_can_gates = false;

	this->conductivity = conductivity;
}

void BG::install_cat_gates(double conductivity) {
	gating_parameter gate_cat_m(3, -36, 1.5);
	gating_parameter gate_cat_h(-5.2, -68, 10);
	solgat = solve_gating(gate_cat_m, gate_cat_h, conductivity, 135, 2, 1);

	inst_cal_gates = false;
	inst_cat_gates = true;
	inst_can_gates = false;

	this->conductivity = conductivity;
}

void BG::install_can_gates(double conductivity) {
	gating_parameter gate_can_m(3.4, -21, 1.5);
	gating_parameter gate_can_h(-2, -40, 75);

	solgat = solve_gating(gate_can_m, gate_can_h, conductivity, 135, 2, 1);

	inst_can_gates = true;
	inst_cal_gates = false;
	inst_cat_gates = false;

	this->conductivity = conductivity;
}

bool BG::installed_cat_gates() {
	return inst_cat_gates;
}

bool BG::installed_cal_gates() {
	return inst_cal_gates;
}

bool BG::installed_can_gates() {
	return inst_can_gates;
}

BG::BG() {
	conductivity = 0;
	permeability = 0;
	Neumann_flux = 0;
	inst_can_gates = false;
	inst_cal_gates = false;
	inst_cat_gates = false;
}

double BG::get_Neumann_Flux() {
	return Neumann_flux; // [C/s]
}

// actually should read 1e-18 (was 1e-12)
double BG::get_Neumann_Flux_as_Concentration(const double dt, const double valency) const { // dt [s]
	return 1e-18 * (dt * Neumann_flux * 6.24e18) / (6.022e23 * valency); // [Mol/s]

}

double BG::F = 9.648 * 10000;
double BG::R = 8.314 * 1000;
double BG::T = 300;
double BG::NA = 6.022e23;
double BG::C = 6.24e18;
