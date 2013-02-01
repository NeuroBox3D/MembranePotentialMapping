/****************************************************************************/
/*                                                                          */
/* File:      globals.cc                                            	*/
/*                                                                          */
/* Purpose:   global variables    */
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

bool BG::installed_cat_gates() {
	return inst_cat_gates;
}

bool BG::installed_cal_gates() {
	return inst_cal_gates;
}

bool BG::installed_can_gates() {
	return inst_can_gates;
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

double BG::get_Neumann_Flux_as_Concentration(const double dt, const double valency) const { // dt [s]
	//return 1e3 * (dt * Neumann_flux * 6.24e18) / (6.022e23 * valency); // Flux [mM] * [s] = [mM * s]
	return 1e-12 * (dt * Neumann_flux * 6.24e18) / (6.022e23 * valency); // [Mol/s]

}

double BG::F = 9.648 * 10000;
double BG::R = 8.314 * 1000;
double BG::T = 300;
double BG::NA = 6.022e23;
double BG::C = 6.24e18;
