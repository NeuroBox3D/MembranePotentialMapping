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
int BG::ap_interval_duration_in_ms = 10;
/////////////////////////////////////////////////////////////////////////////

bool BG::installed_cal_gates() {
   return inst_cal_gates;
}  


void BG::install_cal_gates(double cond) {
   gating_parameter gate_cal_m(3, -36, 1.5);
   gating_parameter gate_cal_h(0, 0, 0);
   solgat = solve_gating(gate_cal_m, gate_cal_h, cond, 135, 2, 0);

   inst_cal_gates = true;
   inst_cat_gates = false;
   inst_can_gates = false;

   conductivity = cond;
}


bool BG::installed_cat_gates() {
   return inst_cat_gates;
}

void BG::install_cat_gates(double cond) {
   gating_parameter gate_cat_m(3, -36, 1.5);
   gating_parameter gate_cat_h(-5.2, -68, 10);
   solgat = solve_gating(gate_cat_m, gate_cat_h, cond, 135, 2, 1);

   inst_cal_gates = false;
   inst_cat_gates = true;
   inst_can_gates = false;

   conductivity = cond;
}


bool BG::installed_can_gates() {
   return inst_can_gates;
}
void BG::install_can_gates(double cond)
{

 // Ca N type m
	gating_parameter gate_can_m( 3.4, -21, 1.5 );
    // Ca N type h
	gating_parameter gate_can_h( -2, - 40, 75  );
    //solgat =  solve_gating( gate_can_m, gate_can_h, 0.06*1, 135, 2, 1 );    
    solgat =  solve_gating( gate_can_m, gate_can_h,  cond, 135, 2, 1 );    

   inst_cal_gates = false;
   inst_cat_gates = false;  
   inst_can_gates = true;

   conductivity = cond;

}

/////////////////////////////////////////////////////////////////////////////

// static objects, values to be choosen for each case seperately

BG::BG() {
   conductivity = 1000 ;
   Neumann_flux = 0;
   output_file_current = string("calculated_values.gdt" );
   inst_can_gates = false;
}

double BG::get_Neumann_Flux() {
   return Neumann_flux; // Mol
}


double BG::get_Neumann_Flux_as_Concentration(const double dt, const double valency) const { // dt [s]
	return 1e-18 * (Neumann_flux * 6.24e18) / (6.022e23 * valency); // [Mol/s]

}

/////////////////////////////////////////////////////////////////////////////
