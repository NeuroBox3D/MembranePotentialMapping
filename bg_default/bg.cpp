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
//int BG::ap_interval_duration_in_ms = 10;
/////////////////////////////////////////////////////////////////////////////

void BG::install_cal_gates_cfp(double permeability) {
	gating_parameter gate_cal_m( 3, -36, 1.5 );
	gating_parameter gate_cal_h( 0, 0, 0 );

	solgat = solve_gating(gate_cal_m, gate_cal_h, permeability, 2, 0, 2.0);

	inst_cal_gates = true;
	inst_cat_gates = false;
	inst_can_gates = false;

	this->permeability = permeability;
}


void BG::install_cat_gates_cfp(double permeability) {

	gating_parameter gate_cat_m( 3, -36, 1.5 );
	gating_parameter gate_cat_h( -5.2, -68, 10 );
	solgat = solve_gating(gate_cat_m, gate_cat_h, permeability, 2, 1, 2.0);

	inst_cal_gates = false;
	inst_cat_gates = true;
	inst_can_gates = false;

	this->permeability = permeability;
}


void BG::install_can_gates_cfp(double permeability) {

	gating_parameter gate_can_m( 3.4, -21, 1.5 );
	gating_parameter gate_can_h( -2, - 40, 75  );
	solgat = solve_gating( gate_can_m, gate_can_h, permeability, 2, 1, 2.0 );

	inst_cal_gates = false;
	inst_cat_gates = false;
	inst_can_gates = true;

	this->permeability = permeability;
}


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
    solgat =  solve_gating( gate_can_m, gate_can_h, cond, 135, 2, 1 );    
      
   inst_can_gates = true;
   inst_cal_gates = false;
   inst_cat_gates = false;

   conductivity = cond;

}

/////////////////////////////////////////////////////////////////////////////

// static objects, values to be choosen for each case seperately

BG::BG() {
   conductivity = 1000;
   permeability = 0;
   Neumann_flux = 0;
  // output_file_current = string("calculated_values.gdt" );
   inst_can_gates = false;
   inst_cal_gates = false;
   inst_cat_gates = false;
}

double BG::get_Neumann_Flux() {
   return Neumann_flux;
}

double BG::get_Flux_As_Concentration(double delta_t, double valency) const {
	return (Neumann_flux * 1000 * 1000 * delta_t * BG::C)/ (1000 * BG::NA* valency);
}

double BG::F = 9.648 * 10000;
double BG::R = 8.314 * 1000;
double BG::T = 300;
double BG::NA = 6.022e23;
double BG::C = 6.24e18;

/////////////////////////////////////////////////////////////////////////////
