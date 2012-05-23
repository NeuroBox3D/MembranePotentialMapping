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
   return Neumann_flux;
}

/////////////////////////////////////////////////////////////////////////////
