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


using namespace bg;
using namespace vug;
int BG::ap_interval_duration_in_ms = 10;
Vm2uG<std::string> mapping = Vm2uG<std::string> ("", ".csv");
double x = 0;
double y = 0;
double z = 0;
std::string TimeStepBaseDir = std::string("");
/////////////////////////////////////////////////////////////////////////////

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
    solgat =  solve_gating( gate_can_m, gate_can_h,  conductivity, 135, 2, 1 );    
      
   inst_can_gates = true;

}

/////////////////////////////////////////////////////////////////////////////

// static objects, values to be choosen for each case seperately

BG::BG() {
   conductivity = 1000 ;
   Neumann_flux = 0;
   output_file_current = string("calculated_values.gdt" );
   inst_can_gates = false;
   x = 0;
   y = 0;
   z = 0;
}

double BG::get_Neumann_Flux() {
   return Neumann_flux;

}

void BG::setXYZ(double xx, double yy, double zz) {
   BG::x = xx;
   BG::y = yy;
   BG::z = zz;
}

void BG::setTimeStepBaseDir(std::string) {
   
}
/////////////////////////////////////////////////////////////////////////////