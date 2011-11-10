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
/////////////////////////////////////////////////////////////////////////////
using namespace bg;

void BG::install_can_gates(double conductivity)
{

 // Ca N type m
	gating_parameter gate_can_m( 3.4, -21, 1.5 );
    // Ca N type h
	gating_parameter gate_can_h( -2, - 40, 75  );

    //solgat =  solve_gating( gate_can_m, gate_can_h, 0.06*1, 135, 2, 1 );    
    solgat.push_back(solve_gating( gate_can_m, gate_can_h, conductivity, 135, 2, 1 ));

}

void BG::install_cat_gates(double conductivity)
{
   gating_parameter gate_cat_m(3, -36, 1.5);
   gating_parameter gate_cat_h(-5.2, -68, 10);
   
   solgat.push_back(solve_gating(gate_cat_m, gate_cat_h, conductivity, 135, 2, 1));
}

void BG::install_cal_gates(double conductivity)
{
   gating_parameter gate_cat_m(3, -36, 1.5);
   gating_parameter gate_cat_h(0, 0, 0);
   
   solgat.push_back(solve_gating(gate_cat_m, gate_cat_h, conductivity, 135, 2, 0));
}

/////////////////////////////////////////////////////////////////////////////

// static objects, values to be choosen for each case seperately

BG::BG() {
   Neumann_flux = 0;
}

double BG::calc_current_at_start( double time, double Vm )
{
    double time_in_ms = time * 1000.;

    Neumann_flux = solgat[0].initial_value_for_current(Vm, time);
	
    return Neumann_flux;
}

/////////////////////////////////////////////////////////////////////////////
