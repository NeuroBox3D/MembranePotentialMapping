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

using namespace bg;
double solve_gating::I( double x, double y, double t )
{
    double curr = - g * pow( x, xp ) * pow( y, yp ) * ( BG::Voltage(t) - V_ret );    

	return curr;
}

/////////////////////////////////////////////////////////////////////////////

double solve_gating::initial_value_for_current( double basic_voltage ) // t = 0
{
  current_x_0 = gate_x.x_infty( basic_voltage );
  current_y_0 = gate_y.x_infty( basic_voltage );

  double time = 0.;

  return I( current_x_0, current_y_0, time );
}

double solve_gating::calculate_current_expli_next_timestep( double time, double delta_t )
{
    
   current_x_0 = gate_x.x_1( time, current_x_0, delta_t );
   current_y_0 = gate_y.x_1( time, current_y_0, delta_t );

   return I( current_x_0, current_y_0, time );
} 

/////////////////////////////////////////////////////////////////////////////
