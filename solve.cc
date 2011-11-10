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
#include <assert.h>

double solve_gating::I( double x, double y, double Vm )
{
    double curr = - g * pow( x, xp ) * pow( y, yp ) * ( Vm - V_ret );    

	return curr;
}

double solve_gating::I( double x, double y, double Vm, double Ca_ext, double Ca_int )
{
    double curr = pow( x, xp ) * pow( y, yp ) * Flux(Ca_ext, Ca_int, Vm);    

	return curr;
}


/////////////////////////////////////////////////////////////////////////////

double solve_gating::initial_value_for_current( double basic_voltage, double time ) // t = 0
{
  current_x_0 = gate_x.x_infty( basic_voltage );
  current_y_0 = gate_y.x_infty( basic_voltage );

  //double time = 0.;

  return I( current_x_0, current_y_0, time );
}

double solve_gating::calculate_current_expli_next_timestep( double time, double delta_t, double Vm )
{
   current_x_0 = gate_x.x_1( time, current_x_0, delta_t, Vm );
   current_y_0 = gate_y.x_1( time, current_y_0, delta_t, Vm );

   return I( current_x_0, current_y_0, time );
} 

double solve_gating::calculate_current_expli_next_timestep( double time, double delta_t, double Vm, double Ca_ext, double Ca_int )
{
    
   current_x_0 = gate_x.x_1( time, current_x_0, delta_t, Vm );
   current_y_0 = gate_y.x_1( time, current_y_0, delta_t, Vm );

   return I( current_x_0, current_y_0, time, Ca_ext, Ca_int );
} 



/////////////////////////////////////////////////////////////////////////////
double solve_gating::Flux( double Ca_ext, double Ca_int, double Vm )
{

   return (permeability * z_flux*z_flux * F*F * Vm / ( R*T ) * ( Ca_int - Ca_ext * exp( - z_flux * F * Vm / ( R*T ) ) ) / ( 1 - exp( - z_flux * Vm / ( R*T ) ) ));

}



double solve_gating::initial_value_for_current( double basic_voltage, double time, double Ca_ext, double Ca_int)
{

  current_x_0 = gate_x.x_infty( basic_voltage );
  current_y_0 = gate_y.x_infty( basic_voltage );

 // double time = 0.;

  return I( current_x_0, current_y_0, time, Ca_ext, Ca_int );
}

