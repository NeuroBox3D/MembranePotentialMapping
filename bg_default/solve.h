/****************************************************************************/
/*                                                                          */
/* File:      solve.h                                           	*/
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


#ifndef __solve_h__
#define __solve_h__

#include "gating.h"
#include "spannung.h"
#include <vector>
#include <fstream>

using namespace std;

class solve_gating
{
public:

    	solve_gating()
			: gate_x( gating_parameter() ), gate_y(gating_parameter() ) 
        {};

		solve_gating( gating_parameter gate_x_, gating_parameter gate_y_, double g_, double V_ret_, int xp_, int yp_ )
			: gate_x( gate_x_ ), gate_y( gate_y_ ), g(g_), V_ret(V_ret_), xp(xp_), yp(yp_)
		{};

        double initial_value_for_current( double basic_voltage = -65, double myVm = -75 ); // t = 0, uses -65 mV as standard value may be changed
	double calculate_current_expli_next_timestep( double time, double delta_t , double myVm); 

        double current_x() { return current_x_0; }
        double current_y() { return current_y_0; }

        double Current_current( double time, double myVm) { return I(  current_x_0, current_y_0, time ,  myVm); }

private:

	double I( double x, double y, double t, double myVm );

	gating_parameter gate_x, gate_y;

	double g, V_ret;

	int xp, yp;

        double current_x_0, current_y_0;


};



#endif
