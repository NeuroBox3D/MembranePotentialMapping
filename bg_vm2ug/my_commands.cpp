/****************************************************************************/
/*                                                                          */
/* File:      my_commands.cc                                            	*/
/*                                                                          */
/* Purpose:   several own functions in order to use the the classes    */ 
/* 
* Author:	  Markus M. Knodel                                              */
/*                Goethe Center for Scientific Computing             */
/*                University of Frankfurt                   */
/*                Kettenhofweg 139                           */
/*                60325 Frankfurt              */
/*                Germany              */
/*                email: markus.knodel@gcsc.uni-frankfurt.de    */              
/*                                                                  */
/* History:   2009 begin            									   
*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/



#include "my_commands.h"
#include "bg.h"
#include "solve.h"
#include <fstream>
#include <iostream>
#include "spannung.h"

using namespace std;

using namespace bg;
/////////////////////////////////////////////////////////////////////////////

double BG::timestepping_of_gates_and_calc_current( double time, double delta_t )
{
    double time_in_ms = time * 1000.;
    double delta_t_in_ms = delta_t * 1000.;

      Neumann_flux = solgat.calculate_current_expli_next_timestep( time_in_ms, delta_t_in_ms ); 

 /*   fstream dat_stream;
    dat_stream.open( output_file_current.c_str(), ios::out | ios::app );

    dat_stream << time << "   " <<  solgat.current_x() << "  " 
                     << solgat.current_y() << "  " 
	       << solgat.Current_current( time_in_ms ) << "   " 
                     << Voltage( time_in_ms ) << "    " << Neumann_flux << endl;

    dat_stream.close(); */
   
 //   return Neumann_flux;
   return solgat.Current_current(time_in_ms);
}

/////////////////////////////////////////////////////////////////////////////

double BG::calc_current_at_start( double time )
{
    double time_in_ms = time * 1000.;

    Neumann_flux = solgat.initial_value_for_current();

    // in case of parallel programming following must be done only by master!

 /*   fstream dat_stream;
    dat_stream.open( output_file_current.c_str(), ios::out );

    dat_stream << time << "   " << solgat.current_x() << "  " 
                     << solgat.current_y() << "  " 
                     << solgat.Current_current( time_in_ms ) << "   " 
                     << Voltage( time_in_ms ) << endl;

    dat_stream.close(); */
	
    return Neumann_flux;
}

/////////////////////////////////////////////////////////////////////////////


