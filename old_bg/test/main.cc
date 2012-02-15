/****************************************************************************/
/*                                                                          */
/* File:      main.cc                                            	*/
/*                                                                          */
/* Purpose:   main program as example in order to get familiar with  */
/*            classes  used      */ 
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


/* standard C library */
#include <stddef.h>
#include <stdio.h>
#include <cmath>

/* standard C library */
#include <stddef.h>
#include <stdio.h>
#include <math.h>

#include "bg.h"
#include "my_commands.h"

using namespace bg;

int main (int argc, char **argv)
{
   BG b;
  double n =  b.Neumann_flux;
   b.install_can_gates(1000.0);

    double delta_t = 0.001;

    double maxtime = 0.5; 

    double time = 0;

    b.calc_current_at_start( time );

    for( int step = 1; step <= 80; step++ )
    {
        time = delta_t * 1;
        b.timestepping_of_gates_and_calc_current( time, delta_t );
    }

    for( int step = 1; step <= 80; step++ )
    {
        time = delta_t * 2;
        b.timestepping_of_gates_and_calc_current( time, delta_t );
    }

    return 0;
}
