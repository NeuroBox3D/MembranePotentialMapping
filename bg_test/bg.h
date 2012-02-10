/****************************************************************************/
/*                                                                          */
/* File:      globals.h                                            	*/
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


#ifndef __globals_h__
#define __globals_h__

/* system includes */
#include <stddef.h>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#include <vector>
#include <string>
#include "solve.h"
#include <fstream>

using namespace std;

namespace bg {

class BG
{
public:
   static int ap_interval_duration_in_ms;
   void install_can_gates(double cond = 1000);
   
 solve_gating solgat;

 double Neumann_flux;

 string output_file_current;
 
 double  conductivity;
 
 static double Voltage( double time_glob );
  
 static double ttrafo_into_ap( double time ); // 100 Hz trafo
   
   double timestepping_of_gates_and_calc_current( double time, double delta_t );
   double calc_current_at_start( double time );

   BG();
};
}	 

#endif

