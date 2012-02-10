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
#include "vm2ug.h"

using namespace std;
using namespace vug;

namespace bg {

class BG
{
public:

   static int ap_interval_duration_in_ms; // which stimulation transformation should be choosen
   void install_can_gates(double cond = 1000);
   bool installed_can_gates();
   
 solve_gating solgat;

 double Neumann_flux;

 string output_file_current;
 
 double  conductivity; // assumed as 1000 mA / mV !
 
 static double Voltage( double time_glob ); // returns the voltage ( time in milliseconds here: if trivial case, i. e. typical AP given)
  
 static double ttrafo_into_ap( double time ); // 100 Hz trafo; time in ms here
   
   double timestepping_of_gates_and_calc_current( double time, double delta_t );
   double calc_current_at_start( double time );
   double get_Neumann_Flux();
   BG();
   static Vm2uG<std::string> mapping;

   static double x, y, z;
   static void setXYZ(double x, double y, double z);
   static void setTimeStepBaseDir(std::string);
   static std::string TimeStepBaseDir;
private:
   bool inst_can_gates;
};
}	 

#endif

