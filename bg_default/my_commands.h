#ifndef __my_commands_h__
#define __my_commands_h__
#include <string>

/****************************************************************************/
/*                                                                          */
/* File:      my_commands.h                                             	*/
/*                                                                          */
/* Purpose:   several own functions in order to use the the classes    */ 
/* 
 Author:	  Markus M. Knodel                                              */
/*                Goethe Center for Scientific Computing             */
/*                University of Frankfurt                   */
/*                Kettenhofweg 139                           */
/*                60325 Frankfurt              */
/*                Germany              */
/*                email: markus.knodel@gcsc.uni-frankfurt.de    */            
/*                                                                */
/* History:   2009 begin            					*/				    
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/


double timestepping_of_gates_and_calc_current( double time, double delta_t, double myVm );

double calc_current_at_start( double time = 0 );

#endif
