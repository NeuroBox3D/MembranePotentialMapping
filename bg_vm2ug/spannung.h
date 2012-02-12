/****************************************************************************/
/*                                                                          */
/* File:      spannung.h                                                	    */
/*                                                                          */
/* Purpose:   voltage model for action potentials                           */
/*			  										                        */
/*                                                                          */
/* Author:	  Markus M. Knodel                                              */
/*            Goethe Center for Scientific Computing GCSC                   */
/*              Kettenhofweg 139                                            */ 
/*              University Frankfurt                                        */
/*              Frankfurt, Germany                                          */
/*              email: markus.knodel@gcsc.uni-frankfurt.de                  */
/*																			*/
/* History:   August 09 begin              									*/
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/


#ifndef __spannung_h__
#define __spannung_h__
#include <string>

double ttrafo_into_ap( double time ); // transforms the time into the interval of one AP

double Voltage( double time, double x, double y, double z, double myVm ); // time: in ms!!!!

#endif
