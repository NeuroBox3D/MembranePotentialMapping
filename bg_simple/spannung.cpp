/****************************************************************************/
/*                                                                          */
/* File:      spannung.cc                                                	    */
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


#include "spannung.h"
#include "bg.h"

using namespace ug::membrane_potential_mapping::bg;
/////////////////////////////////////////////////////////////////////////////

double BG::Voltage( double time_glob ) // time in ms!!!!
{
    double time = ttrafo_into_ap( time_glob );

	if( time < 1 || time > 8 )
        return - 65;

    if( time < 2.5 ) 
		return -65 + 8. * ( time - 1 );

    if( time < 3 )
        return - 53 + 160 * ( time - 2.5 );

    if( time < 4.5 )
        return 27 - 62 * ( time - 3 );

    if( time < 7 )
        return -66 + 0.4 * ( time - 4.5 );

    return - 65;
}

/////////////////////////////////////////////////////////////////////////////

double BG::ttrafo_into_ap( double time ) // 100 Hz trafo
{
        int tt_ganz = static_cast<int>( time );

        double tt_anh = time - tt_ganz;

        //int tt_teil = tt_ganz % 10;
        int tt_teil = tt_ganz % 10; //ap_interval_duration_in_ms;

        double tt = tt_teil + tt_anh; 

        return tt;

}
/////////////////////////////////////////////////////////////////////////////
