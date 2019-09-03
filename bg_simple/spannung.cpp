/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus M. Knodel
 * Creation date: 2009-08
 *
 * This file is part of NeuroBox, which is based on UG4.
 *
 * NeuroBox and UG4 are free software: You can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3
 * (as published by the Free Software Foundation) with the following additional
 * attribution requirements (according to LGPL/GPL v3 §7):
 *
 * (1) The following notice must be displayed in the appropriate legal notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 *
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 *
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating PDE based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * "Stepniewski, M., Breit, M., Hoffer, M. and Queisser, G.
 *   NeuroBox: computational mathematics in multiscale neuroscience.
 *   Computing and visualization in science (2019).
 * "Breit, M. et al. Anatomically detailed and large-scale simulations studying
 *   synapse loss and synchrony using NeuroBox. Front. Neuroanat. 10 (2016), 8"
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include "voltage.h"

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
