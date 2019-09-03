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


#ifndef UG__PLUGINS__MEMBRANE_POTENTIAL_MAPPING__BG_SIMPLE__SOLVE_H
#define UG__PLUGINS__MEMBRANE_POTENTIAL_MAPPING__BG_SIMPLE__SOLVE_H

#include "gating.h"
#include <vector>
#include <fstream>
#include "voltage.h"

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

        double initial_value_for_current( double basic_voltage = -65 ); // t = 0, uses -65 mV as standard value may be changed
        double calculate_current_expli_next_timestep( double time, double delta_t );

        double current_x() { return current_x_0; }
        double current_y() { return current_y_0; }

        double Current_current( double time ) { return I(  current_x_0, current_y_0, time ); }

        double dCa() const {
        		return 0.0; // since no calcium concenctration considered then partial derivative w.r.t zo calcium concentration is zero!
        	}

	private:

		double I( double x, double y, double t );

		gating_parameter gate_x, gate_y;

		double g, V_ret;

		int xp, yp;

		double current_x_0, current_y_0;

};



#endif // UG__PLUGINS__MEMBRANE_POTENTIAL_MAPPING__BG_SIMPLE__SOLVE_H
