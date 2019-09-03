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


#ifndef UG__PLUGINS__MEMBRANE_POTENTIAL_MAPPING__BG_SIMPLE__BG_H
#define UG__PLUGINS__MEMBRANE_POTENTIAL_MAPPING__BG_SIMPLE__BG_H

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

namespace ug {
namespace membrane_potential_mapping {
namespace bg {

class BG
{
public:

   static int ap_interval_duration_in_ms; // which stimulation transformation should be choosen
   void install_can_gates(double cond = 1000);
   void install_cal_gates(double cond = 0.3*0.6);
   void install_cat_gates(double cond = 0.3);
   bool installed_can_gates();
   bool installed_cal_gates();
   bool installed_cat_gates();
   
 solve_gating solgat;

 double Neumann_flux;

 string output_file_current;
 
 double  conductivity; // assumed as 1000 mA / mV !
 
 static double Voltage( double time_glob ); // returns the voltage ( time in milliseconds here: if trivial case, i. e. typical AP given)
  
 static double ttrafo_into_ap( double time ); // 100 Hz trafo; time in ms here
   
   double timestepping_of_gates_and_calc_current( double time, double delta_t );
   double calc_current_at_start( double time );
   double get_Neumann_Flux();
   // returns flux in Mol/s
   double get_Neumann_Flux_as_Concentration(const double dt=1.0, const double valency=2.0) const;
   inline double dCa(double delta_t=1e-4) const {
   		return delta_t * 1e3 * solgat.dCa(); // solgat.dCa() === 0
   	}
   BG();
private:
   bool inst_can_gates;
   bool inst_cal_gates;
   bool inst_cat_gates;
};

} // namespace bg
} // namespace membrane_potential_mapping
} // namespace ug

#endif // UG__PLUGINS__MEMBRANE_POTENTIAL_MAPPING__BG_SIMPLE__BG_H

