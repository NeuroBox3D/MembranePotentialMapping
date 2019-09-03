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


#ifndef UG__PLUGINS__MEMBRANE_POTENTIAL_MAPPING__BG_SIMPLE__GATING_H
#define UG__PLUGINS__MEMBRANE_POTENTIAL_MAPPING__BG_SIMPLE__GATING_H

#include <iostream>

#include <cmath>

#include "voltage.h"

class gating_parameter
{
public:
    
	gating_parameter( double V_12_, double z_, double gamma_, double K_, double tau_0_ )
		: K( K_), z(z_), gamma(gamma_), V_12(V_12_), tau_0( tau_0_),
		F(9.648 * 10000), R(  8.314 * 1000 ), T( 300 ), simple(false)
	{};


    gating_parameter( double z_, double V_12_, double tau_0_ )
                : z(z_), V_12(V_12_), tau_0( tau_0_),
			F(9.648 * 10000), R(  8.314 * 1000 ), T( 300 ), simple(true)
    {};

	gating_parameter() {};
	
	double defect( double x_ip1_k, double x_i, double V,  double Delta_i );
 
	double tau_x ( double V );

	double x_infty( double V );

	double alpha_prime_x( double V );

	double beta_prime_x( double V ); 

        double x_1( double t, double x_0, double tau );

private:

	double K, z, gamma, V_12, tau_0, F, R, T;
	bool simple;
							 
};

#endif // UG__PLUGINS__MEMBRANE_POTENTIAL_MAPPING__BG_SIMPLE__GATING_H

