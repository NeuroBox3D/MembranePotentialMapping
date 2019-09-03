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

#include "gating.h"
#include "bg.h"
#include "voltage.h"

using namespace ug::membrane_potential_mapping::bg;
double gating_parameter::defect(double x_ip1_k, double x_i, double V,
		double Delta_i) {
	return -((x_ip1_k - x_i) * tau_x(V) - Delta_i * (x_infty(V) - x_ip1_k))
			/ (tau_x(V) + Delta_i);
}

double gating_parameter::tau_x(double V) {
	double res = 0;

	if (simple)
		res = tau_0;
	else
		res = 1. / (alpha_prime_x(V) + beta_prime_x(V)) + tau_0;

	return res;
}

double gating_parameter::x_infty(double V) {
	double res = 0;

	if (simple)
		res = 1. / (1. + exp(-z * (V - V_12) * F / (R * T)));
	else
		res = alpha_prime_x(V) / (alpha_prime_x(V) + beta_prime_x(V));

	return res;
}

double gating_parameter::alpha_prime_x(double V) {
	double res = K * exp((z * gamma * (V - V_12) * F) / (R * T));
	return res;
}

double gating_parameter::beta_prime_x(double V) {
	double res = K * exp((-z * (1 - gamma) * (V - V_12) * F) / (R * T));
	return res;
}

double gating_parameter::x_1(double t, double x_0, double tau, double myVm) {
	double x1 = x_0
			+ tau * (x_infty(BG::Voltage(myVm)) - x_0)
					/ tau_x(BG::Voltage(myVm));

	return x1;
}

