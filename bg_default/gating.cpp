/****************************************************************************/
/*                                                                          */
/* File:      gating.cc                                                	    */
/*                                                                          */
/* Purpose:   calculation of gating functions in Borg-Graham model          */
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

#include "gating.h"
#include "spannung.h"
#include "bg.h"

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

