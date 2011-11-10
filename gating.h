/****************************************************************************/
/*                                                                          */
/* File:      gating.h                                                	    */
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


#ifndef __gating_h__
#define __gating_h__

#include <iostream>

#include <cmath>


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
	
        //(  x^{i+1}_{k+1}  - x^{i+1}_{k} )
	double defect( double x_ip1_k, double x_i, double V,  double Delta_i );
 
	double tau_x ( double V );

	double x_infty( double V );

	double alpha_prime_x( double V );

	double beta_prime_x( double V ); 

        double x_1( double t, double x_0, double tau, double Vm );

private:

	double K, z, gamma, V_12, tau_0, F, R, T;
	bool simple;
							 
};

#endif

