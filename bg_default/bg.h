/****************************************************************************/
/*                                                                          */
/* File:      globals.h                                            	*/
/*                                                                          */
/* Purpose:   global variables    */
/*                                                                          */
/* Author:	  Markus M. Knodel                                              */
/*                Goethe Center for Scientific Computing             */
/*                University of Frankfurt                   */
/*                Kettenhofweg 139                           */
/*                60325 Frankfurt              */
/*                Germany              */
/*                email: markus.knodel@gcsc.uni-frankfurt.de    */
/*																			*/
/* History:   2009 begin            									    */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

#ifndef __H__UG__MEMBRANE_POTENTIAL_MAPPING__BG__GLOBALS__
#define __H__UG__MEMBRANE_POTENTIAL_MAPPING__BG__GLOBALS__

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

class BG {
public:
	void install_can_gates(double cond = 1000);
	bool installed_can_gates();
	void install_cal_gates(double cond = 0.3 * 0.6);
	void install_cat_gates(double cond = 0.3);
	bool installed_cal_gates();
	bool installed_cat_gates();
	void install_can_gates_cfp(double permeability = 1e-8);
	void install_cal_gates_cfp(double permeability = 0.3 * 1e-8);
	void install_cat_gates_cfp(double permeability = 0.3 * 1e-8);

	solve_gating solgat;

	double Neumann_flux;

	double conductivity; /* [mA/mV] */
	double permeability;

	static double Voltage(double myVm);

	static double F;
	static double R;
	static double T;
	static double NA;
	static double C;

	double timestepping_of_gates_and_calc_current(double time, double delta_t,
			double myVm);
	double timestepping_of_gates_and_calc_current(double time, double delta_t,
			double myVm, double Ca_i, double Ca_o);
	double calc_current_at_start(double time);
	double calc_current_at_start(double time, double basic_voltage, double myVm,
			double Ca_i, double Ca_o);
	double get_Neumann_Flux();
	// returns flux in Mol/s
	double get_Neumann_Flux_as_Concentration(const double delta_t=1.0,
			const double valency = 2.0) const;

	inline double dCa_dCa_o(double delta_t=1e-4) const {
		return delta_t * 1e3 * solgat.dCadCa_o();
	}
	inline double dCa_dCa_i(double delta_t=1e-4) const {
		return delta_t * 1e3 * solgat.dCadCa_i();
	}
	inline double dCa(double delta_t=1e-4) const {
		return delta_t * 1e3 * solgat.dCa(); // solgat.dCa() === 0
	}

	inline double get_permeability() const {
		return permeability;
	}
	inline double get_conductivity() const {
		return conductivity;
	}
	BG();
private:
	bool inst_can_gates;
	bool inst_cal_gates;
	bool inst_cat_gates;
};
}
}
}

#endif // __H__UG__MEMBRANE_POTENTIAL_MAPPING__BG__GLOBALS__
