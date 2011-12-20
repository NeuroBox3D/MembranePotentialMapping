#include "bg.h"

using namespace bg;

void BG::install_can_gates(double conductivity)
{
	gating_parameter gate_can_m( 3.4, -21, 1.5);
	gating_parameter gate_can_h( -2, - 40, 75);
  solgat.push_back(solve_gating( gate_can_m, gate_can_h, conductivity, 135, 2, 1 ));

}

void BG::install_cat_gates(double conductivity) 
{
   gating_parameter gate_cat_m(3, -36, 1.5);
   gating_parameter gate_cat_h(-5.2, -68, 10);
   solgat.push_back(solve_gating(gate_cat_m, gate_cat_h, conductivity, 135, 2, 1));
}

void BG::install_cal_gates(double conductivity)
{
   gating_parameter gate_cal_m(3, -36, 1.5);
   gating_parameter gate_cal_h(0, 0, 0);
   solgat.push_back(solve_gating(gate_cal_m, gate_cal_h, conductivity, 135, 2, 0));
}

BG::BG()
{
   Neumann_flux = 0;
}

double BG::calc_current_at_start(double Vm)
{
    Neumann_flux = solgat[0].initial_value_for_current(Vm);
    return Neumann_flux;
}
