#ifndef __solve_h__
#define __solve_h__

#include "gating.h"
#include <vector>
#include <fstream>

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

		solve_gating( gating_parameter gate_x_, gating_parameter gate_y_, double g_, double V_ret_, int xp_, int yp_ , double z_flux_, double permeability_)
			: gate_x( gate_x_ ), gate_y( gate_y_ ), g(g_), V_ret(V_ret_), xp(xp_), yp(yp_), z_flux(z_flux_), permeability(permeability_), F(9.548 * 10000), R(8.314 * 1000), T(300)
		{};


        double initial_value_for_current( double basic_voltage = -65, double time = 0.0); // t = 0, uses -65 mV as standard value may be changed
        double initial_value_for_current( double Ca_ext, double Ca_int, double basic_voltage, double time);
	double calculate_current_expli_next_timestep( double time, double delta_t, double Vm, double Ca_ext, double Ca_int ); 
	double calculate_current_expli_next_timestep( double time, double delta_t, double Vm);

        double current_x() { return current_x_0; }
        double current_y() { return current_y_0; }

        double Current_current( double time, double Ca_ext, double Ca_int ) { return I(  current_x_0, current_y_0, time, Ca_ext, Ca_int ); }
      double Current_current (double time) { return I(current_x_0, current_y_0, time); }

private:

	double I( double x, double y, double t );

	double I( double x, double y, double t, double Ca_ext, double Ca_int);
   
  double Flux(double Ca_ext, double Ca_int, double time);

	gating_parameter gate_x, gate_y;

	double g, V_ret;

	int xp, yp;

  double current_x_0, current_y_0;
  
   double z_flux;
   double permeability;
   double F, R, T;


};
#endif
