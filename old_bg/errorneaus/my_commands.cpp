#include "my_commands.h"
#include "bg.h"
#include "solve.h"
#include <fstream>
#include <iostream>

using namespace std;
using namespace bg;
/////////////////////////////////////////////////////////////////////////////

double BG::timestepping_of_gates_and_calc_current( double time, double delta_t, double Vm )
{
    double time_in_ms = time * 1000.;
    double delta_t_in_ms = delta_t * 1000.;

      Neumann_flux = solgat[0].calculate_current_expli_next_timestep( delta_t_in_ms, Vm ); 

    fstream dat_stream;
    dat_stream.open( output_file_current.c_str(), ios::out | ios::app );

    dat_stream << time << "   " <<  solgat[0].current_x() << "  " 
                     << solgat[0].current_y() << "  " 
	       << solgat[0].Current_current( time_in_ms ) << "   " 
                     << Vm << endl;

    dat_stream.close();

    return Neumann_flux;
}

/////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////

double BG::getNeumannFlux(double delta_t, double Vm) {

    double delta_t_in_ms = delta_t * 1000.;

      Neumann_flux  = solgat[0].calculate_current_expli_next_timestep( delta_t_in_ms, Vm ); 

    return Neumann_flux;
}

double BG::getCurrent(double delta_t, double Vm) {
 //  calc_current_at_start(Vm, time);   

   double delta_t_in_ms = delta_t * 1000.;
 
   double flux = 0.;
   for (unsigned int i=0; i < solgat.size(); i++)
      flux+=solgat[i].calculate_current_expli_next_timestep(delta_t_in_ms, Vm);
   
   Neumann_flux = flux;
 //  Neumann_flux = solgat[0].calculate_current_expli_next_timestep( time_in_ms, delta_t_in_ms, Vm ); 

   double current = 0.;

   for (unsigned int i=0; i < solgat.size(); i++)
      current += solgat[i].Current_current(Vm);
 //  return solgat[0].Current_current( time_in_ms );
   return current;

}
