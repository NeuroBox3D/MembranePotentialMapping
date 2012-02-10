#ifndef __my_commands_h__
#define __my_commands_h__

double timestepping_of_gates_and_calc_current( double time, double delta_t , double Vm);

double getNeumannFlux(double delta_t, double Vm);

double getCurrent(double delta_t, double Vm);

#endif
