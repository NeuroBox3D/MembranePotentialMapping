% calculates derivatives of flux with respect to extra- and intracellular calcium concentrations for the BG model (cfp)
%% run as: matlab -nodesktop -nosplash -nojvm -nodesktop -r "bg_cfp, quit" | sed 's/log/std::log/g; s/exp/std::exp/g;'
clear

syms permeability valency
syms F R T
syms Vm C_int C_ext
syms dt

BG = (1e3 * (dt * p * valency * valency * F * F * Vm / R / T * (C_int - C_ext * exp(-valency*F*Vm/(R*T))) / (1-exp(-valency*F*Vm / (R*T))) * 6.24e18)) / (valency * 6.022e23)

simplify(BG)

diff(BG, C_int)

simplify(diff(BG, C_int))
% BG end
