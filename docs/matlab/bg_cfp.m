% calculates derivatives of flux with respect to extra- and intracellular calcium concentrations for the BG model (cfp)
%% run as: matlab -nodesktop -nosplash -nojvm -nodesktop -r "bg, quit" | sed 's/log/std::log/g; s/exp/std::exp/g;'
clear

syms p z
syms F R T
syms Vm C_int C_ext

BG = p * z * z * F * F * Vm / R / T * (C_int - C_ext * exp(-z*F*Vm/(R*T))) / (1-exp(-z*F*Vm / (R*T)))

simplify(BG)

diff(BG, C_int)

simplify(diff(BG, C_int))
% BG end
