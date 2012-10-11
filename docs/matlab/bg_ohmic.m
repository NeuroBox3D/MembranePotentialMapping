% calculates derivatives of flux with respect to extra- and intracellular calcium concentrations for the BG model (ohmic)
%% run as: matlab -nodesktop -nosplash -nojvm -nodesktop -r "bg_ohmic, quit" | sed 's/log/std::log/g; s/exp/std::exp/g;'
clear

syms g x xp y yp Voltage
syms V_ret Ca dt valency

%BG = - g * x^xp * y^yp * Voltage - V_ret
BG = (1e3 * (dt * -g * x^xp * y^yp * Voltage - V_ret * 6.24e18)) / (6.022e23* valency)

simplify(BG)

diff(BG, Ca)

simplify(diff(BG, Ca))
% BG end
