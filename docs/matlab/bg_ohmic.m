% calculates derivatives of flux with respect to extra- and intracellular calcium concentrations for the BG model (ohmic)
%% run as: matlab -nodesktop -nosplash -nojvm -nodesktop -r "bg_ohmic, quit" | sed 's/log/std::log/g; s/exp/std::exp/g;'
clear

syms g x xp y yp Voltage
syms V_ret Ca

BG = - g * x^xp * y^yp * Voltage - V_ret

simplify(BG)

diff(BG, Ca)

simplify(diff(BG, Ca))
% BG end
