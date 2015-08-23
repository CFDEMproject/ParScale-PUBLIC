% Function to calculate the dimensionless conversion time
% for second stage for a given reaction front position and time
% Formulars according Wen, "Noncatalytic Heterogenous Solid Fluid Reaction Models"
% models p. 43 eqn. 40
% Copyright: Stefan Radl, IPPT, TU Graz
function [omega] = Wen_omegaStage2(ThieleMod, N_sh, D_ratio)

    omega = 1 ...
         + ThieleMod.*ThieleMod/6.0 ...
         .*(1+2./N_sh) ...
         + (1-D_ratio) ...
         .* log( sinh(ThieleMod) ./ ThieleMod);

end
