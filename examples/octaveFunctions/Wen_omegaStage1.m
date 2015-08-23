% Function to calculate the dimensionless conversion time
% for first stage for a given reaction front position and time
% Formulars according Wen, "Noncatalytic Heterogenous Solid Fluid Reaction Models"
% models p. 43 eqn. 40
% Copyright: Stefan Radl, IPPT, TU Graz
function [omega] = Wen_omegaStage1(ThieleMod, N_sh)

    omega = 1.0 ...
          +  (1/N_sh) ... 
          .*(ThieleMod.*coth(ThieleMod) -1);

end
