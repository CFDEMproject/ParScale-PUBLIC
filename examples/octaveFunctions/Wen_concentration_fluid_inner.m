% Function to calculate the concentration of the fluid at 0 <= r <= r_m
% for second stage reaction for a given reaction front position and time
% Formulars according Wen, "Noncatalytic Heterogenous Solid Fluid Reaction Models"
% models p. 43 eqn. 36
% Copyright: Stefan Radl / Thomas Forgber, IPPT, TU Graz

function [concentration_fluid_inner] = Wen_concentration_fluid_inner (xi_m, ThieleMod,concentration_reaction_front,xi)

    concentration_fluid_inner =  ...
                                concentration_reaction_front ...
                             .* xi_m./xi ...
                             .* sinh(ThieleMod.*xi) ...
                             ./ sinh(ThieleMod.*xi_m);
end
