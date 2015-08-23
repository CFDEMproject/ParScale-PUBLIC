%%%% Reference solution for simple diffusion with Dirichlet boundary condition
%%%% Target: Calculation of mol-fraction (x_i) along the radius of a spherical
%%%% Based on: 
%%%% Copyright: Forgber/Radl, IPPT 2015, Graz
function [result] = simpleDiffusion1DSpherical(r_starCore, x_core, x_surface, r_dimless)

if (r_dimless>r_starCore) 
result = ...
       x_core ...
     -  (x_core-x_surface) ...
     .* (1-r_starCore./r_dimless) ...
      ./(1-r_starCore);
else
     result = x_core;
end
