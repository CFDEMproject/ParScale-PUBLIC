%%%% Reference solution for 1D Stefan Diffusion with Dirichlet boundary condition
%%%% Target: Calculation of mol-fraction (x_i) along the radius of a spherical
%%%% Based on: 
%%%% Copyright: Forgber/Radl, IPPT 2015, Graz
function [x] = StefanDiffusion1DSpherical(r_starCore, x_core, x_surface, r_dimless) 

 if (r_dimless>r_starCore) 
     x = 1 ...
       - (1-x_core) ...
       .*( ...
              (  (1-x_surface) ...
                /(1-x_core)  ...
              ) ...
              .^ ( ...
                       (1 - r_starCore./r_dimless ) ...
                     / (1 - r_starCore)...
                 ) ...
         );
 else
     x = x_core;
 endif


