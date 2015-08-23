%%%% Reference solution for simple diffusion with Dirichlet boundary condition
%%%% Target: Calculation of mol-fraction (x_i) along the radius of a spherical
%%%% Based on: 
%%%% Copyright: Forgber/Radl, IPPT 2015, Graz
function [flux] = simpleDiffusion1DSphericalFlux(r_dimless,r_starCore, properties, particle, x_surface, x_core)

D_eff     = properties.D_eff;
c_tot     = properties.ctotal;
rParticle = particle.dParticle/2;

 
if(r_dimless < r_starCore)
    flux = 0;
else
    flux = c_tot .* D_eff ./ rParticle                      ...  
           .* 1./(r_dimless.^2*(r_starCore.^(-1)-1)) ...
           .*(x_core-x_surface); 
endif  

end
