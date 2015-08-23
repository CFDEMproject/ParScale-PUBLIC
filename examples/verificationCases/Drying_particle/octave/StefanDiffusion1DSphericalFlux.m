%%%% Calculation of the DIMENSIONAL species flux based on analytical solution in a particle 
%%%% with wet core
%%%% Copyright: Thomas Forgber, IPPT 2015, Graz

function [flux] = StefanDiffusion1DSphericalFlux(r_dimless,r_starCore, properties, particle, x_surface, x_core) 

D_eff     = properties.D_eff;
c_tot     = properties.ctotal;
rParticle = particle.dParticle/2;

 
if(r_dimless < r_starCore)
    flux = 0;
else
    flux = c_tot .* D_eff ./ rParticle                      ...  
           .* 1./(r_dimless.^2*(r_starCore.^(-1)-1)) ...
           .*log((1-x_surface)./(1-x_core));      
endif  

end
