%%%% Reference solution for 1D Stefan Diffusion with Dirichilet boundary condition
%%%% Target: Calculation of the droplet diameter during drying
%%%% particle
%%%% Copyright: Forgber, IPPT 2015, Graz
function [d_evo] = diameter_evo(rho_liquid, d_liq_gas,p_atm,T_iso,MG_liq,Uni_gas_constant,d_particle,x_surface,x_core,time_evap,liquid_phase_frac) 

d_evo = sqrt(d_particle.^2 - ((time_evap .* d_liq_gas.* p_atm.*log((1-x_surface)./(1-x_core)) .* 8.0 .* MG_liq)./(liquid_phase_frac .* rho_liquid.*Uni_gas_constant.*T_iso)));




