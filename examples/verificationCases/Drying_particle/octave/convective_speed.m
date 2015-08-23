%%%% Calculation for species flux based on analytical solution and
%%%% first order discretisation
%%%% Copyright: Thomas Forgber, IPPT 2015, Graz

function [conv_velo] = convective_speed(flux,properties) 

c_tot = properties.ctotal;

conv_velo = flux./c_tot;

