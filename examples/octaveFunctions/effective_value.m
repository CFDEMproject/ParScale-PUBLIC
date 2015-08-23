%%%% Function to calculate an effective property value for a 
%%%% three phase system
%%%% Copyright: Forgber, IPPT 2015, Graz
function [x] = effective_value(void_frac_solid,void_frac_liquid,void_frac_gas,value_solid,value_liquid,value_gas) 

x = void_frac_solid.*value_solid + void_frac_liquid.*value_liquid + void_frac_gas.*value_gas;
