function [pVap] = pVapAntoine(temp, params)
%function will return vapor pressure in Pascal
%standard antoine parameters that use temp in Kelvin and 
%result in a vapor pressure in mmHg (conversion factor mmHg->Pa = 133.3)

    pVap = 133.322368 * 10.^(params.A-params.B./(params.C+temp));

end
