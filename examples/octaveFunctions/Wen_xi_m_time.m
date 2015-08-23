function [xi_m] = Wen_xi_m_time(omega_v,D_ratio,ThieleMod,N_sh)
%Function to calculate the reaction front for the second stage 
%of a heterogeneous reaction for a given time and Thiele modulus
%Formulars according Wen, "Noncatalytic Heterogenous Solid Fluid Reaction Models"
%models p. 43 eqn. 39
% Copyright: Stefan Radl / Thomas Forgber, IPPT, TU Graz
% INPUT
%   TODOs (specify all input)
% OUTPUT
%   TODOs (specify all output)

    xi_m = 0.0 .* omega_v; %initialize output
    [omega_v_full,time_full] = Wen_omega_v(1e-6,ThieleMod,D_ratio,N_sh,1,1);

    %define inline function to be minimized    
    minFunction = inline('@Wen_omega_v(xi,ThieleMod,D_ratio,N_sh,1,1)-omega', ...
                         'xi','omega');

    for iTime=1:length(omega_v)
        if( omega_v(iTime) < Wen_omegaStage1(ThieleMod, N_sh) );
            xi_m(iTime) = 1.0;
%            disp(['xi_m = 1 because omega_v < omegaStage1']);
        elseif( omega_v(iTime) < omega_v_full )
            x0 = 1.0 - ( omega_v(iTime)-1 ) ./ ( omega_v_full-1 ); %guess based on linear function
              [X, FVAL, INFO] = fzero(@(x) minFunction(x,omega_v(iTime)), x0);
            if(INFO~=1)
                error('fzero search in Wen_xi_m_time did not converge!');
            end
            xi_m(iTime) = X;
          end
    end
end
