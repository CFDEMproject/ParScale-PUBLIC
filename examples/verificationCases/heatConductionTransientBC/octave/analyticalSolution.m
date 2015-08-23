%% This is the analytical solution for cumputing the 1-D temperature distribution
%% due to surface evaporation on a single sphere 
%% after Crank - Mathematics of diffusion page 96 eqn. 6.40
%% not suitable for Biot-numbers <10e-3

clear
clc
close all
more off

%add the path to functions for postprocessing
PARSCALE_SRC_DIR = getenv('PASCAL_SRC_DIR');
if(isempty(PARSCALE_SRC_DIR))
    error('The user has not set the environment variable "PASCAL_SRC_DIR". Please do so, e.g., in your .bashrc file, to use this octave script.')
end
addpath([PARSCALE_SRC_DIR,'/../examples/octaveFunctions']);


%format compact;
format long g;

%% System parameters (increase values to get a "better" and better resolved analytical solution
nmax = 10;                   % maximum number of periodic zero values of root function beta_n
numb_grid_points = 25;       % radial grid resolution inside sphere
 
% *** USER INPUT *************************************************
% Physical Parameters
run('parametersHeatConduction.m')

% (Uncomment if sphere internal cooling curve is of interest):
surface_evaporation_solution = 'transient_heat_conduction.dat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (void_frac_gas+ void_frac_liquid+void_frac_solid !=1.0)
    disp('Your void fractions are not summing up to 1 - please check!');
    break;
endif

fid=fopen(surface_evaporation_solution,'w');

%Simulation parameters:
sum_h = 0;                   % sum varialble 
sum_m = 0;                   % sum variable for total transport
n = 1;                       % running variable for sum

%geometric parameters
r =  1E-12;                  % start point inside sphere for temperature function computing
XMIN = 1E-12;                 % minimum distance from pole, necessary for steady solution
r_particle =d_particle/2.0;
dx=(r_particle-r)./(numb_grid_points-1);   %increment between grid points

% Calcuate effective particle properties

lambda_eff = effective_value(void_frac_solid,void_frac_liquid,void_frac_gas,lambda_solid,lambda_liquid,lambda_gas)
rhoP_eff = effective_value(void_frac_solid,void_frac_liquid,void_frac_gas,rhoP_solid,rhoP_liquid,rhoP_gas)
cpP_eff = effective_value(void_frac_solid,void_frac_liquid,void_frac_gas,cpP_solid,cpP_liquid,cpP_gas)
 
% Dimensionless Parameters
a_eff = lambda_eff/(rhoP_eff*cpP_eff)
Fo=((a_eff*t_end)/(d_particle*d_particle))
Bi=(alpha*d_particle)/lambda_eff

%% Analytical solution of surface evaporation %%   

results = zeros(nmax,1); % preallocating the vectorspace of result vector

fun = @(beta_n) beta_n*cot(beta_n)+Bi-1; % root function of beta_n

    while r <=(r_particle+XMIN)                   %Loop over radial distance  

     %fid=fopen(surface_evaporation_solution,'w');   

                    while n <= nmax
                % This function computes the zero values of the function:
                % beta_n *cot(beta_n)+Bi-1=0 , 
                % needed for the analytical solution of surface evaporation on a sphere
                % (n-1)*pi+0.1 is the start value for searching the roots (function is pi-periodic)

                    results = fsolve(fun,(n-1)*pi+0.1);     %finding the roots of function "fun"

                    h= sin(results*r/r_particle) ...
                        /( sin(results) ) ...
                       *  exp(-Fo*results*results) ...
                        /( results*results+Bi*(Bi-1) );

                    sum_h = h + sum_h   ; 

                    m = (6*Bi^2 * exp(-Fo*results*results))/ ...   % total sum of transported scalar
                        (results^2*(results^2+Bi*(Bi-1))) ;

                    sum_m = m + sum_m   ; 
                                %Perform summation for analytical solution
                    n=n+1; 

                    end

        T_rad = (  (T_0_sphere - T_inlet) *2*Bi*r_particle/r*sum_h  )...
               +T_inlet;

        m_total = 1-sum_m   ; % total transported scalar 


        radius = r/r_particle; %relative distance from sphere center [0-1]
        fprintf(fid,['%g \t %g \n'],radius,T_rad);

        r=r+dx;         %increasing increment of radial position
        h = 0;
        sum_h = 0; % resetting the simulation parameters for inner loop
        m =0 ;
        sum_m = 0 ; 
        n = 1;      
        x0 = 1; 

%% If internal temerature profile of sphere is of interest: 
%         subplot(1,2,1), plot (radius, T_rad, 'o');
%         xlabel ('relative distance from sphere center');
%         ylabel ('Temperature [K]');
%         axis ([0 1 T_inlet T_0_sphere]);
%         hold on ;
        
%%    
    end
    T_surf= T_rad;
    
coolingSphere = struct;
coolingSphere.T_surfCranck = T_surf;  
coolingSphere.m_total = m_total;  

%  plot (Fo, T_surf, 'o');
%      xlabel ('Fo');
%      ylabel ('Sphere_Surface_Temperature [K]');
% %     axis ([0 10 T_ambient T_initial]);
%      hold on ;     

%fprintf(fid,[num2str(Fo),'\t',num2str(T_surf),'\n']);    

%end
