R    = 5e-4;            %(outer) radius of the particle, m
k_mA = 10;              %external mass transfer coefficient, m/s
                        %must be set to large if Dirichlet BC is applied!
porosity = 0.5;         %porosity
D_Mol    = 2e-4;        %Molecular rate of diffusion, typical for oxygen in air
D_eA = D_Mol*porosity;  %EFFECTIVE rate of diffusion
DeIA = D_eA;            %rate of diffusion of A (fluid) in reaction zone
k_v = 2e3;              %reaction rate constant (multiplied with Cs0 in analysis)!
CA0 = 2.12e-3;          %initial concentration of gas phase, in kmol/m³_gas, for y=0.20 at 1 bar and 1098 K
Cs0 = 14.1;             %initial concentration of solid, in kmol/m³_tot, for Cu and a volume fraction of 0.1 in the particle
a   = 0.5;              %stoichometric coeff: ratio of gas to solid stoichimetric ratio
radius_interval = 1e-5; % interval for radius output 

%End of user input
k_vIntr = k_v / Cs0;
