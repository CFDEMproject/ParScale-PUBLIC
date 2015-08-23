     d_particle = 5e-3;         %Particle diameter   
     T_0_sphere = 800.0;        %initial temperature sphere
     T_inlet   = 400.0;         %evironment temperature
     t_end = 1;                 %time [sec]
     alpha = 100.0;             %heat transfer coefficient [W/m²K]

% Particle Properties    
    void_frac_solid = 0.8;    %void fraction solid
    lambda_solid =   5.0;     %thermal conductivity solid
    rhoP_solid   =   1000;    %density solid
    cpP_solid    =   300.0;   %thermal capacity solid

    void_frac_liquid = 0.1;   %void fraction liquid
    lambda_liquid =   10.0;   %thermal conductivity liquid
    rhoP_liquid   =   10.0;     %density liquid
    cpP_liquid    =   10.0;   %thermal capacity liquid

    void_frac_gas = 0.1;      %void fraction gas
    lambda_gas =   1.0;       %thermal conductivity gas
    rhoP_gas   =   1.0;       %density gas
    cpP_gas    =   1.0;       %thermal capacity gas
