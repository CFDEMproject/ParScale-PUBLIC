    PARSCALE_SRC_DIR = getenv('PASCAL_SRC_DIR');
    if(isempty(PARSCALE_SRC_DIR))
        error('The user has not set the environment variable "PASCAL_SRC_DIR". Please do so, e.g., in your .bashrc file, to use this octave script.')
    end
    addpath([PARSCALE_SRC_DIR,'/../examples/octaveFunctions']);

%%% PROCESS PARAMETERS %%%%

    c_surface   = 1e-4;                 % vapor concentration at the surface of the particle (Dirichlet boundary)
    D_vapor     = 1.1307e-5             % Diffusion coefficient vapor in the gas [m²/s]
    pTotal      = 1e5;                  % Atmospheric pressure [Pa]
    Temp        = 350;                  % isothermal particle temperature [K]

%%%% PHYSICAL PARTICLE PROPERTIES %%%%

    R_particle        = 1e-3;             % Particle radius [m]
    r_core            = 20/40*R_particle    % Initial wet particle core [m]
    deltaR            = 1/40*R_particle   % Initial wet particle core [m]
    deltaV            = 4/3*pi*((r_core+0.5*deltaR)^3-(r_core-0.5*deltaR)^3)
    liquid_phase_frac = 0.9;               % Liquid phase fraction of wet core - equals gas phase fraction
    gas_phase_frac    = liquid_phase_frac; % Gas phase fraction of dry shell core (must be identical to liquid phase fraction in order for analytical solution to be correct)
    tortuosity        = 1.0;               % Tortuosity of the pores
    Uni_gas_constant  = 8314.462;          % Universal gas constant [J/kmol K]

    %Parameters for Ethanol
    MG_liq      = 46.07;         % Molar weight liquid [kg/kmol]
    rho_liquid  = 772;           % Density of liquid [kg/m³]
    params.A    = 8.20417;
    params.B    = 1642.89;
    params.C    = -42.85;

    %Compute vapor pressure and key reference quantities
    D_eff       = D_vapor * gas_phase_frac / tortuosity;
    p_e_liq     = pVapAntoine(Temp, params)         % vapor pressure [Pa] 
    x_core      = p_e_liq/pTotal                      % Equilibrium mole-fraction (at r_core)
    c_core      = p_e_liq/Uni_gas_constant/Temp     % Equilibrium vapor concentration [kmol/m³] (at r_core)
    c_total     = pTotal/Uni_gas_constant/Temp       % Equilibrium vapor concentration [kmol/m³] (at r_core)
    x_surface   = c_surface / c_total                % Mole-fraction at the surface of the particle (Dirichlet boundary)

    %Put parameters in a struct for easy hand-over to functions
    properties = struct;
    properties.ctotal       = c_total;
    properties.rhoLiquid    = rho_liquid;
    properties.MWLiquid     = MG_liq;
    properties.D_eff        = D_eff;
    properties.pTotal       = pTotal;
    properties.Temp         = Temp ;
    properties

    particle.dParticle      = 2*R_particle;
    particle.rCoreStar      = r_core/R_particle;
    particle.liqPhaseFrac   = liquid_phase_frac;
    particle    
