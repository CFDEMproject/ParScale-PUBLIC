%%%% Main programm for 1D Stefan Diffusion with Dirichilet boundary condition
%%%% Copyright: Forgber, IPPT 2015, Graz

    clear
    clc
    close all
    more off

    PARSCALE_SRC_DIR = getenv('PASCAL_SRC_DIR');
    if(isempty(PARSCALE_SRC_DIR))
        error('The user has not set the environment variable "PASCAL_SRC_DIR". Please do so, e.g., in your .bashrc file, to use this octave script.')
    end
    addpath([PARSCALE_SRC_DIR,'/../examples/octaveFunctions']);

%%%% GRAPHICS/OUTPUT SETTINGS %%%%%
    plotSkip = 3;
    run('formatting.m')

%%%% SET FUNCTIONALLATY

    calc_error_stat                         = true;
    plot_solution                           = true;
    get_analytical_solution                 = true;
    get_numerical_solution_species          = true;
    get_numerical_solution_convective_flux  = true;     % writedebugcontainers in in.file has to be activated!
    calcSpeciesFlux                         = true;
    calc_evaportation_time                  = true;

    run('parametersDrying.m')
#    time_ParScale = ["0.000100"];% Exact folder name of time that will be compared to analytical solution
#    time_ParScale = ["0.000200"];% Exact folder name of time that will be compared to analytical solution
#    time_ParScale = ["0.000300"];% Exact folder name of time that will be compared to analytical solution
#    time_ParScale = ["0.000400"];% Exact folder name of time that will be compared to analytical solution
#    time_ParScale = ["0.700001"];% Exact folder name of time that will be compared to analytical
#    time_ParScale = ["0.020000"];% Exact folder name of time that will be compared to analytical solution
    time_ParScale = ["0.020010"];% Exact folder name of time that will be compared to analytical solution
#    time_ParScale = ["0.000200"];% Exact folder name of time that will be compared to analytical solution
#    time_ParScale = ["0.049500"];% Exact folder name of time that will be compared to analytical solution
#    time_ParScale = ["0.149600"];% Exact folder name of time that will be compared to analytical solution
#    time_ParScale = ["0.010000"];% Exact folder name of time that will be compared to analytical solution
#    time_ParScale = ["0.015000"];% Exact folder name of time that will be compared to analytical solution
#    time_ParScale = ["0.005000"];% Exact folder name of time that will be compared to analytical solution
#    time_ParScale = ["1.000100"];% Exact folder name of time that will be compared to analytical solution
#    time_ParScale = ["2.000100"];% Exact folder name of time that will be compared to analytical solution
#    time_ParScale = ["3.000000"];% Exact folder name of time that will be compared to analytical solution
#    time_ParScale = ["4.000000"];% Exact folder name of time that will be compared to analytical solution
#    time_ParScale = ["5.000000"];% Exact folder name of time that will be compared to analytical solution
    particle_id = 1;             % ID of particle that is compared

% ====================================================================================
%%%%% USER INPUT END %%%%%%%
% ====================================================================================

%%%%% MAIN PROGRAM START %%%%%%%
    time_to_plot = cellstr (time_ParScale);

%%%%%%% GET NUMERICAL SOLUTION SPECIES %%%%%%%%
    prop_name = "species.json"   % file name of property hat will be compared
    d=1;                         % Running variable times

    if(get_numerical_solution_species)

        while d <= size (time_ParScale)
            fileName  = ['../',time_to_plot{d,1},'/',prop_name]
            [xDat, yDat, misc] = jsonGetParScaleData(fileName,'data',particle_id); 
            d=d+1;
        endwhile
    
    endif

    gridpoints  = length(xDat)-1;                   % Number of gridpoints 

%%%%%%% GET NUMERICAL SOLUTION SPECIES %%%%%%%%

    prop_name = "gasConvectiveFlux.json"   % file name of property hat will be compared
    d=1;                                   % Running variable times
    fileName  = ['../',time_to_plot{d,1},'/',prop_name]
    if(~exist(fileName))
        get_numerical_solution_convective_flux = false;
    endif
    if(get_numerical_solution_convective_flux )
    
        while d <= size (time_ParScale)
            
            [xDat, yDat_flux, misc] = jsonGetParScaleData(fileName,'data',particle_id); 
            d=d+1;
        endwhile
    endif



%%%%% SOME INIT CALCULATIONS %%%%%%%

    error_stat = '../verificationStats/error_stat.dat';
    fid=fopen(error_stat,'w');


    if (r_core > R_particle)
        disp('Your Wet-core radius has to be smaller than the particle radius!');
        break;
    endif

    r=1e-16;
    i = 1;
    r_interval = (R_particle/(gridpoints-1));

%%%%%% GET ANALYTICAL SOLUTION %%%%%%%%

    if(get_analytical_solution)

        while (i < (gridpoints+1))

            r_dimless(i)=((i-1)*r_interval)/R_particle;
            [xAnalytical(i)] = StefanDiffusion1DSpherical(r_core / R_particle, x_core,x_surface,r_dimless(i)); 
            [xAnalyticalSimple(i)] = simpleDiffusion1DSpherical(r_core / R_particle, x_core,x_surface,r_dimless(i));

        if(calcSpeciesFlux)
            [flux(i)]               = StefanDiffusion1DSphericalFlux( r_dimless(i),r_core / R_particle, properties, particle, x_surface, x_core);
            integralFlowRate        = flux(i) * 4 * pi * (r_dimless(i)*R_particle)^2;

            [fluxSimple(i)]         = simpleDiffusion1DSphericalFlux( r_dimless(i),r_core / R_particle, properties, particle, x_surface, x_core);
            integralFlowRateSimple  = fluxSimple(i) * 4 * pi * (r_dimless(i)*R_particle)^2;

        endif

        r=r+r_interval;
        i=i+1;

        endwhile 

    endif


%%%%%% CALC EVAPORATION TIME %%%%%%%%
    
    if(calc_evaportation_time)
        [time_evap, time_simple] = evaporation_time(properties, particle, x_surface, x_core,r_core)
    endif


%%%%%%% ERROR STATISTICS %%%%%%

    if(calc_error_stat)

        spartial_index = 1;

        while (spartial_index<=size(xAnalytical,2))
            spartial_error(spartial_index) = abs( ...
                                                  (yDat(spartial_index)...
                                                 -xAnalytical(spartial_index))/yDat(spartial_index)...
                                                );
            %error_sum += spartial_error(spartial_index)
            spartial_index += 1; 
        endwhile

        fprintf(fid,['%g  \n'],spartial_error);
        fprintf(fid,['\n']);
        avg_error = sum(spartial_error)/size(xAnalytical,2)
        fprintf(fid,['%g  \n'],avg_error);

    endif


%%%%% PLOT 1 START %%%%%%%

    if(plot_solution)

        %plot(r_dimless,x,'ro','markersize',markersize,'LineWidth',lineWidth)
        plot(r_dimless,xAnalytical/x_core,'r-','markersize',markersize,'LineWidth',lineWidth)
        hold on   
        plot(r_dimless,xAnalyticalSimple/x_core,'k--','markersize',markersize,'LineWidth',lineWidth)
        plot(xDat,yDat/c_core,'ro','markersize',markersize,'LineWidth',lineWidth);
        hold off;  

        set(gca,'Fontsize',14)
        hLeg= legend('with convection','no convection','ParScale');
        set(hLeg,'box','off')
        set(hLeg,'location','SouthWest','FontSize',stdTextFontSize)
        xlabel('r/R','FontSize',labelFontSize)
        ylabel('x_{i} / x_{eq}','FontSize',labelFontSize)

        %xlim([0,gridpoints]);
        xlim([0,1]);

        if (x_core > x_surface)
            ylim([0,1.05]);
        endif

        if (x_core < x_surface)
            ylim([0,1.05]);
        endif
        print('-dpng ', 'Stefan_diffu_1D_spherical.png')

%%%%% PLOT 2 START %%%%%%%

      if(get_numerical_solution_convective_flux)
        figure
        plot(r_dimless,flux,'r-','markersize',markersize,'LineWidth',lineWidth);
        hold on;
        plot(r_dimless,fluxSimple,'k--','markersize',markersize,'LineWidth',lineWidth);
        hold on;
        plot(xDat,yDat_flux,'ro','markersize',markersize,'LineWidth',lineWidth);
        hold off
        set(gca,'Fontsize',14)
        hLeg= legend('with convection','no convection','ParScale');
        set(hLeg,'box','off')
        set(hLeg,'location','NorthEast','FontSize',stdTextFontSize)
        xlabel('r/R','FontSize',labelFontSize)
        ylabel('n_{g,tot} [kmol/m^2 s]','FontSize',labelFontSize)
        title(['int.flow w. conv.: ', num2str(integralFlowRate),', int.flow no conv.', num2str(integralFlowRateSimple)]);
        xlim([0,1]);
#        ylim([1e-12,1e-3]);
        print('-dpng ', 'Flux.png')
      endif
    endif

%%%%% PLOT 3 START %%%%%%%
     if(plot_solution)

     %plot the time evolution of the liquid phase fraction
      dryingDat = struct;
      fileToLoadName='liquidSpeciesPhaseAv.json';
      myFiles = dir(['../*.*']);
      validData = 0;


        for iDir=1:size(myFiles,1)
            fileToLoad = ['../', myFiles(iDir).name,'/',fileToLoadName];
            if(myFiles(iDir).isdir && exist(fileToLoad))
            validData = validData +1;

            [xDat, yDat, misc] = jsonGetParScaleData(fileToLoad, ...
                                                    'data', ...
                                                     particle_id);
            dryingDat.avLiqContent(validData) = yDat;
            dryingDat.time(validData)         = str2num(myFiles(iDir).name);
            end
        end
    figure
    plot(dryingDat.time(1:plotSkip:end),dryingDat.avLiqContent(1:plotSkip:end), ...
         'bo', ...
         'MarkerFaceColor','b','markersize',markersize,'LineWidth',lineWidth)
    hold on
    set(gca,'Fontsize',14)
    hLeg= legend('av. liq. content');
    set(hLeg,'box','off')
    set(hLeg,'location','NorthEast','FontSize',stdTextFontSize)
    xlabel('t [s]','FontSize',labelFontSize)
    ylabel('<\epsilon_{liq}>','FontSize',labelFontSize)
    print('-dpng ', 'avLiqContent.png')
    


    endif
%%%%% PLOT END %%%%%%%
