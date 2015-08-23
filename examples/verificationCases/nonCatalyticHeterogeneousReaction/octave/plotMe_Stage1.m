%Function to calculate a first-stage reaction
%Formulars according Wen, "Noncatalytic Heterogenous Solid Fluid Reaction Models"
%models p. 42-44

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


%%%% GRAPHICS/OUTPUT SETTINGS %%%%%
run('formatting.m')

%%%%% USER INPUT START %%%%%%%%%%%%%
run('parametersHetReaction.m')
#timeToPlot  = '0.200000';
#timeToPlot  = '0.400000'; 
#timeToPlot  = '1.000100'
timeToPlot  = '2.000100'
#timeToPlot  = '6.000000';  

         
%%%%% USER INPUT END %%%%%%%%%%%%%%%
time_ParScale = str2num(timeToPlot);    %simulation time  
ThieleMod = R * sqrt( a*k_vIntr*Cs0/ D_eA)
Omega_v   = k_vIntr * CA0 * time_ParScale
N_sh      = (k_mA * R)/(D_eA)                %reat of external mass transfer/diffusion 

Omega_vc = Wen_omegaStage1(ThieleMod, N_sh);

disp('Time for complete first stage');
t_complete_first_stage = Omega_vc/(k_vIntr*CA0)


%%%%%%%%%%%%%%ANALYTICAL SOLUTION 
    r = 1e-16;          % small number
    i = 1;

    while (r<=R)
        absolut_radius(i) = r;
        zeta   = r/R;      %dimensionless radius running variable
        dCa(i) = (1/Omega_vc) * ((sinh(zeta*ThieleMod))/(zeta*sinh(ThieleMod)));
        Cs(i)  = (1- ((sinh(zeta*ThieleMod))/(zeta*sinh(ThieleMod))) * (Omega_v / Omega_vc));

        i=i+1;    
        r = r+radius_interval;
    end     

    plot(absolut_radius/R,dCa, ...
           '-','markersize',markersize,'LineWidth',lineWidth);
    hold on;
    plot(absolut_radius/R,Cs, ...
           'r--','markersize',markersize,'LineWidth',lineWidth);
    hold on;
    
% Plot for ParScale
%%%%%%%%%%%%%%%%%% FLUID PLOT %%%%%%%%%%%%%%%%%%%%%%%%
i=1;
fileName                = ['../',timeToPlot,'/speciesA.json'];
[valuesX,particle_data] = jsonGetParScaleData(fileName,'data',i) ;
plot(valuesX,particle_data./CA0,'bo', ...
      'MarkerFaceColor','b','markersize',markersize,'LineWidth',lineWidth)
hold on

%Initial distribution
i=1;
fileName                = ['../0/speciesA.json'];
[valuesX,particle_data] = jsonGetParScaleData(fileName,'data',i) ;
plot(valuesX,particle_data./CA0,'bo','markersize',markersize,'LineWidth',lineWidth)

%%%%%%%%%%%%%%%%%% SOLID PLOT %%%%%%%%%%%%%%%%%%%%%%%%
i=1;
fileName                = ['../',timeToPlot,'/speciesSolid.json'];
[valuesX,particle_data] = jsonGetParScaleData(fileName,'data',i) ;
plot(valuesX,particle_data./Cs0,'ro','markersize',markersize,'LineWidth',lineWidth)
hold on
        
set(gca,'Fontsize',14)
hLeg= legend('Analytical Fluid','Analytical Solid','ParScale Fluid','ParScale Fluid (init.)','ParScale Solid',0);
set(hLeg,'box','off')
set(hLeg,'location','NorthOutside','FontSize',stdTextFontSize)
xlabel('r/R','FontSize',labelFontSize)
ylabel('C_i / C_{i,0} , i=A,S','FontSize',labelFontSize)

xlim([0,1.05]);
ylim([0,1]);
title(['t = ', num2str(time_ParScale,'%.3f'),' [s]']);

print('-dpng ','concProfile_Stage1')

