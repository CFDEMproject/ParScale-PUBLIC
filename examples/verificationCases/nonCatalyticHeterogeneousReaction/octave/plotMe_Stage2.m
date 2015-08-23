%Function to calculate a second stage reaction
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

%%%%% USER INPUT START %%%%%%%
run('parametersHetReaction.m')

#timeToPlot     = '1.100010';% time dir to plot
#timeToPlot     ='1.500000';
#timeToPlot     = '2.000000';
#timeToPlot     = '3.000000';
#timeToPlot     = '3.001000';
#timeToPlot     = '3.002000';
timeToPlot     = '4.000000';
#timeToPlot     = '5.300000';  
#timeToPlot     = '5.000000';
#timeToPlot     = '6.000000'
   
%%%%% USER INPUT END %%%%%%%
ThieleMod    = R * sqrt((a*k_vIntr*Cs0)/(DeIA))
N_sh          = k_mA * R/D_eA                % mass transfer/diffusion
D_ratio       = DeIA/D_eA;
Omega_vc =  Wen_omegaStage2(ThieleMod, N_sh, D_ratio);
disp('Time for complete 1st and 2nd stage');
t_complete_both_stages = Omega_vc/(k_vIntr*CA0)

disp('Diffusion time scale');
diffusionTimeScale = R^2/DeIA

time_ParScale = str2num(timeToPlot);    %simulation time 

omega_vNew    = time_ParScale * k_vIntr * CA0   % dimensionless time 


radius_interval = radius_interval/R
radius_interval_absolut = radius_interval
%%%%% MAIN PROGRAM %%%%%%%

xi_m = Wen_xi_m_time(omega_vNew,D_ratio,ThieleMod,N_sh)
concentration_reaction_front =...
      Wen_concentration_reaction_front (D_ratio,xi_m,N_sh,ThieleMod,CA0)

%Calculate concentration of fluid phase at reaction front;
% not dimensionless!!!
%place loop over xi_wanted here
xi_wanted      = 1e-6;      % dimensionless coordinate at which the concentration profile is wanted
i=1;
while (xi_wanted <=1)
    %Calculate the concentration of fluid phase at wanted positon xi_wanted (=r/R)
    %depending if xi > xi_m_time or xi< xi_m_time different approach is used
    if (xi_wanted > xi_m)
        [concentration_fluid(i)] = ...
        Wen_concentration_fluid_outer (xi_m,N_sh,concentration_reaction_front,xi_wanted,CA0);
    elseif (xi_wanted  <= xi_m)
        [concentration_fluid(i)] = ...
        Wen_concentration_fluid_inner (xi_m,ThieleMod,concentration_reaction_front,xi_wanted);
    else
        disp ('This should not happen - check input parameter xi! Has to be 0 < xi_wanted < 1') 
        return
        close all;
    end

    %Calculate the concentration of solid phase at wanted positon xi_wanted (=r/R)  
    [concentration_solid(i)] = Wen_concentration_solid (xi_m, ThieleMod,Cs0,xi_wanted);
    concentration_solid(i);
    xi_plot(i) =  xi_wanted;
    xi_wanted = xi_wanted + radius_interval; 
    absolut_radius(i) = i*  radius_interval_absolut;
    i=i+1;
end

plot(xi_plot,concentration_fluid./CA0,'-','markersize',markersize,'LineWidth',lineWidth);
hold on;

plot(xi_plot,concentration_solid./Cs0,'r--','markersize',markersize,'LineWidth',lineWidth);
hold on;



%plot(absolut_radius,concentration_solid,'bo');

%%%%%%%%%%%%%%%%%%% ParScale Plot START %%%%%%%%%%%%%%%%%%%
i=1;
%Fluid
fileName                = ['../',timeToPlot,'/speciesA.json'];
[valuesX,particle_data] = jsonGetParScaleData(fileName,'data',i) ;
plot(valuesX,particle_data./CA0,'bo', ...
'MarkerFaceColor','b','markersize',markersize,'LineWidth',lineWidth)
hold on

%Solid
fileName                = ['../',timeToPlot,'/speciesSolid.json'];
[valuesX,particle_data] = jsonGetParScaleData(fileName,'data',i) ;
plot(valuesX,particle_data./Cs0,'ro','markersize',markersize,'LineWidth',lineWidth)
hold on
%%%%%%%%%%%%%%%%%%% ParScale Plot END %%%%%%%%%%%%%%%%%%%
set(gca,'Fontsize',14)
hLeg=legend('Analytical Fluid','Analytical Solid','ParScale Fluid','ParScale Solid',0);
set(hLeg,'box','off');
set(hLeg,'location','NorthOutside','FontSize',stdTextFontSize);
xlabel('r/R','FontSize',labelFontSize)
ylabel('C_i / C_{i,0} , i=A,S','FontSize',labelFontSize)


xlim([0,1.05]);
ylim([0,1]);
title(['t = ', num2str(time_ParScale,'%.3f'),' [s]']);

print('-dpng ','concProfile_Stage2')
