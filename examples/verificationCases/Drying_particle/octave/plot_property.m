clear
clc
close all
more off
run('formatting.m');
format long g;

%%% Read in times %%%
%c=["300.001000";"400.001000";"350.001000";"450.001000"];
c=["400.000100";"370.000100";"360.000100";"350.000100";"410.000100";"420.000100";"390.000100"];
time_to_plot = cellstr (c);

%%% Property to PLot %%%
prop_name = "liquidPhaseFraction.json"
%prop_name = "species.json"

%% Select particle IDs %%%
particle_id = 1;

%%% Loop over time and get values %%%
d=1;
%size(c)
while d <= size (c)
    fileName  = ['../',time_to_plot{d,1},'/',prop_name]
    [xDat, yDat, misc] = jsonGetParScaleData(fileName,'data',particle_id); 
    %sum_ = sum(yDat)/length(yDat);
    sum_ = sum(yDat)
    %plot(xDat,yDat/(sum_),'r-','markersize',markersize,'LineWidth',lineWidth);
    plot(xDat,yDat,'r-','markersize',markersize,'LineWidth',lineWidth);
    hold on;
    d=d+1;
endwhile

%%% Main end %%% 

set(gca,'Fontsize',14)
%hLeg=legend('Analytical Fluid','Analytical Solid','ParScale Fluid','ParScale Solid',0);
%set(hLeg,'box','off');
%set(hLeg,'location','NorthOutside','FontSize',stdTextFontSize);
xlabel('r/R','FontSize',labelFontSize)
ylabel('\phi_{liq}','FontSize',labelFontSize)
xlim ([0.0 1.0])
%ylim ([0 1.0])
print('-dpng ','Property_evolution')

