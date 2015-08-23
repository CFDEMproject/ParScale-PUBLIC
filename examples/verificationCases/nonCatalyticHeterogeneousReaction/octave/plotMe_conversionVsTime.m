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
relativePath  = '../';
fileToLoad    = 'speciesSolidAv.json';
verificationFileName = '../verificationStats/errorConversion.json';
particleIndex = 1;
nPlotPoints   = 100;
plotSkip      = 1;
   
%%%%% USER INPUT END %%%%%%%
ThieleMod     = R * sqrt((a*k_vIntr*Cs0)/(DeIA))
N_sh          = k_mA * R/D_eA;                % mass transfer/diffusion
D_ratio       = DeIA/D_eA;

disp('Diffusion time scale');
diffusionTimeScale = R.^2./DeIA

disp('Reaction time scale');
reactionTimeScale = 1./(k_vIntr*CA0)

%%%%% MAIN PROGRAM %%%%%%%
%1 - Get data from ParScale, this will be unordered
myFiles = dir([relativePath,'*.*']);
validData = 0;

for iDir=1:size(myFiles,1)
    if(myFiles(iDir).isdir)
        validData = validData +1;
        [xDat, yDat, misc] = jsonGetParScaleData([relativePath, myFiles(iDir).name,'/',fileToLoad], ...
                                                'data', ...
                                                 particleIndex);
        parScaleConversion(validData) = 1 - (yDat./Cs0);
        parScaleTime(validData)       = str2num(myFiles(iDir).name);
    end
end

parScaleEndTime         = max(parScaleTime);
parScaleEndTimePosition = parScaleTime==parScaleEndTime;
parScaleFinalConversion = parScaleConversion(parScaleEndTimePosition);
parScaleEndTimeDimLess  = parScaleEndTime./reactionTimeScale ;
parScaleDimLessTime     = parScaleTime./reactionTimeScale ;

%save the key results
results.Cs0             = Cs0;
results.conversion      = parScaleFinalConversion;
results.time            = parScaleEndTime;
disp('final conversion:');
results.conversion 
disp('final solids consumption:');
results.consumption = results.conversion*Cs0
disp('after:');
results.time
jSonDat = savejson('',results, 'result.json');

%Plot the analytical solution
omegaPlot    = linspace(0, parScaleEndTimeDimLess, nPlotPoints);
for iP=1:nPlotPoints
    xi_m(iP)       = Wen_xi_m_time(omegaPlot(iP),D_ratio,ThieleMod,N_sh);
    conversion(iP) = Wen_conversion (ThieleMod, xi_m(iP), N_sh, D_ratio, omegaPlot(iP));
end
plot(omegaPlot,conversion,'b-', ...
'MarkerFaceColor','b','markersize',markersize,'LineWidth',lineWidth)
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute the error
for iP=1:validData
    xi_m_verification(iP)       = Wen_xi_m_time(parScaleDimLessTime(iP),D_ratio,ThieleMod,N_sh);
    conversion_verification(iP) = Wen_conversion (ThieleMod, xi_m_verification(iP), N_sh, D_ratio, parScaleDimLessTime(iP));
end

verification.nSamples      = length(conversion_verification);
verification.ErrorPosition = parScaleDimLessTime';
verification.Error         = (conversion_verification - parScaleConversion)';
verification.ErrorMax      = max(abs(verification.Error));
verification.ErrorMean     = sqrt(mean(verification.Error.^2));
disp('Mean error in conversion:');
verification.ErrorMean
disp('Max error in conversion:');
verification.ErrorMax
jSonDat = savejson('',verification, verificationFileName);

%%%%%%%%%%%%%%%%%%% ParScale Plot START %%%%%%%%%%%%%%%%%%%
plot(parScaleDimLessTime(1:plotSkip:end),parScaleConversion(1:plotSkip:end), ...
     'bo', ...
      'MarkerFaceColor','b','markersize',markersize,'LineWidth',lineWidth)
hold on

%%%%%%%%%%%%%%%%%%% ParScale Plot END %%%%%%%%%%%%%%%%%%%
set(gca,'Fontsize',14)
hLeg=legend('Analytical','ParScale',0);
set(hLeg,'box','off');
set(hLeg,'location','NorthOutside','FontSize',stdTextFontSize);
xlabel('t / t_{react}','FontSize',labelFontSize)
ylabel('X_s','FontSize',labelFontSize)


xlim([0,parScaleEndTimeDimLess]);
ylim([0,1]);

print('-dpng ','conversion')
