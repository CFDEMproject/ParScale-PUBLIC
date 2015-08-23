clear
clc
close all
more off



%%%% GRAPHICS/OUTPUT SETTINGS %%%%%
run('formatting.m')

%%%%% USER INPUT START %%%%%%%
run('parametersHeatConduction.m')
relativePath  = '../';
verificationFileName = '../verificationStats/errorConversion.json';
particleIndex = 1;
nPlotPoints   = 100;
plotSkip      = 2;

iToLoad = 1;
fileToLoad{iToLoad}    = 'liquidSpeciesPhaseAv.json';iToLoad=iToLoad+1;
fileToLoad{iToLoad}    = 'liquidPhaseFraction.json'; 
   
%%%%% USER INPUT END %%%%%%%

%%%%% MAIN PROGRAM %%%%%%%
%1 - Get data from ParScale, this will be unordered
myFiles = dir([relativePath,'*.*']);
validData = 0;

for iDir=1:size(myFiles,1)
    if(myFiles(iDir).isdir)
        validData = validData +1;
        parScaleTime(validData)       = str2num(myFiles(iDir).name);

        for iLoad=1:iToLoad
          [rawX, rawY, misc] = jsonGetParScaleData([relativePath, myFiles(iDir).name,'/',fileToLoad{iLoad}], ...
                                                'data', ...
                                                 particleIndex);

          inData{iLoad}.x(validData,1:length(rawX)) = rawX;
          inData{iLoad}.y(validData,1:length(rawY)) = rawY;
        end

    end
end

%Sort the data
[parScaleTime, iKey] = sort (parScaleTime);

for iLoad=1:iToLoad
    inData{iLoad}.x = inData{iLoad}.x(iKey,:);
    inData{iLoad}.y = inData{iLoad}.y(iKey,:);
end

parScaleEndTime         = max(parScaleTime);
parScaleEndTimePosition = parScaleTime==parScaleEndTime;
parScaleEndTimeDimLess  = parScaleEndTime./dryingTimeScale ;
parScaleDimLessTime     = parScaleTime./dryingTimeScale ;


%2 - Plot mean liquid concentration
iLoad=1
plot(parScaleTime,inData{iLoad}.y,'ro','markersize',markersize,'LineWidth',lineWidth)
hold on
        
set(gca,'Fontsize',14)
hLeg= legend(['\epsilon_{liq}, i = ', num2str(particleIndex)],0);
set(hLeg,'box','off')
set(hLeg,'location','NorthOutside','FontSize',stdTextFontSize)
xlabel('t [s]','FontSize',labelFontSize)
ylabel('\epsilon_{liq}','FontSize',labelFontSize)

#xlim([0,1.05]);
ylim([0,max(inData{iLoad}.y)]);

print('-dpng ', [fileToLoad{iToLoad}(1:end-5),'.png'])

