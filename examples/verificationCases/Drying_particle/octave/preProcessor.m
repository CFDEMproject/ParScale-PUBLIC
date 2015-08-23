% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% Main preprocessing file for a ParScale Simulation
% Copyright: Stefan Radl, TU Graz, 2015
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
clear
clc
more off
close all 
PARSCALE_SRC_DIR = getenv('PASCAL_SRC_DIR');
if(isempty(PARSCALE_SRC_DIR))
    error('The user has not set the environment variable "PASCAL_SRC_DIR". Please do so, e.g., in your .bashrc file, to use this octave script.')
end
addpath([PARSCALE_SRC_DIR,'/../examples/octaveFunctions']);

runDir = '../';

#radialIndex     = 0.50*20;
radialIndex     = 0.50*40;
#radialIndex     = 0.50*60;
porosity        = 0.90;

iVariable=1;
file{iVariable} =           '0/gasPhaseFraction.json';
variableName{iVariable} =   'data.all';
value{iVariable}        =   porosity;
valueBC{iVariable}      =   0;
upToId{iVariable}       =   999;

iVariable=iVariable+1;
file{iVariable} =           '0/gasPhaseFraction.json';
variableName{iVariable} =   'data.all';
value{iVariable}        =   1e-4;
valueBC{iVariable}      =   0;
upToId{iVariable}       =   radialIndex;

iVariable=iVariable+1;
file{iVariable} =           '0/liquidPhaseFraction.json';
variableName{iVariable} =   'data.all';
value{iVariable}        =   0.0;
valueBC{iVariable}      =   0;
upToId{iVariable}       =   999;

iVariable=iVariable+1;
file{iVariable} =           '0/liquidPhaseFraction.json';
variableName{iVariable} =   'data.all';
value{iVariable}        =   porosity;
valueBC{iVariable}      =   0;
upToId{iVariable}       =   radialIndex;

iVariable=iVariable+1;
file{iVariable} =           '0/species.json';
variableName{iVariable} =   'data.all';
value{iVariable}        =   0.00;
valueBC{iVariable}      =   1e-4;
upToId{iVariable}       =   999;

iVariable=iVariable+1;
file{iVariable} =           '0/liquidSpecies.json';
variableName{iVariable} =   'data.all';
value{iVariable}        =   16.757; %1e99; %set to a high value in order to avoid update of liquid phase
valueBC{iVariable}      =   16.757; %1e99; %set to a high value in order to avoid update of liquid phase
upToId{iVariable}       =   999;


iVariable=iVariable+1;
file{iVariable} =           '0/heat.json';
variableName{iVariable} =   'data.all';
value{iVariable}        =   350;
valueBC{iVariable}      =   value{iVariable};
upToId{iVariable}       =   999;

iVariable=iVariable+1;
file{iVariable} =           'settings/liquidToGasEvaporation.json';
variableName{iVariable} =   'vaporPressureModel.evaporationRateConstant';
value{iVariable}        =   3e3;   %This value is critical for the stability of the algorithm! Decrease in case you need more stability
valueBC{iVariable}      =   -1; %not relevant
upToId{iVariable}       =   999;
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% END USER INPUT
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

for iter=1:iVariable
    disp(['applying change to: ', runDir,file{iter}])
    [errorValue] = jsonChangeVariable([runDir,'/',file{iter}], variableName{iter}, value{iter}, valueBC{iter}, upToId{iter})
   
    %on some systems the parsing/writing of JSON files does not work! Must manually change the particle id 
    sysCommand=["sed -i 's/I/1/g' ", runDir,'/',file{iter}];
    system(sysCommand);
end
