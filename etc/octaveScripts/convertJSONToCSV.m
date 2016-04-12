Copyright: Stefan Radl, 2016, TU Graz

clc
clear
more off

################################
#USER INPUT ONLY here
################################


dirToConvert    = '../0';   #specify the directory to convert
fieldToConvert  = 'data';   #by default, all data is in the field 'data'!


################################
#END INPUT
################################

files=dir(dirToConvert);

for iFile=1:size(files,1)
    if(files(iFile).isdir==1)
        disp(['Skipping directory: ', files(iFile).name]);
    elseif (isempty(findstr(files(iFile).name, '.json')))
        disp(['Skipping file: ', files(iFile).name]);
    else #it is now ensure that this is a json file!
        fileFullName    = [dirToConvert,'/',files(iFile).name];
        fileFullNameNew = [fileFullName(1:findstr(fileFullName, '.json')-1),'.csv'];
        disp(['processing File: ', fileFullName, ', will convert to: ', fileFullNameNew]);
        data=loadjson(fileFullName);
        #Analyse input and write each sub-struct to file
        csvwrite(fileFullNameNew, cell2mat(struct2cell(data.(fieldToConvert))));

    end
end
