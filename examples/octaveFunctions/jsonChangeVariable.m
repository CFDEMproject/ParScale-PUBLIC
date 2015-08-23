%Copyright: Stefan Radl, IPPT, 2015
%Function to change a variable (or a sub-variable) in a JSON file
%Can set also intra-particle data (but just for all particles!)
function [errorValue] = jsonChangeVariable(file, variableName, value, valueBC, upToId)

%Read and analyze variable name
inData      = loadjson(file);
dotPosition = findstr(variableName,'.');
setAllPart  = findstr(variableName,'.all');

%in case there is no subvariable, things are easy!
errorValue = 0;
if(isempty(dotPosition))
    currVariable = getfield(inData, variableName);
    outData = setfield(inData,variableName, value);
elseif (~isempty(setAllPart))
    firstName  = variableName(1:dotPosition-1);
    currVariable = getfield(inData, firstName);
    disp(['Setting ',  num2str(length(currVariable)),' particle(s).']);
    for [val, key]  = currVariable
        newValues   = getfield(inData,firstName,key); %copy existing one
        if(~isnan(value))
            if(upToId>length(newValues)-1)
                upToId=length(newValues)-1;
            end
            newValues(1:upToId) = value;
        end
        if(~isnan(valueBC))
            newValues(end) = valueBC;
        end
        outData = setfield(inData,firstName,key,newValues);
    end
else
    firstName  = variableName(1:dotPosition-1);
    secondName = variableName(dotPosition+1:end);
    currVariable = getfield(inData,firstName,secondName);
    outData      = setfield(inData,firstName,secondName,value);
end

%Save to disk
savejson('',outData, file);

end
