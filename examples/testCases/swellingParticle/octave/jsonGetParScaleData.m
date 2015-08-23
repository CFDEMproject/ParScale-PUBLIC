function [xDat, yDat, misc] = jsonGetParScaleData(fileName, fieldName, particleIndex) 

dat         = loadjson(fileName);
field_data  = getfield(dat,fieldName);
cell_data   = struct2cell(field_data);
misc.number_particles = size(cell_data,1);
    
particle_data = cell_data{particleIndex};
yDat          = particle_data;

gridpoints = size(yDat,2);
gridpoint_xDat=1/(gridpoints-2);

m=0;    
while m<gridpoints
 xDat(m+1) = m*gridpoint_xDat;
 m=m+1;
end

end
