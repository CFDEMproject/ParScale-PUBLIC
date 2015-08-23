clear
clc
close all
more off

 

format long g;

dat1 = loadjson('../1.000000/heat.json');
field_data1 =getfield(dat1,'data');
cell_data1 = struct2cell(field_data1);
particle_data1 = cell_data1{1};
gridpoints1=size(particle_data1,2);

refernce_solution = load('transient_heat_conduction.dat');
gridpoints2=size(refernce_solution,1);

error_stat = '../verificationStats/error_stat.dat';

% ========================================================================================   
%%%%% USER INPUT END %%%%%%%
% ========================================================================================


fid=fopen(error_stat,'w');

error_sum=0;
n=1;

while n <=(gridpoints2) 
    
    T_parscale   = particle_data1(n);
    T_analytical = refernce_solution(n,2);
    local_error = abs((T_analytical-T_parscale)/T_analytical);
    error_sum += abs((T_analytical-T_parscale)/T_analytical);
    fprintf(fid,['%g  \n'],local_error);
    n=n+1;
end

avg_error = error_sum/gridpoints2

fprintf(fid,['\n']);
fprintf(fid,['%g  \n'],avg_error);

if avg_error<1e-4
    printf("Average Error is ok \n");
end

if avg_error>1e-4
    printf("Average Error is too high! \n");
end
