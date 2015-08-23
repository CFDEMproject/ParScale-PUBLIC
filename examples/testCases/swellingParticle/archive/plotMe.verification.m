clear
clc
close all
more off

refenrenceDataName = 'referenceSolution';
data2Plot = 'heat.json';
j=1;
gridPoint2Plot{j} = 1;  refDataName{j} = '.center';     symbol{j}='ro'; faceC{j}='r'; j=j+1;
gridPoint2Plot{j} = 15; refDataName{j} = '.inBetween';  symbol{j}='k<'; faceC{j}='k'; j=j+1;
gridPoint2Plot{j} = 99; refDataName{j} = '.surface';    symbol{j}='bv'; faceC{j}='none';j=j+1;


%Scan dir for data
allDirs = dir('.'); isub = [allDirs(:).isdir]; allDirs=allDirs(isub);
allDirs(1)='';allDirs(1)=''; %cut-off current and parent dit

raw = struct;
validTime = 0;
for iFile=1:size(allDirs,1)
    myFile = dir([allDirs(iFile).name,'/',data2Plot]);
    if(size(myFile,1)>0)
        %disp('foundFile!')
        dat = loadjson([allDirs(iFile).name,'/',myFile(1).name]);
        cell_data = struct2cell(getfield(dat,'data'));
        number_particles = size(cell_data,1);
        i=1; %just load first particle!
        particle_data = cell_data{i};
        if(isnan(particle_data)(1)==0)
            validTime = validTime + 1;
            raw.data{validTime}       = particle_data;
            raw.time(validTime)       = str2num(allDirs(iFile).name);
          
            %Extract data to plot
            lastPoint = size(raw.data{validTime},2);
            raw.fluid(validTime)      = raw.data{validTime}(lastPoint);
            for j=1:size(gridPoint2Plot,2)
                if(gridPoint2Plot{j}>=lastPoint)
                    gridPoint2Plot{j} = lastPoint-1;
                end
                raw.plot{j}.dat(validTime) = raw.data{validTime}(gridPoint2Plot{j});
            end
            raw.gridpoints(validTime) = size(particle_data,2);
        end
    end
end


plot(raw.time, raw.fluid,'kd-');
hold on
for j=1:size(gridPoint2Plot,2)
    plot(raw.time, raw.plot{j}.dat,symbol{j},'MarkerFaceColor', faceC{j});
end

%plot reference solution
errorSum  = 0;
nElements = 0;
maxValue = 0;minValue = 9e99;
for j=1:size(gridPoint2Plot,2)
    refData{j} = load([refenrenceDataName,refDataName{j},'.dat']);
    plot(refData{j}(:,1), refData{j}(:,2), [symbol{j}(1),'-']);
    %compute and save error
    for iTime=1:size(raw.time,2)
        refSolutionInterp = interp1(refData{j}(:,1), refData{j}(:,2),raw.time(iTime));
        errorSum=errorSum+(refSolutionInterp-raw.plot{j}.dat(iTime)).^2;
        nElements=nElements+1;
        maxValue=max(maxValue,refSolutionInterp);
        minValue=min(minValue,refSolutionInterp);
    end
end

 

leg=legend('fluid','particle_center','particle_middle','particle_surf');
set(leg,'Location','South');
legend boxoff
xlabel('time (s)','Fontsize',13)
ylabel('Temperature [K]','Fontsize',13)
set(gca,'Fontsize',10);

xlhand = get(gca,'xlabel');ylhand = get(gca,'ylabel');
set(ylhand,'Position',get(ylhand,'Position') - [0.01 0 0])
set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 21 18])


print('-depsc2','plot.verification.eps')


%Error Statistics
referenceValue = maxValue - minValue;
errorSum=sqrt(errorSum/nElements);
disp(['mean absolute error : ', num2str(errorSum)])
disp(['reference value     : ', num2str(referenceValue)])
disp(['mean relative error : ', num2str(errorSum/referenceValue)])
fID=fopen('stats.verification..dat','w')
fprintf(fID,'{ \n');
fprintf(fID,' "mean absolute error" : %g, \n',errorSum);
fprintf(fID,' "referenceValue"      : %g, \n',referenceValue);
fprintf(fID,' "mean relative error" : %g  \n',errorSum/referenceValue);
fprintf(fID,'} \n');
fclose(fID);
