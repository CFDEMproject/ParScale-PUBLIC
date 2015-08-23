clear
clc
close all
more off

data2Plot = 'heat.json';
j=1;
gridPoint2Plot{j} = 1; symbol{j}='ro'; faceC{j}='r'; j=j+1;
gridPoint2Plot{j} = 15; symbol{j}='k<'; faceC{j}='k'; j=j+1;
gridPoint2Plot{j} = 99; symbol{j}='bv'; faceC{j}='none';j=j+1;


%Scan dir for data
allDirs = dir('.'); isub = [allDirs(:).isdir]; allDirs=allDirs(isub);
allDirs(1)='';allDirs(1)=''; %cut-off current and parent dit
myFile = dir(['*/',data2Plot]);

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
  
%axis([1 28 600 650])  
leg=legend('fluid','particle_center','particle_middle','particle_surf');
set(leg,'Location','South');
legend boxoff
xlabel('time (s)','Fontsize',13)
ylabel('Temperature [K]','Fontsize',13)
set(gca,'Fontsize',10);

xlhand = get(gca,'xlabel');ylhand = get(gca,'ylabel');
set(ylhand,'Position',get(ylhand,'Position') - [0.01 0 0])
set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 21 18])

print('-depsc2','plot.eps')
