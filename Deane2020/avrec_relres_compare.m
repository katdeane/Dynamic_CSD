%% README

%  This script is a super simple take on checking the peak amplitude timing
%  difference between relres and avrec. It is only taking the averages as
%  they have already being generated and contains no statistics. Expansion
%  to usable graph/stats would be cool.

%% Run
clear all
cd('C:\Users\Katrina\Documents\CortXplorers\MyCode\Dynamic CSD');
warning('OFF');
dbstop if error

home = pwd;
addpath(genpath(home));

load('plotdataav.mat')
load('plotdatarel.mat')

group = {'Anesthetized','Awake','Muscimol'};
frqz = {'a -3','b -2', 'c -1','d BF', 'e +1', 'f +2', 'g +3'};
plotdata = struct; plotdata.dif = []; 
plotdata.group = []; plotdata.frqz = [];

for igroup = 1:length(group)
    for ifrqz = 1:length(frqz)
        %pull frequency placement out of list in logical 1s and 0s
        currentfrqz = strcmp(frqz(ifrqz),plotdataav.frqz);
        %shorten group list and latency list to match just the current frqz
        frqzlist = plotdataav.peaklate(currentfrqz);
        frqzlistgroup = plotdataav.group(currentfrqz);
        %shorten latency list to match just the current group
        currentgroup = strcmp(group(igroup),frqzlistgroup);
        frqzlist = frqzlist(currentgroup);
        
        %get the mean of that groups latency at that frqz
        Avrec = nanmean(frqzlist);
        
        currentfrqz = strcmp(frqz(ifrqz),plotdatarel.frqz);
        frqzlist = plotdatarel.peaklate(currentfrqz);
        frqzlistgroup = plotdatarel.group(currentfrqz);
        currentgroup = strcmp(group(igroup),frqzlistgroup);
        frqzlist = frqzlist(currentgroup);
        
        Relres = nanmean(frqzlist);
        
        %find the difference and store it
        difference = Relres - Avrec;
        plotdata.dif = horzcat(plotdata.dif, difference);
        plotdata.frqz = horzcat(plotdata.frqz, frqz(ifrqz));
    end
    %add the group list and frqz list to the plot data struct
    plotdata.group = horzcat(plotdata.group, repmat(group(igroup),1,7));
end

%% Grammplot it

cd(home); cd figs; mkdir('AV_REL difference'); cd('AV_REL difference')

% for i1 = 1:length(frqz);
clear g
figure();

g(1,1)=gramm('x',plotdata.frqz,'y',plotdata.dif, 'color', plotdata.group); %,'y',cars.Acceleration,'color',cars.Cylinders,'subset',cars.Cylinders~=3 & cars.Cylinders~=5
g(1,1).stat_summary('type','sem','geom','errorbar'); %mean and sem shown
g(1,1).stat_summary('type','sem','geom','point'); %mean and sem shown
g(1,1).stat_summary('type','sem','geom','area'); %mean and sem shown
g(1,1).set_text_options('base_size',12) 
g(1,1).set_names('x','Tone','y','mV/mm²','color','Group');
g(1,1).axe_property(); %'YLim',[230 250]
g(1,1).set_color_options('map','matlab');
g(1,1).axe_property('XTickLabel',{'-3', '-2', '-1', 'BF', '+1', '+2', '+3'})

g.set_title('Latency Difference (Relres - Avrec)');
g.draw();
g.export('file_name','Latency Difference', 'file_type','pdf');
g.export('file_name','Latency Difference', 'file_type','png');
close all;
% end