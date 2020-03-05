%% Gramm plots for the combination of 3 sinks
% This code requires manual variable input that allows flexibility in which
% sinks from which parameter are combined, whether or not to add muscimol,
% and whether to have GS or ST based tuning

clear
sink = {'IVE','VbE','VIaE'};    %The sinks we would like to compare - currently limited to 3 sinks exactly
para = 1;                       %1 = 'SinkRMS',2 = 'SinkPeakAmp', 3 = 'tempSinkRMS', 4 = 'Sinkonset', 5 = 'SinkPeakLate'
addmusc = 1;                    %1 = include muscimol and 0 = do not include muscimol
which = 2;                      %1 = GS_based and 2 = ST_based

%% Start
cd('D:\MyCode\Dynamic_CSD_Analysis');
warning('OFF');
dbstop if error

home = pwd;
addpath(genpath(home));

%variables for plots
rowofnans = NaN(1,7);
ticks = {'a -3' 'b -2' 'c -1' 'd BF' 'e +1' 'f +2' 'g +3'};
based = {'GS_based','ST_based'};
Parameter = {'SinkRMS','SinkPeakAmp','tempSinkRMS','Sinkonset','SinkPeakLate'};

%% Load in the appropriate files
cd DATA;cd output;
load('AnesthetizedPre_Data.m_Threshold_0.25_Zscore_0_binned_1_mirror_0.mat');
Anesthetized = Data; clear Data;
load('Awake10dB_Data.m_Threshold_0.25_Zscore_0_binned_1_mirror_0.mat')
Awake = Data; clear Data;
load('Muscimol_Data.m_Threshold_0.25_Zscore_0_binned_1_mirror_0.mat')
Muscimol = Data; clear Data;

cd(home);cd figs;
mkdir('Group Plots Gramm Combi'); cd('Group Plots Gramm Combi')

%SINK 1
An_data1 = vertcat(Anesthetized.(based{which}).(Parameter{para}).(sink{1})(:,5:11));
An_datamean1 = nanmean(An_data1,1);
An_datasem1 = nanstd(An_data1,1)/sqrt(sum(~isnan(An_datamean1)));

Aw_data1 = vertcat(Awake.(based{which}).(Parameter{para}).(sink{1})(:,4:10),rowofnans,rowofnans);
Aw_datamean1 = nanmean(Aw_data1,1);
Aw_datasem1 = nanstd(Aw_data1,1)/sqrt(sum(~isnan(Aw_datamean1)));

M_data1 = vertcat(Muscimol.(based{which}).(Parameter{para}).(sink{1})(:,5:11));
M_datamean1 = nanmean(M_data1,1);
M_datasem1 = nanstd(M_data1,1)/sqrt(sum(~isnan(M_datamean1)));

%SINK 2
An_data2 = vertcat(Anesthetized.(based{which}).(Parameter{para}).(sink{2})(:,5:11));
An_datamean2 = nanmean(An_data2,1);
An_datasem2 = nanstd(An_data2,1)/sqrt(sum(~isnan(An_datamean2)));

Aw_data2 = vertcat(Awake.(based{which}).(Parameter{para}).(sink{2})(:,4:10),rowofnans,rowofnans);
Aw_datamean2 = nanmean(Aw_data2,1);
Aw_datasem2 = nanstd(Aw_data2,1)/sqrt(sum(~isnan(Aw_datamean2)));

M_data2 = vertcat(Muscimol.(based{which}).(Parameter{para}).(sink{2})(:,5:11));
M_datamean2 = nanmean(M_data2,1);
M_datasem2 = nanstd(M_data2,1)/sqrt(sum(~isnan(M_datamean2)));

%SINK 3 
An_data3 = vertcat(Anesthetized.(based{which}).(Parameter{para}).(sink{3})(:,5:11));
An_datamean3 = nanmean(An_data3,1);
An_datasem3 = nanstd(An_data3,1)/sqrt(sum(~isnan(An_datamean3)));

Aw_data3 = vertcat(Awake.(based{which}).(Parameter{para}).(sink{3})(:,4:10),rowofnans,rowofnans);
Aw_datamean3 = nanmean(Aw_data3,1);
Aw_datasem3 = nanstd(Aw_data3,1)/sqrt(sum(~isnan(Aw_datamean3)));

M_data3 = vertcat(Muscimol.(based{which}).(Parameter{para}).(sink{3})(:,5:11));
M_datamean3 = nanmean(M_data3,1);
M_datasem3 = nanstd(M_data3,1)/sqrt(sum(~isnan(M_datamean3)));


%% Create Appropriate Structure
plotdata = struct;

tone = ticks';
Angrammdata1 = An_data1(1,:)';
Awgrammdata1 = Aw_data1(1,:)';
Mgrammdata1 = M_data1(1,:)';
Angrammdata2 = An_data2(1,:)';
Awgrammdata2 = Aw_data2(1,:)';
Mgrammdata2 = M_data2(1,:)';
Angrammdata3 = An_data3(1,:)';
Awgrammdata3 = Aw_data3(1,:)';
Mgrammdata3 = M_data3(1,:)';

Angroup = {'Anesthetized' 'Anesthetized' 'Anesthetized' 'Anesthetized' 'Anesthetized' 'Anesthetized' 'Anesthetized'}; %7 to match the amount of ticks
An = Angroup';
Awgroup = {'Awake' 'Awake' 'Awake' 'Awake' 'Awake' 'Awake' 'Awake'};
Aw = Awgroup';
Mgroup = {'Muscimol' 'Muscimol' 'Muscimol' 'Muscimol' 'Muscimol' 'Muscimol' 'Muscimol'};
M = Mgroup';

sinklist1pre = {['a' sink{1}], ['a' sink{1}], ['a' sink{1}], ['a' sink{1}], ['a' sink{1}], ['a' sink{1}], ['a' sink{1}]}; %{'IVE' 'IVE' 'IVE' 'IVE' 'IVE' 'IVE' 'IVE'};
sinklist1 = sinklist1pre';
sinklist2pre = {['b' sink{2}], ['b' sink{2}], ['b' sink{2}], ['b' sink{2}], ['b' sink{2}], ['b' sink{2}], ['b' sink{2}]}; %{'VbE' 'VbE' 'VbE' 'VbE' 'VbE' 'VbE' 'VbE'};
sinklist2 = sinklist2pre';
sinklist3pre = {['c' sink{3}], ['c' sink{3}], ['c' sink{3}], ['c' sink{3}], ['c' sink{3}], ['c' sink{3}], ['c' sink{3}]}; %{'VIaE' 'VIaE' 'VIaE' 'VIaE' 'VIaE' 'VIaE' 'VIaE'};
sinklist3 = sinklist3pre';

for istack = 1:size(An_data1,1)-1
    
    tone = vertcat(tone, ticks');
    
    Angrammdata1 = vertcat(Angrammdata1, An_data1(istack+1,:)');
    Awgrammdata1 = vertcat(Awgrammdata1, Aw_data1(istack+1,:)');
    Mgrammdata1 = vertcat(Mgrammdata1, M_data1(istack+1,:)');
    Angrammdata2 = vertcat(Angrammdata2, An_data2(istack+1,:)');
    Awgrammdata2 = vertcat(Awgrammdata2, Aw_data2(istack+1,:)');
    Mgrammdata2 = vertcat(Mgrammdata2, M_data2(istack+1,:)');
    Angrammdata3 = vertcat(Angrammdata3, An_data3(istack+1,:)');
    Awgrammdata3 = vertcat(Awgrammdata3, Aw_data3(istack+1,:)');
    Mgrammdata3 = vertcat(Mgrammdata3, M_data3(istack+1,:)');
        
    An = vertcat(An, Angroup');
    Aw = vertcat(Aw, Awgroup');
    M = vertcat(M, Mgroup');
    sinklist1 = vertcat(sinklist1, sinklist1pre');
    sinklist2 = vertcat(sinklist2, sinklist2pre');
    sinklist3 = vertcat(sinklist3, sinklist3pre');
    
end

if addmusc == 1 %in this condition, the Muscimol is added to the graph
    plotdata.tone = vertcat(tone, tone, tone, tone, tone, tone, tone, tone, tone);
    plotdata.data = vertcat(Angrammdata1, Awgrammdata1, Mgrammdata1, ...
        Angrammdata2, Awgrammdata2, Mgrammdata2, ...
        Angrammdata3, Awgrammdata3, Mgrammdata3);
    plotdata.group = vertcat(An, Aw, M, An, Aw, M, An, Aw, M);
    plotdata.sink = vertcat(sinklist1, sinklist1, sinklist1, sinklist2, ...
        sinklist2, sinklist2, sinklist3, sinklist3, sinklist3);
    
    fullmean1 = nanmean(vertcat(An_data1,Aw_data1,M_data1),1)';
    fullmean2 = nanmean(vertcat(An_data2,Aw_data2, M_data2),1)';
    fullmean3 = nanmean(vertcat(An_data3,Aw_data3,M_data3),1)';
else %in this condition, only Awake and Anesthetized are compared
    plotdata.tone = vertcat(tone, tone, tone, tone, tone, tone);
    plotdata.data = vertcat(Angrammdata1, Awgrammdata1, Angrammdata2, Awgrammdata2, ...
        Angrammdata3, Awgrammdata3);
    plotdata.group = vertcat(An, Aw, An, Aw, An, Aw);
    plotdata.sink = vertcat(sinklist1, sinklist1, sinklist1, sinklist2, sinklist2, sinklist2, sinklist3, sinklist3, sinklist3);
    
    
    fullmean1 = nanmean(vertcat(An_data1,Aw_data1),1)';
    fullmean2 = nanmean(vertcat(An_data2,Aw_data2),1)';
    fullmean3 = nanmean(vertcat(An_data3,Aw_data3),1)';
end

%for mean graph % this is actually not really in use anymore but there you
%go, easy enough to add back in if wanted (code for it is commented at the
%end of the script)
fullmean = vertcat(fullmean1, fullmean2, fullmean3);
ticksx2 = vertcat(ticks',ticks', ticks');
sink123 = vertcat(sinklist1pre',sinklist2pre', sinklist3pre');


clear g
figure('Position',[100 100 1000 550]);

g(1,1)=gramm('x',plotdata.tone,'y',plotdata.data, 'color', plotdata.group); %,'y',cars.Acceleration,'color',cars.Cylinders,'subset',cars.Cylinders~=3 & cars.Cylinders~=5
g(1,1).facet_grid([],plotdata.sink); %,'scale','free'
g(1,1).stat_summary('type','sem','geom','errorbar'); %mean and sem shown
g(1,1).stat_summary('type','sem','geom','point'); %mean and sem shown
g(1,1).stat_summary('type','sem','geom','area'); %mean and sem shown
g(1,1).set_layout_options('Position',[0 0 0.8 1],...
    'legend_pos',[0.62 0.77 0.15 0.15],... %We detach the legend from the plot and move it to the top right
    'margin_height',[0.1 0.02],...
    'margin_width',[0.1 0.02],...
    'redraw',false);
g(1,1).set_text_options('base_size',11,'legend_title_scaling', 1,'facet_scaling', 1)
g(1,1).axe_property('Ygrid','on'); %,'YLim',[0 0.0015]
if para < 4
    g(1,1).set_names('x','Tone','y','mV/mm²','color','Group');
else
    g(1,1).set_names('x','Tone','y','ms','color','Group');
%     g(1,1).axe_property('Ygrid','on','YLim',[200 350]);
end
g(1,1).set_color_options('map','matlab');
g(1,1).axe_property('XTickLabel',{'-3', '-2', '-1', 'BF', '+1', '+2', '+3'})


g(2,1)=gramm('x',plotdata.data,'color',plotdata.group);
g(2,1).facet_grid(plotdata.sink,[],'row_labels',false); %'scale','free'
g(2,1).set_layout_options('Position',[0.8 0 0.2 1],...
    'legend',false,...
    'margin_height',[0.1 0.02],...
    'margin_width',[0.02 0.05],...
    'redraw',false);
g(2,1).set_text_options('base_size',11);
g(2,1).stat_bin('geom','overlaid_bar','fill','transparent'); %histogram
g(2,1).coord_flip();
g(2,1).axe_property('XTickLabel',''); %,'XLim',[0 0.0015]
g(2,1).set_color_options('map','matlab');


% g.set_title([(based{which}) ' ' (sink{1}) ', ' (sink{2}) ', and ' (sink{3}) ' for ' (Parameter{para})]);
g.draw();
g.export('file_name',[(based{which}) '_' (sink{1}) '_' (sink{2}) '_' (sink{3}) '_' (Parameter{para})], 'file_type','png');
g.export('file_name',[(based{which}) '_' (sink{1}) '_' (sink{2}) '_' (sink{3}) '_' (Parameter{para})], 'file_type','pdf');
close all




% g(1,1)=gramm('x',ticksx2,'y',fullmean);
% g(1,1).facet_grid(sinkIVEVaE,[])
% g(1,1).geom_line();
% g(1,1).set_layout_options('Position',[0 0.8 0.8 0.2],... %Set the position in the figure (as in standard 'Position' axe property)
%     'legend',false,... % No need to display legend for side histograms
%     'margin_height',[0.02 0.05],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
%     'margin_width',[0.1 0.02],...
%     'redraw',false); %We deactivate automatic redrawing/resizing so that the axes stay aligned according to the margin options
% g(1,1).set_names('x','','y','Mean','row','');
% g(1,1).axe_property('XTickLabel',''); % We deactivate the ticks
% g(1,1).set_line_options('base_size',3)
% g(1,1).set_color_options('map',[0 0 0])

