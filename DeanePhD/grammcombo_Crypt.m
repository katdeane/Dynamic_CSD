%% Gramm plots for the combination of 4 sinks
% This code requires manual variable input that allows flexibility in which
% sinks from which parameter are combined. This is specifically for the
% cryptochrome project and should be copied, changed, and renamed for
% different projects

%Input:     D:\MyCode\Dynamic_CSD_Analysis\DATA\Output -> *_Tuning_Avg.mat (manually) 
%Output:    Group tuning curves in figs\Group Plots Gramm Combi\

clear
sink = {'I_II','IV','V','VI'};      %The 4 sinks we would like to compare 
para = 3;                           %1 = 'SinkRMS',2 = 'SinkPeakAmp', 3 = 'SinkPeakLate', 4 = 'Sinkonset'


%% Start
cd('D:\MyCode\Dynamic_CSD');
warning('OFF');
dbstop if error

if ~exist('homedir','var')
    if exist('D:\MyCode\Dynamic_CSD','dir') == 7
        cd('D:\MyCode\Dynamic_CSD');
    elseif exist('D:\Dynamic_CSD_Analysis','dir') == 7
        cd('D:\Dynamic_CSD_Analysis');
    elseif exist('C:\Users\kedea\Documents\Dynamic_CSD_Analysis','dir') == 7
        cd('C:\Users\kedea\Documents\Dynamic_CSD_Analysis')
    end
    
    homedir = pwd;
    addpath(genpath(homedir));
end

%loop through measurements: 
Measurement = {'Pre_4','preAMtono_1','preAMtono_4','preCLtono_1',...
    'preCLtono_4','CLtono_1','AMtono_1'};
%variables for plots
rowofnans = NaN(1,5);
ticks = {'b -2' 'c -1' 'd BF' 'e +1' 'f +2'};
Parameter = {'SinkRMS','SinkPeakAmp','SinkPeakLate','Sinkonset'};

%% Load in the appropriate files
cd DATA; cd Output;
load('KIC_Tuning_Avg.mat');
Control = Tuning; clear Tuning; % no virus
load('KIT_Tuning_Avg.mat');
Treated = Tuning; clear Tuning; % virus w/ cryptochrome
load('KIV_Tuning_Avg.mat');
VControl = Tuning; clear Tuning; % virus control

cd(home); cd figs;
mkdir('Group Plots Gramm Combi'); cd('Group Plots Gramm Combi')

for imeas = 1:length(Measurement)
    %cut out the data associated with each sink at each parameter for each
    %group (BF is center column of the output matrix with 3 columns to each
    %side)
    
    %SINK 1
    C_data1 = vertcat(Control.(sink{1}).(Measurement{imeas}).(Parameter{para})(:,6:10),rowofnans,rowofnans,rowofnans,rowofnans);
    T_data1 = vertcat(Treated.(sink{1}).(Measurement{imeas}).(Parameter{para})(:,6:10));
    V_data1 = vertcat(VControl.(sink{1}).(Measurement{imeas}).(Parameter{para})(:,6:10),rowofnans,rowofnans,rowofnans,rowofnans);
    
    %SINK 2
    C_data2 = vertcat(Control.(sink{2}).(Measurement{imeas}).(Parameter{para})(:,6:10),rowofnans,rowofnans,rowofnans,rowofnans);
    T_data2 = vertcat(Treated.(sink{2}).(Measurement{imeas}).(Parameter{para})(:,6:10));
    V_data2 = vertcat(VControl.(sink{2}).(Measurement{imeas}).(Parameter{para})(:,6:10),rowofnans,rowofnans,rowofnans,rowofnans);
    
    %SINK 3
    C_data3 = vertcat(Control.(sink{3}).(Measurement{imeas}).(Parameter{para})(:,6:10),rowofnans,rowofnans,rowofnans,rowofnans);
    T_data3 = vertcat(Treated.(sink{3}).(Measurement{imeas}).(Parameter{para})(:,6:10));
    V_data3 = vertcat(VControl.(sink{3}).(Measurement{imeas}).(Parameter{para})(:,6:10),rowofnans,rowofnans,rowofnans,rowofnans);
    
    %SINK 4
    C_data4 = vertcat(Control.(sink{4}).(Measurement{imeas}).(Parameter{para})(:,6:10),rowofnans,rowofnans,rowofnans,rowofnans);
    T_data4 = vertcat(Treated.(sink{4}).(Measurement{imeas}).(Parameter{para})(:,6:10));
    V_data4 = vertcat(VControl.(sink{4}).(Measurement{imeas}).(Parameter{para})(:,6:10),rowofnans,rowofnans,rowofnans,rowofnans);
    
    %SINK 1 normalized
    C_data1(C_data1 == 0) = 1; T_data1(T_data1 == 0) = 1; V_data1(V_data1 == 0) = 1;
    C_norm1 = C_data1./(C_data1(:,3));
    T_norm1 = T_data1./T_data1(:,3);
    V_norm1 = V_data1./V_data1(:,3);
    
    %SINK 2 normalized
    C_data2(C_data2 == 0) = 1; T_data2(T_data2 == 0) = 1; V_data2(V_data2 == 0) = 1;
    C_norm2 = C_data2./(C_data2(:,3));
    T_norm2 = T_data2./T_data2(:,3);
    V_norm2 = V_data2./V_data2(:,3);
    
    %SINK 3 normalized
    C_data3(C_data3 == 0) = 1; T_data3(T_data3 == 0) = 1; V_data3(V_data3 == 0) = 1;
    C_norm3 = C_data3./(C_data3(:,3));
    T_norm3 = T_data3./T_data3(:,3);
    V_norm3 = V_data3./V_data3(:,3);
    
    %SINK 4 normalized
    C_data4(C_data4 == 0) = 1; T_data4(T_data4 == 0) = 1; V_data4(V_data4 == 0) = 1;
    C_norm4 = C_data4./(C_data4(:,3));
    T_norm4 = T_data4./T_data4(:,3);
    V_norm4 = V_data4./V_data4(:,3);
    
    
    %% Create appropriate structure
    plotdata = struct;
    
    tone = ticks';
    Cgrammdata1 = C_data1(1,:)';
    Tgrammdata1 = T_data1(1,:)';
    Vgrammdata1 = V_data1(1,:)';
    Cgrammdata2 = C_data2(1,:)';
    Tgrammdata2 = T_data2(1,:)';
    Vgrammdata2 = V_data2(1,:)';
    Cgrammdata3 = C_data3(1,:)';
    Tgrammdata3 = T_data3(1,:)';
    Vgrammdata3 = V_data3(1,:)';
    Cgrammdata4 = C_data4(1,:)';
    Tgrammdata4 = T_data4(1,:)';
    Vgrammdata4 = V_data4(1,:)';
    Cngrammdata1 = C_norm1(1,:)';
    Tngrammdata1 = T_norm1(1,:)';
    Vngrammdata1 = V_norm1(1,:)';
    Cngrammdata2 = C_norm2(1,:)';
    Tngrammdata2 = T_norm2(1,:)';
    Vngrammdata2 = V_norm2(1,:)';
    Cngrammdata3 = C_norm3(1,:)';
    Tngrammdata3 = T_norm3(1,:)';
    Vngrammdata3 = V_norm3(1,:)';
    Cngrammdata4 = C_norm4(1,:)';
    Tngrammdata4 = T_norm4(1,:)';
    Vngrammdata4 = V_norm4(1,:)';
    
    Cgroup = {'Control' 'Control' 'Control' 'Control' 'Control'}; %5 to match the amount of ticks
    C = Cgroup';
    Tgroup = {'Treated' 'Treated' 'Treated' 'Treated' 'Treated'};
    T = Tgroup';
    Vgroup = {'VControl' 'VControl' 'VControl' 'VControl' 'VControl'};
    V = Vgroup';
    
    sinklist1pre = {['a' sink{1}], ['a' sink{1}], ['a' sink{1}], ['a' sink{1}], ['a' sink{1}]};
    sinklist1 = sinklist1pre';
    sinklist2pre = {['b' sink{2}], ['b' sink{2}], ['b' sink{2}], ['b' sink{2}], ['b' sink{2}]};
    sinklist2 = sinklist2pre';
    sinklist3pre = {['c' sink{3}], ['c' sink{3}], ['c' sink{3}], ['c' sink{3}], ['c' sink{3}]};
    sinklist3 = sinklist3pre';
    sinklist4pre = {['d' sink{4}], ['d' sink{4}], ['d' sink{4}], ['d' sink{4}], ['d' sink{4}]};
    sinklist4 = sinklist4pre';
    
    for istack = 1:size(C_data1,1)-1
        
        tone = vertcat(tone, ticks');
        
        Cgrammdata1 = vertcat(Cgrammdata1, C_data1(istack+1,:)');
        Tgrammdata1 = vertcat(Tgrammdata1, T_data1(istack+1,:)');
        Vgrammdata1 = vertcat(Vgrammdata1, V_data1(istack+1,:)');
        Cgrammdata2 = vertcat(Cgrammdata2, C_data2(istack+1,:)');
        Tgrammdata2 = vertcat(Tgrammdata2, T_data2(istack+1,:)');
        Vgrammdata2 = vertcat(Vgrammdata2, V_data2(istack+1,:)');
        Cgrammdata3 = vertcat(Cgrammdata3, C_data3(istack+1,:)');
        Tgrammdata3 = vertcat(Tgrammdata3, T_data3(istack+1,:)');
        Vgrammdata3 = vertcat(Vgrammdata3, V_data3(istack+1,:)');
        Cgrammdata4 = vertcat(Cgrammdata4, C_data4(istack+1,:)');
        Tgrammdata4 = vertcat(Tgrammdata4, T_data4(istack+1,:)');
        Vgrammdata4 = vertcat(Vgrammdata4, V_data4(istack+1,:)');
        Cngrammdata1 = vertcat(Cngrammdata1, C_norm1(istack+1,:)');
        Tngrammdata1 = vertcat(Tngrammdata1, T_norm1(istack+1,:)');
        Vngrammdata1 = vertcat(Vngrammdata1, V_norm1(istack+1,:)');
        Cngrammdata2 = vertcat(Cngrammdata2, C_norm2(istack+1,:)');
        Tngrammdata2 = vertcat(Tngrammdata2, T_norm2(istack+1,:)');
        Vngrammdata2 = vertcat(Vngrammdata2, V_norm2(istack+1,:)');
        Cngrammdata3 = vertcat(Cngrammdata3, C_norm3(istack+1,:)');
        Tngrammdata3 = vertcat(Tngrammdata3, T_norm3(istack+1,:)');
        Vngrammdata3 = vertcat(Vngrammdata3, V_norm3(istack+1,:)');
        Cngrammdata4 = vertcat(Cngrammdata4, C_norm4(istack+1,:)');
        Tngrammdata4 = vertcat(Tngrammdata4, T_norm4(istack+1,:)');
        Vngrammdata4 = vertcat(Vngrammdata4, V_norm4(istack+1,:)');
        
        C = vertcat(C, Cgroup');
        T = vertcat(T, Tgroup');
        V = vertcat(V, Vgroup');
        sinklist1 = vertcat(sinklist1, sinklist1pre');
        sinklist2 = vertcat(sinklist2, sinklist2pre');
        sinklist3 = vertcat(sinklist3, sinklist3pre');
        sinklist4 = vertcat(sinklist4, sinklist4pre');
        
    end
    
    
    plotdata.tone = vertcat(tone, tone, tone, tone, tone, tone, tone, tone, tone, tone, tone, tone);
    plotdata.data = vertcat(Cgrammdata1, Tgrammdata1, Vgrammdata1, ...
        Cgrammdata2, Tgrammdata2, Vgrammdata2, Cgrammdata3, Tgrammdata3, Vgrammdata3,...
        Cgrammdata4, Tgrammdata4, Vgrammdata4);
    plotdata.normdata = vertcat(Cngrammdata1, Tngrammdata1, Vngrammdata1, ...
        Cngrammdata2, Tngrammdata2, Vngrammdata2, Cngrammdata3, Tngrammdata3, Vngrammdata3,...
        Cngrammdata4, Tngrammdata4, Vngrammdata4);
    plotdata.group = vertcat(C, T, V, C, T, V, C, T, V, C, T, V);
    plotdata.sink = vertcat(sinklist1, sinklist1, sinklist1, sinklist2, ...
        sinklist2, sinklist2, sinklist3, sinklist3, sinklist3, sinklist4, sinklist4, sinklist4);
    
    
    
    Tune = struct2table(plotdata); %a smarter person with more time than me would make a table first
    
    %% grammplot for sink tuning curves
    
    clear g
    figure('Position',[100 100 1000 550]);
    
    g(1,1)=gramm('x',Tune.tone,'y',Tune.data, 'color', Tune.group); %,'y',cars.Acceleration,'color',cars.Cylinders,'subset',cars.Cylinders~=3 & cars.Cylinders~=5
    g(1,1).facet_grid([],Tune.sink,'column_labels',true,'row_labels',true); %,'scale','free'
    % g(1,1).geom_point(); %to see each point, removed the YLim (highest points are closer to 2.5
    g(1,1).stat_summary('type','sem','geom','errorbar'); %mean and sem shown
    g(1,1).stat_summary('type','sem','geom','point'); %mean and sem shown
    g(1,1).stat_summary('type','sem','geom','area'); %mean and sem shown
    g(1,1).set_layout_options('Position',[0 0 0.8 1],...
        'legend_pos',[0.62 0.77 0.15 0.15],... %We detach the legend from the plot and move it to the top right
        'margin_height',[0.1 0.02],...
        'margin_width',[0.1 0.02],...
        'redraw',false);
    g(1,1).set_text_options('base_size',11,'legend_title_scaling', 1,'facet_scaling', 1)
    %
    if para == 1
        g(1,1).set_names('x','Tone','y','mV/mm²','color','Group');
        %     g(1,1).axe_property('Ygrid','on','YLim',[0.005 0.04]); %,'YLim',[0 0.0015]
    elseif para == 2
        g(1,1).set_names('x','Tone','y','mV/mm²','color','Group');
        g(1,1).axe_property('Ygrid','on','YLim',[0 0.005]); %,'YLim',[0 0.0015]
    elseif para == 3
        g(1,1).set_names('x','Tone','y','ms','color','Group');
        g(1,1).axe_property('Ygrid','on','YLim',[0 400]);
    else
        g(1,1).set_names('x','Tone','y','ms','color','Group');
    end
    g(1,1).set_color_options('map','matlab');
    g(1,1).axe_property('XTickLabel',{'-2', '-1', 'BF', '+1', '+2'})
    
    
    g.set_title([(sink{1}) ', ' (sink{2}) ', ' (sink{3})  ', and ' (sink{4}) ' for ' (Parameter{para}) '_' (Measurement{imeas})]);
    g.draw();
    g.export('file_name',[(sink{1}) '_' (sink{2}) '_' (sink{3}) '_' (sink{4}) '_' (Parameter{para}) '_' (Measurement{imeas})], 'file_type','png');
    g.export('file_name',[(sink{1}) '_' (sink{2}) '_' (sink{3}) '_' (sink{4}) '_' (Parameter{para}) '_' (Measurement{imeas})], 'file_type','pdf');
    close all
    
    %% grammplot for sink tuning curves NORMALIZED
    
    clear g
    figure('Position',[100 100 1000 550]);
    
    g(1,1)=gramm('x',Tune.tone,'y',Tune.normdata, 'color', Tune.group); %,'y',cars.Acceleration,'color',cars.Cylinders,'subset',cars.Cylinders~=3 & cars.Cylinders~=5
    g(1,1).facet_grid([],Tune.sink); %,'scale','free'
    % g(2,1).geom_point(); %to see each point, removed the YLim (highest points are closer to 2.5
    g(1,1).stat_summary('type','sem','geom','errorbar'); %mean and sem shown
    g(1,1).stat_summary('type','sem','geom','point'); %mean and sem shown
    g(1,1).stat_summary('type','sem','geom','area'); %mean and sem shown
    g(1,1).set_layout_options('Position',[0 0 0.8 1],...
        'legend_pos',[0.62 0.77 0.15 0.15],... %We detach the legend from the plot and move it to the top right
        'margin_height',[0.1 0.02],...
        'margin_width',[0.1 0.02],...
        'redraw',false);
    g(1,1).set_text_options('base_size',11,'legend_title_scaling', 1,'facet_scaling', 1)
    g(1,1).set_names('x','Tone','y','%','color','Group');
    
    g(1,1).set_color_options('map','matlab');
    g(1,1).axe_property('XTickLabel',{'-2', '-1', 'BF', '+1', '+2'})
    
    
    g.set_title([(sink{1}) ', ' (sink{2}) ', ' (sink{3}) ', and ' (sink{4}) ' for ' (Parameter{para}) '_' (Measurement{imeas})]);
    g.draw();
    g.export('file_name',['Normalized ' (sink{1}) '_' (sink{2}) '_' (sink{3}) '_' (sink{4}) '_' (Parameter{para}) '_' (Measurement{imeas})], 'file_type','png');
    g.export('file_name',['Normalized ' (sink{1}) '_' (sink{2}) '_' (sink{3}) '_' (sink{4}) '_' (Parameter{para}) '_' (Measurement{imeas})], 'file_type','pdf');
    close all
    %% got mirror?
    
    % stacks of 7 to sort the lists later (7 animals) - group
    Cgroup = {'Control' 'Control' 'Control' 'Control' 'Control' 'Control' ...
        'Control'};
    Cgroup = Cgroup';
    Tgroup = {'Treated' 'Treated' 'Treated' 'Treated' 'Treated' 'Treated' ...
        'Treated'};
    Tgroup = Tgroup';
    Vgroup = {'VControl' 'VControl' 'VControl' 'VControl' 'VControl' 'VControl' ...
        'VControl'};
    Vgroup = Vgroup';
    % stacks of 7 to sort the lists later (7 animals) - frequency
    BF = {'aBF' 'aBF' 'aBF' 'aBF' 'aBF' 'aBF' 'aBF'};
    BF = BF';
    one = {'b1' 'b1' 'b1' 'b1' 'b1' 'b1' 'b1'};
    one = one';
    two = {'c2' 'c2' 'c2' 'c2' 'c2' 'c2' 'c2'};
    two = two';
    % stacks of 7 to sort the lists later (7 animals) - frequency
    sink1 = {sink{1}, sink{1}, sink{1}, sink{1}, sink{1}, sink{1}, sink{1}};
    sink1 = sink1';
    sink2 = {sink{2}, sink{2}, sink{2}, sink{2}, sink{2}, sink{2}, sink{2}};
    sink2 = sink2';
    sink3 = {sink{3}, sink{3}, sink{3}, sink{3}, sink{3}, sink{3}, sink{3}};
    sink3 = sink3';
    sink4 = {sink{4}, sink{4}, sink{4}, sink{4}, sink{4}, sink{4}, sink{4}};
    sink4 = sink4';
    
    %% sink 1
    %CONTROL
    % pull out the + and - of each side (+-1), average them, normalize to BF
    mirC1_2 = vertcat(C_norm1(:,1),C_norm1(:,5));
    mirC1_1 = vertcat(C_norm1(:,2),C_norm1(:,4));
    mirC1_BF = vertcat(C_norm1(:,3),C_norm1(:,3));
    %concatonate them vertically (BF*gsize, +-1*gsize,...)
    mirC_data1 = vertcat(mirC1_BF,mirC1_1,mirC1_2);
    C1 = vertcat(Cgroup, Cgroup, Cgroup, Cgroup, Cgroup, Cgroup);
    freqC = vertcat(BF, BF, one, one, two, two);
    
    %TREATED
    % pull out the + and - of each side (+-1), average them, normalize to BF
    mirT1_2 = vertcat(T_norm1(:,1),T_norm1(:,5));
    mirT1_1 = vertcat(T_norm1(:,2),T_norm1(:,4));
    mirT1_BF = vertcat(T_norm1(:,3),T_norm1(:,3));
    %concatonate them vertically (BF*gsize, +-1*gsize,...)
    mirT_data1 = vertcat(mirT1_BF,mirT1_1,mirT1_2);
    T1 = vertcat(Tgroup, Tgroup, Tgroup, Tgroup, Tgroup, Tgroup);
    freqT = vertcat(BF, BF, one, one, two, two);
    
    %VIRUS CONTROL
    % pull out the + and - of each side (+-1), average them, normalize to BF
    mirV1_2 = vertcat(V_norm1(:,1),V_norm1(:,5));
    mirV1_1 = vertcat(V_norm1(:,2),V_norm1(:,4));
    mirV1_BF = vertcat(V_norm1(:,3),V_norm1(:,3));
    %concatonate them vertically (BF*gsize, +-1*gsize,...)
    mirV_data1 = vertcat(mirV1_BF,mirV1_1,mirV1_2);
    V1 = vertcat(Vgroup, Vgroup, Vgroup, Vgroup, Vgroup, Vgroup);
    freqV = vertcat(BF, BF, one, one, two, two);
    
    %sanity check: all following variables should be equal in length
    mirData1 = vertcat(mirC_data1,mirT_data1,mirV_data1);
    Groups1 = vertcat(C1,T1,V1);
    Freq1 = vertcat(freqC,freqT,freqV);
    Sinks1 = vertcat(sink1,sink1,sink1,sink1,sink1,sink1,sink1,sink1,sink1,sink1,sink1,sink1,sink1,sink1,sink1,sink1,sink1,sink1);
    %% sink 2
    %CONTROL
    % pull out the + and - of each side (+-1), average them, normalize to BF
    mirC2_2 = vertcat(C_norm2(:,1),C_norm2(:,5));
    mirC2_1 = vertcat(C_norm2(:,2),C_norm2(:,4));
    mirC2_BF = vertcat(C_norm2(:,3),C_norm2(:,3));
    %concatonate them vertically (BF*gsize, +-1*gsize,...)
    mirC_data2 = vertcat(mirC2_BF,mirC2_1,mirC2_2);
    C2 = vertcat(Cgroup, Cgroup, Cgroup, Cgroup, Cgroup, Cgroup);
    freqC = vertcat(BF, BF, one, one, two, two);
    
    %TREATED
    % pull out the + and - of each side (+-1), average them, normalize to BF
    mirT2_2 = vertcat(T_norm2(:,1),T_norm2(:,5));
    mirT2_1 = vertcat(T_norm2(:,2),T_norm2(:,4));
    mirT2_BF = vertcat(T_norm2(:,3),T_norm2(:,3));
    %concatonate them vertically (BF*gsize, +-1*gsize,...)
    mirT_data2 = vertcat(mirT2_BF,mirT2_1,mirT2_2);
    T2 = vertcat(Tgroup, Tgroup, Tgroup, Tgroup, Tgroup, Tgroup);
    freqT = vertcat(BF, BF, one, one, two, two);
    
    %VIRUS CONTROL
    % pull out the + and - of each side (+-1), average them, normalize to BF
    mirV2_2 = vertcat(V_norm2(:,1),V_norm2(:,5));
    mirV2_1 = vertcat(V_norm2(:,2),V_norm2(:,4));
    mirV2_BF = vertcat(V_norm2(:,3),V_norm2(:,3));
    %concatonate them vertically (BF*gsize, +-1*gsize,...)
    mirV_data2 = vertcat(mirV2_BF,mirV2_1,mirV2_2);
    V2 = vertcat(Vgroup, Vgroup, Vgroup, Vgroup, Vgroup, Vgroup);
    freqV = vertcat(BF, BF, one, one, two, two);
    
    %sanity check: all following variables should be equal in length
    mirData2 = vertcat(mirC_data2,mirT_data2,mirV_data2);
    Groups2 = vertcat(C2,T2,V2);
    Freq2 = vertcat(freqC,freqT,freqV);
    Sinks2 = vertcat(sink2,sink2,sink2,sink2,sink2,sink2,sink2,sink2,sink2,sink2,sink2,sink2,sink2,sink2,sink2,sink2,sink2,sink2);
    %% sink 3
    %CONTORL
    % pull out the + and - of each side (+-1), average them, normalize to BF
    mirC3_2 = vertcat(C_norm3(:,1),C_norm3(:,5));
    mirC3_1 = vertcat(C_norm3(:,2),C_norm3(:,4));
    mirC3_BF = vertcat(C_norm3(:,3),C_norm3(:,3));
    %concatonate them vertically (BF*gsize, +-1*gsize,...)
    mirC_data3 = vertcat(mirC3_BF,mirC3_1,mirC3_2);
    C3 = vertcat(Cgroup, Cgroup, Cgroup, Cgroup, Cgroup, Cgroup);
    freqC = vertcat(BF, BF, one, one, two, two);
    
    %TREATED
    % pull out the + and - of each side (+-1), average them, normalize to BF
    mirT3_2 = vertcat(T_norm3(:,1),T_norm3(:,5));
    mirT3_1 = vertcat(T_norm3(:,2),T_norm3(:,4));
    mirT3_BF = vertcat(T_norm3(:,3),T_norm3(:,3));
    %concatonate them vertically (BF*gsize, +-1*gsize,...)
    mirT_data3 = vertcat(mirT3_BF,mirT3_1,mirT3_2);
    T3 = vertcat(Tgroup, Tgroup, Tgroup, Tgroup, Tgroup, Tgroup);
    freqT = vertcat(BF, BF, one, one, two, two);
    
    %VIRUS CONTROL
    % pull out the + and - of each side (+-1), average them, normalize to BF
    mirV3_2 = vertcat(V_norm3(:,1),V_norm3(:,5));
    mirV3_1 = vertcat(V_norm3(:,2),V_norm3(:,4));
    mirV3_BF = vertcat(V_norm3(:,3),V_norm3(:,3));
    %concatonate them vertically (BF*gsize, +-1*gsize,...)
    mirV_data3 = vertcat(mirV3_BF,mirV3_1,mirV3_2);
    V3 = vertcat(Vgroup, Vgroup, Vgroup, Vgroup, Vgroup, Vgroup);
    freqV = vertcat(BF, BF, one, one, two, two);
    
    %sanity check: all following variables should be equal in length
    mirData3 = vertcat(mirC_data3,mirT_data3,mirV_data3);
    Groups3 = vertcat(C3,T3,V3);
    Freq3 = vertcat(freqC,freqT,freqV);
    Sinks3 = vertcat(sink3,sink3,sink3,sink3,sink3,sink3,sink3,sink3,sink3,sink3,sink3,sink3,sink3,sink3,sink3,sink3,sink3,sink3);
    %% sink 4
    %CONTROL
    % pull out the + and - of each side (+-1), average them, normalize to BF
    mirC4_2 = vertcat(C_norm4(:,1),C_norm4(:,5));
    mirC4_1 = vertcat(C_norm4(:,2),C_norm4(:,4));
    mirC4_BF = vertcat(C_norm4(:,3),C_norm4(:,3));
    %concatonate them vertically (BF*gsize, +-1*gsize,...)
    mirC_data4 = vertcat(mirC4_BF,mirC4_1,mirC4_2);
    C4 = vertcat(Cgroup, Cgroup, Cgroup, Cgroup, Cgroup, Cgroup);
    freqC = vertcat(BF, BF, one, one, two, two);
    
    %TREATED
    % pull out the + and - of each side (+-1), average them, normalize to BF
    mirT4_2 = vertcat(T_norm4(:,1),T_norm4(:,5));
    mirT4_1 = vertcat(T_norm4(:,2),T_norm4(:,4));
    mirT4_BF = vertcat(T_norm4(:,3),T_norm4(:,3));
    %concatonate them vertically (BF*gsize, +-1*gsize,...)
    mirT_data4 = vertcat(mirT4_BF,mirT4_1,mirT4_2);
    T4 = vertcat(Tgroup, Tgroup, Tgroup, Tgroup, Tgroup, Tgroup);
    freqT = vertcat(BF, BF, one, one, two, two);
    
    %VIRUS CONTROL
    % pull out the + and - of each side (+-1), average them, normalize to BF
    mirV4_2 = vertcat(V_norm4(:,1),V_norm4(:,5));
    mirV4_1 = vertcat(V_norm4(:,2),V_norm4(:,4));
    mirV4_BF = vertcat(V_norm4(:,3),V_norm4(:,3));
    %concatonate them vertically (BF*gsize, +-1*gsize,...)
    mirV_data4 = vertcat(mirV4_BF,mirV4_1,mirV4_2);
    V4 = vertcat(Vgroup, Vgroup, Vgroup, Vgroup, Vgroup, Vgroup);
    freqV = vertcat(BF, BF, one, one, two, two);
    
    %sanity check: all following variables should be equal in length
    mirData4 = vertcat(mirC_data4,mirT_data4,mirV_data4);
    Groups4 = vertcat(C4,T4,V4);
    Freq4 = vertcat(freqC,freqT,freqV);
    Sinks4 = vertcat(sink4,sink4,sink4,sink4,sink4,sink4,sink4,sink4,sink4,sink4,sink4,sink4,sink4,sink4,sink4,sink4,sink4,sink4);
    
    %%
    mirData = vertcat(mirData1,mirData2,mirData3,mirData4);
    Groups = vertcat(Groups1,Groups2,Groups3,Groups4);
    Freq = vertcat(Freq1,Freq2,Freq3,Freq4);
    Sinks = vertcat(Sinks1,Sinks2,Sinks3,Sinks4);
    
    clear g
    figure('Position',[100 100 1000 550]);
    
    g(1,1)=gramm('x',Freq,'y',mirData, 'color', Groups ); %,'y',cars.Acceleration,'color',cars.Cylinders,'subset',cars.Cylinders~=3 & cars.Cylinders~=5
    g(1,1).facet_grid([],Sinks,'row_labels',false); %,'scale','free'
    % g(2,1).geom_point(); %to see each point, removed the YLim (highest points are closer to 2.5
    g(1,1).stat_summary('type','sem','geom','errorbar'); %mean and sem shown
    g(1,1).stat_summary('type','sem','geom','point'); %mean and sem shown
    g(1,1).stat_summary('type','sem','geom','area'); %mean and sem shown
    g(1,1).set_layout_options('Position',[0 0 0.8 1],...
        'legend_pos',[0.62 0.77 0.15 0.15],... %We detach the legend from the plot and move it to the top right
        'margin_height',[0.1 0.02],...
        'margin_width',[0.1 0.02],...
        'redraw',false);
    g(1,1).set_text_options('base_size',11,'legend_title_scaling', 1,'facet_scaling', 1)
    g(1,1).set_names('x','Tone','y','%','color','Group');
    % g(1,1).axe_property('Ygrid','on','YLim',[0 1.5]);
    
    g(1,1).set_color_options('map','matlab');
    g(1,1).axe_property('XTickLabel',{'BF', '+-1', '+-2', '+-3'})
    
    
    g.set_title([(sink{1}) ', ' (sink{2}) ', ' (sink{3}) ', and ' (sink{4}) ' for ' (Parameter{para}) '_' (Measurement{imeas})]);
    g.draw();
    g.export('file_name',['Normalized mirror ' (sink{1}) '_' ...
        (sink{2}) '_' (sink{3}) '_' (sink{4}) '_' (Parameter{para}) '_' (Measurement{imeas})], 'file_type','png');
    g.export('file_name',['Normalized mirror ' (sink{1}) '_' ...
        (sink{2}) '_' (sink{3}) '_' (sink{4}) '_' (Parameter{para}) '_' (Measurement{imeas})], 'file_type','pdf');
    close all
    
end


