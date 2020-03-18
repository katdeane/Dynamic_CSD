%% Gramm plots for the combination of 3 sinks
% This code requires manual variable input that allows flexibility in which
% sinks from which parameter are combined, whether or not to add muscimol,
% and whether or not to have GS or ST based tuning

clear
sink = {'I_IIE','IVE','VbE','VIaE'};%The sinks we would like to compare - currently limited to 4 sinks exactly
para = 3;                           %1 = 'SinkRMS',2 = 'SinkPeakAmp', 3 = 'SinkPeakLate', 4 = 'Sinkonset'
addmusc = 1;                        %1 = include muscimol and 0 = do not include muscimol
which = 2;                          %1 = GS_based and 2 = ST_based

%% Start
cd('C:\Users\kedea\Documents\Work Stuff\Dynamic_CSD');
warning('OFF');
dbstop if error

home = pwd;
addpath(genpath(home));

%variables for plots
rowofnans = NaN(1,5);
ticks = {'b -2' 'c -1' 'd BF' 'e +1' 'f +2'};
based = {'GS_based','ST_based'};
Parameter = {'SinkRMS','SinkPeakAmp','SinkPeakLate','Sinkonset'};

%% Load in the appropriate files
cd DATA;cd Output;
load('AnesthetizedPre_Data.m_Threshold_0.25_Zscore_0_binned_1_mirror_0.mat');
Ketamine = Data; clear Data;
load('Awake10dB_Data.m_Threshold_0.25_Zscore_0_binned_1_mirror_0.mat')
Awake = Data; clear Data;
load('Muscimol_Data.m_Threshold_0.25_Zscore_0_binned_1_mirror_0.mat')
Muscimol = Data; clear Data;

cd(home);cd figs;
mkdir('Group Plots Gramm Combi'); cd('Group Plots Gramm Combi')

%cut out the data associated with each sink at each parameter for each
%group (BF is center column of the output matrix with 3 columns to each
%side)

%SINK 1
K_data1 = vertcat(Ketamine.(based{which}).(Parameter{para}).(sink{1})(:,6:10));
A_data1 = vertcat(Awake.(based{which}).(Parameter{para}).(sink{1})(:,5:9),rowofnans,rowofnans);
M_data1 = vertcat(Muscimol.(based{which}).(Parameter{para}).(sink{1})(:,6:10));

%SINK 2 
K_data2 = vertcat(Ketamine.(based{which}).(Parameter{para}).(sink{2})(:,6:10));
A_data2 = vertcat(Awake.(based{which}).(Parameter{para}).(sink{2})(:,5:9),rowofnans,rowofnans);
M_data2 = vertcat(Muscimol.(based{which}).(Parameter{para}).(sink{2})(:,6:10));

%SINK 3 
K_data3 = vertcat(Ketamine.(based{which}).(Parameter{para}).(sink{3})(:,6:10));
A_data3 = vertcat(Awake.(based{which}).(Parameter{para}).(sink{3})(:,5:9),rowofnans,rowofnans);
M_data3 = vertcat(Muscimol.(based{which}).(Parameter{para}).(sink{3})(:,6:10));

%SINK 4 
K_data4 = vertcat(Ketamine.(based{which}).(Parameter{para}).(sink{4})(:,6:10));
A_data4 = vertcat(Awake.(based{which}).(Parameter{para}).(sink{4})(:,5:9),rowofnans,rowofnans);
M_data4 = vertcat(Muscimol.(based{which}).(Parameter{para}).(sink{4})(:,6:10));

%SINK 1 normalized
K_data1(K_data1 == 0) = 1; A_data1(A_data1 == 0) = 1; M_data1(M_data1 == 0) = 1;
K_norm1 = K_data1./(K_data1(:,3));
A_norm1 = A_data1./A_data1(:,3);
M_norm1 = M_data1./M_data1(:,3);

%SINK 2 normalized
K_data2(K_data2 == 0) = 1; A_data2(A_data2 == 0) = 1; M_data2(M_data2 == 0) = 1;
K_norm2 = K_data2./(K_data2(:,3));
A_norm2 = A_data2./A_data2(:,3);
M_norm2 = M_data2./M_data2(:,3);

%SINK 3 normalized 
K_data3(K_data3 == 0) = 1; A_data3(A_data3 == 0) = 1; M_data3(M_data3 == 0) = 1;
K_norm3 = K_data3./(K_data3(:,3));
A_norm3 = A_data3./A_data3(:,3);
M_norm3 = M_data3./M_data3(:,3);

%SINK 4 normalized 
K_data4(K_data4 == 0) = 1; A_data4(A_data4 == 0) = 1; M_data4(M_data4 == 0) = 1;
K_norm4 = K_data4./(K_data4(:,3));
A_norm4 = A_data4./A_data4(:,3);
M_norm4 = M_data4./M_data4(:,3);


%% Create appropriate structure
plotdata = struct;

tone = ticks';
Kgrammdata1 = K_data1(1,:)';
Agrammdata1 = A_data1(1,:)';
Mgrammdata1 = M_data1(1,:)';
Kgrammdata2 = K_data2(1,:)';
Agrammdata2 = A_data2(1,:)';
Mgrammdata2 = M_data2(1,:)';
Kgrammdata3 = K_data3(1,:)';
Agrammdata3 = A_data3(1,:)';
Mgrammdata3 = M_data3(1,:)';
Kgrammdata4 = K_data4(1,:)';
Agrammdata4 = A_data4(1,:)';
Mgrammdata4 = M_data4(1,:)';
Kngrammdata1 = K_norm1(1,:)';
Angrammdata1 = A_norm1(1,:)';
Mngrammdata1 = M_norm1(1,:)';
Kngrammdata2 = K_norm2(1,:)';
Angrammdata2 = A_norm2(1,:)';
Mngrammdata2 = M_norm2(1,:)';
Kngrammdata3 = K_norm3(1,:)';
Angrammdata3 = A_norm3(1,:)';
Mngrammdata3 = M_norm3(1,:)';
Kngrammdata4 = K_norm4(1,:)';
Angrammdata4 = A_norm4(1,:)';
Mngrammdata4 = M_norm4(1,:)';

Kgroup = {'Ketamine' 'Ketamine' 'Ketamine' 'Ketamine' 'Ketamine'}; %5 to match the amount of ticks
K = Kgroup';
Agroup = {'Awake' 'Awake' 'Awake' 'Awake' 'Awake'};
A = Agroup';
Mgroup = {'Muscimol' 'Muscimol' 'Muscimol' 'Muscimol' 'Muscimol'};
M = Mgroup';

sinklist1pre = {['a' sink{1}], ['a' sink{1}], ['a' sink{1}], ['a' sink{1}], ['a' sink{1}]}; %{'IVE' 'IVE' 'IVE' 'IVE' 'IVE' 'IVE' 'IVE'};
sinklist1 = sinklist1pre';
sinklist2pre = {['b' sink{2}], ['b' sink{2}], ['b' sink{2}], ['b' sink{2}], ['b' sink{2}]}; %{'VbE' 'VbE' 'VbE' 'VbE' 'VbE' 'VbE' 'VbE'};
sinklist2 = sinklist2pre';
sinklist3pre = {['c' sink{3}], ['c' sink{3}], ['c' sink{3}], ['c' sink{3}], ['c' sink{3}]}; %{'VIaE' 'VIaE' 'VIaE' 'VIaE' 'VIaE' 'VIaE' 'VIaE'};
sinklist3 = sinklist3pre';
sinklist4pre = {['d' sink{4}], ['d' sink{4}], ['d' sink{4}], ['d' sink{4}], ['d' sink{4}]}; %{'VIaE' 'VIaE' 'VIaE' 'VIaE' 'VIaE' 'VIaE' 'VIaE'};
sinklist4 = sinklist4pre';

for istack = 1:size(K_data1,1)-1
    
    tone = vertcat(tone, ticks');
    
    Kgrammdata1 = vertcat(Kgrammdata1, K_data1(istack+1,:)');
    Agrammdata1 = vertcat(Agrammdata1, A_data1(istack+1,:)');
    Mgrammdata1 = vertcat(Mgrammdata1, M_data1(istack+1,:)');
    Kgrammdata2 = vertcat(Kgrammdata2, K_data2(istack+1,:)');
    Agrammdata2 = vertcat(Agrammdata2, A_data2(istack+1,:)');
    Mgrammdata2 = vertcat(Mgrammdata2, M_data2(istack+1,:)');
    Kgrammdata3 = vertcat(Kgrammdata3, K_data3(istack+1,:)');
    Agrammdata3 = vertcat(Agrammdata3, A_data3(istack+1,:)');
    Mgrammdata3 = vertcat(Mgrammdata3, M_data3(istack+1,:)');
    Kgrammdata4 = vertcat(Kgrammdata4, K_data4(istack+1,:)');
    Agrammdata4 = vertcat(Agrammdata4, A_data4(istack+1,:)');
    Mgrammdata4 = vertcat(Mgrammdata4, M_data4(istack+1,:)');
    Kngrammdata1 = vertcat(Kngrammdata1, K_norm1(istack+1,:)');
    Angrammdata1 = vertcat(Angrammdata1, A_norm1(istack+1,:)');
    Mngrammdata1 = vertcat(Mngrammdata1, M_norm1(istack+1,:)');
    Kngrammdata2 = vertcat(Kngrammdata2, K_norm2(istack+1,:)');
    Angrammdata2 = vertcat(Angrammdata2, A_norm2(istack+1,:)');
    Mngrammdata2 = vertcat(Mngrammdata2, M_norm2(istack+1,:)');
    Kngrammdata3 = vertcat(Kngrammdata3, K_norm3(istack+1,:)');
    Angrammdata3 = vertcat(Angrammdata3, A_norm3(istack+1,:)');
    Mngrammdata3 = vertcat(Mngrammdata3, M_norm3(istack+1,:)');
    Kngrammdata4 = vertcat(Kngrammdata4, K_norm4(istack+1,:)');
    Angrammdata4 = vertcat(Angrammdata4, A_norm4(istack+1,:)');
    Mngrammdata4 = vertcat(Mngrammdata4, M_norm4(istack+1,:)');
    
    K = vertcat(K, Kgroup');
    A = vertcat(A, Agroup');
    M = vertcat(M, Mgroup');
    sinklist1 = vertcat(sinklist1, sinklist1pre');
    sinklist2 = vertcat(sinklist2, sinklist2pre');
    sinklist3 = vertcat(sinklist3, sinklist3pre');
    sinklist4 = vertcat(sinklist4, sinklist4pre');
    
end

if addmusc == 1 %in this condition, the Muscimol is added to the graph
    plotdata.tone = vertcat(tone, tone, tone, tone, tone, tone, tone, tone, tone, tone, tone, tone);
    plotdata.data = vertcat(Kgrammdata1, Agrammdata1, Mgrammdata1, ...
        Kgrammdata2, Agrammdata2, Mgrammdata2, Kgrammdata3, Agrammdata3, Mgrammdata3,...
        Kgrammdata4, Agrammdata4, Mgrammdata4);
    plotdata.normdata = vertcat(Kngrammdata1, Angrammdata1, Mngrammdata1, ...
        Kngrammdata2, Angrammdata2, Mngrammdata2, Kngrammdata3, Angrammdata3, Mngrammdata3,...
        Kngrammdata4, Angrammdata4, Mngrammdata4);
    plotdata.group = vertcat(K, A, M, K, A, M, K, A, M, K, A, M);
    plotdata.sink = vertcat(sinklist1, sinklist1, sinklist1, sinklist2, ...
        sinklist2, sinklist2, sinklist3, sinklist3, sinklist3, sinklist4, sinklist4, sinklist4);
    
else %in this condition, only Awake and Ketamine are compared
    plotdata.tone = vertcat(tone, tone, tone, tone, tone, tone, tone, tone);
    plotdata.data = vertcat(Kgrammdata1, Agrammdata1,...
        Kgrammdata2, Agrammdata2, Kgrammdata3, Agrammdata3, Kgrammdata4, Agrammdata4);
    plotdata.normdata = vertcat(Kngrammdata1, Angrammdata1,...
        Kngrammdata2, Angrammdata2, Kngrammdata3, Angrammdata3, Kngrammdata4, Angrammdata4);
    plotdata.group = vertcat(K, A, K, A, K, A, K, A);
    plotdata.sink = vertcat(sinklist1, sinklist1, sinklist2, sinklist2, ...
        sinklist3, sinklist3, sinklist4, sinklist4);

end

T = struct2table(plotdata); %a smarter person with more time than me would make a table first

%% grammplot for sink tuning curves

clear g
figure('Position',[100 100 1000 550]);

g(1,1)=gramm('x',T.tone,'y',T.data, 'color', T.group); %,'y',cars.Acceleration,'color',cars.Cylinders,'subset',cars.Cylinders~=3 & cars.Cylinders~=5
g(1,1).facet_grid([],T.sink,'column_labels',true,'row_labels',true); %,'scale','free'
% g(1,1).geom_point(); %to see each point, removed the YLim (highest points are closer to 2.5
g(1,1).stat_summary('type','std','geom','errorbar'); %mean and std shown
g(1,1).stat_summary('type','std','geom','point'); %mean and std shown
g(1,1).stat_summary('type','std','geom','area'); %mean and std shown
g(1,1).set_layout_options('Position',[0 0 0.8 1],...
    'legend_pos',[0.62 0.77 0.15 0.15],... %We detach the legend from the plot and move it to the top right
    'margin_height',[0.1 0.02],...
    'margin_width',[0.1 0.02],...
    'redraw',false);
g(1,1).set_text_options('base_size',11,'legend_title_scaling', 1,'facet_scaling', 1)
% 
if para == 1
    g(1,1).set_names('x','Tone','y','mV/mm?','color','Group');
%     g(1,1).axe_property('Ygrid','on','YLim',[0.005 0.04]); %,'YLim',[0 0.0015]
elseif para == 2
    g(1,1).set_names('x','Tone','y','mV/mm?','color','Group');
    g(1,1).axe_property('Ygrid','on','YLim',[0 0.005]); %,'YLim',[0 0.0015]
elseif para == 3
    g(1,1).set_names('x','Tone','y','ms','color','Group');
    g(1,1).axe_property('Ygrid','on','YLim',[200 350]);
else
    g(1,1).set_names('x','Tone','y','ms','color','Group');
end
g(1,1).set_color_options('map','matlab');
g(1,1).axe_property('XTickLabel',{'-2', '-1', 'BF', '+1', '+2'})


g.set_title([(based{which}) ' ' (sink{1}) ', ' (sink{2}) ', ' (sink{3})  ', and ' (sink{4}) ' for ' (Parameter{para}) 'std']);
g.draw();
g.export('file_name',[(based{which}) '_' (sink{1}) '_' (sink{2}) '_' (sink{3}) '_' (sink{4}) '_' (Parameter{para}) 'std'], 'file_type','png');
g.export('file_name',[(based{which}) '_' (sink{1}) '_' (sink{2}) '_' (sink{3}) '_' (sink{4}) '_' (Parameter{para}) 'std'], 'file_type','pdf');
close all

%% grammplot for sink tuning curves NORMALIZED

clear g
figure('Position',[100 100 1000 550]);

g(1,1)=gramm('x',T.tone,'y',T.normdata, 'color', T.group); %,'y',cars.Acceleration,'color',cars.Cylinders,'subset',cars.Cylinders~=3 & cars.Cylinders~=5
g(1,1).facet_grid([],T.sink); %,'scale','free'
% g(2,1).geom_point(); %to see each point, removed the YLim (highest points are closer to 2.5
g(1,1).stat_summary('type','std','geom','errorbar'); %mean and std shown
g(1,1).stat_summary('type','std','geom','point'); %mean and std shown
g(1,1).stat_summary('type','std','geom','area'); %mean and std shown
g(1,1).set_layout_options('Position',[0 0 0.8 1],...
    'legend_pos',[0.62 0.77 0.15 0.15],... %We detach the legend from the plot and move it to the top right
    'margin_height',[0.1 0.02],...
    'margin_width',[0.1 0.02],...
    'redraw',false);
g(1,1).set_text_options('base_size',11,'legend_title_scaling', 1,'facet_scaling', 1)
g(1,1).set_names('x','Tone','y','%','color','Group');

g(1,1).set_color_options('map','matlab');
g(1,1).axe_property('XTickLabel',{'-2', '-1', 'BF', '+1', '+2'})


g.set_title([(based{which}) ' ' (sink{1}) ', ' (sink{2}) ', ' (sink{3}) ', and ' (sink{4}) ' for ' (Parameter{para}) 'std']);
g.draw();
g.export('file_name',['Normalized ' (based{which}) '_' (sink{1}) '_' (sink{2}) '_' (sink{3}) '_' (sink{4}) '_' (Parameter{para}) 'std'], 'file_type','png');
g.export('file_name',['Normalized ' (based{which}) '_' (sink{1}) '_' (sink{2}) '_' (sink{3}) '_' (sink{4}) '_' (Parameter{para}) 'std'], 'file_type','pdf');
close all
%% got mirror?

%largest group size:
if length(Ketamine.names) > length(Awake.names)
    gsize = length(Ketamine.names);
else
    gsize = length(Awake.names);
end

% stacks of 11 to sort the lists later (11 animals) - group
Kgroup = {'Ketamine' 'Ketamine' 'Ketamine' 'Ketamine' 'Ketamine' 'Ketamine' ...
    'Ketamine' 'Ketamine' 'Ketamine' 'Ketamine' 'Ketamine'};
Kgroup = Kgroup';
Agroup = {'Awake' 'Awake' 'Awake' 'Awake' 'Awake' 'Awake' 'Awake' 'Awake' 'Awake' 'Awake' 'Awake'};
Agroup = Agroup';
Mgroup = {'Muscimol' 'Muscimol' 'Muscimol' 'Muscimol' 'Muscimol' 'Muscimol' ...
    'Muscimol' 'Muscimol' 'Muscimol' 'Muscimol' 'Muscimol'};
Mgroup = Mgroup';
% stacks of 11 to sort the lists later (11 animals) - frequency
BF = {'aBF' 'aBF' 'aBF' 'aBF' 'aBF' 'aBF' 'aBF' 'aBF' 'aBF' 'aBF' 'aBF'};
BF = BF';
one = {'b1' 'b1' 'b1' 'b1' 'b1' 'b1' 'b1' 'b1' 'b1' 'b1' 'b1'};
one = one';
two = {'c2' 'c2' 'c2' 'c2' 'c2' 'c2' 'c2' 'c2' 'c2' 'c2' 'c2'};
two = two';
% stacks of 11 to sort the lists later (11 animals) - frequency
sink1 = {sink{1}, sink{1}, sink{1}, sink{1}, sink{1}, sink{1}, sink{1}, sink{1}, sink{1}, sink{1}, sink{1}};
sink1 = sink1';
sink2 = {sink{2}, sink{2}, sink{2}, sink{2}, sink{2}, sink{2}, sink{2}, sink{2}, sink{2}, sink{2}, sink{2}};
sink2 = sink2';
sink3 = {sink{3}, sink{3}, sink{3}, sink{3}, sink{3}, sink{3}, sink{3}, sink{3}, sink{3}, sink{3}, sink{3}}; 
sink3 = sink3';
sink4 = {sink{4}, sink{4}, sink{4}, sink{4}, sink{4}, sink{4}, sink{4}, sink{4}, sink{4}, sink{4}, sink{4}}; 
sink4 = sink4';
 
%% sink 1
%KETAMINE
% pull out the + and - of each side (+-1), average them, normalize to BF
mirK1_2 = vertcat(K_norm1(:,1),K_norm1(:,5));
mirK1_1 = vertcat(K_norm1(:,2),K_norm1(:,4));
mirK1_BF = vertcat(K_norm1(:,3),K_norm1(:,3));
%concatonate them vertically (BF*gsize, +-1*gsize,...)
mirK_data1 = vertcat(mirK1_BF,mirK1_1,mirK1_2);
K1 = vertcat(Kgroup, Kgroup, Kgroup, Kgroup, Kgroup, Kgroup);
freqK = vertcat(BF, BF, one, one, two, two);


%AWAKE
% pull out the + and - of each side (+-1), average them, normalize to BF
mirA1_2 = vertcat(A_norm1(:,1),A_norm1(:,5));
mirA1_1 = vertcat(A_norm1(:,2),A_norm1(:,4));
mirA1_BF = vertcat(A_norm1(:,3),A_norm1(:,3));
%concatonate them vertically (BF*gsize, +-1*gsize,...)
mirA_data1 = vertcat(mirA1_BF,mirA1_1,mirA1_2);
A1 = vertcat(Agroup, Agroup, Agroup, Agroup, Agroup, Agroup);
freqA = vertcat(BF, BF, one, one, two, two);

%MUSCIMOL
% pull out the + and - of each side (+-1), average them, normalize to BF
mirM1_2 = vertcat(M_norm1(:,1),M_norm1(:,5));
mirM1_1 = vertcat(M_norm1(:,2),M_norm1(:,4));
mirM1_BF = vertcat(M_norm1(:,3),M_norm1(:,3));
%concatonate them vertically (BF*gsize, +-1*gsize,...)
mirM_data1 = vertcat(mirM1_BF,mirM1_1,mirM1_2);
M1 = vertcat(Mgroup, Mgroup, Mgroup, Mgroup, Mgroup, Mgroup);
freqM = vertcat(BF, BF, one, one, two, two);

%sanity check: all following variables should be equal in length
mirData1 = vertcat(mirK_data1,mirA_data1,mirM_data1);
Groups1 = vertcat(K1,A1,M1);
Freq1 = vertcat(freqK,freqA,freqM);
Sinks1 = vertcat(sink1,sink1,sink1,sink1,sink1,sink1,sink1,sink1,sink1,sink1,sink1,sink1,sink1,sink1,sink1,sink1,sink1,sink1);
%% sink 2
%KETAMINE
% pull out the + and - of each side (+-1), average them, normalize to BF
mirK2_2 = vertcat(K_norm2(:,1),K_norm2(:,5));
mirK2_1 = vertcat(K_norm2(:,2),K_norm2(:,4));
mirK2_BF = vertcat(K_norm2(:,3),K_norm2(:,3));
%concatonate them vertically (BF*gsize, +-1*gsize,...)
mirK_data2 = vertcat(mirK2_BF,mirK2_1,mirK2_2);
K2 = vertcat(Kgroup, Kgroup, Kgroup, Kgroup, Kgroup, Kgroup);
freqK = vertcat(BF, BF, one, one, two, two);


%AWAKE
% pull out the + and - of each side (+-1), average them, normalize to BF
mirA2_2 = vertcat(A_norm2(:,1),A_norm2(:,5));
mirA2_1 = vertcat(A_norm2(:,2),A_norm2(:,4));
mirA2_BF = vertcat(A_norm2(:,3),A_norm2(:,3));
%concatonate them vertically (BF*gsize, +-1*gsize,...)
mirA_data2 = vertcat(mirA2_BF,mirA2_1,mirA2_2);
A2 = vertcat(Agroup, Agroup, Agroup, Agroup, Agroup, Agroup);
freqA = vertcat(BF, BF, one, one, two, two);

%MUSCIMOL
% pull out the + and - of each side (+-1), average them, normalize to BF
mirM2_2 = vertcat(M_norm2(:,1),M_norm2(:,5));
mirM2_1 = vertcat(M_norm2(:,2),M_norm2(:,4));
mirM2_BF = vertcat(M_norm2(:,3),M_norm2(:,3));
%concatonate them vertically (BF*gsize, +-1*gsize,...)
mirM_data2 = vertcat(mirM2_BF,mirM2_1,mirM2_2);
M2 = vertcat(Mgroup, Mgroup, Mgroup, Mgroup, Mgroup, Mgroup);
freqM = vertcat(BF, BF, one, one, two, two);

%sanity check: all following variables should be equal in length
mirData2 = vertcat(mirK_data2,mirA_data2,mirM_data2);
Groups2 = vertcat(K2,A2,M2);
Freq2 = vertcat(freqK,freqA,freqM);
Sinks2 = vertcat(sink2,sink2,sink2,sink2,sink2,sink2,sink2,sink2,sink2,sink2,sink2,sink2,sink2,sink2,sink2,sink2,sink2,sink2);
%% sink 3
%KETAMINE
% pull out the + and - of each side (+-1), average them, normalize to BF
mirK3_2 = vertcat(K_norm3(:,1),K_norm3(:,5));
mirK3_1 = vertcat(K_norm3(:,2),K_norm3(:,4));
mirK3_BF = vertcat(K_norm3(:,3),K_norm3(:,3));
%concatonate them vertically (BF*gsize, +-1*gsize,...)
mirK_data3 = vertcat(mirK3_BF,mirK3_1,mirK3_2);
K3 = vertcat(Kgroup, Kgroup, Kgroup, Kgroup, Kgroup, Kgroup);
freqK = vertcat(BF, BF, one, one, two, two);


%AWAKE
% pull out the + and - of each side (+-1), average them, normalize to BF
mirA3_2 = vertcat(A_norm3(:,1),A_norm3(:,5));
mirA3_1 = vertcat(A_norm3(:,2),A_norm3(:,4));
mirA3_BF = vertcat(A_norm3(:,3),A_norm3(:,3));
%concatonate them vertically (BF*gsize, +-1*gsize,...)
mirA_data3 = vertcat(mirA3_BF,mirA3_1,mirA3_2);
A3 = vertcat(Agroup, Agroup, Agroup, Agroup, Agroup, Agroup);
freqA = vertcat(BF, BF, one, one, two, two);

%MUSCIMOL
% pull out the + and - of each side (+-1), average them, normalize to BF
mirM3_2 = vertcat(M_norm3(:,1),M_norm3(:,5));
mirM3_1 = vertcat(M_norm3(:,2),M_norm3(:,4));
mirM3_BF = vertcat(M_norm3(:,3),M_norm3(:,3));
%concatonate them vertically (BF*gsize, +-1*gsize,...)
mirM_data3 = vertcat(mirM3_BF,mirM3_1,mirM3_2);
M3 = vertcat(Mgroup, Mgroup, Mgroup, Mgroup, Mgroup, Mgroup);
freqM = vertcat(BF, BF, one, one, two, two);

%sanity check: all following variables should be equal in length
mirData3 = vertcat(mirK_data3,mirA_data3,mirM_data3);
Groups3 = vertcat(K3,A3,M3);
Freq3 = vertcat(freqK,freqA,freqM);
Sinks3 = vertcat(sink3,sink3,sink3,sink3,sink3,sink3,sink3,sink3,sink3,sink3,sink3,sink3,sink3,sink3,sink3,sink3,sink3,sink3);
%% sink 4
%KETAMINE
% pull out the + and - of each side (+-1), average them, normalize to BF
mirK4_2 = vertcat(K_norm4(:,1),K_norm4(:,5));
mirK4_1 = vertcat(K_norm4(:,2),K_norm4(:,4));
mirK4_BF = vertcat(K_norm4(:,3),K_norm4(:,3));
%concatonate them vertically (BF*gsize, +-1*gsize,...)
mirK_data4 = vertcat(mirK4_BF,mirK4_1,mirK4_2);
K4 = vertcat(Kgroup, Kgroup, Kgroup, Kgroup, Kgroup, Kgroup);
freqK = vertcat(BF, BF, one, one, two, two);


%AWAKE
% pull out the + and - of each side (+-1), average them, normalize to BF
mirA4_2 = vertcat(A_norm4(:,1),A_norm4(:,5));
mirA4_1 = vertcat(A_norm4(:,2),A_norm4(:,4));
mirA4_BF = vertcat(A_norm4(:,3),A_norm4(:,3));
%concatonate them vertically (BF*gsize, +-1*gsize,...)
mirA_data4 = vertcat(mirA4_BF,mirA4_1,mirA4_2);
A4 = vertcat(Agroup, Agroup, Agroup, Agroup, Agroup, Agroup);
freqA = vertcat(BF, BF, one, one, two, two);

%MUSCIMOL
% pull out the + and - of each side (+-1), average them, normalize to BF
mirM4_2 = vertcat(M_norm4(:,1),M_norm4(:,5));
mirM4_1 = vertcat(M_norm4(:,2),M_norm4(:,4));
mirM4_BF = vertcat(M_norm4(:,3),M_norm4(:,3));
%concatonate them vertically (BF*gsize, +-1*gsize,...)
mirM_data4 = vertcat(mirM4_BF,mirM4_1,mirM4_2);
M4 = vertcat(Mgroup, Mgroup, Mgroup, Mgroup, Mgroup, Mgroup);
freqM = vertcat(BF, BF, one, one, two, two);

%sanity check: all following variables should be equal in length
mirData4 = vertcat(mirK_data4,mirA_data4,mirM_data4);
Groups4 = vertcat(K4,A4,M4);
Freq4 = vertcat(freqK,freqA,freqM);
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
g(1,1).stat_summary('type','std','geom','errorbar'); %mean and std shown
g(1,1).stat_summary('type','std','geom','point'); %mean and std shown
g(1,1).stat_summary('type','std','geom','area'); %mean and std shown
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


g.set_title([(based{which}) ' ' (sink{1}) ', ' (sink{2}) ', ' (sink{3}) ', and ' (sink{4}) ' for ' (Parameter{para}) 'std']);
g.draw();
g.export('file_name',['Normalized mirror ' (based{which}) '_' (sink{1}) '_' ...
    (sink{2}) '_' (sink{3}) '_' (sink{4}) '_' (Parameter{para}) 'std'], 'file_type','png');
g.export('file_name',['Normalized mirror ' (based{which}) '_' (sink{1}) '_' ...
    (sink{2}) '_' (sink{3}) '_' (sink{4}) '_' (Parameter{para}) 'std'], 'file_type','pdf');
close all




