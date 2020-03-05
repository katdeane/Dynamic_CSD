%% Stats Loop ST
% copy of stats loop but with ST data!
%This loop is an interface for the teg_repeated_measures_Anova code. It
%runs through each parameter (pulled from group analysis Data structure) in
%each specified sink. It generates simple errorbar line plots in a group 
%comparison centered around the BF. Variables to change, just
%below, are varinames (in our case groups and frequencies), levels (how many
%groups and how many frequencies), Order, and Parameter based on what is
%needed for comparisons. Changing the specific input will require renaming
%variable within the code itself. 

%% Start
clear all
cd('D:\MyCode\Dynamic_CSD_Analysis');
warning('OFF');
dbstop if error

home = pwd;
addpath(genpath(home));

%for teg_repeated measures_ANOVA
levels = [2,7];
varinames = {'Groups','Frequencies'};
rowofnans = NaN(1,7);
srowofnans = [{NaN(1,50)} {NaN(1,50)} {NaN(1,50)} {NaN(1,50)} {NaN(1,50)} {NaN(1,50)} {NaN(1,50)}];

ticks = {'-3' '-2' '-1' 'BF' '+1' '+2' '+3'};

% Order = {'IVE','IVL','I_IIE','I_IIL', 'VaE','VaL','VbE','VbL','VIE','VIL'};
Order = {'IVE','IVL','I_IIE','I_IIL', 'VaE','VaL','VbE','VbL','VIaE','VIaL','VIbE','VIbL'};
Parameter = {'SinkRMS','SinkPeakAmp','SinkPeakLate','Sinkonset'};
SParameter = {'SingleSinkRMS','SingleSinkPeakAmp','SingleSinkPeakLat'};


%% Load in the appropriate files
cd DATA;cd output;
load('AnesthetizedPre_Data.m_Threshold_0.25_Zscore_0_binned_1_mirror_0.mat')
Anesthetized = Data; clear Data; 
load('Awake10dB_Data.m_Threshold_0.25_Zscore_0_binned_1_mirror_0.mat')
Awake = Data; clear Data;
load('ANChronic_Data.m_Threshold_0.25_Zscore_0_binned_1_mirror_0.mat')
AnChronic = Data; clear Data;
load('Muscimol_Data.m_Threshold_0.25_Zscore_0_binned_1_mirror_0.mat')
Muscimol = Data; clear Data;


cd(home);cd figs;
mkdir('Group Tuning Curves ST'); % change this name if you ever need to run more stats so as not to over-write
cd('Group Tuning Curves ST');mkdir('AVG');mkdir('SINGLE');addpath(genpath(home));

%% Sink Loop
for isink = 1:length(Order)
    
    disp(['***********************ST STATS FOR ' (Order{isink}) '***********************'])
    %% Avg Trial Loop
    cd(home); cd figs; cd('Group Tuning Curves ST');cd('AVG');
    for ipara = 1:length(Parameter)
        
        disp(['********' (Parameter{ipara}) '********']) %these headers allow a copy and paste of all post-hoc tests which aren't saved in the data at the end
        
        % Below are the data, means, and sems of the current parameter and
        % sink
        An_data = vertcat(Anesthetized.ST_based.(Parameter{ipara}).(Order{isink})(:,5:11)); %, rowofnans
        An_datamean = nanmean(An_data,1);
        An_datasem = nanstd(An_data,1)/sqrt(sum(~isnan(An_datamean)));
        
        Aw_data = vertcat(Awake.ST_based.(Parameter{ipara}).(Order{isink})(:,4:10), rowofnans,rowofnans);
        Aw_datamean = nanmean(Aw_data,1);
        Aw_datasem = nanstd(Aw_data,1)/sqrt(sum(~isnan(Aw_datamean)));
        
        AnC_data = vertcat(AnChronic.GS_based.(Parameter{ipara}).(Order{isink})(:,4:10), rowofnans,rowofnans);
        AnC_datamean = nanmean(AnC_data,1);
        AnC_datasem = nanstd(AnC_data,1)/sqrt(sum(~isnan(AnC_datamean)));
        
        M_data = vertcat(Muscimol.ST_based.(Parameter{ipara}).(Order{isink})(:,5:11));
        M_datamean = nanmean(M_data,1);
        M_datasem = nanstd(M_data,1)/sqrt(sum(~isnan(M_datamean)));
        
        %now we generate a figure based on these three groups
        h=figure;
        
        errorbar(1:7, An_datamean, An_datasem, 'Linewidth', 2); hold on;
        errorbar(1:7, Aw_datamean, Aw_datasem, 'Linewidth', 2);
        errorbar(1:7, AnC_datamean, AnC_datasem, 'Linewidth', 2);
        errorbar(1:7, M_datamean, M_datasem, 'Linewidth', 2);
        legend('Anesthetized', 'Awake', 'Anesthetized Chronic', 'Muscimol')
        title(['Sink' (Order{isink}) (Parameter{ipara})])
        set(gca,'XTick',1:7); set(gca,'XTickLabel',ticks,'FontSize',10);
        
        cd([home '\figs\' 'Group Tuning Curves ST']);
        savefig(h,['Sink' (Order{isink}) (Parameter{ipara})]); close all;
        
        %finally, we run the statistics. The stats will only appear in the
        %command window if they are significant
        com_AnvsAw = horzcat(Aw_data, An_data);
        com_AnCvsAw = horzcat(Aw_data, AnC_data);
        com_AnvsM = horzcat(An_data, M_data);
        
        disp('**Anesthetized vs Awake**')
        AnesthetizedvsAwake = teg_repeated_measures_ANOVA(com_AnvsAw, levels, varinames);
        
        disp('**Anesthetized Chronic vs Awake**')
        AnChronicvsAwake = teg_repeated_measures_ANOVA(com_AnCvsAw, levels, varinames);
        
        disp('**Anesthetized vs Muscimol**')
        AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(com_AnvsM, levels, varinames);
        
        %significant or not, the F = #(#,#) and p= # will be saves in this file
        cd(home); cd DATA
        mkdir('StatsST'); cd StatsST
        savefile = [(Order{isink}) (Parameter{ipara}) 'SinkStats.mat'];
        save(savefile, 'AnesthetizedvsAwake', 'AnChronicvsAwake', 'AnesthetizedvsMuscimol')
    end
    
end

%% Normalized Sink Group
for isink = 1:length(Order)
    
    disp(['***********************Normalized ST STATS FOR ' (Order{isink}) '***********************'])
    %% Avg Trial Loop
    cd(home); cd figs; cd('Group Tuning Curves ST');cd('AVG');
    for ipara = 1:length(Parameter)
        
        disp(['********Normalized ' (Parameter{ipara}) '********']) %these headers allow a copy and paste of all post-hoc tests which aren't saved in the data at the end
        
        % Below are the data, means, and sems of the current parameter and
        % sink
        An_data = vertcat(Anesthetized.ST_based.(Parameter{ipara}).(Order{isink})(:,5:11)); %, rowofnans
        An_data = An_data./An_data(:,4);
        An_datamean = nanmean(An_data,1);
        An_datasem = nanstd(An_data,1)/sqrt(sum(~isnan(An_datamean)));
        
        Aw_data = vertcat(Awake.ST_based.(Parameter{ipara}).(Order{isink})(:,4:10), rowofnans,rowofnans);
        Aw_data = Aw_data./Aw_data(:,4);
        Aw_datamean = nanmean(Aw_data,1);
        Aw_datasem = nanstd(Aw_data,1)/sqrt(sum(~isnan(Aw_datamean)));
        
        AnC_data = vertcat(AnChronic.GS_based.(Parameter{ipara}).(Order{isink})(:,4:10), rowofnans,rowofnans);
        AnC_data = AnC_data./AnC_data(:,4);
        AnC_datamean = nanmean(AnC_data,1);
        AnC_datasem = nanstd(AnC_data,1)/sqrt(sum(~isnan(AnC_datamean)));
        
        M_data = vertcat(Muscimol.ST_based.(Parameter{ipara}).(Order{isink})(:,5:11));
        M_data = M_data./M_data(:,4);
        M_datamean = nanmean(M_data,1);
        M_datasem = nanstd(M_data,1)/sqrt(sum(~isnan(M_datamean)));
        
        %now we generate a figure based on these three groups
        h=figure;
        
        errorbar(1:7, An_datamean, An_datasem, 'Linewidth', 2); hold on;
        errorbar(1:7, Aw_datamean, Aw_datasem, 'Linewidth', 2);
        errorbar(1:7, AnC_datamean, AnC_datasem, 'Linewidth', 2);
        errorbar(1:7, M_datamean, M_datasem, 'Linewidth', 2);
        legend('Anesthetized', 'Awake', 'Anesthetized Chronic', 'Muscimol')
        title(['Sink' (Order{isink}) (Parameter{ipara})])
        set(gca,'XTick',1:7); set(gca,'XTickLabel',ticks,'FontSize',10);
        
        cd([home '\figs\' 'Group Tuning Curves ST']);
        savefig(h,['Sink' (Order{isink}) (Parameter{ipara})]); close all;
        
        %finally, we run the statistics. The stats will only appear in the
        %command window if they are significant
        com_AnvsAw = horzcat(Aw_data, An_data);
        com_AnCvsAw = horzcat(Aw_data, AnC_data);
        com_AnvsM = horzcat(An_data, M_data);
        
        disp('**Norm Anesthetized vs Awake**')
        AnesthetizedvsAwakeN = teg_repeated_measures_ANOVA(com_AnvsAw, levels, varinames);
        
        disp('**Norm Anesthetized Chronic vs Awake**')
        AnChronicvsAwakeN = teg_repeated_measures_ANOVA(com_AnCvsAw, levels, varinames);
        
        disp('**Norm Anesthetized vs Muscimol**')
        AnesthetizedvsMuscimolN = teg_repeated_measures_ANOVA(com_AnvsM, levels, varinames);
        
        %significant or not, the F = #(#,#) and p= # will be saves in this file
        cd(home); cd DATA
        mkdir('NormStatsST'); cd NormStatsST
        savefile = [(Order{isink}) (Parameter{ipara}) 'NormSinkStats.mat'];
        save(savefile, 'AnesthetizedvsAwakeN', 'AnChronicvsAwakeN', 'AnesthetizedvsMuscimolN')
    end
    
end
