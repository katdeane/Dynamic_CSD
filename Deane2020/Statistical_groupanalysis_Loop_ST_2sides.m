function Statistical_groupanalysis_Loop_ST_2sides(homedir)
%% Stats Loop FOR 5 FREQUENCIES self-tuned

%This loop is an interface for the teg_repeated_measures_Anova code. It
%runs through each parameter (pulled from group analysis Data structure) in
%each specified sink. It generates simple errorbar line plots in a group 
%comparison centered around the BF. Variables to change, just
%below, are varinames (in our case groups and frequencies), levels (how many
%groups and how many frequencies), Order, and Parameter based on what is
%needed for comparisons. Changing the specific input will require renaming
%variable within the code itself. 

%% Start
warning('OFF');
dbstop if error

% Change directory to your working folder
if ~exist('homedir','var')
    if exist('D:\MyCode\Dynamic_CSD','dir') == 7
        cd('D:\MyCode\Dynamic_CSD');
    elseif exist('C:\Users\kedea\Documents\Work Stuff\Dynamic_CSD','dir') == 7
        cd('C:\Users\kedea\Documents\Work Stuff\Dynamic_CSD')
    end
    
    homedir = pwd;
    addpath(genpath(homedir));
end

%for teg_repeated measures_ANOVA
levels = [2,5];
levelsMir = [2,3];
varinames = {'Groups','Frequencies'};
rowofnans = NaN(1,5);
srowofnans = [{NaN(1,50)} {NaN(1,50)} {NaN(1,50)} {NaN(1,50)} {NaN(1,50)}];

ticks = {'-2' '-1' 'BF' '+1' '+2'};

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


cd(homedir);cd figs;
mkdir('Group Tuning Curves ST'); % change this name if you ever need to run more stats so as not to over-write
cd('Group Tuning Curves ST');mkdir('AVG');mkdir('SINGLE');

%% Sink Loop
for isink = 1:length(Order)
    
    disp(['***********************ST STATS FOR ' (Order{isink}) '***********************'])
    %% Avg Trial Loop
    cd(homedir); cd figs; cd('Group Tuning Curves ST');cd('AVG');
    for ipara = 1:length(Parameter)
        
        disp(['********' (Parameter{ipara}) '********']) %these headers allow a copy and paste of all post-hoc tests which aren't saved in the data at the end
        
        % Below are the data, means, and sems of the current parameter and
        % sink
        An_data = vertcat(Anesthetized.ST_based.(Parameter{ipara}).(Order{isink})(:,6:10)); %, rowofnans
        An_datamean = nanmean(An_data,1);
        An_datasem = nanstd(An_data,1)/sqrt(sum(~isnan(An_datamean)));
        
        Aw_data = vertcat(Awake.ST_based.(Parameter{ipara}).(Order{isink})(:,5:9), rowofnans,rowofnans);
        Aw_datamean = nanmean(Aw_data,1);
        Aw_datasem = nanstd(Aw_data,1)/sqrt(sum(~isnan(Aw_datamean)));
        
        AnC_data = vertcat(AnChronic.GS_based.(Parameter{ipara}).(Order{isink})(:,5:9), rowofnans,rowofnans);
        AnC_datamean = nanmean(AnC_data,1);
        AnC_datasem = nanstd(AnC_data,1)/sqrt(sum(~isnan(AnC_datamean)));
        
        M_data = vertcat(Muscimol.ST_based.(Parameter{ipara}).(Order{isink})(:,6:10));
        M_datamean = nanmean(M_data,1);
        M_datasem = nanstd(M_data,1)/sqrt(sum(~isnan(M_datamean)));
        
        %now we generate a figure based on these three groups
        h=figure;
        
        errorbar(1:5, An_datamean, An_datasem, 'Linewidth', 2); hold on;
        errorbar(1:5, Aw_datamean, Aw_datasem, 'Linewidth', 2);
        errorbar(1:5, AnC_datamean, AnC_datasem, 'Linewidth', 2);
        errorbar(1:5, M_datamean, M_datasem, 'Linewidth', 2);
        legend('Anesthetized', 'Awake', 'Anesthetized Chronic', 'Muscimol')
        title(['Sink' (Order{isink}) (Parameter{ipara})])
        set(gca,'XTick',1:5); set(gca,'XTickLabel',ticks,'FontSize',10);
        
        cd([homedir '\figs\' 'Group Tuning Curves ST']);
        savefig(h,['Sink' (Order{isink}) (Parameter{ipara})]); close all;
        
        %finally, we run the statistics. The stats will only appear in the
        %command window if they are significant
        com_AnvsAw = horzcat(Aw_data, An_data);
        com_AnCvsAw = horzcat(Aw_data, AnC_data);
        com_AnvsM = horzcat(An_data, M_data);
        
        disp('**Anesthetized vs Awake**')
        AnesthetizedvsAwake = teg_repeated_measures_ANOVA(com_AnvsAw, levels, varinames);
        AvK_CD = iMakeCohensD(Aw_data, An_data)
            
        disp('**Anesthetized Chronic vs Awake**')
        AnChronicvsAwake = teg_repeated_measures_ANOVA(com_AnCvsAw, levels, varinames);
        AvN_CD = iMakeCohensD(Aw_data, AnC_data)
        
        disp('**Anesthetized vs Muscimol**')
        AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(com_AnvsM, levels, varinames);
        MvK_CD = iMakeCohensD(M_data, An_data)
        
        %significant or not, the F = #(#,#) and p= # will be saves in this file
        cd(homedir); cd DATA
        mkdir('StatsST'); cd StatsST
        savefile = [(Order{isink}) (Parameter{ipara}) 'SinkStats.mat'];
        save(savefile, 'AnesthetizedvsAwake', 'AvK_CD', 'AnChronicvsAwake', 'AvN_CD', 'AnesthetizedvsMuscimol', 'MvK_CD')
    end
    
end

%% Normalized Sink Group
for isink = 1:length(Order)
    
    disp(['***********************Normalized ST STATS FOR ' (Order{isink}) '***********************'])
    %% Avg Trial Loop
    cd(homedir); cd figs; cd('Group Tuning Curves ST');cd('AVG');
    for ipara = 1:length(Parameter)
        
        disp(['********Normalized ' (Parameter{ipara}) '********']) %these headers allow a copy and paste of all post-hoc tests which aren't saved in the data at the end
        
        % Below are the data, means, and sems of the current parameter and
        % sink
        An_data = vertcat(Anesthetized.ST_based.(Parameter{ipara}).(Order{isink})(:,6:10)); %, rowofnans
        An_data(An_data == 0) = 1;
        An_data = An_data./An_data(:,3);
        An_datamean = nanmean(An_data,1);
        An_datasem = nanstd(An_data,1)/sqrt(sum(~isnan(An_datamean)));
        
        Aw_data = vertcat(Awake.ST_based.(Parameter{ipara}).(Order{isink})(:,5:9), rowofnans,rowofnans);
        Aw_data(Aw_data == 0) = 1;
        Aw_data = Aw_data./Aw_data(:,3);
        Aw_datamean = nanmean(Aw_data,1);
        Aw_datasem = nanstd(Aw_data,1)/sqrt(sum(~isnan(Aw_datamean)));
        
        AnC_data = vertcat(AnChronic.GS_based.(Parameter{ipara}).(Order{isink})(:,5:9), rowofnans,rowofnans);
        AnC_data(AnC_data == 0) = 1;
        AnC_data = AnC_data./AnC_data(:,3);
        AnC_datamean = nanmean(AnC_data,1);
        AnC_datasem = nanstd(AnC_data,1)/sqrt(sum(~isnan(AnC_datamean)));
        
        M_data = vertcat(Muscimol.ST_based.(Parameter{ipara}).(Order{isink})(:,6:10));
        M_data(M_data == 0) = 1;
        M_data = M_data./M_data(:,3);
        M_datamean = nanmean(M_data,1);
        M_datasem = nanstd(M_data,1)/sqrt(sum(~isnan(M_datamean)));
        
        %now we generate a figure based on these three groups
        h=figure;
        
        errorbar(1:5, An_datamean, An_datasem, 'Linewidth', 2); hold on;
        errorbar(1:5, Aw_datamean, Aw_datasem, 'Linewidth', 2);
        errorbar(1:5, AnC_datamean, AnC_datasem, 'Linewidth', 2);
        errorbar(1:5, M_datamean, M_datasem, 'Linewidth', 2);
        legend('Anesthetized', 'Awake', 'Anesthetized Chronic', 'Muscimol')
        title(['Sink' (Order{isink}) (Parameter{ipara})])
        set(gca,'XTick',1:5); set(gca,'XTickLabel',ticks,'FontSize',10);
        
        cd([homedir '\figs\' 'Group Tuning Curves ST']);
        savefig(h,['Sink' (Order{isink}) (Parameter{ipara})]); close all;
        
        %finally, we run the statistics. The stats will only appear in the
        %command window if they are significant
       
        com_AnvsAw = horzcat(Aw_data, An_data);
        com_AnCvsAw = horzcat(Aw_data, AnC_data);
        com_AnvsM = horzcat(An_data, M_data);
        
        disp('**Anesthetized vs Awake**')
        AnesthetizedvsAwakeN = teg_repeated_measures_ANOVA(com_AnvsAw, levels, varinames);
        AvK_CDN = iMakeCohensD(Aw_data, An_data)
            
        disp('**Anesthetized Chronic vs Awake**')
        AnChronicvsAwakeN = teg_repeated_measures_ANOVA(com_AnCvsAw, levels, varinames);
        AvN_CDN = iMakeCohensD(Aw_data, AnC_data)
        
        disp('**Anesthetized vs Muscimol**')
        AnesthetizedvsMuscimolN = teg_repeated_measures_ANOVA(com_AnvsM, levels, varinames);
        MvK_CDN = iMakeCohensD(M_data, An_data)
        
        %significant or not, the F = #(#,#) and p= # will be saves in this file
        cd(homedir); cd DATA
        mkdir('StatsST'); cd StatsST
        savefile = [(Order{isink}) (Parameter{ipara}) 'SinkStats_Norm.mat'];
        save(savefile, 'AnesthetizedvsAwakeN', 'AvK_CDN', 'AnChronicvsAwakeN', 'AvN_CDN', 'AnesthetizedvsMuscimolN', 'MvK_CDN')
    end
    
end

%% Normalized Mirror Sink Group
for isink = 1:length(Order)
    
    disp(['***********************Mirror Normalized ST STATS FOR ' (Order{isink}) '***********************'])
    %% Avg Trial Loop
    cd(homedir); cd figs; cd('Group Tuning Curves ST');cd('AVG');
    for ipara = 1:length(Parameter)
        
        disp(['********Mirror Normalized ' (Parameter{ipara}) '********']) %these headers allow a copy and paste of all post-hoc tests which aren't saved in the data at the end
        
        % Below are the data, means, and sems of the current parameter and
        % sink
        An_data = vertcat(Anesthetized.ST_based.(Parameter{ipara}).(Order{isink})(:,6:10)); %, rowofnans
        An_data(An_data == 0) = 1;
        An_data = An_data./An_data(:,3);
        % stack the left under the right and double the central colume
        An_Mir_data = vertcat(An_data(:,3:5),(horzcat(An_data(:,3),An_data(:,2),An_data(:,1))));
        An_datamean = nanmean(An_Mir_data,1);
        An_datasem = nanstd(An_Mir_data,1)/sqrt(sum(~isnan(An_datamean)));
        
        Aw_data = vertcat(Awake.ST_based.(Parameter{ipara}).(Order{isink})(:,5:9), rowofnans,rowofnans);
        Aw_data(Aw_data == 0) = 1;
        Aw_data = Aw_data./Aw_data(:,3);
        % stack the left under the right and double the central colume
        Aw_Mir_data = vertcat(Aw_data(:,3:5),(horzcat(Aw_data(:,3),Aw_data(:,2),Aw_data(:,1))));
        Aw_datamean = nanmean(Aw_Mir_data,1);
        Aw_datasem = nanstd(Aw_Mir_data,1)/sqrt(sum(~isnan(Aw_datamean)));
        
        AnC_data = vertcat(AnChronic.GS_based.(Parameter{ipara}).(Order{isink})(:,5:9), rowofnans,rowofnans);
        AnC_data(AnC_data == 0) = 1;
        AnC_data = AnC_data./AnC_data(:,3);
        % stack the left under the right and double the central colume
        AnC_Mir_data = vertcat(AnC_data(:,3:5),(horzcat(AnC_data(:,3),AnC_data(:,2),AnC_data(:,1))));
        AnC_datamean = nanmean(AnC_Mir_data,1);
        AnC_datasem = nanstd(AnC_Mir_data,1)/sqrt(sum(~isnan(AnC_datamean)));
        
        M_data = vertcat(Muscimol.ST_based.(Parameter{ipara}).(Order{isink})(:,6:10));
        M_data(M_data == 0) = 1;
        M_data = M_data./M_data(:,3);
        % stack the left under the right and double the central colume
        M_Mir_data = vertcat(M_data(:,3:5),(horzcat(M_data(:,3),M_data(:,2),M_data(:,1))));
        M_datamean = nanmean(M_Mir_data,1);
        M_datasem = nanstd(M_Mir_data,1)/sqrt(sum(~isnan(M_datamean)));
        
        %now we generate a figure based on these three groups
        h=figure;
        
        errorbar(1:3, An_datamean, An_datasem, 'Linewidth', 2); hold on;
        errorbar(1:3, Aw_datamean, Aw_datasem, 'Linewidth', 2);
        errorbar(1:3, AnC_datamean, AnC_datasem, 'Linewidth', 2);
        errorbar(1:3, M_datamean, M_datasem, 'Linewidth', 2);
        legend('Anesthetized', 'Awake', 'Anesthetized Chronic', 'Muscimol')
        title(['Sink' (Order{isink}) (Parameter{ipara})])
        set(gca,'XTick',1:5); set(gca,'XTickLabel',ticks,'FontSize',10);
        
        cd([homedir '\figs\' 'Group Tuning Curves ST']);
        savefig(h,['Sink' (Order{isink}) (Parameter{ipara})]); close all;
        
        %finally, we run the statistics. The stats will only appear in the
        %command window if they are significant
       
        com_AnvsAw = horzcat(Aw_Mir_data, An_Mir_data);
        com_AnCvsAw = horzcat(Aw_Mir_data, AnC_Mir_data);
        com_AnvsM = horzcat(An_Mir_data, M_Mir_data);
        
        disp('**Anesthetized vs Awake**')
        AnesthetizedvsAwakeNM = teg_repeated_measures_ANOVA(com_AnvsAw, levelsMir, varinames);
        AvK_CDNM = iMakeCohensD(Aw_Mir_data, An_Mir_data)
            
        disp('**Anesthetized Chronic vs Awake**')
        AnChronicvsAwakeNM = teg_repeated_measures_ANOVA(com_AnCvsAw, levelsMir, varinames);
        AvN_CDNM = iMakeCohensD(Aw_Mir_data, AnC_Mir_data)
        
        disp('**Anesthetized vs Muscimol**')
        AnesthetizedvsMuscimolNM = teg_repeated_measures_ANOVA(com_AnvsM, levelsMir, varinames);
        MvK_CDNM = iMakeCohensD(M_Mir_data, An_Mir_data)
        
        %significant or not, the F = #(#,#) and p= # will be saves in this file
        cd(homedir); cd DATA
        mkdir('StatsST'); cd StatsST
        savefile = [(Order{isink}) (Parameter{ipara}) 'SinkStats_Mir_Norm.mat'];
        save(savefile, 'AnesthetizedvsAwakeNM', 'AvK_CDNM', 'AnChronicvsAwakeNM', 'AvN_CDNM', 'AnesthetizedvsMuscimolNM', 'MvK_CDNM')
    end
    
end
