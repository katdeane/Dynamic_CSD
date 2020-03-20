function Stats_Loop_5Freq(homedir)
% just the important stats
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

ticks = {'-2' '-1' 'BF' '+1' '+2'};

Order = {'I_IIE','IVE','VbE','VIaE'};
Parameter = {'SinkRMS','SinkPeakAmp','SinkPeakLate'};


%% Load in the appropriate files

cd(homedir);cd DATA;cd Output;
load('AnesthetizedPre_Data.m_Threshold_0.25_Zscore_0_binned_1_mirror_0.mat','Data')
Anesthetized = Data; clear Data; 
load('Awake10dB_Data.m_Threshold_0.25_Zscore_0_binned_1_mirror_0.mat','Data')
Awake = Data; clear Data;
load('Muscimol_Data.m_Threshold_0.25_Zscore_0_binned_1_mirror_0.mat','Data')
Muscimol = Data; clear Data;

cd(homedir);cd figs;
mkdir('Group Tuning Curves ST'); % change this name if you ever need to run more stats so as not to over-write
cd('Group Tuning Curves ST');
%% Sink Loop
for isink = 1:length(Order)
    
    disp(['***********************ST STATS FOR ' (Order{isink}) '***********************'])
    %% Avg Trial Loop
    for ipara = 1:length(Parameter)
        cd(homedir); cd figs; cd('Group Tuning Curves ST');
        disp(['********' (Parameter{ipara}) '********']) %these headers allow a copy and paste of all post-hoc tests which aren't saved in the data at the end
        
        % Below are the data, means, and sems of the current parameter and
        % sink
        An_data = vertcat(Anesthetized.ST_based.(Parameter{ipara}).(Order{isink})(:,6:10)); %, rowofnans
        An_datamean = nanmean(An_data,1);
        An_datastd = nanstd(An_data,1);
        
        Aw_data = vertcat(Awake.ST_based.(Parameter{ipara}).(Order{isink})(:,5:9), rowofnans,rowofnans);
        Aw_datamean = nanmean(Aw_data,1);
        Aw_datastd = nanstd(Aw_data,1);
        
        M_data = vertcat(Muscimol.ST_based.(Parameter{ipara}).(Order{isink})(:,6:10));
        M_datamean = nanmean(M_data,1);
        M_datastd = nanstd(M_data,1);
        
        %now we generate a figure based on these three groups
        h=figure;
        
        errorbar(1:5, An_datamean, An_datastd, 'Linewidth', 2); hold on;
        errorbar(1:5, Aw_datamean, Aw_datastd, 'Linewidth', 2);
        errorbar(1:5, M_datamean, M_datastd, 'Linewidth', 2);
        legend('Anesthetized', 'Awake', 'Muscimol')
        title(['Sink' (Order{isink}) (Parameter{ipara})])
        set(gca,'XTick',1:5); set(gca,'XTickLabel',ticks,'FontSize',10);
        
        savefig(h,['Sink' (Order{isink}) (Parameter{ipara})]); close all;
        
        %finally, we run the statistics. The stats will only appear in the
        %command window if they are significant
        com_AnvsAw = horzcat(Aw_data, An_data);
        com_AnvsM = horzcat(An_data, M_data);
        
        disp('**Anesthetized vs Awake**')
        AnesthetizedvsAwake = teg_repeated_measures_ANOVA(com_AnvsAw, levels, varinames);
        AvK_CD = iMakeCohensD(Aw_data, An_data)  %#ok<*NOPTS>
                    
        disp('**Anesthetized vs Muscimol**')
        AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(com_AnvsM, levels, varinames);
        MvK_CD = iMakeCohensD(M_data, An_data)
        
        %significant or not, the F = #(#,#) and p= # will be saves in this file
        cd(homedir); cd DATA
        mkdir('StatsST'); cd StatsST
        savefile = [(Order{isink}) (Parameter{ipara}) 'SinkStats.mat'];
        save(savefile, 'AnesthetizedvsAwake', 'AvK_CD', 'AnesthetizedvsMuscimol', 'MvK_CD')
    end
    
end

%% Normalized Sink Group
for isink = 1:length(Order)
    
    disp(['***********************Normalized ST STATS FOR ' (Order{isink}) '***********************'])
    %% Avg Trial Loop
    for ipara = 1:length(Parameter)
        cd(homedir); cd figs; cd('Group Tuning Curves ST');
        disp(['********Normalized ' (Parameter{ipara}) '********']) %these headers allow a copy and paste of all post-hoc tests which aren't saved in the data at the end
        
        % Below are the data, means, and sems of the current parameter and
        % sink
        An_data = vertcat(Anesthetized.ST_based.(Parameter{ipara}).(Order{isink})(:,6:10)); %, rowofnans
        An_data(An_data == 0) = 1;
        An_data = An_data./An_data(:,3);
        An_datamean = nanmean(An_data,1);
        An_datastd = nanstd(An_data,1);
        
        Aw_data = vertcat(Awake.ST_based.(Parameter{ipara}).(Order{isink})(:,5:9), rowofnans,rowofnans);
        Aw_data(Aw_data == 0) = 1;
        Aw_data = Aw_data./Aw_data(:,3);
        Aw_datamean = nanmean(Aw_data,1);
        Aw_datastd = nanstd(Aw_data,1);
        
        M_data = vertcat(Muscimol.ST_based.(Parameter{ipara}).(Order{isink})(:,6:10));
        M_data(M_data == 0) = 1;
        M_data = M_data./M_data(:,3);
        M_datamean = nanmean(M_data,1);
        M_datastd = nanstd(M_data,1);
        
        %now we generate a figure based on these three groups
        h=figure;
        
        errorbar(1:5, An_datamean, An_datastd, 'Linewidth', 2); hold on;
        errorbar(1:5, Aw_datamean, Aw_datastd, 'Linewidth', 2);
        errorbar(1:5, M_datamean, M_datastd, 'Linewidth', 2);
        legend('Anesthetized', 'Awake', 'Muscimol')
        title(['Sink' (Order{isink}) (Parameter{ipara})])
        set(gca,'XTick',1:5); set(gca,'XTickLabel',ticks,'FontSize',10);
        
        savefig(h,['Sink' (Order{isink}) (Parameter{ipara}) '_Norm']); close all;
        
        %finally, we run the statistics. The stats will only appear in the
        %command window if they are significant
       
        com_AnvsAw = horzcat(Aw_data, An_data);
        com_AnvsM = horzcat(An_data, M_data);
        
        disp('**Anesthetized vs Awake**')
        AnesthetizedvsAwakeN = teg_repeated_measures_ANOVA(com_AnvsAw, levels, varinames);
        AvK_CDN = iMakeCohensD(Aw_data, An_data)
                    
        disp('**Anesthetized vs Muscimol**')
        AnesthetizedvsMuscimolN = teg_repeated_measures_ANOVA(com_AnvsM, levels, varinames);
        MvK_CDN = iMakeCohensD(M_data, An_data)
        
        %significant or not, the F = #(#,#) and p= # will be saves in this file
        cd(homedir); cd DATA
        mkdir('StatsST'); cd StatsST
        savefile = [(Order{isink}) (Parameter{ipara}) 'SinkStats_Norm.mat'];
        save(savefile, 'AnesthetizedvsAwakeN', 'AvK_CDN', 'AnesthetizedvsMuscimolN', 'MvK_CDN')
    end
    
end

%% Normalized Mirror Sink Group
ParameterMir = {'SinkRMS','SinkPeakAmp'};
for isink = 1:length(Order)
    
    disp(['***********************Mirror Normalized ST STATS FOR ' (Order{isink}) '***********************'])
    %% Avg Trial Loop
    for ipara = 1:length(ParameterMir)
        cd(homedir); cd figs; cd('Group Tuning Curves ST');
        disp(['********Mirror Normalized ' (ParameterMir{ipara}) '********']) %these headers allow a copy and paste of all post-hoc tests which aren't saved in the data at the end
        
        % Below are the data, means, and sems of the current parameter and
        % sink
        An_data = vertcat(Anesthetized.ST_based.(ParameterMir{ipara}).(Order{isink})(:,6:10)); %, rowofnans
        An_data(An_data == 0) = 1;
        An_data = An_data./An_data(:,3);
        % stack the left under the right and double the central colume
        An_Mir_data = vertcat(An_data(:,3:5),(horzcat(An_data(:,3),An_data(:,2),An_data(:,1))));
        An_datamean = nanmean(An_Mir_data,1);
        An_datastd = nanstd(An_Mir_data,1);
        
        Aw_data = vertcat(Awake.ST_based.(ParameterMir{ipara}).(Order{isink})(:,5:9), rowofnans,rowofnans);
        Aw_data(Aw_data == 0) = 1;
        Aw_data = Aw_data./Aw_data(:,3);
        % stack the left under the right and double the central colume
        Aw_Mir_data = vertcat(Aw_data(:,3:5),(horzcat(Aw_data(:,3),Aw_data(:,2),Aw_data(:,1))));
        Aw_datamean = nanmean(Aw_Mir_data,1);
        Aw_datastd = nanstd(Aw_Mir_data,1);
        
        %now we generate a figure based on these three groups
        h=figure;
        
        errorbar(1:3, An_datamean, An_datastd, 'Linewidth', 2); hold on;
        errorbar(1:3, Aw_datamean, Aw_datastd, 'Linewidth', 2);
        legend('Anesthetized', 'Awake', 'Muscimol')
        title(['Sink' (Order{isink}) (ParameterMir{ipara})])
        set(gca,'XTick',1:5); set(gca,'XTickLabel',ticks,'FontSize',10);
       
        savefig(h,['Sink' (Order{isink}) (ParameterMir{ipara}) '_Mir_Norm']); close all;
        
        %finally, we run the statistics. The stats will only appear in the
        %command window if they are significant
       
        com_AnvsAw = horzcat(Aw_Mir_data, An_Mir_data);
        
        disp('**Anesthetized vs Awake**')
        AnesthetizedvsAwakeNM = teg_repeated_measures_ANOVA(com_AnvsAw, levelsMir, varinames);
        AvK_CDNM = iMakeCohensD(Aw_Mir_data, An_Mir_data)
        
        %significant or not, the F = #(#,#) and p= # will be saved in this file
        cd(homedir); cd DATA
        mkdir('StatsST'); cd StatsST
        savefile = [(Order{isink}) (ParameterMir{ipara}) 'SinkStats_Mir_Norm.mat'];
        save(savefile, 'AnesthetizedvsAwakeNM', 'AvK_CDNM')
    end
    
end

%% Normalized Mirror Sink Group no BF
levelsMir2bin = [2,2];
for isink = 1:length(Order)
    
    disp(['***********************Mirror Normalized ST STATS FOR ' (Order{isink}) ' 2 BINS ***********************'])
    %% Avg Trial Loop
    for ipara = 1:length(ParameterMir)
        cd(homedir); cd figs; cd('Group Tuning Curves ST');
        disp(['********Mirror Normalized ' (ParameterMir{ipara}) ' 2 BINS********']) %these headers allow a copy and paste of all post-hoc tests which aren't saved in the data at the end
        
        % Below are the data, means, and sems of the current parameter and
        % sink
        An_data = vertcat(Anesthetized.ST_based.(ParameterMir{ipara}).(Order{isink})(:,6:10)); %, rowofnans
        An_data(An_data == 0) = 1;
        An_data = An_data./An_data(:,3);
        % stack the left under the right and double the central colume
        An_Mir_data = vertcat(An_data(:,4:5),(horzcat(An_data(:,2),An_data(:,1))));
        
        Aw_data = vertcat(Awake.ST_based.(ParameterMir{ipara}).(Order{isink})(:,5:9), rowofnans,rowofnans);
        Aw_data(Aw_data == 0) = 1;
        Aw_data = Aw_data./Aw_data(:,3);
        % stack the left under the right and double the central colume
        Aw_Mir_data = vertcat(Aw_data(:,4:5),(horzcat(Aw_data(:,2),Aw_data(:,1))));
         
        %finally, we run the statistics. The stats will only appear in the
        %command window if they are significant
       
        com_AnvsAw = horzcat(Aw_Mir_data, An_Mir_data);
        
        disp('**Anesthetized vs Awake**')
        AnesthetizedvsAwakeNM = teg_repeated_measures_ANOVA(com_AnvsAw, levelsMir2bin, varinames);
        AvK_CDNM = iMakeCohensD(Aw_Mir_data, An_Mir_data)
                    
        %significant or not, the F = #(#,#) and p= # will be saved in this file
        cd(homedir); cd DATA
        mkdir('StatsST'); cd StatsST
        savefile = [(Order{isink}) (ParameterMir{ipara}) 'SinkStats_Mir_Norm_2bin.mat'];
        save(savefile, 'AnesthetizedvsAwakeNM', 'AvK_CDNM')
    end
    
end
