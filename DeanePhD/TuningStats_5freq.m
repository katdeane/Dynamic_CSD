%% Stats Loop FOR 5 FREQUENCIES
%This loop is an interface for the teg_repeated_measures_Anova code. It
%runs through each parameter (pulled from group analysis Data structure) in
%each specified sink. It generates simple errorbar line plots in a group 
%comparison centered around the BF. Variables to change, just
%below, are varinames (in our case groups and frequencies), levels (how many
%groups and how many frequencies), sink, and Parameter based on what is
%needed for comparisons. Changing the specific input will require renaming
%variable within the code itself. 

%% Start
clear
cd('D:\MyCode\Dynamic_CSD_Analysis');
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

%for teg_repeated measures_ANOVA
levels = [2,5];
levelsMir = [2,3];
varinames = {'Groups','Frequencies'};
rowofnans = NaN(1,5);
% srowofnans = [{NaN(1,50)} {NaN(1,50)} {NaN(1,50)} {NaN(1,50)} {NaN(1,50)}];

ticks = {'-2' '-1' 'BF' '+1' '+2'};

sink = {'I_II','IV','V','V'};
Parameter = {'SinkRMS','SinkPeakAmp','SinkPeakLate','Sinkonset'};
% SParameter = {'SingleSinkRMS','SingleSinkPeakAmp','SingleSinkPeakLat'};
Measurement = {'Pre_4','preAMtono_1','preAMtono_4','preCLtono_1',...
    'preCLtono_4','CLtono_1','AMtono_1'};


%% Load in the appropriate files
cd DATA;cd Output;
load('KIC_Tuning_Avg.mat');
Control = Tuning; clear Tuning; 
load('KIT_Tuning_Avg.mat');
Treated = Tuning; clear Tuning;
load('KIV_Tuning_Avg.mat');
VControl = Tuning; clear Tuning;

cd(homedir); cd DATA
mkdir('Stats_Crypto'); cd Stats_Crypto

%% Measurement Loop
for imeas = 1:length(Measurement)
    %% Sink Loop
    for isink = 1:length(sink)
        
        disp(['***********************ST STATS FOR ' (sink{isink}) ' ' (Measurement{imeas}) '***********************'])
        %% Avg Trial Loop
        for ipara = 1:length(Parameter)
            
            disp(['********' (Parameter{ipara}) '********']) %these headers allow a copy and paste of all post-hoc tests which aren't saved in the data at the end
            
            % Below are the data, means, and sems of the current parameter and
            % sink
            C_data = vertcat(Control.(sink{1}).(Measurement{imeas}).(Parameter{ipara})(:,6:10),...
                rowofnans,rowofnans,rowofnans,rowofnans); %, rowofnans
            C_datamean = nanmean(C_data,1);
            C_datasem = nanstd(C_data,1)/sqrt(sum(~isnan(C_datamean)));
            
            T_data = vertcat(Treated.(sink{1}).(Measurement{imeas}).(Parameter{ipara})(:,6:10));
            T_datamean = nanmean(T_data,1);
            T_datasem = nanstd(T_data,1)/sqrt(sum(~isnan(T_datamean)));
            
            V_data = vertcat(VControl.(sink{1}).(Measurement{imeas}).(Parameter{ipara})(:,6:10),...
                rowofnans,rowofnans,rowofnans,rowofnans);
            V_datamean = nanmean(V_data,1);
            V_datasem = nanstd(V_data,1)/sqrt(sum(~isnan(V_datamean)));
            
            
            %Run the statistics. The stats will only appear in the
            %command window if they are significant
            com_CvsT = horzcat(T_data, C_data);
            com_VvsT = horzcat(T_data, V_data);
            
            disp('**Treated vs Naive Control**')
            TreatedvsControl = teg_repeated_measures_ANOVA(com_CvsT, levels, varinames);
            TvC_CD = iMakeCohensD(T_data, C_data)
            
            disp('**Treated vs Virus Control**')
            TreatedvsVControl = teg_repeated_measures_ANOVA(com_VvsT, levels, varinames);
            TvV_CD = iMakeCohensD(T_data, V_data)
            
            %significant or not, the F = #(#,#) and p= # will be saves in this file
            savefile = [(sink{isink}) '_' (Parameter{ipara}) '_' (Measurement{imeas}) '_' 'SinkStats.mat'];
            save(savefile, 'TreatedvsControl','TvC_CD','TreatedvsVControl','TvV_CD')
        end %parameter
    end %sink
end %measurement

%% Normalized Sink Group

for imeas = 1:length(Measurement)
    for isink = 1:length(sink)
        
        disp(['***********************Normalized ST STATS FOR ' (sink{isink}) (Measurement{imeas}) '***********************'])
        for ipara = 1:length(Parameter)
            
            disp(['********Normalized ' (Parameter{ipara}) '********']) %these headers allow a copy and paste of all post-hoc tests which aren't saved in the data at the end
            
            % Below are the data, means, and sems of the current parameter and
            % sink
            C_data = vertcat(Control.(sink{1}).(Measurement{imeas}).(Parameter{ipara})(:,6:10),...
                rowofnans,rowofnans,rowofnans,rowofnans); %, rowofnans
            C_data(C_data == 0) = 1;
            C_data = C_data./C_data(:,3);
            C_datamean = nanmean(C_data,1);
            C_datasem = nanstd(C_data,1)/sqrt(sum(~isnan(C_datamean)));
            
            T_data = vertcat(Treated.(sink{1}).(Measurement{imeas}).(Parameter{ipara})(:,6:10));
            T_data(T_data == 0) = 1;
            T_data = T_data./T_data(:,3);
            T_datamean = nanmean(T_data,1);
            T_datasem = nanstd(T_data,1)/sqrt(sum(~isnan(T_datamean)));
                        
            V_data = vertcat(VControl.(sink{1}).(Measurement{imeas}).(Parameter{ipara})(:,6:10),...
                rowofnans,rowofnans,rowofnans,rowofnans);
            V_data(V_data == 0) = 1;
            V_data = V_data./V_data(:,3);
            V_datamean = nanmean(V_data,1);
            V_datasem = nanstd(V_data,1)/sqrt(sum(~isnan(V_datamean)));
            
            %Run the statistics. The stats will only appear in the
            %command window if they are significant
            
            com_CvsT = horzcat(T_data, C_data);
            com_VvsT = horzcat(T_data, V_data);
            
            disp('**Treated vs Naive Control**')
            TreatedvsControlN = teg_repeated_measures_ANOVA(com_CvsT, levels, varinames);
            TvC_CD = iMakeCohensD(T_data, C_data)
                        
            disp('**Treated vs Virus Control**')
            TreatedvsVControlN = teg_repeated_measures_ANOVA(com_VvsT, levels, varinames);
            TvV_CD = iMakeCohensD(T_data, V_data)
            
            %significant or not, the F = #(#,#) and p= # will be saves in this file
            savefile = [(sink{isink}) '_' (Parameter{ipara}) '_' (Measurement{imeas}) '_' 'SinkStats_Norm.mat'];
            save(savefile, 'TreatedvsControlN','TvC_CD','TreatedvsVControlN','TvV_CD')
        end
    end
end

%% Normalized Mirror Sink Group
for imeas = 1:length(Measurement)
    for isink = 1:length(sink)
        
        disp(['***********************Mirror Normalized ST STATS FOR ' (sink{isink}) (Measurement{imeas}) '***********************'])
        for ipara = 1:length(Parameter)
            
            disp(['********Mirror Normalized ' (Parameter{ipara}) '********']) %these headers allow a copy and paste of all post-hoc tests which aren't saved in the data at the end
            
            % Below are the data, means, and sems of the current parameter and
            % sink
            C_data = vertcat(Control.(sink{1}).(Measurement{imeas}).(Parameter{ipara})(:,6:10),...
                rowofnans,rowofnans,rowofnans,rowofnans); 
            C_data(C_data == 0) = 1;
            C_data = C_data./C_data(:,3);
            % stack the left under the right and double the central colume
            C_Mir_data = vertcat(C_data(:,3:5),(horzcat(C_data(:,3),C_data(:,2),C_data(:,1))));
            C_datamean = nanmean(C_Mir_data,1);
            C_datasem = nanstd(C_Mir_data,1)/sqrt(sum(~isnan(C_datamean)));
            
            T_data = vertcat(Treated.(sink{1}).(Measurement{imeas}).(Parameter{ipara})(:,6:10));
            T_data(T_data == 0) = 1;
            T_data = T_data./T_data(:,3);
            % stack the left under the right and double the central colume
            T_Mir_data = vertcat(T_data(:,3:5),(horzcat(T_data(:,3),T_data(:,2),T_data(:,1))));
            T_datamean = nanmean(T_Mir_data,1);
            T_datasem = nanstd(T_Mir_data,1)/sqrt(sum(~isnan(T_datamean)));
            
            V_data = vertcat(VControl.(sink{1}).(Measurement{imeas}).(Parameter{ipara})(:,6:10),...
                rowofnans,rowofnans,rowofnans,rowofnans);
            V_data(V_data == 0) = 1;
            V_data = V_data./V_data(:,3);
            % stack the left under the right and double the central colume
            V_Mir_data = vertcat(V_data(:,3:5),(horzcat(V_data(:,3),V_data(:,2),V_data(:,1))));
            V_datamean = nanmean(V_Mir_data,1);
            V_datasem = nanstd(V_Mir_data,1)/sqrt(sum(~isnan(V_datamean)));
            
            %Run the statistics. The stats will only appear in the
            %command window if they are significant
            
            com_CvsT = horzcat(T_Mir_data, C_Mir_data);
            com_VvsT = horzcat(T_Mir_data, V_Mir_data);
            
            disp('**Treated vs Naive Control**')
            TreatedvsControlNM = teg_repeated_measures_ANOVA(com_CvsT, levelsMir, varinames);
            TvC_CD = iMakeCohensD(T_Mir_data, C_Mir_data)
            
            disp('**Treated vs Virus Control**')
            TreatedvsVControlNM = teg_repeated_measures_ANOVA(com_VvsT, levelsMir, varinames);
            TvV_CD = iMakeCohensD(T_Mir_data, V_Mir_data)
            
            %significant or not, the F = #(#,#) and p= # will be saves in this file
            savefile = [(sink{isink}) '_' (Parameter{ipara}) '_' (Measurement{imeas}) '_' 'SinkStats_Norm.mat'];
            save(savefile, 'TreatedvsControlNM','TvC_CD','TreatedvsVControlNM','TvV_CD')
            
        end
    end
end
