function ChangeInSpikes(homedir)

% This script takes *.mat files out of the DATA/ folder. It checks the
% condition names and finds the measurements associated with both clicks
% and AMs (seperately). It then plots the Avrecs for each animal over each
% frequency presentation. 

%Input:     D:\MyCode\Dynamic_CSD_Analysis\DATA -> Spikes_*_Data.mat 
%Output:    Figures of in "ChangeIn_Spikes"

clear
%% standard operations
warning('OFF');
dbstop if error

% Change directory to your working folder
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
cd (homedir),cd DATA;

%% Load in
input = dir('*.mat');
entries = length(input);
layers = {'All', 'I_II', 'IV', 'V', 'VI'};
stimlist = [2,5,10,20,40];

% set up simple cell sheets to hold all data: avrec of total/layers and
% peaks of pre conditions
SpikesAll = cell(length(stimlist),length(layers),entries);
PeakofPre = zeros(length(stimlist),length(layers),entries);
AnCount = 1; % total number of animals counted

SP_SpikesAll = cell(length(layers),entries);
SP_PeakofPre = zeros(length(layers),entries);
SP_AnCount = 1;

% loop through number of Data mats in folder
for i1 = 1:entries
    
    if ~contains(input(i1).name,'Spikes')
        continue
    end
    
    clear Data
    load(input(i1).name);
    names = fieldnames(Data); 
    thisG = [input(i1).name(8:10) '.m'];
    
    cd (homedir),cd groups;
    run(thisG);

    cd (homedir),cd figs;
    mkdir Single_Spikes; cd Single_Spikes;
    
    %% Clicks    
    for iAn = 1:length(names)
        % 5 figures per animal
        
        for iLay = 1:length(layers)
            
            SpikeCurves = figure('Name',['Spikes_Clicks_' names{iAn}],'Position',[-1000 100 800 1200]);
            
            for iStim = 1:length(stimlist)
                % create container for lables
                CondN = cell(1,size(Data,2));
                
                
                % 1 subplot per stimulus
                subplot(length(stimlist),1,iStim);
                allmeas = [];
                
                for iMeas = 1:size(Data,2)
                    if isempty(Data(iMeas).(names{iAn}))
                        continue
                    end
                    
                    if ~contains((Data(iMeas).(names{iAn}).Condition),'CL_')
                        continue
                    end
                    % take an average of all channels (already averaged across trials)
                    if contains(layers{iLay}, 'All')
                        avgchan = mean(Data(iMeas).(names{iAn}).AvgSpikes{1, iStim});
                    else
                        avgchan = mean(Data(iMeas).(names{iAn}).AvgSpikes{1, iStim}(str2num(Layer.(layers{iLay}){iAn}),:));
                    end
                    % smooth the data by gaussian window:
                    g_avgall = smoothdata(avgchan,'gaussian',5);
                    
                    % pull out condition
                    CondN{1,iMeas} = Data(iMeas).(names{iAn}).Condition;
                    % plot it
                    plot(g_avgall, 'LineWidth', 1)
                    hold on
                    % store this avg temporarily with buddies
                    allmeas = vertcat(allmeas,g_avgall);
                    
                    % store peak if preCL condition
                    if contains(CondN{1,iMeas},'preCL')
                        PeakofPre(iStim,iLay,AnCount) = max(avgchan);
                    end
                end
                %
                % and store the lot
                SpikesAll{iStim,iLay,AnCount} = allmeas;
                                
                CondN = CondN(~cellfun('isempty',CondN));
                legend(CondN)
                title([num2str(stimlist(iStim)) ' Hz'])
                hold off
                
            end
            h = gcf;
            savefig(h,['Spikes_Clicks_' layers{iLay} '_' names{iAn}],'compact')
            saveas(h,['Spikes_Clicks_' layers{iLay} '_' names{iAn} '.pdf'])
            close (h)
        end
        
        AnCount = AnCount + 1;
    end
    
    %% Amplitude Modulation
%     for iAn = 1:length(names)
%         % 1 figure per animal
%         SpikeCurves = figure('Name',['Spikes_AMs_' names{iAn}],'Position',[-1000 100 800 1200]);
%         
%         stimlist = [2,5,10,20,40];
%         
%         for iStim = 1:length(stimlist)
%             % create container for lables
%             CondN = cell(1,size(Data,2));
%             
%             % 1 subplot per stimulus
%             subplot(length(stimlist),1,iStim);
%             
%             for iMeas = 1:size(Data,2)
%                 if isempty(Data(iMeas).(names{iAn}))
%                     continue
%                 end
%                 
%                 if ~contains((Data(iMeas).(names{iAn}).Condition),'AM_')
%                     continue
%                 end
%                 % take an average of all channels (already averaged across trials)
%                 
%                 avgchan = mean(Data(iMeas).(names{iAn}).AvgSpikes{1, iStim});
%                 % smooth the data by gaussian window:
%                 g_avgall = smoothdata(avgchan,'gaussian',5);
%                 % pull out condition
%                 CondN{1,iMeas} = Data(iMeas).(names{iAn}).Condition;
%                 % plot it
%                 plot(g_avgall, 'LineWidth', 1)
%                 hold on
%             end
%             
%             CondN = CondN(~cellfun('isempty',CondN));
%             legend(CondN)
%             title([num2str(stimlist(iStim)) ' Hz'])
%             hold off
%             
%             
%         end
%         h = gcf;
%         savefig(h,['Spikes_AMs_' names{iAn}],'compact')
%         saveas(h,['Spikes_AMs_' names{iAn} '.pdf'])
%         close (h)
%     end
%% Spontaneous

for iAn = 1:length(names)
    % 5 figures per animal
    
    for iLay = 1:length(layers)
        
        SpikeCurves = figure('Name',['Spikes_Spontaneous_Clicks_' names{iAn}],'Position',[-1000 100 800 400]);
        
        CondN = cell(1,size(Data,2));
        
        % 1 subplot per stimulus
        allmeas = [];
        
        for iMeas = 1:size(Data,2)
            if isempty(Data(iMeas).(names{iAn}))
                continue
            end
            
            if ~contains((Data(iMeas).(names{iAn}).Condition),'sp')
                continue
            end
            % take an average of all channels (already averaged across trials)
            if contains(layers{iLay}, 'All')
                avgchan = mean(vertcat(Data(iMeas).(names{iAn}).AvgSpikes{1, :}));
            else
                chan1 = Data(iMeas).(names{iAn}).AvgSpikes{1, 1}(str2num(Layer.(layers{iLay}){iAn}),:);
                chan2 = Data(iMeas).(names{iAn}).AvgSpikes{1, 2}(str2num(Layer.(layers{iLay}){iAn}),:);
                chan3 = Data(iMeas).(names{iAn}).AvgSpikes{1, 3}(str2num(Layer.(layers{iLay}){iAn}),:);
                chan4 = Data(iMeas).(names{iAn}).AvgSpikes{1, 4}(str2num(Layer.(layers{iLay}){iAn}),:);
                chan5 = Data(iMeas).(names{iAn}).AvgSpikes{1, 5}(str2num(Layer.(layers{iLay}){iAn}),:);
                avgchan = mean(vertcat(chan1,chan2,chan3,chan4,chan5));
            end
            % smooth the data by gaussian window:
            g_avgall = smoothdata(avgchan,'gaussian',5);
            
            % pull out condition
            CondN{1,iMeas} = Data(iMeas).(names{iAn}).Condition;
            % plot it
            plot(g_avgall, 'LineWidth', 1)
            hold on
            % store this avg temporarily with buddies
            allmeas = vertcat(allmeas,g_avgall);
            
            % store peak if preCL condition
            if contains(CondN{1,iMeas},'spPre')
                SP_PeakofPre(iLay,SP_AnCount) = max(avgchan);
            end
        end
        
        if isempty(allmeas)
            h = gcf; close(h)
            continue
        end
        
        % and store the lot
        SP_SpikesAll{iLay,SP_AnCount} = allmeas;
        
        CondN = CondN(~cellfun('isempty',CondN));
        legend(CondN)
        title('Spontaneous')
        hold off
        
        
        h = gcf;
        savefig(h,['Spikes_Spontaneous_Clicks_' layers{iLay} '_' names{iAn}],'compact')
        saveas(h,['Spikes_Spontaneous_Clicks_' layers{iLay} '_' names{iAn} '.pdf'])
        close (h)
    end
    
    SP_AnCount = SP_AnCount + 1;
end

    cd(homedir); cd DATA;
end

% save it out
cd (homedir),cd figs;
mkdir Group_Spikes; cd Group_Spikes;
% save('SpikesAll','SpikesAll','PeakofPre');
save('Spont_SpikesAll','SP_SpikesAll','SP_PeakofPre');