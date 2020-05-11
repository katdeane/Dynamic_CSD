function Dynamic_spike(homedir)

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
cd (homedir),cd groups;

%% Load in
input = dir('*.m');
entries = length(input);

% loop through number of Data mats in folder
for i1 = 1:entries
    
    % get group file for animal details
    run(input(i1).name);
    
    %% Choose Condition
    Condition = {'Pre' 'preAM' 'preAMtono' 'preCL' 'preCLtono' 'spPre1' 'spPost1' ...
        'CL' 'CLtono' 'spPre2' 'spPost2' 'AM' 'AMtono' 'spEnd'}; 
%     Condition = {'preCL' 'CL'};
%     Condition = {'Pre'};
    
    %% Condition and Indexer
    Data = struct;
    cd (homedir); cd figs;
    mkdir(['Single_Spike_' input(i1).name(1:end-2)]);
    
    Indexer = imakeIndexer(Condition,animals,Cond);
    
    % loop through animal subjects
    for iAn = 1:length(animals)
        cd (homedir); cd raw;
        name = animals{iAn};
        
        % loop through conditions from condition list
        for iCond = 1:length(Condition)
            
            % loop through measurements in each condition
            for iMeas = 1:length(Cond.(Condition{iCond}){iAn})
                if iMeas == 1
                    CondIDX = Indexer(2).(Condition{iCond});
                else
                    CondIDX = Indexer(2).(Condition{iCond})+iMeas-1;
                end
                tic
                
                % load in measurement
                measurement = Cond.(Condition{iCond}){iAn}{iMeas};
                if ~isempty(measurement)
                    disp(['Analyzing animal: ' name '_' measurement])
                    try
                        load ([name '_' measurement]);
                    catch
                        fprintf('the name or measurement does not exist/n')
                    end
                    
                    % generate some important features
                    Fs          = P.Fs_AD(1); %ms  sampling rate
                    BL          = Header.t_pre*Fs; %ms  baseline
                    last_click  = Header.t_sig.*Fs-2; %ms  time of last click
                    frqz        = Header.stimlist(:,3)'; %Hz  set of click frequencies
                    chan_order  = str2num(channels{iAn});
                    StimID      = unique(cat(SWEEP.Header));
                    numTrials   = max(Header.n_rep); 
                    layers      = {'I_II','IV','V','VI','All'};
                    
                    if contains(Condition{iCond},'tono') 
                        StimLabel   = {'1 kHz', '2 kHz', '4 kHz', '8 kHz', '16 kHz', '32 kHz', 'pause', 'click'};
                    elseif contains(Condition{iCond},'CL') || contains(Condition{iCond},'AM')
                        StimLabel   = {'2 Hz', '5 Hz', '10 Hz', '20 Hz', '40 Hz'};
                    elseif contains(Condition{iCond},'sp')
                        StimLabel   = {'pause', 'pause', 'pause', 'pause', 'pause'};
                    else
                        StimLabel   = {'1 kHz', '2 kHz', '4 kHz', '8 kHz', '16 kHz', '32 kHz', 'pause', 'click'};
                    end
                    
                       
                    % sort sweep spikes list into stimuli groups
                    SinglSpikes = cell(numTrials,length(StimID));
                    SinglTimes = cell(numTrials,length(StimID));
                    AvgSpikes = cell(1,length(StimID));
                    SumSpikes = cell(1,length(StimID));
                    for iStim = 1:length(StimID)
                        ii = 0;
                        for iTrial = 1:size(SWEEP,1)
                            if SWEEP(iTrial).Header == iStim
                                
                                [trial_spikes, trial_times] = SweepSpikeChan(SWEEP(iTrial).Spikes,chan_order,Fs);
                                
                                ii = ii + 1;
                                SinglSpikes{ii,iStim} = trial_spikes;
                                SinglTimes{ii,iStim} = trial_times;
                            end
                        end
                       AvgSpikes{1,iStim} = mean(cat(3,SinglSpikes{:,iStim}),3);
                       SumSpikes{1,iStim} = sum(cat(3,SinglSpikes{:,iStim}),3);  
                    end
                    
                    
                    
                    %% Visualize Spikes 
%                     
%                     cd(homedir); cd figs;
%                     cd(['Single_Spike_' input(i1).name(1:end-2)]);
%                     
%                     % run through stimuli
%                     for iStim = 1:length(StimID)
%                         Stim_times = horzcat(SinglTimes{:,iStim});
%                         
%                         % run through layers
%                         for iLay = 1:length(layers)
%                             
%                             % condition to skip if the layer is empty
%                             if ~strcmp(layers{iLay},'All')
%                                 if strcmp('[]',Layer.(layers{iLay}){1,iAn})
%                                     continue
%                                 end
%                             end
%                             
%                             % function to concatonate spike data
%                             layer_spikes = StimLayerCat(Stim_times, Layer, ...
%                                 layers{iLay}, iAn, length(chan_order));
%                             
%                             % function to display spike as psth
%                             ph = psth(layer_spikes, 5, Fs, 49, 1400);
%                             h = gcf;
%                             set(h,'Name',[name '_' measurement ': ' Condition{iCond} ...
%                                 '_' num2str(iMeas) ' ' layers{iLay}])
%                             set(h, 'PaperType', 'A4');
%                             set(h, 'PaperOrientation', 'landscape');
%                             set(h, 'PaperUnits', 'centimeters');
%                             savefig(h,['Spike_PSTH ' name '_' measurement '_' ...
%                                 Condition{iCond} '_' num2str(iMeas) ' ' ...
%                                 StimLabel{iStim} ' ' layers{iLay}],'compact')
%                             saveas(h,['Spike_PSTH ' name '_' measurement '_' ...
%                                 Condition{iCond} '_' num2str(iMeas) ' ' ...
%                                 StimLabel{iStim} ' ' layers{iLay} '.pdf'])
%                             close (h)
%                         end
%                     end
                    
                    %% save and quit
                    Data(CondIDX).(name).measurement    = [name '_' measurement];
                    Data(CondIDX).(name).Condition      = [Condition{iCond} '_' num2str(iMeas)];
                    Data(CondIDX).(name).BL             = BL;
                    Data(CondIDX).(name).StimDur        = last_click;
                    Data(CondIDX).(name).Frqz           = frqz;
                    Data(CondIDX).(name).SinglTimes     = SinglTimes;
                    Data(CondIDX).(name).SinglSpikes    = SinglSpikes;
                    Data(CondIDX).(name).AvgSpikes      = AvgSpikes;
                    Data(CondIDX).(name).SumSpikes      = SumSpikes;
                    
                end
            end
        end
    end
    
    cd ([homedir '/DATA'])
    save(['Spikes_' input(i1).name(1:end-2) '_Data'],'Data');
    clear Data
                    
end