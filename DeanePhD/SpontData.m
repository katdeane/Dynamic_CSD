%% Spontaneous Activity

% What do I want out of this code? 
% I want a measure of how many spikes occur in a given time period (2
% minutes) so that I can compare that spiking activity to my stimulated
% periods. 

% And what do I have to put into this code?
% Some of the plx files converted and some didn't, although I don't know
% why. I have all of those files and just need to see what the difference
% is and which is better to start from. 

% one of the converted ones:
% load('D:\MyCode\Dynamic_CSD_Analysis\raw\KIC01_spont-pre01.mat')
% SWEEP contains no Header
% spike list is very tiny... could be this animal which had a terrible
% signal anyway. 

%% Let's start by pretending this is a CSD analysis. 

if exist('D:\MyCode\Dynamic_CSD_Analysis','dir') == 7
    cd('D:\MyCode\Dynamic_CSD_Analysis');
elseif exist('D:\Dynamic_CSD_Analysis','dir') == 7
    cd('D:\Dynamic_CSD_Analysis');
elseif exist('C:\Users\kedea\Documents\Dynamic_CSD_Analysis','dir') == 7
    cd('C:\Users\kedea\Documents\Dynamic_CSD_Analysis')
end

home = pwd; 
addpath(genpath(home));
cd (home),cd groups;

%% Load in
input = dir('*.m');
entries = length(input);


for i1 = 1:entries
    
    run(input(i1).name);
    Condition = {'spPre1' 'spPost1' 'spPre2' 'spPost2' 'spEnd'};
    
    Data = struct;
    cd (home); cd figs;
    mkdir(['Spont_' input(i1).name(1:end-2)]);
    
    Indexer = imakeIndexer(Condition,animals,Cond);
    %% Spont CSD 
    
    for iA = 1:length(animals)
        cd (home); cd raw;
        name = animals{iA};
        
        for iC = 1:length(Condition)
            for i4 = 1:length(Cond.(Condition{iC}){iA})
                if i4 == 1
                    CondIDX = Indexer(2).(Condition{iC});
                else
                    CondIDX = Indexer(2).(Condition{iC})+i4-1;
                end
                
                measurement = Cond.(Condition{iC}){iA}{i4};
                if ~isempty(measurement)
                    disp(['Analyzing animal: ' name '_' measurement])
                    clear SWEEP
                    try
                        load ([name '_' measurement]);clear avgFP;
                    catch
                        fprintf('the name or measurement does not exist/n')
                    end
                    
                    BL = 0;
                    Fs = P.Fs_AD(1); %sampling rate
                    chan_order  = str2num(channels{iA}); 
                    
                    [SingleLayer_AVREC,~,~, SingleTrialFP, ~,...
                        SingleTrialCSD, ~, SingleTrialAvgRecCSD,...
                        SingleTrialRelResCSD, ~,~,...
                        SingleTrialAbsResCSD, LayerRelRes, ~] =...
                        SingleTrialCSD_full(SWEEP, chan_order,1:length(chan_order),BL);

                    %delete empty columns to have the correct amount of stimuli
                    %present (needed for attenuation 30)
                    SingleLayer_AVREC = SingleLayer_AVREC(~cellfun('isempty',SingleLayer_AVREC));
                    SingleTrialFP = SingleTrialFP(~cellfun('isempty',SingleTrialFP));
                    SingleTrialCSD = SingleTrialCSD(~cellfun('isempty',SingleTrialCSD));
                    SingleTrialRelResCSD = SingleTrialRelResCSD(~cellfun('isempty',SingleTrialRelResCSD));
                    SingleTrialAbsResCSD = SingleTrialAbsResCSD(~cellfun('isempty',SingleTrialAbsResCSD));
                    LayerRelRes = LayerRelRes(~cellfun('isempty',LayerRelRes));
                    SingleTrialAvgRecCSD = SingleTrialAvgRecCSD(~cellfun('isempty',SingleTrialAvgRecCSD));
                   
                    
                    % pull out 1 second bins and stack them in cells
                    len = 1000; % 1s
                    STCSDps = cell(1,120000/len);
                    count = 0; 
                    for icsd = 1:len:size(SingleTrialCSD{1},2)
                        count = count + 1;
                        STCSDps{1,count} = SingleTrialCSD{1}(:,icsd:icsd+len-1);
                    end
                    
                    % averaging over 1 second windows
                    AvgCSDps = mean(cat(3,STCSDps{:}),3);
                    
                    % ploting the average
                    figure('Name',[name ' ' measurement ])
                    disp('Plotting CSD')
                    
                    imagesc(AvgCSDps)
                    title ([name ' ' measurement ': ' Condition{iC} ' Hz'])
                    colormap (jet)
                    caxis([-0.0005 0.0005])
                    

                    cd([home '\figs\' 'Spont_' input(i1).name(1:end-2)])
                    h = gcf;
                    set(h, 'PaperType', 'A4');
                    set(h, 'PaperOrientation', 'landscape');
                    set(h, 'PaperUnits', 'centimeters');
                    savefig(h,[name '_' measurement '_CSDs' ],'compact')
                    close (h)
                    
                    %% Spont Spikes
                                        
                    % going through these layers to collect spike data
                    layers      = {'I_IIE','IVE','VaE','VbE','VIaE','VIbE','All'};
                                
                    [SinglSpikes, SinglTimes] = SweepSpikeChan(SWEEP.Spikes,chan_order,Fs);
                    
                    disp('Plotting Spikes')
                    % run through layers
                    for iLay = 1:length(layers)
                        
                        % condition to skip if the layer is empty
                        if ~strcmp(layers{iLay},'All')
                            if strcmp('[]',Layer.(layers{iLay}){1,iA})
                                continue
                            end
                        end
                        
                        % function to concatonate spike data
                        layer_spikes = StimLayerCat(SinglTimes, Layer, ...
                            layers{iLay}, iA, length(chan_order));
                        
                        % function to display spike as psth
                        ph = psth(layer_spikes, 5, Fs, 120, 1000);
                        h = gcf;
                        set(h,'Name',[name '_' measurement ': ' Condition{iC} ...
                            '_' layers{iLay}])
                        set(h, 'PaperType', 'A4');
                        set(h, 'PaperOrientation', 'landscape');
                        set(h, 'PaperUnits', 'centimeters');
                        savefig(h,['Spike_PSTH ' name '_' measurement '_' ...
                            Condition{iC} '_' layers{iLay}],'compact')
                        saveas(h,['Spike_PSTH ' name '_' measurement '_' ...
                            Condition{iC} '_' layers{iLay} '.pdf'])
                        close (h)
                    end
                    
                    
                    % main
                    Data(CondIDX).(name).measurement =[name '_' measurement];
                    Data(CondIDX).(name).Condition = [Condition{iC} '_' num2str(i4)];
                    
                    % csd etc. 
                    Data(CondIDX).(name).singletrialLFP = SingleTrialFP;
%                     Data(CondIDX).(name).LFP = AvgFP;
                    Data(CondIDX).(name).singletrialCSD =SingleTrialCSD;
                    Data(CondIDX).(name).ST_CSD_persecond = STCSDps;
                    Data(CondIDX).(name).AvgCSD_persecond = AvgCSDps;
%                     Data(CondIDX).(name).LayerRelRes = AvgLayerRelRes;
%                     Data(CondIDX).(name).AVREC_raw = AvgRecCSD;
                    Data(CondIDX).(name).SingleTrial_AVREC_raw = SingleTrialAvgRecCSD;
                    Data(CondIDX).(name).SingleTrial_RelRes_raw = SingleTrialRelResCSD;
                    Data(CondIDX).(name).SingleTrial_AbsRes_raw = SingleTrialAbsResCSD;
%                     Data(CondIDX).(name).RELRES_raw = AvgRelResCSD;
%                     Data(CondIDX).(name).ABSRES_raw = AvgAbsResCSD;

                    % spikes
                    Data(CondIDX).(name).SinglTimes     = SinglTimes;
                    Data(CondIDX).(name).SinglSpikes    = SinglSpikes;
%                     Data(CondIDX).(name).AvgSpikes      = AvgSpikes;
%                     Data(CondIDX).(name).SumSpikes      = SumSpikes;
                    
                end
            end
        end
    end
    cd ([home '/DATA'])
    save(['Spont_' input(i1).name(1:end-2) '_Data'],'Data');
    clear Data
end



















