% 03-05-2019 AC: first pass. set up basic analysis loop with figure output
% and accounting for Cone of Influence
% 05-05-2019 AC: Important updates:
%   - repackaging data to be gramm compatible
%   - adding option to params to suppress figure output
% 28-05-2019 AC: changing to store one scalogram for each animal altogether
%   - rel 2 sink center removed: average across 3 middle channels for each
%   layer taken
% 03-06-2019 AC: Removed restriction for doing layer 4 only
% 16-09-2019 AC: branched from original CSD_allLayers_scalogram to compute
% csd for single trials.

%% INIT PARAMETERS
clear
%close all
params.sampleRate = 1000; % Hz
params.startTime = -0.2; % seconds
params.timeLimits = [-0.2 0.399]; % seconds
params.frequencyLimits = [5 params.sampleRate/2]; % Hz
params.voicesPerOctave = 8;
params.timeBandWidth = 54;
% where do I find the group input scripts:
params.groupFold = 'C:\Users\kedea\Documents\Dynamic_CSD_Analysis\groups\';
% where do I find the split up CSD data files
params.datFold = 'C:\Users\kedea\Documents\Dynamic_CSD_Analysis\AndrewSpectralData\Data\';
params.doFig = false;

% Change these for specific runs
layerNames = {'I_IIE'}; 
whichtone = {'offBF'};
savename = 'STscalograms_IIoffBF';


%% ANEST MEASURE

[awake] = runCwtCsd('Awake10dB',params,layerNames,whichtone);
[anest] = runCwtCsd('AnesthetizedPre',params,layerNames,whichtone);
[musc]  = runCwtCsd('Muscimol',params,layerNames,whichtone);

%% Reorganize data into Gramm-compatible structure

wtTable = [struct2table(awake); struct2table(anest); struct2table(musc)];
save(['C:\Users\kedea\Documents\Dynamic_CSD_Analysis\AndrewSpectralData\Data\' savename],'wtTable')


%% SUBFUNCTIONS

function [datStruct] = runCwtCsd(groupFile,params,layerNames,whichtone)

run([params.groupFold groupFile '.m'])
% Generates variables called animals and Layer

% Init datastruct to pass out
datStruct = struct();
datStruct.scalogram     = [];
datStruct.animal        = [];
datStruct.freq          = [];
datStruct.condition     = [];
%datStruct.rel2Bf        = [];
datStruct.layer         = [];
datStruct.trial         = [];
%datStruct.channel       = [];

idx = 1;
for hh = 1:length(layerNames)
    for ii = 1:length(animals)
        
        thisAnimal = animals{ii};
        thisDatFold = [params.datFold thisAnimal '\' groupFile '\'];
        d = dir(thisDatFold);
        d = d(~[d.isdir]);
        bfFile = d(contains({d.name},whichtone)).name;
        
        % Format here forces the use of eval to get channel locs
        thisLayerChans = eval(Layer.(layerNames{hh}){ii});
        centerChan = thisLayerChans(ceil(length(thisLayerChans)/2));
        theseChans = thisLayerChans-centerChan;
        % Select only center 3 channels
        thisLayerChans = thisLayerChans(theseChans >=-1 & theseChans <=1);        
        
        load([thisDatFold bfFile]) 
        try 
            singleCSD = singleCSD_offBF;
        catch
            singleCSD = singleCSD_BF;
        end
        
        % Repeat over trials
        for ll = 1:size(singleCSD,3)
            
            hold_cwts = zeros(54,size(singleCSD,2));
            % Repeat over channels
            for kk = 1:length(thisLayerChans)
                
                csdChan = squeeze(singleCSD(thisLayerChans(kk),:,:));
                ROI = csdChan(:,ll);
                
                if params.doFig; sub = figure; subplot(2,ceil(size(csdChan,1)/2),1); end
                
                timeValues = params.startTime + (0:length(ROI)-1).'/params.sampleRate;
                % Limit the cwt frequency limits
                params.frequencyLimits(1) = max(params.frequencyLimits(1),...
                    cwtfreqbounds(numel(ROI),params.sampleRate,...
                    'TimeBandWidth',params.timeBandWidth));
                % Get the cone of influence in a very round-about way....
                %                 thisInd = bfToneInd(jj)-bfInd;
                %                 tit = sprintf('ID_%s_%dfromBF_Layer_%s_condition_%s',thisAnimal,thisInd, layerNames{hh},groupFile);
                %                 cwtfig = figure('Name',tit); cwt(ROI,params.sampleRate, ...
                %                     'VoicesPerOctave',params.voicesPerOctave, ...
                %                     'TimeBandWidth',params.timeBandWidth, ...
                %                     'FrequencyLimits',params.frequencyLimits);
                %                 ax = gca;
                %                 xdat = ax.Children(2).XData-200;
                %                 ydat = ax.Children(2).YData;
                
                [WT,F] = cwt(ROI,params.sampleRate, ...
                    'VoicesPerOctave',params.voicesPerOctave, ...
                    'TimeBandWidth',params.timeBandWidth, ...
                    'FrequencyLimits',params.frequencyLimits);
                
                % to get average of channels: add each and divide after loop
                hold_cwts = hold_cwts + WT;
            end
                
                WT = hold_cwts/3;
                
                datStruct(idx).scalogram            = WT;
                datStruct(idx).condition            = {groupFile};
                datStruct(idx).animal               = animals(ii);
                datStruct(idx).layer                = layerNames{hh};
                %datStruct(idx).rel2Bf               = thisInd;
                datStruct(idx).freq                 = F;
                datStruct(idx).trial                = ll;
                %datStruct(idx).channel              = thisLayerChans(kk) - median(thisLayerChans);
                
                idx = idx + 1;
                clear WT hold_cwts
                close all
            
        end
        
    end
end
end


