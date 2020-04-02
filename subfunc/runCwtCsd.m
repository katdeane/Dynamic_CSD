function [datStruct] = runCwtCsd(groupFile,params,homedir)

% Input:    group name to match *Data.mat in \Dynamic_CSD\DATA, parameters
%           set for CWT analysis, home director
% Output:   Runs CWT analysis using the Wavelet toolbox. figures of
%           animal-wise scalograms -> Dynamic_CSD\figs\Spectral_MagPlots
%           and group table output to be formatted to full scalograms.mat 

cd(homedir);

% Generates variables called animals and Layer
run([groupFile '.m'])
load([groupFile '_Data.mat'],'Data')
names = fieldnames(Data);

cd figs; mkdir('Spectral_MagPlot'); cd('Spectral_MagPlot');
% Init datastruct to pass out
datStruct = struct();
datStruct.scalogram     = [];
datStruct.animal        = [];
datStruct.freq          = [];
datStruct.condition     = [];
datStruct.rel2Bf        = [];
datStruct.layer         = [];
% layerNames = fieldnames(Layer);
%layerNames = {'IVE'}; % REMEMBER TO REMOVE THIS AFTER TESTING!
count = 1;
for iAn = 1:length(animals) %#ok<*USENS>
    
    BF_Frqz = find(Data.(names{iAn}).Frqz == Data.(names{iAn}).GS_BF);
    rel2BFlist = [0 -2];
    CSD_0 = Data.(names{iAn}).singletrialCSD{1,BF_Frqz};
    if BF_Frqz > 2
        CSD_min2 = Data.(names{iAn}).singletrialCSD{1,BF_Frqz-2};
    else
        CSD_min2 = 'skip me';
    end
        
    for iLay = 1:length(params.layers)
        
        % For constructing position of recording relative to layer, whole
        % numbers are best, but if you have an even number of channels, you get
        % halves. If even number of channels present, for now, using ceil() to
        % assign the next number and base all other numbers on that. This
        % means that odd numbers (e.g. 7) will use exact middle (e.g. 4), while
        % even (e.g. 6) returns middle/slightly lower (e.g. 3). 8 would return
        % 4. etc etc.
        % May switch to an averaging strategy OR simply take best values in the
        % future.
        
        % Format here forces the use of eval to get channel locs
        curChan = str2num(Layer.(params.layers{iLay}){iAn}); %#ok<*ST2NM>
        centerChan = curChan(ceil(length(curChan)/2));
        theseChans = curChan-centerChan;
        % Select only center 3 channels
        curChan = curChan(theseChans >=-1 & theseChans <=1);
        
        for ifrqz = 1:length(rel2BFlist) % just take BF and offBF
        
            rel2BF = rel2BFlist(ifrqz);
            
            % in case the off-BF condition does not exist
            if BF_Frqz + rel2BF <= 0
                continue
            end
            
            %average over trials
            if rel2BF == 0
                meanLayCSD = mean(CSD_0(curChan,:,:),3);
            elseif rel2BF == -2
                meanLayCSD = mean(CSD_min2(curChan,:,:),3);
            end
            %average over channels
            ROI = mean(meanLayCSD,1);
            
            % Limit the cwt frequency limits
            params.frequencyLimits(1) = max(params.frequencyLimits(1),...
                cwtfreqbounds(numel(ROI),params.sampleRate,...
                'TimeBandWidth',params.timeBandWidth));
            % Get the cone of influence in a very round-about way....
            tit = sprintf('ID_%s_%dfromBF_Layer_%s_condition_%s',animals{iAn},rel2BF, ...
                params.layers{iLay},groupFile);
            cwtfig = figure('Name',tit); cwt(ROI,params.sampleRate, ...
                'VoicesPerOctave',params.voicesPerOctave, ...
                'TimeBandWidth',params.timeBandWidth, ...
                'FrequencyLimits',params.frequencyLimits);
            saveas(cwtfig,tit)
            close all
            
            [WT,F] = cwt(ROI,params.sampleRate, ...
                'VoicesPerOctave',params.voicesPerOctave, ...
                'TimeBandWidth',params.timeBandWidth, ...
                'FrequencyLimits',params.frequencyLimits);
            datStruct(count).scalogram            = WT;
            datStruct(count).condition            = {groupFile};
            datStruct(count).animal               = animals(iAn);
            datStruct(count).layer                = params.layers{iLay};
            datStruct(count).rel2Bf               = rel2BF;
            datStruct(count).freq                 = F;
            
            count = count + 1;
            clear WT
        end
    end
end
