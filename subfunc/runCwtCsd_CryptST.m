function [datTune,datCL,datAM,datSP] = runCwtCsd_CryptST(groupFile,params,homedir)

% Input:    group name to match *Data.mat in \Dynamic_CSD\DATA, parameters
%           set for CWT analysis, home director
% Output:   Runs CWT analysis using the Wavelet toolbox. figures of
%           animal-wise scalograms -> Dynamic_CSD\figs\Spectral_MagPlots
%           and group table output to be formatted to full scalograms.mat 

% Note:     Output of single trials!

cd(homedir);

% Generates variables called animals and Layer
run([groupFile '.m'])

% Init datastruct to pass out
% structure for tonotopic measurements
datTune               = struct();
datTune.scalogram     = []; % Continuous wavelet transform
datTune.animal        = []; % animal name in group
datTune.freq          = []; % oscillatory frequencies
datTune.group         = []; % KIC, KIT, KIV
datTune.layer         = []; % as determined by group page
datTune.stimulus      = []; % 0 = BF 
datTune.measurement   = []; % measurement 
% structure for click measurements
datCL               = struct();
datCL.scalogram     = []; 
datCL.animal        = []; 
datCL.group         = []; 
datCL.layer         = []; 
datCL.stimulus      = []; % 2,5,10,20,40 
datCL.measurement   = []; 
% structure for amplitude modulation measurements
datAM               = struct();
datAM.scalogram     = []; 
datAM.animal        = []; 
datAM.freq          = []; 
datAM.group         = []; 
datAM.layer         = []; 
datAM.stimulus      = []; 
datAM.measurement   = []; 
% structure for spontaneous measurements
datSP               = struct();
datSP.scalogram     = []; 
datSP.animal        = []; 
datSP.freq          = []; 
datSP.group         = []; 
datSP.layer         = []; 
datSP.measurement   = []; 

tunecount   = 1;
amcount     = 1;
clcount     = 1; 
spcount     = 1;
clickfrq    = [2,5,10,20,40];

for iAn = 1:length(animals)  %#ok<*USENS>
    
    clear Data
    load([animals{iAn} '_Data.mat'],'Data') %#ok<*IDISVAR>
    cond = {Data.Condition};
    
   % loop through conditions
    for iCond = 1:length(cond)
        
        if isempty(Data(iCond).Condition)
            continue
        end
        
        if contains(Data(iCond).Condition,'Pre_') ...
                || contains(Data(iCond).Condition,'tono')
            
            BF_Frqz = find(Data(iCond).Frqz == Data(iCond).BF_IV);
            
            if isempty(BF_Frqz)
                BF_Frqz = find(Data(iCond).Frqz == median([Data.BF_IV]));
            end
            
            % loop through stimuli
            for istim = 1:length(Data(iCond).Frqz)
                
                CSD = Data(iCond).SglTrl_CSD{1,istim};
                stimulus = istim - BF_Frqz;
                
                for iLay = 1:length(params.layers)
                    
                    if isempty(str2num(Layer.(params.layers{iLay}){iAn}))
                        continue
                    end
                    
                    % For constructing position of recording relative to layer, use ceil() to
                    % assign the next number and base all other numbers on that. This
                    % means that odd numbers (e.g. 7) will use exact middle (e.g. 4), while
                    % even (e.g. 6) returns middle/slightly lower (e.g. 3).
                    
                    % Format here forces the use of eval to get channel locs
                    curChan = str2num(Layer.(params.layers{iLay}){iAn}); %#ok<*ST2NM>
                    centerChan = curChan(ceil(length(curChan)/2));
                    theseChans = curChan-centerChan;
                    % Select only center 3 channels
                    curChan = curChan(theseChans >=-1 & theseChans <=1);
                    
                    %repeat over trials
                    for itrial = 1:size(CSD,3)
                        meanLayCSD = mean(CSD(curChan,:,:),3);
                        
                        %average over channels
                        ROI = mean(meanLayCSD,1);
                        
                        % Limit the cwt frequency limits
                        params.frequencyLimits(1) = max(params.frequencyLimits(1),...
                            cwtfreqbounds(numel(ROI),params.sampleRate,...
                            'TimeBandWidth',params.timeBandWidth));
                        
                        [WT,F] = cwt(ROI,params.sampleRate, ...
                            'VoicesPerOctave',params.voicesPerOctave, ...
                            'TimeBandWidth',params.timeBandWidth, ...
                            'FrequencyLimits',params.frequencyLimits);
                        datTune(tunecount).scalogram            = WT;
                        datTune(tunecount).group                = groupFile;
                        datTune(tunecount).animal               = animals{iAn};
                        datTune(tunecount).layer                = params.layers{iLay};
                        datTune(tunecount).stimulus             = stimulus;
                        datTune(tunecount).freq                 = F;
                        datTune(tunecount).measurement          = Data(iCond).Condition;
                        
                        tunecount = tunecount + 1;
                        clear WT
                    end
                end %layer
            end %stimuli
        elseif contains(Data(iCond).Condition,'CL_')
            for istim = 1:length(Data(iCond).Frqz)
                
                CSD = Data(iCond).SglTrl_CSD{1,istim};
                stimulus = clickfrq(istim);

                for iLay = 1:length(params.layers)
                    
                    if isempty(str2num(Layer.(params.layers{iLay}){iAn}))
                        continue
                    end
                    
                    % Format here forces the use of eval to get channel locs
                    curChan = str2num(Layer.(params.layers{iLay}){iAn}); %#ok<*ST2NM>
                    centerChan = curChan(ceil(length(curChan)/2));
                    theseChans = curChan-centerChan;
                    % Select only center 3 channels
                    curChan = curChan(theseChans >=-1 & theseChans <=1);
                    
                    %average over trials
                    meanLayCSD = mean(CSD(curChan,:,:),3);
                    
                    %average over channels
                    ROI = mean(meanLayCSD,1);
                    
                    % Limit the cwt frequency limits
                    params.frequencyLimits(1) = max(params.frequencyLimits(1),...
                        cwtfreqbounds(numel(ROI),params.sampleRate,...
                        'TimeBandWidth',params.timeBandWidth));
                                        
                    [WT,F] = cwt(ROI,params.sampleRate, ...
                        'VoicesPerOctave',params.voicesPerOctave, ...
                        'TimeBandWidth',params.timeBandWidth, ...
                        'FrequencyLimits',params.frequencyLimits);
                    datCL(clcount).scalogram            = WT;
                    datCL(clcount).group                = groupFile;
                    datCL(clcount).animal               = animals{iAn};
                    datCL(clcount).layer                = params.layers{iLay};
                    datCL(clcount).stimulus             = stimulus;
                    datCL(clcount).freq                 = F;
                    datCL(clcount).measurement          = Data(iCond).Condition;
                    
                    clcount = clcount + 1;
                    clear WT
                end %layer
            end %stimulus
        elseif contains(Data(iCond).Condition,'AM_')
            for istim = 1:length(Data(iCond).Frqz)
                
                CSD = Data(iCond).SglTrl_CSD{1,istim};
                stimulus = clickfrq(istim);

                for iLay = 1:length(params.layers)
                    
                    if isempty(str2num(Layer.(params.layers{iLay}){iAn}))
                        continue
                    end
                    
                    % Format here forces the use of eval to get channel locs
                    curChan = str2num(Layer.(params.layers{iLay}){iAn}); %#ok<*ST2NM>
                    centerChan = curChan(ceil(length(curChan)/2));
                    theseChans = curChan-centerChan;
                    % Select only center 3 channels
                    curChan = curChan(theseChans >=-1 & theseChans <=1);
                    
                    %average over trials
                    meanLayCSD = mean(CSD(curChan,:,:),3);
                    
                    %average over channels
                    ROI = mean(meanLayCSD,1);
                    
                    % Limit the cwt frequency limits
                    params.frequencyLimits(1) = max(params.frequencyLimits(1),...
                        cwtfreqbounds(numel(ROI),params.sampleRate,...
                        'TimeBandWidth',params.timeBandWidth));
                                        
                    [WT,F] = cwt(ROI,params.sampleRate, ...
                        'VoicesPerOctave',params.voicesPerOctave, ...
                        'TimeBandWidth',params.timeBandWidth, ...
                        'FrequencyLimits',params.frequencyLimits);
                    datAM(amcount).scalogram            = WT;
                    datAM(amcount).group                = groupFile;
                    datAM(amcount).animal               = animals{iAn};
                    datAM(amcount).layer                = params.layers{iLay};
                    datAM(amcount).stimulus             = stimulus;
                    datAM(amcount).freq                 = F;
                    datAM(amcount).measurement          = Data(iCond).Condition;
                    
                    amcount = amcount + 1;
                    clear WT
                end %layer
            end %stimulus
        elseif contains(Data(iCond).Condition,'sp')
                
            % all stim are the same, concatonate along trials (dim 3)
            CSD = cat(3,Data(iCond).SglTrl_CSD{1,:});
            
            for iLay = 1:length(params.layers)
                
                if isempty(str2num(Layer.(params.layers{iLay}){iAn}))
                    continue
                end
                
                % Format here forces the use of eval to get channel locs
                curChan = str2num(Layer.(params.layers{iLay}){iAn}); %#ok<*ST2NM>
                centerChan = curChan(ceil(length(curChan)/2));
                theseChans = curChan-centerChan;
                % Select only center 3 channels
                curChan = curChan(theseChans >=-1 & theseChans <=1);
                
                %average over trials
                meanLayCSD = mean(CSD(curChan,:,:),3);
                
                %average over channels
                ROI = mean(meanLayCSD,1);
                
                % Limit the cwt frequency limits
                params.frequencyLimits(1) = max(params.frequencyLimits(1),...
                    cwtfreqbounds(numel(ROI),params.sampleRate,...
                    'TimeBandWidth',params.timeBandWidth));
                
                [WT,F] = cwt(ROI,params.sampleRate, ...
                    'VoicesPerOctave',params.voicesPerOctave, ...
                    'TimeBandWidth',params.timeBandWidth, ...
                    'FrequencyLimits',params.frequencyLimits);
                datSP(spcount).scalogram            = WT;
                datSP(spcount).group                = groupFile;
                datSP(spcount).animal               = animals{iAn};
                datSP(spcount).layer                = params.layers{iLay};
                datSP(spcount).freq                 = F;
                datSP(spcount).measurement          = Data(iCond).Condition;
                
                spcount = spcount + 1;
                clear WT
            end %layer
        end %tono, click, or amplitude modulated
    end %conditions
end
