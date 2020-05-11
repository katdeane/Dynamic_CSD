function [dat] = runCwtCsd_CryptST(groupFile,params,homedir,meas)

% Input:    group name to match *Data.mat in \Dynamic_CSD\DATA, parameters
%           set for CWT analysis, home director
% Output:   Runs CWT analysis using the Wavelet toolbox. figures of
%           animal-wise scalograms -> Dynamic_CSD\figs\Spectral_MagPlots
%           and group table output to be formatted to full scalograms.mat

% Note:     Output of single trials!

warning('OFF');
dbstop if error

% Change directory to your working folder
if ~exist('homedir','var')
    if exist('D:\MyCode\Dynamic_CSD','dir') == 7
        cd('D:\MyCode\Dynamic_CSD');
    elseif exist('C:\Users\kedea\Documents\Work Stuff\Dynamic_CSD','dir') == 7
        cd('C:\Users\kedea\Documents\Work Stuff\Dynamic_CSD')
    elseif exist('D:\Dynamic_CSD','dir') == 7
        cd('D:\Dynamic_CSD')
    end
    
    homedir = pwd;
    addpath(genpath(homedir));
end

cd(homedir)

% Generates variables called animals and Layer
run([groupFile '.m'])

% Init datastruct to pass out
% structure in common
dat               = struct();
dat.scalogram     = []; % Continuous wavelet transform
dat.animal        = []; % animal name in group
dat.group         = []; % KIC, KIT, KIV
dat.layer         = []; % as determined by group page
dat.measurement   = []; % measurement
dat.freq          = []; % oscillatory frequencies
if ~contains(meas,'SP')
    dat.stimulus      = []; % 0 = BF or 2,5,10,20,40
end

count   = 1;
clickfrq    = [2,5,10,20,40];

for iAn = 1:length(animals)  %#ok<*USENS>
    
    clear curData
    load([animals{iAn} '_Data.mat'],'Data') %#ok<*IDISVAR>
    curData = Data((strcmp({Data.Condition},meas)));
    clear Data
     
    % Tonotopy measurements
    if (contains(curData.Condition,'Pre_') ...
            || contains(curData.Condition,'tono')) ...
        
        BF_Frqz = find(curData.Frqz == curData.BF_IV);
        
        if isempty(BF_Frqz)
            BF_Frqz = 8000;
        end
        
        % loop through stimuli
        for istim = 1:length(curData.Frqz)
            
            CSD = curData.SglTrl_CSD{1,istim};
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
                    
                    hold_cwts = zeros(54,size(CSD,2));
                    
                    %repeat over channels
                    for iChan = 1:length(curChan)
                        
                        csdChan = squeeze(CSD(curChan(iChan),:,:));
                        ROI = csdChan(:,itrial);
                        
                        % Limit the cwt frequency limits
                        params.frequencyLimits(1) = max(params.frequencyLimits(1),...
                            cwtfreqbounds(numel(ROI),params.sampleRate,...
                            'TimeBandWidth',params.timeBandWidth));
                        
                        [WT,F] = cwt(ROI,params.sampleRate, ...
                            'VoicesPerOctave',params.voicesPerOctave, ...
                            'TimeBandWidth',params.timeBandWidth, ...
                            'FrequencyLimits',params.frequencyLimits);
                        % to get average of channels: add each and divide after loop
                        hold_cwts = hold_cwts + WT;
                    end %channel
                    
                    WT = hold_cwts/length(curChan);
                    
                    dat(count).scalogram            = WT;
                    dat(count).group                = groupFile;
                    dat(count).animal               = animals{iAn};
                    dat(count).layer                = params.layers{iLay};
                    dat(count).stimulus             = stimulus;
                    dat(count).freq                 = F;
                    dat(count).measurement          = curData.Condition;
                    
                    count = count + 1;
                    clear WT hold_cwts
                end %trial
            end %layer
        end %stimuli
    % Click or Amplitude Modulation Measurements
    elseif contains(curData.Condition,'CL_') || contains(curData.Condition,'AM_')
        for istim = 1:length(curData.Frqz)
            
            CSD = curData.SglTrl_CSD{1,istim};
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
                
                %repeat over trials
                for itrial = 1:size(CSD,3)
                    
                    hold_cwts = zeros(54,size(CSD,2));
                    
                    %repeat over channels
                    for iChan = 1:length(curChan)
                        
                        csdChan = squeeze(CSD(curChan(iChan),:,:));
                        ROI = csdChan(:,itrial);
                        
                        % Limit the cwt frequency limits
                        params.frequencyLimits(1) = max(params.frequencyLimits(1),...
                            cwtfreqbounds(numel(ROI),params.sampleRate,...
                            'TimeBandWidth',params.timeBandWidth));
                        
                        [WT,F] = cwt(ROI,params.sampleRate, ...
                            'VoicesPerOctave',params.voicesPerOctave, ...
                            'TimeBandWidth',params.timeBandWidth, ...
                            'FrequencyLimits',params.frequencyLimits);
                        % to get average of channels: add each and divide after loop
                        hold_cwts = hold_cwts + WT;
                    end %channel
                    
                    WT = hold_cwts/length(curChan);
                    dat(count).scalogram            = WT;
                    dat(count).group                = groupFile;
                    dat(count).animal               = animals{iAn};
                    dat(count).layer                = params.layers{iLay};
                    dat(count).stimulus             = stimulus;
                    dat(count).freq                 = F;
                    dat(count).measurement          = curData.Condition;
                    
                    count = count + 1;
                    clear WT hold_cwts
                end %trial
            end %layer
        end %stimulus
    % Spontaneous Measurement
    elseif contains(curData.Condition,'sp') 
        
        % all stim are the same, concatonate along trials (dim 3)
        CSD = cat(3,curData.SglTrl_CSD{1,:});
        
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
            
            %repeat over trials
            for itrial = 1:size(CSD,3)
                
                hold_cwts = zeros(54,size(CSD,2));
                
                %repeat over channels
                for iChan = 1:length(curChan)
                    
                    csdChan = squeeze(CSD(curChan(iChan),:,:));
                    ROI = csdChan(:,itrial);
                    
                    % Limit the cwt frequency limits
                    params.frequencyLimits(1) = max(params.frequencyLimits(1),...
                        cwtfreqbounds(numel(ROI),params.sampleRate,...
                        'TimeBandWidth',params.timeBandWidth));
                    
                    [WT,F] = cwt(ROI,params.sampleRate, ...
                        'VoicesPerOctave',params.voicesPerOctave, ...
                        'TimeBandWidth',params.timeBandWidth, ...
                        'FrequencyLimits',params.frequencyLimits);
                    % to get average of channels: add each and divide after loop
                    hold_cwts = hold_cwts + WT;
                end %channel
                
                WT = hold_cwts/length(curChan);
                dat(count).scalogram            = WT;
                dat(count).group                = groupFile;
                dat(count).animal               = animals{iAn};
                dat(count).layer                = params.layers{iLay};
                dat(count).freq                 = F;
                dat(count).measurement          = curData.Condition;
                
                count = count + 1;
                clear WT hold_cwts
            end %trial
        end %layer
    end %tono, click, or amplitude modulated
end
