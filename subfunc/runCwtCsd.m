function [datStruct] = runCwtCsd(groupFile,params,homedir,runcase)

% Input:    group name to match *Data.mat in \Dynamic_CSD\DATA, parameters
%           set for CWT analysis, home director
% Output:   Runs CWT analysis using the Wavelet toolbox. figures of
%           animal-wise scalograms -> Dynamic_CSD\figs\Spectral_MagPlots
%           and group table output to be formatted to full scalograms.mat 

if ~contains(runcase,'average') && ~contains(runcase,'single')
    error('The run case needs to be either "average" or "single"')
end

cd(homedir);

% Generates variables called animals and Layer
run([groupFile '.m'])
load([groupFile '_Data.mat'],'Data')
names = fieldnames(Data);

% Init datastruct to pass out
datStruct = struct();
datStruct.scalogram     = [];
datStruct.animal        = [];
datStruct.freq          = [];
datStruct.condition     = [];
datStruct.layer         = [];
if contains(runcase,'average')
    datStruct.rel2Bf    = [];
    cd figs; mkdir('Spectral_MagPlot'); cd('Spectral_MagPlot');
elseif contains(runcase,'single')
    datStruct.trial    = [];
    cd figs; mkdir('Spectral_AngPlot'); cd('Spectral_AngPlot');
end
    

count = 1;
for iAn = 1:length(animals) %#ok<*USENS>
    
    BF_Frqz = find(Data.(names{iAn}).Frqz == Data.(names{iAn}).GS_BF);
    CSD_0 = Data.(names{iAn}).singletrialCSD{1,BF_Frqz};
    if BF_Frqz > 2
        CSD_min2 = Data.(names{iAn}).singletrialCSD{1,BF_Frqz-2};
    else
        CSD_min2 = 'skip me';
    end
        
    for iLay = 1:length(params.layers)
        
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
        
        if contains(runcase,'average')
            for ifrqz = 1:length(params.rel2BFlist) % just take BF and offBF
                
                rel2BF = params.rel2BFlist(ifrqz);
                
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
        elseif contains(runcase,'single')
            % Repeat over trials
            for itrial = 1:size(CSD_0,3)
                
                hold_cwts = zeros(54,size(CSD_0,2));
                
                if BF_Frqz + params.rel2BFlist <= 0
                    continue
                end
                
                % Repeat over channels
                for iChan = 1:length(curChan)
                    
                    if params.rel2BFlist == 0
                        csdChan = squeeze(CSD_0(curChan(iChan),:,:));
                    elseif params.rel2BFlist == -2
                        csdChan = squeeze(CSD_min2(curChan(iChan),:,:));
                    end
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
                end
                
                WT = hold_cwts/3;
                %Get the cone of influence in a very round-about way....
                if itrial == 1
                    tit = sprintf('ID_%s_%dfromBF_Layer_%s_condition_%s_trialone',animals{iAn},...
                        params.rel2BFlist, params.layers{iLay},groupFile);
                    [X,Y]=meshgrid(F,1:size(WT,2));
                    cwtfig = figure('Name',tit);
                    surf(Y,X,angle(WT)','EdgeColor','None'); view(2);
                    saveas(cwtfig,tit)
                    close(cwtfig)
                end
                    
                datStruct(count).scalogram            = WT;
                datStruct(count).condition            = {groupFile};
                datStruct(count).animal               = animals(iAn);
                datStruct(count).layer                = params.layers{iLay};
                datStruct(count).freq                 = F;
                datStruct(count).trial                = itrial;
                                
                count = count + 1;
                clear WT hold_cwts
            end
        end
    end
end
