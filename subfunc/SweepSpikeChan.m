function [trial_spikes, trial_times] = SweepSpikeChan(inSWEEP,channels,Fs)

%% Sweep Spike Channel sorting 
% author: Katrina Deane 8/30/12

% input:   a single trial's SWEEP.Spikes
%           the current animal's channel list
%           frequency sampling rate
%  ouput:   that trial's spikes sorted into given channels for the animal
%           each channel is represented by a single vector of ms time
%           points; 0 = no spikes; 1 = 1 spike; etc. (shape of CSD)
% Used in code "dynamic_spike" to pre-sort spikes. 

% NOTE: the spike times are *rounded* to the nearest ms and are therefore 
% no longer raw data beyond this point!


%%
    % set up output matrix (channels x measurement length)
    len_meas = 1400; %(maybe make function input)
    len_chan = length(channels);
    trial_spikes = false(len_chan,len_meas);
    trial_times = nan(len_chan,55); %not the perfect solution...
   
for i = 1:len_chan
    
    % generate 2 boolean vectors the length of the signal:
    boo_v = false(1,len_meas);
    boo_v2 = false(1,len_meas);
    
    % NOTE: currently the stimulus onset is set to 0, meaning spikes before 
    % are negative numbers -> Pull out spike time points and add .2 
    trialSWEEP = inSWEEP(inSWEEP(:,2) == channels(i))+0.2;
    % time points to the closes ms of each spike in this channel:
    trialSWEEP = round(trialSWEEP*Fs)';
    
    % if a spike rounded to 0, make it 1 ms
    if ~isempty(trialSWEEP) && sum(trialSWEEP(:) == 0) >= 1
        trialSWEEP(trialSWEEP(:) == 0) = 1;
    end
    
    % for 2 spikes in 1 ms:
    for ii = 1:length(trialSWEEP)-1
        if trialSWEEP(ii) == trialSWEEP(ii+1)
            boo_v2(trialSWEEP(ii)) = true;
        end
    end
    
    % make true each spike time point in the boolean:
    boo_v(trialSWEEP) = true; 
    
    chan_spikes = boo_v + boo_v2;
    
    % place vector into output matrix
    trial_spikes(i,:) = chan_spikes;
    if ~isempty(trialSWEEP)
        trial_times(i,1:length(trialSWEEP)) = trialSWEEP;
    end
    
end

