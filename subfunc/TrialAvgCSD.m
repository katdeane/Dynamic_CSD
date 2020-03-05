function [AvgFP, TrialAvgCSD, TrialAvgAvgRecCSD, TrialAvgRelResCSD,AvgLayerCSD] = TrialAvgCSD(SWEEP,chan_order,chan,n_SamplesPreevent,avgFP,kernel,chan_dist)

% check inputs
if ~exist('kernel','var') % filter parameters based on kernel width
    kernel   = 300; % micrometer kernel width
end
if ~exist('chan_dist','var')
    chan_dist= 50; % % NeuronexusProbes distance
end
if ~exist('chan_order','var')
    chan_order = [1:size(SWEEP(1,1).AD,2)];
    disp('Warning: channel order not given; taking all channels as is!')
end
if ~exist('n_SamplesPreevent','var')
    n_SamplesPreevent = 0; % % NeuronexusProbes distance
end

kernel   = kernel./chan_dist; % kernel size in terms of number of channels
hammsiz  = kernel+(rem(kernel,2)-1)*-1; % linear extrapolation and running avg
paddsiz  = floor(hammsiz/2)+1;

stim_list = [SWEEP.Header];
un_stim = unique(stim_list);
un_stim = un_stim(find(~isnan(un_stim)));
n_stim  = size(un_stim,2);

AvgFP = {};
SingleTrialCSD = {};
SingleTrialAvgRecCSD = {};
SingleTrialRelResCSD = {};
for i_stim = un_stim
    i_TrialCurStim = find(stim_list == i_stim);
    
    % create matrix that mimicks avgFP in its dimensions (chan, time,
    % stim) - chan, time, trial
    
    dummy_AD = [SWEEP(i_TrialCurStim).AD];
    dummy_AD = reshape([dummy_AD],size(SWEEP(i_TrialCurStim(1)).AD,1),size(SWEEP(i_TrialCurStim(1)).AD,2),[]);
    dummy_AD = dummy_AD(:,chan_order,:);
    %dummy_AD = dummy_AD(:,chan,:);

    x1  = dummy_AD;
    
    % now average the field potential for a stimulus
    x1 = mean(x1,3)';
    if n_SamplesPreevent > 0 
        x1 = base_corr(x1,n_SamplesPreevent,2); % correct the baseline based on the first n_SamplesPrevent number of samples in the 2nd dimension
    end
    
    % RelResCSD based on unfiltered data
    x4 = (get_csd_00(x1,1,chan_dist,1))*10^3; % correct data to be in V/mm^2 and not mV/mm^2
    x5 = get_Harding(x4,1); 
    
    dummy_RelResCSD = x5.var3';
    
    
    % CSD based on spatially filtered data
    
    x2 = padd_linex(x1,paddsiz,1);
    x3 = filt_Hamm(x2,hammsiz,1);
    x4 = get_csd(x3,1,chan_dist,1)*10^3; % correct data to be in V/mm^2 and not mV/mm^2
    
    x5 = get_Harding(x4,1);
    x6 = nanmean(x4(chan,:));

    % handover to new SWEEP structure
    AvgFP{i_stim} = x1';
    TrialAvgCSD{i_stim} = x4;
    AvgLayerCSD{i_stim}=x6';
    TrialAvgAvgRecCSD{i_stim} = x5.var4';
    TrialAvgRelResCSD{i_stim} = dummy_RelResCSD;

    clear x1 x2 x3 x4 x5 dummy_AD
end
clear SWEEP % not needed anymore
