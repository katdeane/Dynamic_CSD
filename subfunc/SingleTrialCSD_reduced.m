function [AvgFP, AvgCSD, AvgRecCSD, SingleTrialAvgRecCSD, AvgRelResCSD, SingleTrialAvgRelResCSD, AvgAbsResCSD,...
    AvgLayerRelRes] =SingleTrialCSD_reduced(SWEEP,chan_order,n_SamplesPreevent,kernel,chan_dist)
%currenty all single trial output is silenced 

%% check inputs
if ~exist('kernel','var') % filter parameters based on kernel width
     kernel   = 300;
end
if ~exist('chan_dist','var')
    chan_dist= 50; % NeuronexusProbes distance
end
if ~exist('chan_order','var')
    chan_order = 1:size(SWEEP(1,1).AD,2);
    disp('Warning: channel order not given; taking all channels as is!')
end
if ~exist('n_SamplesPreevent','var')
    n_SamplesPreevent = 0; 
end

%% Calculate parameters for filtering
kernelxchannel   = kernel./chan_dist; % kernel size in terms of number of channels
hammsiz  = kernelxchannel+(rem(kernelxchannel,2)-1)*-1; % linear extrapolation and running avg
paddsiz  = floor(hammsiz/2)+1;

%% Generate output
stim_list = [SWEEP.Header];
for i_trial = 1:size(stim_list,2)

    % create matrix that mimicks avgFP in its dimensions (chan, time,
    % stim) - chan, time, trial
    dummy_AD = SWEEP(i_trial).AD(:,chan_order,:);
    
    x1  = dummy_AD';
    if n_SamplesPreevent > 0 
        x1 = base_corr(x1,n_SamplesPreevent,2); % correct the baseline based on the first n_SamplesPrevent number of samples in the 2nd dimension
    end
    % AvgRecCSD and RelResCSD %based on unfiltered data
    x2 = padd_linex(x1,paddsiz,1);%MB 27032018 addition of X2 and X3 before X4 to have 32 Channels for Harding
    x3 = filt_Hamm(x2,hammsiz,1);
    x4 = (get_csd(x3,1,chan_dist,1))*10^3; % correct data to be in V/mm^2 and not mV/mm^2 % X1 changed in X3
    x5 = get_Harding(x4,1); 
    
    % CSD based on spatially filtered data
    x2 = padd_linex(x1,paddsiz,1);
    x3 = filt_Hamm(x2,hammsiz,1);
    x4 = (get_csd(x3,1,chan_dist,1))*10^3; % correct data to be in V/mm^2 and not mV/mm^2
 
    % handover to new SWEEP structure     
    SWEEP2(i_trial).AD  = x1';
    SWEEP2(i_trial).CSD = x4;
    SWEEP2(i_trial).AvgRecCSD = x5.var4';
    SWEEP2(i_trial).RelResCSD = x5.var3';
    SWEEP2(i_trial).AbsResCSD = x5.var1';% X5.var7 = layer specific RelRes
    SWEEP2(i_trial).LayerRelRes = x5.var6;% MB 27032018
    
    clear x1 x2 x3 x4 x5 dummy_AD
end

un_stim = unique(stim_list);
un_stim = un_stim(find(~isnan(un_stim)));

AvgFP = {};
SingleTrialAvgRecCSD = {};
    
  for i_stim = un_stim
    i_TrialCurStim = find(stim_list == i_stim);
    AvgFP{i_stim} = mean(reshape([SWEEP2(i_TrialCurStim).AD],size(SWEEP2(i_TrialCurStim(1)).AD,1),size(SWEEP2(i_TrialCurStim(1)).AD,2),[]),3);
    AvgLayerRelRes{i_stim} = mean(reshape([SWEEP2(i_TrialCurStim).LayerRelRes],size(SWEEP2(i_TrialCurStim(1)).LayerRelRes,1),size(SWEEP2(i_TrialCurStim(1)).LayerRelRes,2),[]),3);
    AvgCSD{i_stim} = mean(reshape([SWEEP2(i_TrialCurStim).CSD],size(SWEEP2(i_TrialCurStim(1)).CSD,1),size(SWEEP2(i_TrialCurStim(1)).CSD,2),[]),3);   
    AvgRecCSD{i_stim} = mean(reshape([SWEEP2(i_TrialCurStim).AvgRecCSD],size(SWEEP2(i_TrialCurStim(1)).AvgRecCSD,1),size(SWEEP2(i_TrialCurStim(1)).AvgRecCSD,2),[]),3);
    SingleTrialAvgRecCSD{i_stim} = reshape([SWEEP2(i_TrialCurStim).AvgRecCSD],size(SWEEP2(i_TrialCurStim(1)).AvgRecCSD,1),size(SWEEP2(i_TrialCurStim(1)).AvgRecCSD,2),[]);
    AvgRelResCSD{i_stim} = mean(reshape([SWEEP2(i_TrialCurStim).RelResCSD],size(SWEEP2(i_TrialCurStim(1)).RelResCSD,1),size(SWEEP2(i_TrialCurStim(1)).RelResCSD,2),[]),3);
    SingleTrialAvgRelResCSD{i_stim} = reshape([SWEEP2(i_TrialCurStim).RelResCSD],size(SWEEP2(i_TrialCurStim(1)).RelResCSD,1),size(SWEEP2(i_TrialCurStim(1)).RelResCSD,2),[]);
    AvgAbsResCSD{i_stim} = mean(reshape([SWEEP2(i_TrialCurStim).AbsResCSD],size(SWEEP2(i_TrialCurStim(1)).AbsResCSD,1),size(SWEEP2(i_TrialCurStim(1)).AbsResCSD,2),[]),3);
  end
