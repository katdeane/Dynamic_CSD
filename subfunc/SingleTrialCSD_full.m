function [SingleRec_preAVREC,Rec_preAVREC,AvgFP, SingleTrialFP, AvgCSD,...
    SingleTrialCSD, AvgRecCSD, SingleTrialAvgRecCSD,...
    SingleTrialRelResCSD, AvgRelResCSD,AvgAbsResCSD,...
    SingleTrialAbsResCSD, SingleLayerRelRes, AvgLayerRelRes] =...
    SingleTrialCSD_full(SWEEP,chan_order,chan,n_SamplesPreevent,kernel,chan_dist)
%% check inputs
if ~exist('kernel','var') % filter parameters based on kernel width
     kernel   = 300; %450 also preferred by some
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

    dummy_AD = SWEEP(i_trial).AD(:,chan_order,:);
    x1  = dummy_AD';

    
    if n_SamplesPreevent > 0 
        x1 = base_corr(x1,n_SamplesPreevent,2); % correct the baseline based on the first n_SamplesPrevent number of samples in the 2nd dimension
    end
    % CSD
    x2 = padd_linex(x1,paddsiz,1);%MB 27032018 addition of X2 and X3 before X4 to have 32 Channels for Harding
    x3 = filt_Hamm(x2,hammsiz,1);
    x4 = (get_csd(x3,1,chan_dist,1))*10^3; % correct data to be in V/mm^2 and not mV/mm^2 % X1 changed in X3
    % AvgRecCSD and RelResCSD 
    x5 = get_Harding(x4,1); 
    x6 = abs(x4);
    
    % handover to new SWEEP structure     
    SWEEP2(i_trial).AD  = x1';
    SWEEP2(i_trial).CSD = x4;
    SWEEP2(i_trial).AvgRecCSD = x5.var4';
    SWEEP2(i_trial).RelResCSD = x5.var3';
    SWEEP2(i_trial).AbsResCSD = x5.var1';% X5.var7 = layer specific RelRes
    SWEEP2(i_trial).LayerRelRes = x5.var6;% MB 27032018
    SWEEP2(i_trial).Layer_AVREC = x6;
    clear x1 x2 x3 x4 x5 x6 dummy_AD
end

un_stim = unique(stim_list);
un_stim = un_stim(find(~isnan(un_stim)));

% LFP and CSD
AvgFP =                 cell(1,length(un_stim));
SingleTrialFP =         cell(1,length(un_stim));
AvgCSD =                cell(1,length(un_stim));
SingleTrialCSD =        cell(1,length(un_stim));
% calculations on top of CSD
AvgRecCSD =             cell(1,length(un_stim));
SingleTrialAvgRecCSD =  cell(1,length(un_stim));
AvgRelResCSD =          cell(1,length(un_stim));
SingleTrialRelResCSD =  cell(1,length(un_stim));
AvgAbsResCSD =          cell(1,length(un_stim));
SingleTrialAbsResCSD =  cell(1,length(un_stim));
% channelwise
AvgLayerRelRes =        cell(1,length(un_stim));
SingleLayerRelRes =     cell(1,length(un_stim));
SingleRec_preAVREC =    cell(1,length(un_stim)); % all channels rectified per trial for use in layerwise avrec analysis
Rec_preAVREC =          cell(1,length(un_stim)); % all channels rectified per trial for use in layerwise avrec analysis

for i_stim = un_stim
    i_TrialCurStim = find(stim_list == i_stim);
    dummy_AvgFP = reshape([SWEEP2(i_TrialCurStim).AD],size(SWEEP2(i_TrialCurStim(1)).AD,1),size(SWEEP2(i_TrialCurStim(1)).AD,2),[]);
    SingleTrialFP{i_stim} = dummy_AvgFP;
    AvgFP{i_stim} = mean(dummy_AvgFP,3);
    clear dummy_AvgFP
    
    dummy_LayerRelRes = reshape([SWEEP2(i_TrialCurStim).LayerRelRes],size(SWEEP2(i_TrialCurStim(1)).LayerRelRes,1),size(SWEEP2(i_TrialCurStim(1)).LayerRelRes,2),[]);
    SingleLayerRelRes{i_stim} = dummy_LayerRelRes;
    AvgLayerRelRes{i_stim} = mean(dummy_LayerRelRes,3);
    clear dummy_LayerRelRes
    
    dummy_SingleTrialCSD = reshape([SWEEP2(i_TrialCurStim).CSD],size(SWEEP2(i_TrialCurStim(1)).CSD,1),size(SWEEP2(i_TrialCurStim(1)).CSD,2),[]);
    SingleTrialCSD{i_stim} = dummy_SingleTrialCSD;
    AvgCSD{i_stim} = mean(dummy_SingleTrialCSD,3);
    clear dummy_SingleTrialCSD
    
    dummy_SingleTrialAvgRecCSD = reshape([SWEEP2(i_TrialCurStim).AvgRecCSD],size(SWEEP2(i_TrialCurStim(1)).AvgRecCSD,1),size(SWEEP2(i_TrialCurStim(1)).AvgRecCSD,2),[]);
    SingleTrialAvgRecCSD{i_stim} = dummy_SingleTrialAvgRecCSD;
    AvgRecCSD{i_stim} = mean(dummy_SingleTrialAvgRecCSD,3);
    clear dummy_SingleTrialAvgRecCSD
    
    dummy_SingleTrialRelResCSD = reshape([SWEEP2(i_TrialCurStim).RelResCSD],size(SWEEP2(i_TrialCurStim(1)).RelResCSD,1),size(SWEEP2(i_TrialCurStim(1)).RelResCSD,2),[]);
    SingleTrialRelResCSD{i_stim} = dummy_SingleTrialRelResCSD;
    AvgRelResCSD{i_stim} = mean(dummy_SingleTrialRelResCSD,3);
    clear dummy_SingleTrialRelResCSD
    
    dummy_SingleTrialAbsResCSD = reshape([SWEEP2(i_TrialCurStim).AbsResCSD],size(SWEEP2(i_TrialCurStim(1)).AbsResCSD,1),size(SWEEP2(i_TrialCurStim(1)).AbsResCSD,2),[]);
    SingleTrialAbsResCSD{i_stim} = dummy_SingleTrialAbsResCSD;
    AvgAbsResCSD{i_stim} = mean(dummy_SingleTrialAbsResCSD,3);
    clear dummy_SingleTrialAbsResCSD
    
    dummy_SingleTrialLayerAVREC = reshape([SWEEP2(i_TrialCurStim).Layer_AVREC],size(SWEEP2(i_TrialCurStim(1)).Layer_AVREC,1),size(SWEEP2(i_TrialCurStim(1)).Layer_AVREC,2),[]);
    SingleRec_preAVREC{i_stim} = dummy_SingleTrialLayerAVREC ;
    Rec_preAVREC{i_stim} = mean(dummy_SingleTrialLayerAVREC,3 );
    clear dummy_SingleTrialAbsResCSD
end
