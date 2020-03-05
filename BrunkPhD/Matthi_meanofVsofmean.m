%%% Input
warning('OFF');
dbstop if error

G = {'Input_Control_7post_n7.m' 'Input_HighP_7post.m' 'Input_YFP_7post.m'};

Kernel = 300;
chan_dist= 50;


home = cd;
addpath(genpath(cd))
cd groups


for i1 = 1:length(G)
   run(G{i1}); cd (home);cd raw;
   keyboard
   X ;% animals measurements Para(AVREC,AbsRes,RelRes),(mean of/of mean) 
   
    
   for i2 = 1:length(animals)
%        keyboard  
       Index = 1;      
       Measure =fieldnames(Cond);
       for i3 = 1:size(fieldnames(Cond),1)
           
           for i4 = 1:length(Cond.(Measure{i3}){i2})
               
                try
               load([animals{i2} '_' Cond.(Measure{i3}){i2}{i4}]);
               
               ChanOrder =str2num(channels{i1});
               BL = Header.t_pre*P.Fs_AD(1);
               tone = Header.t_sig(1)*P.Fs_AD(1);
               frqz = Header.stimlist(:,1); 
               frqz(find(frqz == inf))=[];
               frqz(find(frqz == 0))=[];
               
               kernel   = Kernel./chan_dist; % kernel size in terms of number of channels
               hammsiz  = kernel+(rem(kernel,2)-1)*-1; % linear extrapolation and running avg
               paddsiz  = floor(hammsiz/2)+1;
               
               stim_list = [SWEEP.Header];
               All = nan(size(SWEEP(1).AD(:,ChanOrder,:),2),size(SWEEP(1).AD(:,ChanOrder,:),1),length(SWEEP));
               for i_trial = 1:size(stim_list,2)

                % create matrix that mimicks avgFP in its dimensions (chan, time,
                % stim) - chan, time, trial
                dummy_AD = SWEEP(i_trial).AD(:,ChanOrder,:);

                x1  = dummy_AD';
                if BL > 0 
                    x1 = base_corr_00(x1,BL,2); % correct the baseline based on the first n_SamplesPrevent number of samples in the 2nd dimension
                end
                % AvgRecCSD and RelResCSD %based on unfiltered data
                x2 = padd_linex_01(x1,paddsiz,1);%MB 27032018 addition of X2 and X3 before X4 to have 32 Channels for Harding
                x3 = filt_Hamm_00(x2,hammsiz,1);
                x4 = (get_csd_00(x3,1,chan_dist,1))*10^3; % correct data to be in V/mm^2 and not mV/mm^2 % X1 changed in X3
 
                All(:,:,i_trial)= x4;
                x5 = get_Harding_03(x4,1); %x4 = 32?
            %     x6 = nanmean(x4(chan,:));

                % CSD based on spatially filtered data
                x2 = padd_linex_01(x1,paddsiz,1);
                x3 = filt_Hamm_00(x2,hammsiz,1);
                x4 = (get_csd_00(x3,1,chan_dist,1))*10^3; % correct data to be in V/mm^2 and not mV/mm^2
            %     x6 = nanmean(x4(chan,:));

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
                n_stim  = size(un_stim,2);
                % computer ran out of memory so now do this for every stimulus
                % separately
                AVRECofMean = [];
                MeanofAVREC = [];
                AbsResofMean = [];
                MeanofAbsRes = [];
                RelResofMean = [];
                MeanofRelRes = [];
%                                keyboard
                  for i_stim = un_stim
                    i_TrialCurStim = find(stim_list == i_stim);
%                     keyboard
                    dummy_SingleTrialAvgRecCSD = reshape([SWEEP2(i_TrialCurStim).AvgRecCSD],size(SWEEP2(i_TrialCurStim(1)).AvgRecCSD,1),size(SWEEP2(i_TrialCurStim(1)).AvgRecCSD,2),[]);
                    AVRECofMean(i_stim,:) = mean(abs(mean(All(:,:,i_TrialCurStim ),3)));                    
                    MeanofAVREC(i_stim,:) = mean(mean(abs(All(:,:,(i_TrialCurStim))),1),3);
                    clear dummy_SingleTrialAvgRecCSD
                    
                    dummy_SingleTrialAbsResCSD = reshape([SWEEP2(i_TrialCurStim).AbsResCSD],size(SWEEP2(i_TrialCurStim(1)).AbsResCSD,1),size(SWEEP2(i_TrialCurStim(1)).AbsResCSD,2),[]);
                    AbsResofMean(i_stim,:) = sum(mean(All(:,:,(i_TrialCurStim)),3));
                    MeanofAbsRes(i_stim,:) = mean(mean(sum(All(:,:,(i_TrialCurStim))),1),3);
                    clear dummy_SingleTrialAbsResCSD

                    dummy_SingleTrialRelResCSD = reshape([SWEEP2(i_TrialCurStim).RelResCSD],size(SWEEP2(i_TrialCurStim(1)).RelResCSD,1),size(SWEEP2(i_TrialCurStim(1)).RelResCSD,2),[]);
                    RelResofMean(i_stim,:) = sum(mean(All(:,:,(i_TrialCurStim)),3))./mean(abs(mean(All(:,:,i_TrialCurStim ),3)));
                    MeanofRelRes(i_stim,:) = mean(sum(All(:,:,(i_TrialCurStim)))./mean(abs(All(:,:,(i_TrialCurStim)))),3);
                    clear dummy_SingleTrialRelResCSD
% 

                  end
 
               catch
                   
               end
%                 keyboard
               BF = max(sqrt(mean(AVRECofMean(:,200:250).^2,2)));
               BF_Pos = find (sqrt(mean(AVRECofMean(:,200:250).^2,2)) == BF);
               matrix(i2,1:400,Index,1,1) = AVRECofMean(BF_Pos,201:600);
               matrix(i2,1:400,Index,2,1) = AbsResofMean(BF_Pos,201:600);
               matrix(i2,1:400,Index,3,1) = RelResofMean(BF_Pos,201:600);
               
               matrix(i2,1:400,Index,1,2) = MeanofAVREC(BF_Pos,201:600);
               matrix(i2,1:400,Index,2,2) = MeanofAbsRes(BF_Pos,201:600);
               matrix(i2,1:400,Index,3,2) = MeanofRelRes(BF_Pos,201:600);
               Index = Index+1;
           end
       end
   end 
end