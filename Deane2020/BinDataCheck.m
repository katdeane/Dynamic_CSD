% quick test junk code

% load('AvrecChunkData.mat')
load('RelresChunkData.mat')

%split the groups
KetRMS = T.RMS(strcmp('Ket', T.Group));%T.RMS_Avg
AwaRMS = T.RMS(strcmp('Aw', T.Group));
MusRMS = T.RMS(strcmp('Mus', T.Group));
%% Analysis 1: pre, tone, post (1:20, 20:40, 40:60)

Ketpre = cellfun(@(x) mean(x(1:20)),KetRMS);
Kettone = cellfun(@(x) mean(x(21:40)),KetRMS);
Ketpost = cellfun(@(x) mean(x(41:60)),KetRMS);

%manually tried ttest2 comparison with the following:
Awapre = cellfun(@(x) mean(x(1:20)),AwaRMS); %ns
[Hpre, Ppre, CIpre, Tpre] = ttest2(Ketpre,Awapre);
Awatone = cellfun(@(x) mean(x(21:40)),AwaRMS); %ns
[Htone, Ptone, CItone, Ttone] = ttest2(Kettone,Awatone);
Awapost = cellfun(@(x) mean(x(41:60)),AwaRMS); %ns
[Hpost, Ppost, CIpost, Tpost] = ttest2(Ketpost,Awapost);

Muspre = cellfun(@(x) mean(x(1:20)),MusRMS); %***
Mustone = cellfun(@(x) mean(x(21:40)),MusRMS); %***
Muspost = cellfun(@(x) mean(x(41:60)),MusRMS); %**

%% Analysis 2: pre20, tone10,tone20, post10 (10:20, 20:30, 30:40, 40:50)

Ketpre20 = cellfun(@(x) mean(x(10:20)),KetRMS);
Kettone10 = cellfun(@(x) mean(x(21:30)),KetRMS);
Kettone20 = cellfun(@(x) mean(x(31:40)),KetRMS);
Ketpost10 = cellfun(@(x) mean(x(41:50)),KetRMS);

%manually tried ttest2 comparison with the following:
Awapre20 = cellfun(@(x) mean(x(10:20)),AwaRMS);%ns
[Hpre20, Ppre20, CIpre20, Tpre20] = ttest2(Ketpre20,Awapre20);
Awatone10 = cellfun(@(x) mean(x(21:30)),AwaRMS);%ns
[Htone10, Ptone10, CItone10, Ttone10] = ttest2(Kettone10,Awatone10);
Awatone20 = cellfun(@(x) mean(x(31:40)),AwaRMS);%ns
[Htone20, Ptone20, CItone20, Ttone20] = ttest2(Kettone20,Awatone20);
Awapost10 = cellfun(@(x) mean(x(41:50)),AwaRMS);%ns
[Hpost10, Ppost10, CIpost10, Tpost10] = ttest2(Ketpost10,Awapost10);

Muspre20 = cellfun(@(x) mean(x(10:20)),MusRMS);%***
Mustone10 = cellfun(@(x) mean(x(21:30)),MusRMS);%***
Mustone20 = cellfun(@(x) mean(x(31:40)),MusRMS);%***
Muspost10 = cellfun(@(x) mean(x(41:50)),MusRMS);%***

% super conservative gets us a little closer but stil ns
Kettone5 = cellfun(@(x) mean(x(21:25)),KetRMS);
Awatone5 = cellfun(@(x) mean(x(21:25)),AwaRMS); %ns
[Htone5, Ptone5, CItone5, Ttone5] = ttest2(Kettone5,Awatone5);

%% Analysis 3 permutation test?

KetRMSall = horzcat(KetRMS{:,1});
AwaRMSall = horzcat(AwaRMS{:,1}); %ns
[Hall, Pall, CIall, Statsall] = ttest2(KetRMSall,AwaRMSall);
MusRMSall = horzcat(MusRMS{:,1});

% Highest T
Obs_T = Statsall.tstat;

allRMS = vertcat(KetRMS, AwaRMS);
Perm_T = [];

for nperm = 1:1000
    permlist = randperm(length(allRMS));
    
    for ii = 1:length(KetRMS)
        Perm1{ii,1}=allRMS{permlist(ii),1};
    end
    
    for ii = length(KetRMS)+1:length(KetRMS)+length(AwaRMS)
        Perm2{ii-length(KetRMS),1}=allRMS{permlist(ii),1};
    end
    
    Perm1all = horzcat(Perm1{:,1});
    Perm2all = horzcat(Perm2{:,1});
    [Hp, Pp, CIp, Statsp] = ttest2(Perm1all,Perm2all);
    
    Perm_T = [Perm_T Statsp.tstat];
end

pthresh = 2.101;
sig_T = sum(Perm_T>Obs_T,2); 
sig_T(sig_T <= pthresh) = true; %if it's not more than 50, it's siginificant for that frequency
sig_T(sig_T > pthresh) = false; 
sig_T
%ns!! 
