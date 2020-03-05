%% Cluster Permutation Test

% Author: Katrina Deane

% following: https://benediktehinger.de/blog/science/statistics-cluster-permutation-test/
% This script is for the purpose of varifying first that there is a
% significant change between conditions by comparing them to chance through
% permutations and then to check at which layers and sinks that
% significance holds


%Input:     is from folder DATA; specifically (not automatically) named per 
%           Kat's MT groups
%Output:    is placed in figs - PermuteCSDs; observed difference matrix,
%           observed t matrix, observed cluster matrix in .fig / the
%           observed vs permustation curves across frequency bins in .pdf
%           and .png (grammplots), ALSO the stats output calculated in .mat

% NOTE: currently everything is set at the beginning of the code except the
% layer selection (e.g. LIV early is (8:11,200:300) and needs to be
% manually changed if there are other specifications required

clear

%% standard operations

% Change directory to your working folder
if exist('D:\MyCode\Dynamic_CSD_Analysis','dir') == 7
    cd('D:\MyCode\Dynamic_CSD_Analysis');
elseif exist('D:\Dynamic_CSD_Analysis','dir') == 7
    cd('D:\Dynamic_CSD_Analysis');
elseif exist('C:\Users\kedea\Documents\Dynamic_CSD_Analysis','dir') == 7
    cd('C:\Users\kedea\Documents\Dynamic_CSD_Analysis')
else
    error('Please add your working folder to the list in the header')
end

home = pwd; 
addpath(genpath(home));
cd (home),cd DATA;

Ticknames = {'BF - 3','BF - 2','BF - 1','BF','BF + 1', 'BF + 2', 'BF + 3'};
ticklen = length(Ticknames);
tickorder = {'a-3', 'b-2', 'c-1', 'dBF', 'e+1', 'f+2', 'g+3'};
nperms = 1000; %number of times we want the permutation

%% Which Comparison

which = 3;                  % 1 = Anesthetized Chronic vs Awake
                            % 2 = Anesthetized vs Muscimol
                            % 3 = Anesthetized vs Awake
                            
%% Load in required data

%load in data and call up the list of animal names
if which == 1
    load('ANChronic_Data.mat')
    First_Group = Data; clear Data; first_names = fieldnames(First_Group);
    load('Awake10dB_Data.mat')
    Second_Group = Data; clear Data; second_names = fieldnames(Second_Group);
    
    compName = 'Chronic';
elseif which == 2
    load('AnesthetizedPre_Data.mat')
    First_Group = Data; clear Data; first_names = fieldnames(First_Group);
    load('Muscimol_Data.mat')
    Second_Group = Data; clear Data; second_names = fieldnames(Second_Group);
    
    compName = 'Acute';
elseif which == 3
    load('AnesthetizedPre_Data.mat')
    First_Group = Data; clear Data; first_names = fieldnames(First_Group);
    load('Awake10dB_Data.mat')
    Second_Group = Data; clear Data; second_names = fieldnames(Second_Group);
    
    compName = 'Mixed';
end
%% CHANGE THIS IF GROUP SIZE CHANGES

%clusters defined by t value associated with with p=0.05 for n subjects
if which == 1 % current t value threshold ONLY good for this df
    df = length(second_names) - 1; % get true df to make sure it matches
    if df == 8
        t_thresh = 2.306; % http://www.ttable.org/
    else
        error('You need to change your tvalue threshold! Check this link: http://www.ttable.org/ (df is group size - 1)')
    end
elseif which == 2
    df = length(second_names) - 1; % get true df to make sure it matches
    if df == 10
        t_thresh = 2.228;
    else
        error('You need to change your tvalue threshold! Check this link: http://www.ttable.org/ (df is group size - 1)')
    end
elseif which == 3
    df = length(first_names)+length(second_names) - 2; % get true df to make sure it matches
    if df == 18
        t_thresh = 2.101;
    else
        error('You need to change your tvalue threshold! Check this link: http://www.ttable.org/ (df is group size - 1)')
    end
end

%% line up CSD's with +-3 around the BF for each group
first_cont = {}; second_cont = {}; 

for iA = 1:length(first_names) %for each animal
    BF = find((First_Group.(first_names{iA}).GS_BF) == (First_Group.(first_names{iA}).Frqz')); %get the BF Position
    
    for itick = -3:ticklen-4 %through each freqency around charted BF
        if BF+itick > 0 && BF+itick <= length(First_Group.(first_names{iA}).CSD) 
            first_cont{iA,itick+4} = First_Group.(first_names{iA}).CSD{BF+itick}(1:28,1:600); %had to cut to smallest form of each matrix (28x600)
        else
            first_cont{iA,itick+4} = NaN(28,600); %in case the frequency around the BF doesn't exist
        end
    end
end

for iA = 1:length(second_names) %for each animal
    BF = find((Second_Group.(second_names{iA}).GS_BF) == (Second_Group.(second_names{iA}).Frqz')); %get the BF Position
    
    for itick = -3:ticklen-4 %through each freqency around charted BF
        if BF+itick > 0 && BF+itick <= length(Second_Group.(second_names{iA}).CSD)
            second_cont{iA,itick+4} = Second_Group.(second_names{iA}).CSD{BF+itick}(1:28,1:600);
        else
            second_cont{iA,itick+4} = NaN(28,600);
        end
    end
end

% check containers for size and shape verification (should be number of
% animals by number of defined ticks with each cell having 28x600)

%% Conditions and differences - Permutation step 1
if which == 1 || which == 2 % these are within group and we can find difference between each animal
    obs_difmat = {};
    
    % for every animal, find the difference between condition 1 and condition 2
    % rows = animals, columns = frequency (4 = BF)
    
    for itick = 1:ticklen
        for iAnimal = 1:length(second_names)
            Cond1 = first_cont{iAnimal,itick};
            Cond2 = second_cont{iAnimal,itick};
            obs_difmat{iAnimal,itick} = Cond1 - Cond2;
        end
    end
    
    obs_difmeans = {};
    obs_difstd = {};
    
    for itick = 1:ticklen %creates 1x7 CSDs with 4 being the BF
        
        obs_difcat = cat(3,obs_difmat{1,itick},obs_difmat{2,itick},obs_difmat{3,itick},obs_difmat{4,itick}...
            ,obs_difmat{5,itick},obs_difmat{6,itick},obs_difmat{7,itick},obs_difmat{8,itick},obs_difmat{9,itick}); %9 animals
        
        obs_difmeans{itick} = nanmean(obs_difcat,3);
        obs_difstd{itick} = nanstd(obs_difcat,0,3); % currently using 0 for W but don't know what it should be KD
        
    end
    
elseif which == 3 % between group comparison
    obs1_means = {};
    obs1_std = {};
    obs2_means = {};
    obs2_std = {};
    obs_difmeans = {}; % fot the figure
    
    for itick = 1:length(Ticknames) %creates 1x7 CSDs with 4 being the BF
        
        cond1cat = cat(3,first_cont{1,itick},first_cont{2,itick},first_cont{3,itick},first_cont{4,itick}...
            ,first_cont{5,itick},first_cont{6,itick},first_cont{7,itick},first_cont{8,itick},first_cont{9,itick}...
            ,first_cont{10,itick},first_cont{11,itick});
        
        obs1_means{itick} = nanmean(cond1cat,3);
        
        cond2cat = cat(3,second_cont{1,itick},second_cont{2,itick},second_cont{3,itick},second_cont{4,itick}...
            ,second_cont{5,itick},second_cont{6,itick},second_cont{7,itick},second_cont{8,itick},second_cont{9,itick}); %9 animals
        
        obs2_means{itick} = nanmean(cond2cat,3);
      
        obs1_std{itick} = nanstd(cond1cat,0,3); % currently using 0 for W but don't know what it should be KD
        obs2_std{itick} = nanstd(cond2cat,0,3);
        
        obs_difmeans{itick} = obs1_means{itick} - obs2_means{itick};
    end
end

% check mean for size and shape verification (should be 1 by number of
% defined ticks with each cell still having 28x600)

%% t test - step 2
%find the t values along all data points for each frequency bin
obs_t = {}; %observed conditions t values
obs_maxt = {}; % to compare max t rather than clustermass
obs_clusters = {}; %clusters = 1, all else is NaN
obs_clustermass = {}; %total number of 1's (not currently done: 8400/mass = ratio)

obs_tstat = {}; %to hold the t scores of the inbuilt ttest function (single vector output)
obs_H = {}; %to hold the H scores (1 or 0)


for itick = 1:ticklen
    if which == 1 || which == 2 % within group ttest
        obs_t{itick} = obs_difmeans{itick}./(obs_difstd{itick}/sqrt(length(second_names)));
        obs_clusters{itick} = abs(obs_t{itick});
        obs_clusters{itick}(obs_clusters{itick} < t_thresh) = NaN;
        obs_clusters{itick}(obs_clusters{itick} >= t_thresh) = 1;
        
        % check cluster mass for 300 ms from tone onset
        obs_clustermass{itick} = nansum(nansum(obs_clusters{itick}(:,200:500)));
        
        % check max t value of matrix
        obs_maxt{itick} = nanmax(nanmax(abs(obs_t{itick})));
        
        % ttest for line
        [H,~,~,STATS] = ttest(obs_difmeans{itick});
        obs_tstat{itick} = STATS.tstat;
        obs_H{itick} = H;
        
    elseif which == 3 % between group ttest2
        obs_t{itick} = (obs1_means{itick} - obs2_means{itick})./...
            sqrt((obs1_std{itick}.^2/length(first_names)) + (obs2_std{itick}.^2/length(second_names)));
        obs_clusters{itick} = abs(obs_t{itick});
        obs_clusters{itick}(obs_clusters{itick} < t_thresh) = NaN;
        obs_clusters{itick}(obs_clusters{itick} >= t_thresh) = 1;
        
        % check cluster mass for 300 ms from tone onset
        obs_clustermass{itick} = nansum(nansum(obs_clusters{itick}(:,200:500)));
        
        % check max t value of matrix
        obs_maxt{itick} = nanmax(nanmax(abs(obs_t{itick})));
        
        % ttest for line
        [H,~,~,STATS] = ttest2(obs1_means{itick},obs2_means{itick},'vartype','unequal');
        obs_tstat{itick} = STATS.tstat;
        obs_H{itick} = H;
    end
end

%% for layer specific: 
% % pull out clusters
obs_layer = struct;
obs_layer.II_early = cellfun(@(x) x(1:4,200:300),obs_clusters, 'UniformOutput', false);% LIVearly = (8:11,200:300)
obs_layer.IV_early = cellfun(@(x) x(8:11,200:300),obs_clusters, 'UniformOutput', false);% LIVearly = (8:11,200:300)
obs_layer.Vb_early = cellfun(@(x) x(19:22,200:300),obs_clusters, 'UniformOutput', false);% LVbearly = (19:22,200:300)
obs_layer.VIa_early = cellfun(@(x) x(23:26,200:300),obs_clusters, 'UniformOutput', false);% LVIaearly = (23:26,200:300)
obs_layer.II_late = cellfun(@(x) x(1:4,300:500),obs_clusters, 'UniformOutput', false); % late = 300:500 ; a bit long...
obs_layer.IV_late = cellfun(@(x) x(8:11,300:500),obs_clusters, 'UniformOutput', false); % late = 300:500 ; a bit long...
obs_layer.Vb_late = cellfun(@(x) x(19:22,300:500),obs_clusters, 'UniformOutput', false);
obs_layer.VIa_late = cellfun(@(x) x(23:26,300:500),obs_clusters, 'UniformOutput', false);

% % sum clusters (twice to get final value)
for i = 1:2
    obs_layer.II_early = cellfun(@nansum, obs_layer.II_early, 'UniformOutput', false);
    obs_layer.IV_early = cellfun(@nansum, obs_layer.IV_early, 'UniformOutput', false);
    obs_layer.Vb_early = cellfun(@nansum, obs_layer.Vb_early, 'UniformOutput', false);
    obs_layer.VIa_early = cellfun(@nansum, obs_layer.VIa_early, 'UniformOutput', false);
    obs_layer.II_late = cellfun(@nansum, obs_layer.II_late, 'UniformOutput', false);
    obs_layer.IV_late = cellfun(@nansum, obs_layer.IV_late, 'UniformOutput', false);
    obs_layer.Vb_late = cellfun(@nansum, obs_layer.Vb_late, 'UniformOutput', false);
    obs_layer.VIa_late = cellfun(@nansum, obs_layer.VIa_late, 'UniformOutput', false);
end

% pull out areas of t mat
obs_layermax = struct;
obs_layermax.II_early = cellfun(@(x) x(1:4,200:300),obs_t, 'UniformOutput', false);% LIIearly = (1:4,200:300)
obs_layermax.IV_early = cellfun(@(x) x(8:11,200:300),obs_t, 'UniformOutput', false);% LIVearly = (8:11,200:300)
obs_layermax.Vb_early = cellfun(@(x) x(19:22,200:300),obs_t, 'UniformOutput', false);% LVbearly = (19:22,200:300)
obs_layermax.VIa_early = cellfun(@(x) x(23:26,200:300),obs_t, 'UniformOutput', false);% LVIaearly = (23:26,200:300)
obs_layermax.II_late = cellfun(@(x) x(1:4,300:500),obs_t, 'UniformOutput', false); % late = 300:500 ; a bit long...
obs_layermax.IV_late = cellfun(@(x) x(8:11,300:500),obs_t, 'UniformOutput', false); % late = 300:500 ; a bit long...
obs_layermax.Vb_late = cellfun(@(x) x(19:22,300:500),obs_t, 'UniformOutput', false);
obs_layermax.VIa_late = cellfun(@(x) x(23:26,300:500),obs_t, 'UniformOutput', false);

%make it absolute values
obs_layermax.II_early = cellfun(@abs,obs_layermax.II_early, 'UniformOutput', false);% LIIearly = (1:4,200:300)
obs_layermax.IV_early = cellfun(@abs,obs_layermax.IV_early, 'UniformOutput', false);% LIVearly = (8:11,200:300)
obs_layermax.Vb_early = cellfun(@abs,obs_layermax.Vb_early, 'UniformOutput', false);% LVbearly = (19:22,200:300)
obs_layermax.VIa_early = cellfun(@abs,obs_layermax.VIa_early, 'UniformOutput', false);% LVIaearly = (23:26,200:300)
obs_layermax.II_late = cellfun(@abs,obs_layermax.II_late, 'UniformOutput', false); % late = 300:500 ; a bit long...
obs_layermax.IV_late = cellfun(@abs,obs_layermax.IV_late, 'UniformOutput', false); % late = 300:500 ; a bit long...
obs_layermax.Vb_late = cellfun(@abs,obs_layermax.Vb_late, 'UniformOutput', false);
obs_layermax.VIa_late = cellfun(@abs,obs_layermax.VIa_late, 'UniformOutput', false);

% find the very maximum
for i = 1:2
    obs_layermax.II_early = cellfun(@nanmax, obs_layermax.II_early, 'UniformOutput', false);
    obs_layermax.IV_early = cellfun(@nanmax, obs_layermax.IV_early, 'UniformOutput', false);
    obs_layermax.Vb_early = cellfun(@nanmax, obs_layermax.Vb_early, 'UniformOutput', false);
    obs_layermax.VIa_early = cellfun(@nanmax, obs_layermax.VIa_early, 'UniformOutput', false);
    obs_layermax.II_late = cellfun(@nanmax, obs_layermax.II_late, 'UniformOutput', false);
    obs_layermax.IV_late = cellfun(@nanmax, obs_layermax.IV_late, 'UniformOutput', false);
    obs_layermax.Vb_late = cellfun(@nanmax, obs_layermax.Vb_late, 'UniformOutput', false);
    obs_layermax.VIa_late = cellfun(@nanmax, obs_layermax.VIa_late, 'UniformOutput', false);
end


%% permute and calculate the above - step 3
permstruct = struct;
permstruct.difmat = {};
permstruct.difmeans = {};
permstruct.maxt = {};
permstruct.clusters = {};
permstruct.clustermass = [];
permstruct.t = [];

perm_layer = struct;
perm_layermassII_early = [];perm_layermassII_late = [];
perm_layermassIV_early = [];perm_layermassIV_late = [];
perm_layermassVb_early = [];perm_layermassVb_late = [];
perm_layermassVIa_early = [];perm_layermassVIa_late = [];

perm_layermax = struct;
perm_layermaxII_early = [];perm_layermaxII_late = [];
perm_layermaxIV_early = [];perm_layermaxIV_late = [];
perm_layermaxVb_early = [];perm_layermaxVb_late = [];
perm_layermaxVIa_early = [];perm_layermaxVIa_late = [];

for iperm = 1:nperms
    
    if which == 1 || which == 2 % within groups
        perm_difmat = {};
        
        % for every animal, find the difference between condition 1 and condition 2
        % rows = animals, columns = frequency (4 = BF)
        
        for itick = 1:ticklen
            for iAnimal = 1:length(second_names)
                P = randi([1,2],1); % get a random P either 1 or 2 for each animal
                
                Cond1 = first_cont{iAnimal,itick};
                Cond2 = second_cont{iAnimal,itick};
                
                if P == 1 %Cond1 leading
                    perm_difmat{iAnimal,itick} = Cond1 - Cond2;
                    
                elseif P ==2 %Cond2 leading
                    perm_difmat{iAnimal,itick} = Cond2 - Cond1;
                end
                
            end
        end
        
        perm_difmeans = {};
        for itick = 1:ticklen %creates 1x7 CSDs with 4 being the BF
            
            perm_difcat = cat(3,perm_difmat{1,itick},perm_difmat{2,itick},perm_difmat{3,itick},perm_difmat{4,itick}...
                ,perm_difmat{5,itick},perm_difmat{6,itick},perm_difmat{7,itick},perm_difmat{8,itick},perm_difmat{9,itick}); %9 animals
            
            perm_difmeans{itick} = nanmean(perm_difcat,3);
            perm_difstd{itick} = nanstd(perm_difcat,0,3); % currently using 0 for W but don't know what it should be KD
            
        end
    elseif which == 3 % between groups
        perm1_means = {};
        perm1_std = {};
        perm2_means = {};
        perm2_std = {};
        
        for itick = 1:ticklen %creates 1x7 CSDs with 4 being the BF
            
            allcondcat = cat(3,first_cont{1,itick},first_cont{2,itick},first_cont{3,itick},first_cont{4,itick}...
                ,first_cont{5,itick},first_cont{6,itick},first_cont{7,itick},first_cont{8,itick},first_cont{9,itick}...
                ,first_cont{10,itick},first_cont{11,itick},second_cont{1,itick},second_cont{2,itick},...
                second_cont{3,itick},second_cont{4,itick},second_cont{5,itick},second_cont{6,itick},...
                second_cont{7,itick},second_cont{8,itick},second_cont{9,itick});
            
            Pv = randperm(df+2);
            
            cond1perm = cat(3,allcondcat(:,:,Pv(1)),allcondcat(:,:,Pv(2)),allcondcat(:,:,Pv(3)),allcondcat(:,:,Pv(4))...
                ,allcondcat(:,:,Pv(5)),allcondcat(:,:,Pv(6)),allcondcat(:,:,Pv(7)),allcondcat(:,:,Pv(8))...
                ,allcondcat(:,:,Pv(9)),allcondcat(:,:,Pv(10)),allcondcat(:,:,Pv(11)));
            
            perm1_means{itick} = nanmean(cond1perm,3);
            
            cond2perm = cat(3,allcondcat(:,:,Pv(12)),allcondcat(:,:,Pv(13)),allcondcat(:,:,Pv(14)),allcondcat(:,:,Pv(15))...
                ,allcondcat(:,:,Pv(16)),allcondcat(:,:,Pv(17)),allcondcat(:,:,Pv(18)),allcondcat(:,:,Pv(19))...
                ,allcondcat(:,:,Pv(20)));
            
            perm2_means{itick} = nanmean(cond2perm,3);
            
            perm1_std{itick} = nanstd(cond1perm,0,3); % currently using 0 for W but don't know what it should be KD
            perm2_std{itick} = nanstd(cond2perm,0,3);

        end
    end
    
    % t test
   
    % check mean for size and shape verification (should be 1 by number of
    % defined ticks with each cell still having 28x600)
    
    %find the t values along all data points for each frequency bin
    perm_t = {}; %observed conditions t values
    perm_clusters = {}; %clusters = 1, all else is NaN
    perm_clustermass = {}; %total number of 1's (not currently done: 8400/mass = ratio)
    % t = mean/(std/sqrt(n))
    for itick = 1:ticklen
        if which == 1 || which == 2
            perm_t{itick} = perm_difmeans{itick}./(perm_difstd{itick}/sqrt(length(second_names)));
            perm_clusters{itick} = abs(perm_t{itick});
            perm_clusters{itick}(perm_clusters{itick} < t_thresh) = NaN;
            perm_clusters{itick}(perm_clusters{itick} >= t_thresh) = 1;
            
            % check cluster mass for 300 ms from tone onset
            perm_clustermass{itick} = nansum(nansum(perm_clusters{itick}(:,200:500)));
            
            % check max t value of matrix
            perm_maxt{itick} = nanmax(nanmax(abs(perm_t{itick})));
            
        elseif which == 3 % between group ttest2
            perm_t{itick} = (perm1_means{itick} - perm2_means{itick})./...
                sqrt((perm1_std{itick}.^2/length(first_names)) + (perm2_std{itick}.^2/length(second_names)));
            perm_clusters{itick} = abs(perm_t{itick});
            perm_clusters{itick}(perm_clusters{itick} < t_thresh) = NaN;
            perm_clusters{itick}(perm_clusters{itick} >= t_thresh) = 1;
            
            % check cluster mass for 300 ms from tone onset
            perm_clustermass{itick} = nansum(nansum(perm_clusters{itick}(:,200:500)));
            
            % check max t value of matrix
            perm_maxt{itick} = nanmax(nanmax(abs(perm_t{itick})));
        end
    end
    
%     permstruct.difmat{iperm} = perm_difmat;
%     permstruct.difmeans{iperm} = perm_difmeans;
    permstruct.t{iperm} = perm_t;
    permstruct.maxt = horzcat(permstruct.maxt,perm_maxt);
    permstruct.clusters{iperm} = perm_clusters;
    permstruct.clustermass = horzcat(permstruct.clustermass,perm_clustermass);
    
    % for layer specific:
    % % pull out clusters
    perm_layer.II_early{iperm} = cellfun(@(x) x(1:4,200:300),permstruct.clusters{iperm}, 'UniformOutput', false);% LIIearly = (1:4,200:300)
    perm_layer.IV_early{iperm} = cellfun(@(x) x(8:11,200:300),permstruct.clusters{iperm}, 'UniformOutput', false);% LIVearly = (8:11,200:300)
    perm_layer.Vb_early{iperm} = cellfun(@(x) x(19:22,200:300),permstruct.clusters{iperm}, 'UniformOutput', false);% LVbearly = (19:22,200:300)
    perm_layer.VIa_early{iperm} = cellfun(@(x) x(23:26,200:300),permstruct.clusters{iperm}, 'UniformOutput', false);% LVIaearly = (23:26,200:300)
    perm_layer.II_late{iperm} = cellfun(@(x) x(1:4,300:500),permstruct.clusters{iperm}, 'UniformOutput', false); % late = 300:500 ; a bit long...
    perm_layer.IV_late{iperm} = cellfun(@(x) x(8:11,300:500),permstruct.clusters{iperm}, 'UniformOutput', false); % late = 300:500 ; a bit long...
    perm_layer.Vb_late{iperm} = cellfun(@(x) x(19:22,300:500),permstruct.clusters{iperm}, 'UniformOutput', false);
    perm_layer.VIa_late{iperm} = cellfun(@(x) x(23:26,300:500),permstruct.clusters{iperm}, 'UniformOutput', false);
    
    % % sum clusters (twice to get final value)
    layer_clustermassII_early = {};layer_clustermassII_late = {};
    layer_clustermassIV_early = {};layer_clustermassIV_late = {};
    layer_clustermassVb_early = {};layer_clustermassVb_late = {};
    layer_clustermassVIa_early = {};layer_clustermassVIa_late = {};
    for itick = 1:ticklen
        layer_clustermassII_early{itick} = nansum(nansum(perm_layer.II_early{1,iperm}{1,itick}));
        layer_clustermassIV_early{itick} = nansum(nansum(perm_layer.IV_early{1,iperm}{1,itick}));
        layer_clustermassVb_early{itick} = nansum(nansum(perm_layer.Vb_early{1,iperm}{1,itick}));
        layer_clustermassVIa_early{itick} = nansum(nansum(perm_layer.VIa_early{1,iperm}{1,itick}));
        
        layer_clustermassII_late{itick} = nansum(nansum(perm_layer.II_late{1,iperm}{1,itick}));
        layer_clustermassIV_late{itick} = nansum(nansum(perm_layer.IV_late{1,iperm}{1,itick}));
        layer_clustermassVb_late{itick} = nansum(nansum(perm_layer.Vb_late{1,iperm}{1,itick}));
        layer_clustermassVIa_late{itick} = nansum(nansum(perm_layer.VIa_late{1,iperm}{1,itick}));
    end
    % % add them to the list
    perm_layermassII_early = horzcat(perm_layermassII_early,layer_clustermassII_early);
    perm_layermassIV_early = horzcat(perm_layermassIV_early,layer_clustermassIV_early);
    perm_layermassVb_early = horzcat(perm_layermassVb_early,layer_clustermassVb_early);
    perm_layermassVIa_early = horzcat(perm_layermassVIa_early,layer_clustermassVIa_early);
    
    perm_layermassII_late = horzcat(perm_layermassII_late,layer_clustermassII_late);
    perm_layermassIV_late = horzcat(perm_layermassIV_late,layer_clustermassIV_late);
    perm_layermassVb_late= horzcat(perm_layermassVb_late,layer_clustermassVb_late);
    perm_layermassVIa_late = horzcat(perm_layermassVIa_late,layer_clustermassVIa_late);
    

    % % pull out t mat
    perm_layermax.II_early{iperm} = cellfun(@(x) x(1:4,200:300),permstruct.t{iperm}, 'UniformOutput', false);% LIIearly = (1:4,200:300)
    perm_layermax.IV_early{iperm} = cellfun(@(x) x(8:11,200:300),permstruct.t{iperm}, 'UniformOutput', false);% LIVearly = (8:11,200:300)
    perm_layermax.Vb_early{iperm} = cellfun(@(x) x(19:22,200:300),permstruct.t{iperm}, 'UniformOutput', false);% LVbearly = (19:22,200:300)
    perm_layermax.VIa_early{iperm} = cellfun(@(x) x(23:26,200:300),permstruct.t{iperm}, 'UniformOutput', false);% LVIaearly = (23:26,200:300)
    perm_layermax.II_late{iperm} = cellfun(@(x) x(1:4,300:500),permstruct.t{iperm}, 'UniformOutput', false); % late = 300:500 ; a bit long...
    perm_layermax.IV_late{iperm} = cellfun(@(x) x(8:11,300:500),permstruct.t{iperm}, 'UniformOutput', false); % late = 300:500 ; a bit long...
    perm_layermax.Vb_late{iperm} = cellfun(@(x) x(19:22,300:500),permstruct.t{iperm}, 'UniformOutput', false);
    perm_layermax.VIa_late{iperm} = cellfun(@(x) x(23:26,300:500),permstruct.t{iperm}, 'UniformOutput', false);
    
    % get absolute values
    perm_layermax.II_early{iperm} = cellfun(@abs,perm_layermax.II_early{iperm}, 'UniformOutput', false);% LIIearly = (1:4,200:300)
    perm_layermax.IV_early{iperm} = cellfun(@abs,perm_layermax.IV_early{iperm}, 'UniformOutput', false);% LIVearly = (8:11,200:300)
    perm_layermax.Vb_early{iperm} = cellfun(@abs,perm_layermax.Vb_early{iperm}, 'UniformOutput', false);% LVbearly = (19:22,200:300)
    perm_layermax.VIa_early{iperm} = cellfun(@abs,perm_layermax.VIa_early{iperm}, 'UniformOutput', false);% LVIaearly = (23:26,200:300)
    perm_layermax.II_late{iperm} = cellfun(@abs,perm_layermax.II_late{iperm}, 'UniformOutput', false); % late = 300:500 ; a bit long...
    perm_layermax.IV_late{iperm} = cellfun(@abs,perm_layermax.IV_late{iperm}, 'UniformOutput', false); % late = 300:500 ; a bit long...
    perm_layermax.Vb_late{iperm} = cellfun(@abs,perm_layermax.Vb_late{iperm}, 'UniformOutput', false);
    perm_layermax.VIa_late{iperm} = cellfun(@abs,perm_layermax.VIa_late{iperm}, 'UniformOutput', false);
    
    % % sum clusters (twice to get final value)
    layer_clustermaxII_early = {};layer_clustermaxII_late = {};
    layer_clustermaxIV_early = {};layer_clustermaxIV_late = {};
    layer_clustermaxVb_early = {};layer_clustermaxVb_late = {};
    layer_clustermaxVIa_early = {};layer_clustermaxVIa_late = {};
    for itick = 1:ticklen
        layer_clustermaxII_early{itick} = nanmax(nanmax(perm_layermax.II_early{1,iperm}{1,itick}));
        layer_clustermaxIV_early{itick} = nanmax(nanmax(perm_layermax.IV_early{1,iperm}{1,itick}));
        layer_clustermaxVb_early{itick} = nanmax(nanmax(perm_layermax.Vb_early{1,iperm}{1,itick}));
        layer_clustermaxVIa_early{itick} = nanmax(nanmax(perm_layermax.VIa_early{1,iperm}{1,itick}));
        
        layer_clustermaxII_late{itick} = nanmax(nanmax(perm_layermax.II_late{1,iperm}{1,itick}));
        layer_clustermaxIV_late{itick} = nanmax(nanmax(perm_layermax.IV_late{1,iperm}{1,itick}));
        layer_clustermaxVb_late{itick} = nanmax(nanmax(perm_layermax.Vb_late{1,iperm}{1,itick}));
        layer_clustermaxVIa_late{itick} = nanmax(nanmax(perm_layermax.VIa_late{1,iperm}{1,itick}));
    end
    % % add them to the list
    perm_layermaxII_early = horzcat(perm_layermaxII_early,layer_clustermaxII_early);
    perm_layermaxIV_early = horzcat(perm_layermaxIV_early,layer_clustermaxIV_early);
    perm_layermaxVb_early = horzcat(perm_layermaxVb_early,layer_clustermaxVb_early);
    perm_layermaxVIa_early = horzcat(perm_layermaxVIa_early,layer_clustermaxVIa_early);
    
    perm_layermaxII_late = horzcat(perm_layermaxII_late,layer_clustermaxII_late);
    perm_layermaxIV_late = horzcat(perm_layermaxIV_late,layer_clustermaxIV_late);
    perm_layermaxVb_late= horzcat(perm_layermaxVb_late,layer_clustermaxVb_late);
    perm_layermaxVIa_late = horzcat(perm_layermaxVIa_late,layer_clustermaxVIa_late);
  
end

%% Visualize observed phenomonon
cd(home);cd figs; mkdir('PermuteCSDs'); cd('PermuteCSDs')

% pull out all visuals associated
figure('Name',[compName ' Observed Difference Matrix'])
for itick = 1:ticklen
    subplot(2,round(ticklen/2),itick)
    imagesc(obs_difmeans{1,itick}(:,200:500))
    caxis([-0.0005 0.0005])
%     colormap('jet')
    colorbar
    title(Ticknames{itick})
end

h = gcf;
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(h,[compName ' Observed Difference Matrix'],'compact')
close (h)

figure('Name',[compName ' Observed T Matrix'])
for itick = 1:ticklen
    subplot(2,round(ticklen/2),itick)
    imagesc(obs_t{1,itick}(:,200:500))
    caxis([-5 5])
    colormap('jet')
    title(Ticknames{itick})
end

h = gcf;
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(h,[compName ' Observed T Matrix'],'compact')
close (h)

figure('Name',[compName ' Observed Clusters Matrix'])
for itick = 1:ticklen
    subplot(2,round(ticklen/2),itick)
    imagesc(obs_clusters{1,itick}(:,200:500))
    caxis([-1 1])
    colormap('winter')
    title(obs_clustermass{itick})
end

h = gcf;
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(h,[compName ' Observed Cluster Matrix'],'compact')
close (h)


%% Get some lines for the luls

figure('Name',[compName ' Observed T Curves'])
for itick = 1:ticklen
    
%     for the full thing, we can run a ttest also
%         full = nanmean(obs_t{1,itick}(:,200:500));
%         fullp = nansum(obs_clusters{1,itick}(:,200:500));
%         sigfull = full; sigfull(fullp < 14) = NaN;
 
    fullttest = obs_tstat{itick}(:,200:500);
    sigH = fullttest;
    sigH(obs_H{itick}(:,200:500) == 0) = NaN;
    
    layerII = nanmean(obs_t{1,itick}(1:4,200:500)); %the average t's
    layer_IIp = nansum(obs_clusters{1,itick}(1:4,200:500)); %the sum of sig values in each column
    layerIV = nanmean(obs_t{1,itick}(8:11,200:500)); %the average t's
    layer_IVp = nansum(obs_clusters{1,itick}(8:11,200:500)); %the sum of sig values in each column
    layerVb = nanmean(obs_t{1,itick}(19:22,200:500));
    layer_Vbp = nansum(obs_clusters{1,itick}(19:22,200:500));
    layerVIa = nanmean(obs_t{1,itick}(23:26, 200:500));
    layer_VIap = nansum(obs_clusters{1,itick}(23:26,200:500));
    
    siglayerII = layerII; siglayerII(layer_IIp < 2) = NaN;  %if the sum is less than half, don't show it
    siglayerIV = layerIV; siglayerIV(layer_IVp < 2) = NaN;  %if the sum is less than half, don't show it
    siglayerVb = layerVb; siglayerVb(layer_Vbp < 2) = NaN;
    siglayerVIa = layerVIa; siglayerVIa(layer_VIap < 2) = NaN;
    
    subplot(2,round(ticklen/2),itick)
    plot(fullttest); hold on
    plot(layerII); 
    plot(layerIV); 
    plot(layerVb)
    plot(layerVIa)
    plot(sigH, 'LineWidth',2)
    plot(siglayerII, 'LineWidth',2)
    plot(siglayerIV, 'LineWidth',2)
    plot(siglayerVb, 'LineWidth',2)
    plot(siglayerVIa, 'LineWidth',2)
    legend('Full','LII','LIV','LVb','LVIa','full*','IV*','II*','Vb*','VIa*')
    title(Ticknames{itick})
end

h = gcf;
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(h,[compName ' Observed T Curves'],'compact')
close (h)

%% Visualize cluster masses from permutations

%function saves grammplots and gives out the final significance calculation
[sig_mass, pVal, permMean, permSTD] = Perm_makegramm(nperms,obs_clustermass, permstruct.clustermass, ...
    obs_maxt, permstruct.maxt, tickorder, ticklen, Ticknames, ...
    [compName ' Observed vs Permutation']);


%% Visualize cluster masses in layer bins and check significance

%II Early
[IIE_sig_mass, IIE_pVal, IIE_permMean, IIE_permSTD] = Perm_makegramm(nperms,obs_layer.II_early, perm_layermassII_early, ...
    obs_layermax.II_early,perm_layermaxII_early, tickorder, ticklen, Ticknames, ...
    [compName ' LII early Observed vs Permutation']);

%IV Early
[IVE_sig_mass, IVE_pVal, IVE_permMean, IVE_permSTD] = Perm_makegramm(nperms,obs_layer.IV_early, perm_layermassIV_early, ...
    obs_layermax.IV_early,perm_layermaxIV_early, tickorder, ticklen, Ticknames, ...
    [compName ' LIV early Observed vs Permutation']);


%Vb Early
[VbE_sig_mass, VbE_pVal, VbE_permMean, VbE_permSTD] = Perm_makegramm(nperms,obs_layer.Vb_early, perm_layermassVb_early, ...
    obs_layermax.Vb_early,perm_layermaxVb_early,tickorder, ticklen, Ticknames, ...
    [compName ' LVb early Observed vs Permutation']);


%VIa Early
[VIaE_sig_mass, VIaE_pVal, VIaE_permMean, VIaE_permSTD] = Perm_makegramm(nperms,obs_layer.VIa_early, perm_layermassVIa_early, ...
    obs_layermax.VIa_early,perm_layermaxVIa_early,tickorder, ticklen, Ticknames, ...
    [compName ' LVIa early Observed vs Permutation']);


%II Late
[IIL_sig_mass, IIL_pVal, IIL_permMean, IIL_permSTD] = Perm_makegramm(nperms,obs_layer.II_late, perm_layermassII_late, ...
    obs_layermax.II_late,perm_layermaxII_late, tickorder, ticklen, Ticknames, ...
    [compName ' LII late Observed vs Permutation']);

%IV Late
[IVL_sig_mass, IVL_pVal, IVL_permMean, IVL_permSTD] = Perm_makegramm(nperms,obs_layer.IV_late, perm_layermassIV_late, ...
    obs_layermax.IV_late,perm_layermaxIV_late, tickorder, ticklen, Ticknames, ...
    [compName ' LIV late Observed vs Permutation']);


%Vb Early
[VbL_sig_mass, VbL_pVal, VbL_permMean, VbL_permSTD] = Perm_makegramm(nperms,obs_layer.Vb_late, perm_layermassVb_late, ...
    obs_layermax.Vb_late,perm_layermaxVb_late, tickorder, ticklen, Ticknames,...
    [compName ' LVb late Observed vs Permutation']);


%VIa Early
[VIaL_sig_mass, VIaL_pVal, VIaL_permMean, VIaL_permSTD] = Perm_makegramm(nperms,obs_layer.VIa_late, perm_layermassVIa_late, ...
    obs_layermax.VIa_late,perm_layermaxVIa_late, tickorder, ticklen,...
    Ticknames, [compName ' LVIa late Observed vs Permutation']);


% LIVearly = (8:11,200:300)
% LVbearly = (19:22,200:300)
% LVIaearly = (23:26,200:300)

% LIVlate = (8:11,300:500)
% LVblate = (19:22,300:500)
% LVIalate = (23:26,300:500)

%% Save all the stats to a .mat


if which == 1
    save Chronic_permstats.mat sig_mass sig_max IIE_sig_mass IIE_sig_max IVE_sig_mass IVE_sig_max ...
        VbE_sig_mass VbE_sig_max VIaE_sig_mass VIaE_sig_max IIL_sig_mass IIL_sig_max IVL_sig_mass ...
        IVL_sig_max VbL_sig_mass VbL_sig_max VIaL_sig_mass VIaL_sig_max obs_tstat obs_H ...
        IIE_pVal IIL_pVal IVE_pVal IVL_pVal VbL_pVal VbE_pVal VIaE_pVal VIaL_pVal pVal
elseif which == 2
    save Acute_permstats.mat sig_mass sig_max IIE_sig_mass IIE_sig_max IVE_sig_mass IVE_sig_max ...
        VbE_sig_mass VbE_sig_max VIaE_sig_mass VIaE_sig_max IIL_sig_mass IIL_sig_max IVL_sig_mass ...
        IVL_sig_max VbL_sig_mass VbL_sig_max VIaL_sig_mass VIaL_sig_max obs_tstat obs_H ...
        IIE_pVal IIL_pVal IVE_pVal IVL_pVal VbL_pVal VbE_pVal VIaE_pVal VIaL_pVal pVal
elseif which == 3
    save Mixed_permstats.mat sig_mass IIE_sig_mass IVE_sig_mass ...
        VbE_sig_mass VIaE_sig_mass IIL_sig_mass IVL_sig_mass ...
        VbL_sig_mass VIaL_sig_mass obs_tstat obs_H...
        IIE_pVal IIL_pVal IVE_pVal IVL_pVal VbL_pVal VbE_pVal VIaE_pVal VIaL_pVal pVal...
        permMean permSTD IIE_permMean IIE_permSTD IVE_permMean IVE_permSTD ...
        VbE_permMean VbE_permSTD VIaE_permMean VIaE_permSTD IIL_permMean ...
        IIL_permSTD IVL_permMean IVL_permSTD ...
        VbL_permMean VbL_permSTD VIaL_permMean VIaL_permSTD
end








