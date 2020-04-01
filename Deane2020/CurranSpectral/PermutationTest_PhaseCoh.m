function PermutationTest_PhaseCoh(inputmat,inputname)

% Input:    Layer to analyze, (possible input: relative to BF)
%           Needs scalogramsfull.mat from Andrew Curran's wavelet analysis
% Output:   Figures for means and observed difference of awake/ketamine
%           comparison; figures for observed t values, clusters, ttest line
%           output; boxplot and significance of permutation test -> Pictures folder

%% standard operations

warning('OFF');
dbstop if error

% Change directory to your working folder
if ~exist('homedir','var')
    if exist('D:\MyCode\Dynamic_CSD','dir') == 7
        cd('D:\MyCode\Dynamic_CSD');
    elseif exist('C:\Users\kedea\Documents\Work Stuff\Dynamic_CSD','dir') == 7
        cd('C:\Users\kedea\Documents\Work Stuff\Dynamic_CSD')
    end
    
    homedir = pwd;
    addpath(genpath(homedir));
end

cd (homedir),cd DATA;

nperms = 1000;
pthresh = nperms*(0.05/7); % Bonferroni corrected for 7 tests
grpsizeA = 9; grpsizeK = 11;

% frequencies can be found in wtTable.freq{1} to clarify the following
% rows choices; actual intended rows commented
theta = (49:54);        %(4:7);
alpha = (44:48);        %(8:12);
beta_low = (39:43);     %(13:18);
beta_high = (34:38);    %(19:30);
gamma_low = (26:33);    %(31:60);
gamma_high = (19:25);   %(61:100);

osciName = {'theta' 'alpha' 'beta_low' 'beta_high' 'gamma_low' 'gamma_high'};
osciRows = {theta alpha beta_low beta_high gamma_low gamma_high};

%% Load in and seperate Data
load(inputmat, 'wtTable')
name = inputname;

% varNames = unique(wtTable.layer);
params.startTime = -0.2; % seconds
params.limit = 600;
allnames = unique(wtTable.animal);

% Get phase coherence for groups
AwakeAll = nan(54,params.limit,grpsizeA);
KetAll   = nan(54,params.limit,grpsizeK);
MuscAll  = nan(54,params.limit,grpsizeK);

% through each animal 
for iAn = 1: grpsizeA + grpsizeK
    
    if iAn <= grpsizeA
        % pull out the data for this group and this animal
        awake = table2cell(wtTable(contains(wtTable.condition,'Awake')& contains(wtTable.animal,allnames{iAn}),1)); %&wtTable.rel2Bf==rel2BFin
        awake =  cellfun(@(x) x(:,1:params.limit),awake,'UniformOutput',false);
        % set up output cells for transformed data
        transAwake = cell(size(awake));
        
        % take the phase z/abs(z) for each single trial
        for iTrial = 1:length(awake)
            curTrial = awake{iTrial};
            transAwake{iTrial} = curTrial./abs(curTrial);
        end
        
        % for the subject:
        % get the mean of single trials and the absolute for that mean
        curAn = abs(mean(cat(3,transAwake{:}),3));
        AwakeAll(:,:,iAn) = curAn;
        
    else
        anest = table2cell(wtTable(contains(wtTable.condition,'Anesth')& contains(wtTable.animal,allnames{iAn}),1));
        anest =  cellfun(@(x) x(:,1:params.limit),anest,'UniformOutput',false);
        % set up output cells for transformed data
        transKet = cell(size(anest));
        
        % take the phase z/abs(z) for each single trial
        for iTrial = 1:length(anest)
            curTrial = anest{iTrial};
            transKet{iTrial} = curTrial./abs(curTrial);
        end
        
        % get the mean of single trials and the absolute for that mean
        curAn = abs(mean(cat(3,transKet{:}),3));
        KetAll(:,:,iAn-grpsizeA) = curAn;
        
        musc = table2cell(wtTable(contains(wtTable.condition,'Muscimol')& contains(wtTable.animal,allnames{iAn}),1));
        musc =  cellfun(@(x) x(:,1:params.limit),musc,'UniformOutput',false);
        % set up output cells for transformed data
        transMusc = cell(size(musc));
        
        % take the phase z/abs(z) for each single trial
        for iTrial = 1:length(musc)
            curTrial = musc{iTrial};
            transMusc{iTrial} = curTrial./abs(curTrial);
        end
        
        % get the mean of single trials and the absolute for that mean
        curAn = abs(mean(cat(3,transMusc{:}),3));
        MuscAll(:,:,iAn-grpsizeA) = curAn;
    end
end

% for the group: 
% take the mean across subjects
grpA = mean(AwakeAll,3);
grpK = mean(KetAll,3);
% grpM = mean(MuscAll,3);

diffK_A = grpK - grpA;

%% dif fig
cd(homedir); cd figs; mkdir('Spectral_PhCPerm'); cd('Spectral_PhCPerm');

[X,Y]=meshgrid(wtTable.freq{1},params.startTime*1000:(params.limit-201));
figure('Name',['Observed Phase Difference' name],'Position',[100 100 1065 400]); 
ketFig = subplot(131);
surf(Y',X',grpK,'EdgeColor','None'); view(2);
set(gca,'YScale','log'); title('Ketamine')
yticks([0 10 20 30 40 50 60 80 100 200 300 500])
colorbar
% clim = get(gca,'clim');

awakeFig = subplot(132);
surf(Y',X',grpA,'EdgeColor','None'); view(2);
set(gca,'YScale','log'); title('Awake')
yticks([4 7 12 18 30 60 100])
colorbar
clim = get(gca,'clim');
% clim = [clim; get(gca,'clim')];

diffFig = subplot(133);
surf(Y',X',diffK_A,'EdgeColor','None'); view(2);
set(gca,'YScale','log'); title('Observed Diff')
yticks([4 7 12 18 30 60 100])
colorbar
clim = [clim; get(gca,'clim')];

newC = [min(clim(:)) max(clim(:))];

%figure(awakeFig); 
set(awakeFig,'Clim',newC);colorbar;
%figure(ketFig); 
set(ketFig,'Clim',newC);colorbar;
%figure(muscFig); 
set(diffFig,'Clim',newC);colorbar;

h = gcf;
h.Renderer = 'Painters';
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(['Observed Phase Difference' name])
% saveas(gcf, ['Observed Phase Difference' name '.pdf'])

%% Mann Whitney U Test (ranksum)

shiftedA = shiftdim(AwakeAll,1);
shiftedK = shiftdim(KetAll,1);
obs_Esize = zeros(size(grpA));
obs_Cmass = zeros(size(grpA));

% this test works on vectors. so we're doing it pointwise. Awyis.
for row = 1:size(AwakeAll,1)
    for col = 1:size(AwakeAll,2)
        
        % take out the point of interest
        pointA = shiftedA(col,:,row);
        pointK = shiftedK(col,:,row);
        
        % two-tailed test (stats.p(2)) like with the ttest2 used in mag plots, the
        % samples are large enough to use method 'normal approximation'
        stats = mwwtest(pointA,pointK);
        if stats.p(2) <= 0.05
            sigtest = 1;
        else
            sigtest = 0;
        end
        
        % effect size of mwwtest is r = abs(z/sqrt(n1+n2)) / 0.1 is small, 0.3 is
        % medium, 0.5 is large
        Esize = abs(stats.Z/sqrt(grpsizeA+grpsizeK));
        if Esize <= 0.1
            Esize = 1;
        elseif Esize >= 0.5
            Esize = 3;
        else
            Esize = 2;
        end
        
        obs_Esize(row,col) = Esize;
        obs_Cmass(row,col) = sigtest;
        
    end
end

% full matrix
obs_clustermass = nansum(nansum(obs_Cmass));

% layer specific
obs_layer = struct;

for ispec = 1:length(osciName)
    obs_layer.(osciName{ispec}) = obs_Cmass(osciRows{ispec},:);
    
    % % sum clusters (twice to get final value)
    for i = 1:2
        obs_layer.(osciName{ispec}) = nansum(obs_layer.(osciName{ispec}));
    end
end


%% cluster fig 
[X,Y]=meshgrid(wtTable.freq{1},params.startTime*1000:(params.limit-201));
figure('Name',['Observed Effect Size and Clustermass' name],'Position',[-1070 900 1065 400]); 
subplot(121);
surf(Y,X,squeeze(obs_Esize)','EdgeColor','None'); view(2);
set(gca,'YScale','log'); title('Effect size mat')
yticks([0 10 20 30 40 50 60 80 100 200 300 500])
colorbar

subplot(122);
surf(Y,X,squeeze(obs_Cmass)','EdgeColor','None'); view(2);
set(gca,'YScale','log'); title('Clusters where p>0.05')
yticks([0 10 20 30 40 50 60 80 100 200 300 500])

h = gcf;
h.Renderer = 'Painters';
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(['Observed Effect Size and Clustermass' name])
% saveas(gcf, ['Observed Effect Size and Clustermass' name '.pdf'])

close all

%% Permutation Step 3 - do the permute

% set up some containers
mass_clustermass = NaN([1 nperms]);
tic
perm_layer = struct;
for ispec = 1:length(osciName)
    perm_layer.(osciName{ispec}) = NaN([1 nperms]);
end

% echo loop to find observed values
for iperm = 1:nperms
    
    % randomize the order of animal selection
    order = randperm(grpsizeA+grpsizeK);
    
    % Get phase coherence for groups
    permbox1 = nan(54,params.limit,grpsizeA);
    permbox2   = nan(54,params.limit,grpsizeK);
    
    % through each animal
    for iAn = 1: grpsizeA + grpsizeK
        
        if iAn <= grpsizeA
            % pull out the data for the random animal but make sure not to
            % take Muscimol trials
            perm1 = table2cell(wtTable(~contains(wtTable.condition,'Musc')& contains(wtTable.animal,allnames{order(iAn)}),1));
            perm1 =  cellfun(@(x) x(:,1:params.limit),perm1,'UniformOutput',false);
            % set up output cells for transformed data
            transPerm1 = cell(size(perm1));
            
            % take the phase z/abs(z) for each single trial
            for iTrial = 1:length(perm1)
                curTrial = perm1{iTrial};
                transPerm1{iTrial} = curTrial./abs(curTrial);
            end
            
            % for the subject:
            % get the mean of single trials and the absolute for that mean
            curAn = abs(mean(cat(3,transPerm1{:}),3));
            permbox1(:,:,iAn) = curAn;
            
        else
            perm2 = table2cell(wtTable(contains(wtTable.condition,'Anesth')& contains(wtTable.animal,allnames{iAn}),1));
            perm2 =  cellfun(@(x) x(:,1:params.limit),perm2,'UniformOutput',false);
            % set up output cells for transformed data
            transPerm2 = cell(size(perm2));
            
            % take the phase z/abs(z) for each single trial
            for iTrial = 1:length(perm2)
                curTrial = perm2{iTrial};
                transPerm2{iTrial} = curTrial./abs(curTrial);
            end
            
            % get the mean of single trials and the absolute for that mean
            curAn = abs(mean(cat(3,transPerm2{:}),3));
            permbox2(:,:,iAn-grpsizeA) = curAn;
            
        end
    end
    
    % Mann Whitney U Test (ranksum) %%
    
    shifted1 = shiftdim(permbox1,1);
    shifted2 = shiftdim(permbox2,1);
    perm_Cmass = zeros(size(grpA));
    
    % pointwise
    for row = 1:size(permbox1,1)
        for col = 1:size(permbox1,2)
            
            % take out the point of interest
            point1 = shifted1(col,:,row);
            point2 = shifted2(col,:,row);
            
            % two-tailed test (stats.p(2)) 
            stats = mwwtest(point1,point2);
            if stats.p(2) <= 0.05
                sigtest = 1;
            else
                sigtest = NaN;
            end
            perm_Cmass(row,col) = sigtest;
        end
    end
    
     % check cluster mass for 300 ms from tone onset
    per_clustermass = nansum(nansum(perm_Cmass));
    mass_clustermass(iperm) = per_clustermass;
    
    % for layer specific: %%%
    % % pull out clusters
    
    for ispec = 1:length(osciName)
        hold_permlayer = perm_Cmass(osciRows{ispec},:);
        
        % % sum clusters (twice to get final value)
        for i = 1:2
            hold_permlayer = nansum(hold_permlayer);
        end
        perm_layer.(osciName{ispec})(iperm) = hold_permlayer;
    end
    
end
toc

cd(homedir); cd DATA; cd Spectral; mkdir('Spectral_PhCPerm'); cd('Spectral_PhCPerm');
%% Check Significance of full clustermass

% In how many instances is the clustermass of the permutation LESS than
% the observed clustermass
sig_massLess = sum(mass_clustermass<obs_clustermass,2); 
pValL = sig_massLess/nperms;
sig_massLess(sig_massLess <= pthresh) = true; %if it's not more than 50, it's siginificant for that frequency
sig_massLess(sig_massLess > pthresh) = false;

% In how many instances is the clustermass of the permutation MORE than
% the observed clustermass
sig_massMore = sum(mass_clustermass>obs_clustermass,2); %compare how many instances the sum of the permutation is greater than observed
pValM = sig_massMore/nperms;
permMean = mean(mass_clustermass);
permSTD = std(mass_clustermass);
sig_massMore(sig_massMore <= pthresh) = true; %if it's not more than the threshold, it's siginificant for that frequency
sig_massMore(sig_massMore > pthresh) = false;

if pthresh <= 1
    disp('      pthresh is less than 1')
end

figure('Name','Observed cluster vs Permutation'); 
boxplot(mass_clustermass); hold on;

if sig_massLess == 0 && sig_massMore == 0
    
    plot(1,obs_clustermass,'ro','LineWidth',4)
    legend('ns')
    axis 'auto y'
else
    
    plot(1,obs_clustermass,'go','LineWidth',4)
    legend('p<0.007')
    axis 'auto y'
end

h = gcf;
h.Renderer = 'Painters';
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(['Full Permutation' name])
saveas(gcf, ['Full Permutation' name '.pdf'])

save(['Permutation' name ' Full.mat'],'pValL','pValM','permMean','permSTD')

%% Check Significance of layers' clustermass

for ispec = 1:length(osciName)
% In how many instances is the clustermass of the permutation LESS than
% the observed clustermass
sig_massLess = sum(perm_layer.(osciName{ispec})<obs_layer.(osciName{ispec}),2); 
pValL = sig_massLess/nperms;
sig_massLess(sig_massLess <= pthresh) = true; %if it's not more than 50, it's siginificant for that frequency
sig_massLess(sig_massLess > pthresh) = false;

% In how many instances is the clustermass of the permutation MORE than
% the observed clustermass
sig_massMore = sum(perm_layer.(osciName{ispec})>obs_layer.(osciName{ispec}),2); %compare how many instances the sum of the permutation is greater than observed
pValM = sig_massMore/nperms;
permMean = mean(perm_layer.(osciName{ispec}));
permSTD = std(perm_layer.(osciName{ispec}));
sig_massMore(sig_massMore <= pthresh) = true; %if it's not more than 50, it's siginificant for that frequency
sig_massMore(sig_massMore > pthresh) = false;

figure('Name',['Observed cluster vs Permutation' name ' ' osciName{ispec}]); 
boxplot(perm_layer.(osciName{ispec})); hold on;

if sig_massLess == 0 && sig_massMore == 0
    
    plot(1,obs_layer.(osciName{ispec}),'ro','LineWidth',4)
    legend('ns')
    axis 'auto y'
else
    
    plot(1,obs_layer.(osciName{ispec}),'go','LineWidth',4)
    legend('p<0.007')
    axis 'auto y'
end

h = gcf;
h.Renderer = 'Painters';
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(['Permutation' name ' ' osciName{ispec}])
% saveas(gcf, ['Permutation' name ' ' osciName{ispec} '.pdf'])

save(['Permutation' name ' ' osciName{ispec} '.mat'],'pValL','pValM','permMean','permSTD')
% fprintf(['The P value after clustermass permutation for ' layer num2str(rel2BFin) ' at ' osciName{ispec} ' is ' num2str(PMore) ' \n'])

end


