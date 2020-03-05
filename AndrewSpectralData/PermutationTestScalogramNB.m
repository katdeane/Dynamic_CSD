function PermutationTestScalogramNB(layer,rel2BFin)
% Input:    Layer to analyze, (possible input: relative to BF)
%           Needs scalogramsfull.mat from Andrew Curran's wavelet analysis
% Output:   Figures for means and observed difference of awake/ketamine
%           comparison; figures for observed t values, clusters, ttest line
%           output; boxplot and significance of permutation test -> Pictures folder

%% INIT 

% Change directory to your working folder
if exist('D:\MyCode\Dynamic_CSD_Analysis','dir') == 7
    cd('D:\MyCode\Dynamic_CSD_Analysis');
elseif exist('D:\Dynamic_CSD_Analysis','dir') == 7
    cd('D:\Dynamic_CSD_Analysis');
elseif exist('C:\Users\kedea\Documents\Dynamic_CSD_Analysis','dir') == 7
    cd('C:\Users\kedea\Documents\Dynamic_CSD_Analysis')
end

home = pwd;
addpath(genpath(home));

cd AndrewSpectralData; cd Data;

% function in
% rel2BFin = 0; % 0 = BF, 1 = BF+1, etc. (stay close to BF, awake animals...
                                             % don't have a full range)
% layer = 'IVE';

nperms = 1000;
pthresh = nperms*(0.05/7);

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
load('scalograms.mat', 'wtTable')

% varNames = unique(wtTable.layer);
params.startTime = -0.2; % seconds
params.limit = 600;
                                        
% check if layer or full mat needed
if ~strcmp(layer, 'ALL')
    %Pull out layer 
    wt2 = wtTable(contains(wtTable.layer,layer),:);
    wtTable = wt2;
else
    % if condition is all, it needs to still pull out only early sinks
    % because there is no time distinction here (i.e. = VIE == VIL)
    wt2 = wtTable(contains(wtTable.layer,'E'),:);
    wtTable = wt2;
end

%Pull out conditions and limit them to 600ms AT BF!
awake = table2cell(wtTable(contains(wtTable.condition,'Awake')&wtTable.rel2Bf==rel2BFin,1));
awake =  cellfun(@(x) x(:,1:params.limit),awake,'UniformOutput',false);

anest = table2cell(wtTable(contains(wtTable.condition,'Anesth')&wtTable.rel2Bf==rel2BFin,1));
anest =  cellfun(@(x) x(:,1:params.limit),anest,'UniformOutput',false);

musc = table2cell(wtTable(contains(wtTable.condition,'Muscimol')&wtTable.rel2Bf==rel2BFin,1));
musc =  cellfun(@(x) x(:,1:params.limit),musc,'UniformOutput',false);

%Stack the individual animals' data (animal#x54x600)
for ii = 1:length(awake)
    awake2(ii,:,:)=awake{ii};
end
awake2 = abs(awake2);

for ii = 1:length(anest)
    anest2(ii,:,:)=anest{ii};
end
anest2 = abs(anest2);

% for ii = 1:length(musc)
%     musc2(ii,:,:)=musc{ii};
% end
% musc2 = abs(musc2);

if ~strcmp(layer, 'ALL')
    grpsizeA = length(awake); %groups have 1 point per data (1 sink)
    grpsizeK = length(anest);
else
    grpsizeA = length(awake)/7; %groups have 7 points per data (7 sinks)
    grpsizeK = length(anest)/7;
end

%% Degrees of Freedom and t Threshold
% 
% df = grpsizeA+grpsizeK-2;
% 
% if df == 18
%     t_thresh = 2.101; %two tailed: 2.101, one tailed: 1.734
% elseif df == 16
%     t_thresh = 2.120; %two tailed: 2.120, one tailed: 1.746
% else
%     error('You need to change your tvalue threshold! Check this link: http://www.ttable.org/')
% end

%% Permutation Step 1 - Observed Differences

obs1_mean3 = nanmean(anest2,1);
obs1_std3 = nanstd(anest2,0,1);

obs2_mean3 = nanmean(awake2,1);
obs2_std3 = nanstd(awake2,0,1);

obs_difmeans = obs1_mean3 - obs2_mean3;
% absolute value of difference in means! - take note for future KD
% obs_difmeans = abs(obs_difmeans);

%% Permutation Step 2 - t test
%find the t values along all data points for each frequency bin

obs_t = (obs1_mean3 - obs2_mean3)./...
    sqrt((obs1_std3.^2/grpsizeK) + (obs2_std3.^2/grpsizeA));
obs_clusters = abs(obs_t);
% obs_clusters(obs_clusters < t_thresh) = NaN;
% obs_clusters(obs_clusters >= t_thresh) = 1;

% check cluster mass for 300 ms from tone onset
obs_clustermass = nansum(nansum(obs_clusters));

% ttest for line
[H,~,~,STATS] = ttest2(obs1_mean3,obs2_mean3,'vartype','unequal');
obs_tstat = STATS.tstat;
obs_H = H;

%% for layer specific: 
% % pull out clusters

obs_layer = struct;

for ispec = 1:length(osciName)
    obs_layer.(osciName{ispec}) = obs_clusters(:,osciRows{ispec},:);
    
    % % sum clusters (twice to get final value)
    for i = 1:2
        obs_layer.(osciName{ispec}) = nansum(obs_layer.(osciName{ispec}));
    end
end

cd(home); cd AndrewSpectralData; cd Pictures; cd PermutationsNB
%% dif fig
[X,Y]=meshgrid(wtTable.freq{1},params.startTime*1000:(params.limit-201));
ObservedDiff = figure('Name','Observed Difference Values BF','Position',[-1070 500 1065 400]); 
ketFig = subplot(131);
surf(Y,X,squeeze(obs1_mean3)','EdgeColor','None'); view(2);
set(gca,'YScale','log'); title('Ketamine')
yticks([0 10 20 30 40 50 60 80 100 200 300 500])
colorbar
% clim = get(gca,'clim');

awakeFig = subplot(132);
surf(Y,X,squeeze(obs2_mean3)','EdgeColor','None'); view(2);
set(gca,'YScale','log'); title('Awake')
yticks([0 10 20 30 40 50 60 80 100 200 300 500])
colorbar
clim = get(gca,'clim');
% clim = [clim; get(gca,'clim')];

diffFig = subplot(133);
surf(Y,X,squeeze(obs_difmeans)','EdgeColor','None'); view(2);
set(gca,'YScale','log'); title('Observed Diff')
yticks([0 10 20 30 40 50 60 80 100 200 300 500])
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
savefig(['Observed Difference ' layer ' ' num2str(rel2BFin)])
saveas(gcf, ['Observed Difference ' layer ' ' num2str(rel2BFin) '.pdf'])
% print('Observed Difference', '-dpdf') - not saving as vector, leaving for
% now KD

%% t fig 
[X,Y]=meshgrid(wtTable.freq{1},params.startTime*1000:(params.limit-201));
tValues = figure('Name','Observed t Values BF','Position',[-1070 900 1065 400]); 
tFig = subplot(131);
surf(Y,X,squeeze(obs_t)','EdgeColor','None'); view(2);
set(gca,'YScale','log'); title('t Value Mat')
yticks([0 10 20 30 40 50 60 80 100 200 300 500])
colorbar

clusFig = subplot(132);
surf(Y,X,squeeze(obs_clusters)','EdgeColor','None'); view(2);
set(gca,'YScale','log'); title('Clusters where p>0.05')
colormap('winter')
yticks([0 10 20 30 40 50 60 80 100 200 300 500])


tcurveFig = subplot(133);
ttest = obs_tstat;
sigH = ttest;
sigH(obs_H == 0) = NaN;
plot(squeeze(ttest)); hold on
plot(squeeze(sigH), 'LineWidth',2)
legend('t-values','significant')

h = gcf;
h.Renderer = 'Painters';
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(['Observed t and p ' layer ' ' num2str(rel2BFin)])
saveas(gcf, ['Observed t and p ' layer ' ' num2str(rel2BFin) '.pdf'])

%% Permutation Step 3 - do the permute
mass_clustermass = NaN([1 nperms]);
% put the whole group in one container
contAll = vertcat(awake,anest);
perm_layer = struct;

for ispec = 1:length(osciName)
    perm_layer.(osciName{ispec}) = NaN([1 nperms]);
end

for iperm = 1:nperms
    % determine random list order to pull
    order = randperm(length(contAll));
    % pull based on random list order
    for ii = 1:grpsizeA
        awakePerm(ii,:,:)=contAll{order(ii)};
    end
    awakePerm = abs(awakePerm);
    
    for ii = grpsizeA+1:grpsizeK+grpsizeA
        anestPerm(ii-grpsizeA,:,:)=contAll{order(ii)};
    end
    anestPerm = abs(anestPerm);
    % did the permute.
    
    
    per1_mean3 = nanmean(anestPerm,1);
    per1_std3 = nanstd(anestPerm,0,1);
    
    per2_mean3 = nanmean(awakePerm,1);
    per2_std3 = nanstd(awakePerm,0,1);
    
    per_difmeans = per1_mean3 - per2_mean3;
    
    % t test %%%
    per_t = (per1_mean3 - per2_mean3)./...
        sqrt((per1_std3.^2/grpsizeK) + (per2_std3.^2/grpsizeA));
    per_clusters = abs(per_t);
%     per_clusters(per_clusters < t_thresh) = NaN;
%     per_clusters(per_clusters >= t_thresh) = 1;
    
    % check cluster mass for 300 ms from tone onset
    per_clustermass = nansum(nansum(per_clusters));
    mass_clustermass(iperm) = per_clustermass;
    
    %sanity check for permutations (visualization of diff and clustermass)
%     if mod(iperm,100) == 0
%         figure('Position',[-1070 20 900 400]);
%         subplot(121)
%         surf(Y,X,squeeze(per_difmeans)','EdgeColor','None'); view(2);
%         caxis(newC); colorbar
%         set(gca,'YScale','log');
%         yticks([0 10 20 30 40 50 60 80 100 200 300 500])
%         subplot(122)
%         surf(Y,X,squeeze(per_clusters)','EdgeColor','None'); view(2);
%         set(gca,'YScale','log');
%         yticks([0 10 20 30 40 50 60 80 100 200 300 500])
%     end
    
    % for layer specific: %%%
    % % pull out clusters
    
    for ispec = 1:length(osciName)
        hold_permlayer = per_clusters(:,osciRows{ispec},:);
        
        % % sum clusters (twice to get final value)
        for i = 1:2
            hold_permlayer = nansum(hold_permlayer);
        end
        perm_layer.(osciName{ispec})(iperm) = hold_permlayer;
    end
    
end


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
sig_massMore(sig_massMore <= pthresh) = true; %if it's not more than 50, it's siginificant for that frequency
sig_massMore(sig_massMore > pthresh) = false;

PermutSpread = figure('Name',['Observed cluster vs Permutation' layer]); 
boxplot(mass_clustermass); hold on;

if sig_massLess == 0 && sig_massMore == 0
    
    plot(1,obs_clustermass,'ro','LineWidth',4)
    legend('ns')
else
    
    plot(1,obs_clustermass,'go','LineWidth',4)
    legend('p<0.007')
end

h = gcf;
h.Renderer = 'Painters';
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(['Full Permutation ' layer ' ' num2str(rel2BFin)])
saveas(gcf, ['Full Permutation ' layer '.pdf'])

save(['Permutation ' layer ' ' num2str(rel2BFin) '.mat'],'pValL','pValM')

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
sig_massMore(sig_massMore <= pthresh) = true; %if it's not more than 50, it's siginificant for that frequency
sig_massMore(sig_massMore > pthresh) = false;

PermutSpreadLayer = figure('Name',['Observed cluster vs Permutation ' layer ' ' osciName{ispec}]); 
boxplot(perm_layer.(osciName{ispec})); hold on;

if sig_massLess == 0 && sig_massMore == 0
    
    plot(1,obs_layer.(osciName{ispec}),'ro','LineWidth',4)
    legend('ns')
else
    
    plot(1,obs_layer.(osciName{ispec}),'go','LineWidth',4)
    legend('p<0.007')
end

h = gcf;
h.Renderer = 'Painters';
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(['Permutation ' layer ' ' osciName{ispec} ' ' num2str(rel2BFin)])
saveas(gcf, ['Permutation ' layer ' ' osciName{ispec} ' ' num2str(rel2BFin) '.pdf'])

save(['Permutation ' layer ' ' osciName{ispec} ' ' num2str(rel2BFin) '.mat'],'pValL','pValM')
fprintf(['The P value after clustermass permutation for ' layer num2str(rel2BFin) ' at ' osciName{ispec} ' is ' num2str(pValM) ' \n'])

end
