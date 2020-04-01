function cohensDScalogram(layer, rel2BFin, homedir)
% Input:    Layer to analyze, (possible input: relative to BF)
%           Needs scalogramsfull.mat from Andrew Curran's wavelet analysis
% Output:   Figures for means and observed difference of awake/ketamine
%           comparison; figures for observed t values, clusters, ttest line
%           output; boxplot and significance of permutation test -> Pictures folder
%   2019-06-12 AC: added some cohen's d stuff to output a discritized and
%   normal version of plot, changed from PermutationTestScalogram
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
if ~exist('rel2BFin','var')
    rel2BFin = 0; % run the BF if not specified
end
if ~exist('layer','var')
    layer = 'IVE'; % run the granular layer if not specified
end 

cd (homedir),cd DATA;

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
awake = cellfun(@(x) x(:,1:params.limit),awake,'UniformOutput',false);

anest = table2cell(wtTable(contains(wtTable.condition,'Anesth')&wtTable.rel2Bf==rel2BFin,1));
anest = cellfun(@(x) x(:,1:params.limit),anest,'UniformOutput',false);

% musc = table2cell(wtTable(contains(wtTable.condition,'Muscimol')&wtTable.rel2Bf==rel2BFin,1));
% musc = cellfun(@(x) x(:,1:params.limit),musc,'UniformOutput',false);

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


%% Cohen's D

obs1_mean3 = nanmean(anest2,1);
obs1_std3 = nanstd(anest2,0,1);

obs2_mean3 = nanmean(awake2,1);
obs2_std3 = nanstd(awake2,0,1);

%find the t values along all data points for each frequency bin
S = sqrt((((grpsizeK - 1).*obs1_std3.^2)+((grpsizeA-1).*obs2_std3.^2))/(grpsizeK + grpsizeA - 2));
cohensD = (obs1_mean3-obs2_mean3)./S;
cohensDLine = nanmean(squeeze(abs(cohensD)),1);

cd(homedir); cd figs; mkdir('Spectral_MagPerm'); cd('Spectral_MagPerm');

%% Cohen's D fig 
newD = cohensD;
newD(newD<0.5   & newD>-0.5)    = 0; % small or very small
newD(newD>0.5   & newD<0.8)     = 0.5; % medium
newD(newD<-0.5  & newD>-0.8)    = 0.5; % medium
newD(newD>0.8   & newD<1.2)     = 1; % large
newD(newD<-0.8  & newD>-1.2)    = 1; % large
newD(newD>1.2)                  = 1.5; % very large
newD(newD<-1.2)                 = 1.5; % very large

[X,Y]=meshgrid(wtTable.freq{1},params.startTime*1000:(params.limit-201));
figure('Name','Cohens D','Position',[-1070 300 1065 1000]);
subplot(221);
surf(Y,X,squeeze(newD)','EdgeColor','None'); view(2);
set(gca,'YScale','log'); title('Cohens D Surface')
yticks([0 10 20 30 40 50 60 80 100 200 300 500])
colorbar

subplot(222);
contour(Y,X,squeeze(newD)')
set(gca,'YScale','log'); title('Cohens D Contour')
yticks([0 10 20 30 40 50 60 80 100 200 300 500])
colorbar

subplot(223);
plot(cohensDLine);
title('Cohens D Line')
yticks([-0.5 -0.2 0 0.2 0.5 0.8 1 1.2])

h = gcf;
h.Renderer = 'Painters';
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(['Cohens D Observed ' layer ' ' num2str(rel2BFin)])
saveas(gcf, ['Cohens D Observed ' layer ' ' num2str(rel2BFin) '.pdf'])
close(h)
