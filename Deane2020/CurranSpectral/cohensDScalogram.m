function cohensDScalogram(layer, rel2BFin)
% Input:    Layer to analyze, (possible input: relative to BF)
%           Needs scalogramsfull.mat from Andrew Curran's wavelet analysis
% Output:   Figures for means and observed difference of awake/ketamine
%           comparison; figures for observed t values, clusters, ttest line
%           output; boxplot and significance of permutation test -> Pictures folder
%   2019-06-12 AC: added some cohen's d stuff to output a discritized and
%   normal version of plot, changed from PermutationTestScalogram
%% INIT 

cd('D:\MyCode\Dynamic_CSD_Analysis');
home = pwd;
addpath(genpath(home));

cd AndrewSpectralData; cd Data;

% rel2BFin = 0; % 0 = BF, 1 = BF+1, etc. (stay close to BF, awake animals...
                                        % don't have a full range)
% layer = 'IVE';

% frequencies can be found in wtTable.freq{1} to clarify the following
% rows choices; actual intended rows commented
% theta = (49:54);        %(4:7);
% alpha = (44:48);        %(8:12);
% beta_low = (39:43);     %(13:18);
% beta_high = (34:38);    %(19:30);
% gamma_low = (26:33);    %(31:60);
% gamma_high = (19:25);   %(61:100);

thetaLine = ones([1 400])*7;
alphaLine = ones([1 400])*12;
beta_lowLine = ones([1 400])*18;
beta_highLine = ones([1 400])*30;
gamma_lowLine = ones([1 400])*60;
gamma_highLine = ones([1 400])*100;

%% Load in and seperate Data
load('scalogramsfull.mat', 'wtTable')

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

% musc = table2cell(wtTable(contains(wtTable.condition,'Muscimol')&wtTable.rel2Bf==rel2BFin,1));
% musc =  cellfun(@(x) x(:,1:params.limit),musc,'UniformOutput',false);

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

cd(home); cd AndrewSpectralData; cd Pictures; cd Permutations

%% Cohen's D fig 
newD = cohensD;
newD(newD<0.5   & newD>-0.5)    = 0;
newD(newD>0.5   & newD<0.8)     = 0.5;
newD(newD<-0.5  & newD>-0.8)    = 0.5;
newD(newD>0.8   & newD<1.2)     = 1;
newD(newD<-0.8  & newD>-1.2)    = 1;
newD(newD>1.2)                  = 1.5;
newD(newD<-1.2)                 = 1.5;

[X,Y]=meshgrid(wtTable.freq{1},params.startTime*1000:(params.limit-201));
CohensDFig = figure('Name','Cohens D','Position',[-1070 300 1065 1000]);
surfCDFig = subplot(221);
surf(Y,X,squeeze(newD)','EdgeColor','None'); view(2);
set(gca,'YScale','log'); title('Cohens D Surface')
yticks([0 10 20 30 40 50 60 80 100 200 300 500])
colorbar
% hold on;
% plot(thetaLine,'m')
% plot(alphaLine, 'm')
% plot(beta_lowLine, 'm')
% plot(beta_highLine, 'm')
% plot(gamma_lowLine, 'm')
% plot(gamma_highLine, 'm')
% hold off;

contCDFig = subplot(222);
contour(Y,X,squeeze(newD)')
set(gca,'YScale','log'); title('Cohens D Contour')
yticks([0 10 20 30 40 50 60 80 100 200 300 500])
colorbar
% hold on;
% plot(thetaLine,'m')
% plot(alphaLine, 'm')
% plot(beta_lowLine, 'm')
% plot(beta_highLine, 'm')
% plot(gamma_lowLine, 'm')
% plot(gamma_highLine, 'm')
% hold off;

LineCDFig = subplot(223);
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
