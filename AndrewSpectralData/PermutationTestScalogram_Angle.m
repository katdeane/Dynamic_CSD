function PermutationTestScalogram_Angle(layer,rel2BFin)
% Input:    Layer to analyze, (possible input: relative to BF)
%           Needs scalogramsfull.mat from Andrew Curran's wavelet analysis
% Output:   Figures for means and observed difference of awake/ketamine
%           comparison; figures for observed t values, clusters, ttest line
%           output; boxplot and significance of permutation test -> Pictures folder

%% INIT 

cd('D:\MyCode\Dynamic_CSD_Analysis');
home = pwd;
addpath(genpath(home));

cd AndrewSpectralData; cd Data;

% rel2BFin = 0; % 0 = BF, 1 = BF+1, etc. (stay close to BF, awake animals...
                                        % don't have a full range)
% layer = 'IVE';
nperms = 1000;
pthresh = nperms*0.05;

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

%Stack the individual animals' data (animal#x54x600)
for ii = 1:length(awake)
    awake_complex(ii,:,:)=awake{ii};
end
awake_angle = angle(awake_complex);
awake_complex = shiftdim(awake_complex, 1);

for ii = 1:length(anest)
    anest_complex(ii,:,:)=anest{ii};
end
anest_angle = angle(anest_complex);
anest_complex = shiftdim(anest_complex, 1);


if ~strcmp(layer, 'ALL')
    grpsizeA = length(awake); %groups have 1 point per data (1 sink)
    grpsizeK = length(anest);
else
    grpsizeA = length(awake)/7; %groups have 7 points per data (7 sinks)
    grpsizeK = length(anest)/7;
end


%% Permutation Step 1 - Observed Differences

obs1_mean3 = circ_mean(anest_angle);
obs2_mean3 = circ_mean(awake_angle);


%% Permutation Step 2 - t test -- F test in this case
%find the t values along all data points for each frequency bin

[obs_clusters, obs_Fmatrix] = matrix_circ_wwtest(anest_complex,awake_complex);

obs_clusters(obs_clusters > 0.05) = NaN;
obs_clusters(obs_clusters <= 0.05) = 1;

% check cluster mass for 300 ms from tone onset
obs_clustermass = nansum(nansum(obs_clusters));


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

cd(home); cd AndrewSpectralData; cd Pictures; cd Permutations
%% dif fig
[X,Y]=meshgrid(wtTable.freq{1},params.startTime*1000:(params.limit-201));
Observed = figure('Name','AP Observed Values BF','Position',[-1070 500 1065 400]); 
ketFig = subplot(121);
contourf(Y,X,squeeze(obs1_mean3)','EdgeColor','None'); view(2);
set(gca,'YScale','log'); title('Ketamine')
yticks([0 10 20 30 40 50 60 80 100 200 300 500])
colorbar
clim = get(gca,'clim');

awakeFig = subplot(122);
contourf(Y,X,squeeze(obs2_mean3)','EdgeColor','None'); view(2);
set(gca,'YScale','log'); title('Awake')
yticks([0 10 20 30 40 50 60 80 100 200 300 500])
colorbar
clim = [clim; get(gca,'clim')];


newC = [min(clim(:)) max(clim(:))];
%figure(awakeFig); 
set(awakeFig,'Clim',newC);colorbar;
%figure(ketFig); 
set(ketFig,'Clim',newC);colorbar;

h = gcf;
h.Renderer = 'Painters';
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(['AP Observed ' layer ' ' num2str(rel2BFin)])
% saveas(gcf, ['AP Observed ' layer ' ' num2str(rel2BFin) '.pdf'])


 %% ww fig
[X,Y]=meshgrid(wtTable.freq{1},params.startTime*1000:(params.limit-201));
FValues = figure('Name','AP Observed F Values BF','Position',[-1070 900 1065 400]); 
fFig = subplot(121);
surf(Y,X,squeeze(obs_Fmatrix)','EdgeColor','None'); view(2);
set(gca,'YScale','log'); title('F Value Mat')
yticks([0 10 20 30 40 50 60 80 100 200 300 500])
colorbar

clusFig = subplot(122);
surf(Y,X,squeeze(obs_clusters)','EdgeColor','None'); view(2);
set(gca,'YScale','log'); title('Clusters where p>0.05')
colormap('winter')
yticks([0 10 20 30 40 50 60 80 100 200 300 500])


h = gcf;
h.Renderer = 'Painters';
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(['AP Observed F and p ' layer ' ' num2str(rel2BFin)])
% saveas(gcf, ['AP Observed F and p ' layer ' ' num2str(rel2BFin) '.pdf'])

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
        awakePerm(:,:,ii)=contAll{order(ii)};
    end
    
    for ii = grpsizeA+1:grpsizeK+grpsizeA
        anestPerm(:,:,ii-grpsizeA)=contAll{order(ii)};
    end

    % did the permute.
    
    
    % F test %%%
    [per_clusters, ~] = matrix_circ_wwtest(anestPerm,awakePerm);
    
    per_clusters(per_clusters > 0.05) = NaN;
    per_clusters(per_clusters <= 0.05) = 1;
    
    % check cluster mass for 300 ms from tone onset
    per_clustermass = nansum(nansum(per_clusters));
    mass_clustermass(iperm) = per_clustermass;
    
    
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
sig_massLess = sum(mass_clustermass>obs_clustermass,2); 
sig_massLess(sig_massLess <= pthresh) = true; %if it's not more than 50, it's siginificant for that frequency
sig_massLess(sig_massLess > pthresh) = false;

% In how many instances is the clustermass of the permutation MORE than
% the observed clustermass
sig_massMore = sum(mass_clustermass<obs_clustermass,2); %compare how many instances the sum of the permutation is greater than observed
sig_massMore(sig_massMore <= pthresh) = true; %if it's not more than 50, it's siginificant for that frequency
sig_massMore(sig_massMore > pthresh) = false;

PermutSpread = figure('Name',['AP Observed cluster vs Permutation' layer]); 
boxplot(mass_clustermass); hold on;

if sig_massLess == 0 && sig_massMore == 0
    
    plot(1,obs_clustermass,'ro','LineWidth',4)
    legend('ns')
else
    
    plot(1,obs_clustermass,'go','LineWidth',4)
    legend('p<0.05')
end

h = gcf;
h.Renderer = 'Painters';
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(['AP Full Permutation ' layer ' ' num2str(rel2BFin)])
% saveas(gcf, ['AP Full Permutation ' layer '.pdf'])

%% Check Significance of layers' clustermass

for ispec = 1:length(osciName)
% In how many instances is the clustermass of the permutation LESS than
% the observed clustermass
sig_massLess = sum(perm_layer.(osciName{ispec})>obs_layer.(osciName{ispec}),2); 
sig_massLess(sig_massLess <= pthresh) = true; %if it's not more than 50, it's siginificant for that frequency
sig_massLess(sig_massLess > pthresh) = false;

% In how many instances is the clustermass of the permutation MORE than
% the observed clustermass
sig_massMore = sum(perm_layer.(osciName{ispec})<obs_layer.(osciName{ispec}),2); %compare how many instances the sum of the permutation is greater than observed
sig_massMore(sig_massMore <= pthresh) = true; %if it's not more than 50, it's siginificant for that frequency
sig_massMore(sig_massMore > pthresh) = false;

PermutSpreadLayer = figure('Name',['AP Observed cluster vs Permutation ' layer ' ' osciName{ispec}]); 
boxplot(perm_layer.(osciName{ispec})); hold on;

if sig_massLess == 0 && sig_massMore == 0
    
    plot(1,obs_layer.(osciName{ispec}),'ro','LineWidth',4)
    legend('ns')
    axis 'auto y'
else
    
    plot(1,obs_layer.(osciName{ispec}),'go','LineWidth',4)
    legend('p<0.05')
    axis 'auto y'
end

h = gcf;
h.Renderer = 'Painters';
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(['AP Permutation ' layer ' ' osciName{ispec} ' ' num2str(rel2BFin)])
% saveas(gcf, ['AP Permutation ' layer ' ' osciName{ispec} ' ' num2str(rel2BFin) '.pdf'])

end












