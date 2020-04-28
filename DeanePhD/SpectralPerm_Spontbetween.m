function SpectralPerm_Spontbetween(layer,homedir)

% Input:    Layer to analyze and home directory,
%           Loads wavelet transform (WT) mat file from Data/Spectral
% Output:   Figures for means and observed difference of groups KIT vs KIC,
%           KIT vs KIV, and KIC vs KIV comparison and
%           figures for observed t values, and clusters in figs/Crypt_MagPerm_SP;
%           boxplot and significance of permutation test in figs/Crypt_MagPerm_SP

% note:     Permutation done between groups on each measurement type

%% standard operations

warning('OFF');
dbstop if error

% Change directory to your working folder
if ~exist('homedir','var')
    if exist('D:\MyCode\Dynamic_CSD','dir') == 7
        cd('D:\MyCode\Dynamic_CSD');
    elseif exist('C:\Users\kedea\Documents\Work Stuff\Dynamic_CSD','dir') == 7
        cd('C:\Users\kedea\Documents\Work Stuff\Dynamic_CSD')
    else
        error('Please add your directory to the list of detected directories')
    end
    homedir = pwd;
    addpath(genpath(homedir));
end
if ~exist('layer','var')
    layer = 'IV'; % run the granular layer if not specified
end

cd (homedir),cd DATA;

nperms = 1000;
pthresh = nperms*(0.05); %if bonferroni wanted, *(0.05/7)

% frequencies can be found in WT_SP.freq{1} to clarify the following
% rows choices; actual intended rows commented
theta = (49:54);        %(4:7);
alpha = (44:48);        %(8:12);
beta_low = (39:43);     %(13:18);
beta_high = (34:38);    %(19:30);
gamma_low = (26:33);    %(31:60);
gamma_high = (19:25);   %(61:100);

osciName    = {'theta' 'alpha' 'beta_low' 'beta_high' 'gamma_low' 'gamma_high'};
osciRows    = {theta alpha beta_low beta_high gamma_low gamma_high};
% sp pre1 = before first laser, pre2 = before second laser, end = after AMs
measurement = {'spPre1_1','spPost1_1','spPre2_1','spPost2_1','spEnd_1'};

%% Load in and seperate Data
disp('loading mat file WT of spontaneous recordings')
tic
load('WT_SP.mat', 'WT_SP')
WT_SP = WT_SP(startsWith(WT_SP.layer,layer) & endsWith(WT_SP.layer,layer),:);
toc

% varNames = unique(wtTable.layer);
params.startTime = -0.2; % seconds
params.limit = 1300;

for iMea = 1:length(measurement)
    %Pull out conditions and limit them to 1300ms
    KIC = table2cell(WT_SP(contains(WT_SP.group,'KIC') & ...
        startsWith(WT_SP.measurement,measurement{iMea}),1));
    KIC =  cellfun(@(x) x(:,1:params.limit),KIC,'UniformOutput',false);
    
    KIT = table2cell(WT_SP(contains(WT_SP.group,'KIT') & ...
        startsWith(WT_SP.measurement,measurement{iMea}),1));
    KIT =  cellfun(@(x) x(:,1:params.limit),KIT,'UniformOutput',false);
    
    KIV = table2cell(WT_SP(contains(WT_SP.group,'KIV') & ...
        startsWith(WT_SP.measurement,measurement{iMea}),1));
    KIV =  cellfun(@(x) x(:,1:params.limit),KIV,'UniformOutput',false);
    
    %Stack the individual animals' data (animal#x54x600)
    KIC2 = nan(size(KIC,1),size(KIC{1,1},1),size(KIC{1,1},2));
    KIT2 = nan(size(KIT,1),size(KIT{1,1},1),size(KIT{1,1},2));
    KIV2 = nan(size(KIV,1),size(KIV{1,1},1),size(KIV{1,1},2));
    
    for ii = 1:length(KIC)
        KIC2(ii,:,:)= KIC{ii};
    end
    KIC2 = abs(KIC2);
    
    for ii = 1:length(KIT)
        KIT2(ii,:,:)= KIT{ii};
    end
    KIT2 = abs(KIT2);
    
    for ii = 1:length(KIV)
        KIV2(ii,:,:)= KIV{ii};
    end
    KIV2 = abs(KIV2);
    
    grpsizeC = length(KIC);
    grpsizeT = length(KIT);
    grpsizeV = length(KIV);
    
    
    %% Degrees of Freedom and t Threshold
    
    df_TvC = grpsizeC+grpsizeT-2;
    df_TvV = grpsizeV+grpsizeT-2;
    df_CvV = grpsizeV+grpsizeC-2;
    
    Tchart = [12.71,4.303,3.182,2.776,2.571,2.447,2.365,2.306,2.262,2.228,...
        2.201,2.179,2.160,2.145,2.131,2.120,2.110,2.101,2.093,2.086];
    if df_TvC > 20 || df_TvV > 20 || df_CvV > 20
        error('Tchart only goes up to df == 20, to expand check this link: http://www.ttable.org/')
    end
    
    Tthresh_TvC = Tchart(df_TvC);
    Tthresh_TvV = Tchart(df_TvV);
    Tthresh_CvV = Tchart(df_CvV);
    
    %% Permutation Step 1 - Observed Differences
    
    obsT_mean3 = nanmean(KIT2,1);
    obsT_std3 = nanstd(KIT2,0,1);
    
    obsC_mean3 = nanmean(KIC2,1);
    obsC_std3 = nanstd(KIC2,0,1);
    
    obsV_mean3 = nanmean(KIV2,1);
    obsV_std3 = nanstd(KIV2,0,1);
    
    obs_difmeans_TvC = obsT_mean3 - obsC_mean3;
    obs_difmeans_TvV = obsT_mean3 - obsV_mean3;
    obs_difmeans_CvV = obsC_mean3 - obsV_mean3;
    
    %% Permutation Step 2 - t test
    %find the t values along all data points for each frequency bin
    % Treated vs control
    obs_t_TvC = (obsT_mean3 - obsC_mean3)./...
        sqrt((obsT_std3.^2/grpsizeT) + (obsC_std3.^2/grpsizeC));
    obs_clusters_TvC = abs(obs_t_TvC);
    obs_clusters_TvC(obs_clusters_TvC < Tthresh_TvC) = NaN;
    obs_clusters_TvC(obs_clusters_TvC >= Tthresh_TvC) = 1;
    
    % check cluster mass KD: this will need to be improved to include a shorter
    % time window for detection (after each stimulus click)
    obs_clustermass_TvC = nansum(nansum(obs_clusters_TvC));
    
    % Treated vs virus control
    obs_t_TvV = (obsT_mean3 - obsV_mean3)./...
        sqrt((obsT_std3.^2/grpsizeT) + (obsV_std3.^2/grpsizeV));
    obs_clusters_TvV = abs(obs_t_TvV);
    obs_clusters_TvV(obs_clusters_TvV < Tthresh_TvV) = NaN;
    obs_clusters_TvV(obs_clusters_TvV >= Tthresh_TvV) = 1;
    
    % check cluster mass KD: this will need to be improved to include a shorter
    % time window for detection (after each stimulus click)
    obs_clustermass_TvV = nansum(nansum(obs_clusters_TvV));
    
    % control vs virus control
    obs_t_CvV = (obsC_mean3 - obsV_mean3)./...
        sqrt((obsC_std3.^2/grpsizeC) + (obsV_std3.^2/grpsizeV));
    obs_clusters_CvV = abs(obs_t_CvV);
    obs_clusters_CvV(obs_clusters_CvV < Tthresh_CvV) = NaN;
    obs_clusters_CvV(obs_clusters_CvV >= Tthresh_CvV) = 1;
    
    % check cluster mass KD: this will need to be improved to include a shorter
    % time window for detection (after each stimulus click)
    obs_clustermass_CvV = nansum(nansum(obs_clusters_CvV));
    
    
    %% Oscillation frequency bands:
    % % pull out clusters
    
    obs_layer_TvC = struct;
    for iOsc = 1:length(osciName)
        obs_layer_TvC.(osciName{iOsc}) = obs_clusters_TvC(:,osciRows{iOsc},:);
        
        % % sum clusters (twice to get final value)
        for i = 1:2
            obs_layer_TvC.(osciName{iOsc}) = nansum(obs_layer_TvC.(osciName{iOsc}));
        end
    end
    
    obs_layer_TvV = struct;
    for iOsc = 1:length(osciName)
        obs_layer_TvV.(osciName{iOsc}) = obs_clusters_TvV(:,osciRows{iOsc},:);
        
        % % sum clusters (twice to get final value)
        for i = 1:2
            obs_layer_TvV.(osciName{iOsc}) = nansum(obs_layer_TvV.(osciName{iOsc}));
        end
    end
    
    obs_layer_CvV = struct;
    for iOsc = 1:length(osciName)
        obs_layer_CvV.(osciName{iOsc}) = obs_clusters_CvV(:,osciRows{iOsc},:);
        
        % % sum clusters (twice to get final value)
        for i = 1:2
            obs_layer_CvV.(osciName{iOsc}) = nansum(obs_layer_CvV.(osciName{iOsc}));
        end
    end
    
    cd(homedir); cd figs; mkdir('Crypt_MagPerm_SP'); cd('Crypt_MagPerm_SP');
    %% dif figs
    [X,Y]=meshgrid(WT_SP.freq{1}(19:54),params.startTime*1000:(params.limit-201));
    figure('Name','Observed Spectral Power and Differences',...
        'units','normalized','outerposition',[0 0 1 1]);
    % observed values
    kitFig = subplot(231);
    surf(Y,X,squeeze(obsT_mean3(:,19:54,:))','EdgeColor','None'); view(2);
    set(gca,'YScale','log'); title('KIT')
    yticks([0 10 20 30 40 50 60 80 100])
    clim = get(gca,'clim');
    
    kicFig = subplot(232);
    surf(Y,X,squeeze(obsC_mean3(:,19:54,:))','EdgeColor','None'); view(2);
    set(gca,'YScale','log'); title('KIC')
    yticks([0 10 20 30 40 50 60 80 100])
    clim = [clim; get(gca,'clim')]; %#ok<*AGROW>
    
    kivFig = subplot(233);
    surf(Y,X,squeeze(obsV_mean3(:,19:54,:))','EdgeColor','None'); view(2);
    set(gca,'YScale','log'); title('KIV')
    yticks([0 10 20 30 40 50 60 80 100])
    clim = [clim; get(gca,'clim')];
    
    % observed differences
    tvcFig = subplot(234);
    surf(Y,X,squeeze(obs_difmeans_TvC(:,19:54,:))','EdgeColor','None'); view(2);
    set(gca,'YScale','log'); title('KIT-KIC')
    yticks([0 10 20 30 40 50 60 80 100 200 300 500])
    clim = [clim; get(gca,'clim')];
    
    tvvFig = subplot(235);
    surf(Y,X,squeeze(obs_difmeans_TvV(:,19:54,:))','EdgeColor','None'); view(2);
    set(gca,'YScale','log'); title('KIT-KIV')
    yticks([0 10 20 30 40 50 60 80 100 200 300 500])
    clim = [clim; get(gca,'clim')];
    
    cvvFig = subplot(236);
    surf(Y,X,squeeze(obs_difmeans_CvV(:,19:54,:))','EdgeColor','None'); view(2);
    set(gca,'YScale','log'); title('KIC-KIV')
    yticks([0 10 20 30 40 50 60 80 100 200 300 500])
    clim = [clim; get(gca,'clim')];
    
    newC = [min(clim(:)) max(clim(:))];
    
    set(kicFig,'Clim',newC);colorbar;
    set(kitFig,'Clim',newC);colorbar;
    set(kivFig,'Clim',newC);colorbar;
    set(tvcFig,'Clim',newC);colorbar;
    set(tvvFig,'Clim',newC);colorbar;
    set(cvvFig,'Clim',newC);colorbar;
    
    h = gcf;
    h.Renderer = 'Painters';
    set(h, 'PaperType', 'A4');
    set(h, 'PaperOrientation', 'landscape');
    set(h, 'PaperUnits', 'centimeters');
    h_title = ['Spont between Observed Spectral Power of ' layer ' ' measurement{iMea}];
    sgtitle(h_title)
    savefig(h_title)
    saveas(gcf, [h_title '.png'])
    close(h)
    
    %% t fig
    [X,Y]=meshgrid(WT_SP.freq{1}(19:54),params.startTime*1000:(params.limit-201));
    figure('Name','Observed t Values BF', ...
        'units','normalized','outerposition',[0 0 1 1]);
    colormap('winter')
    
    subplot(231);
    surf(Y,X,squeeze(obs_t_TvC(:,19:54,:))','EdgeColor','None'); view(2);
    set(gca,'YScale','log'); title('t Value Mat TvC')
    yticks([0 10 20 30 40 50 60 80 100])
    colorbar
    
    subplot(232);
    surf(Y,X,squeeze(obs_t_TvV(:,19:54,:))','EdgeColor','None'); view(2);
    set(gca,'YScale','log'); title('t Value Mat TvV')
    yticks([0 10 20 30 40 50 60 80 100])
    colorbar
    
    subplot(233);
    surf(Y,X,squeeze(obs_t_CvV(:,19:54,:))','EdgeColor','None'); view(2);
    set(gca,'YScale','log'); title('t Value Mat CvV')
    yticks([0 10 20 30 40 50 60 80 100])
    colorbar
    
    subplot(234);
    surf(Y,X,squeeze(obs_clusters_TvC(:,19:54,:))','EdgeColor','None'); view(2);
    set(gca,'YScale','log'); title('Clusters where p<0.05 TvC')
    yticks([0 10 20 30 40 50 60 80 100])
    colorbar
    
    subplot(235);
    surf(Y,X,squeeze(obs_clusters_TvV(:,19:54,:))','EdgeColor','None'); view(2);
    set(gca,'YScale','log'); title('Clusters where p<0.05 TvV')
    yticks([0 10 20 30 40 50 60 80 100])
    colorbar
    
    subplot(236);
    surf(Y,X,squeeze(obs_clusters_CvV(:,19:54,:))','EdgeColor','None'); view(2);
    set(gca,'YScale','log'); title('Clusters where p<0.05 CvV')
    yticks([0 10 20 30 40 50 60 80 100])
    colorbar
    
    h = gcf;
    h.Renderer = 'Painters';
    set(h, 'PaperType', 'A4');
    set(h, 'PaperOrientation', 'landscape');
    set(h, 'PaperUnits', 'centimeters');
    h_title = ['Observed t and p ' layer ' ' measurement{iMea} ];
    sgtitle(h_title)
    savefig(h_title)
    saveas(gcf, [h_title '.png'])
    close(h)
    
    %% Permutation Step 3 - do the permute
    % preallocate
    mass_clustermass_TvC = NaN([1 nperms]);
    mass_clustermass_TvV = NaN([1 nperms]);
    mass_clustermass_CvV = NaN([1 nperms]);
    
    perm_layer_TvC = struct;
    for iOsc = 1:length(osciName)
        perm_layer_TvC.(osciName{iOsc}) = NaN([1 nperms]);
    end
    perm_layer_TvV = struct;
    for iOsc = 1:length(osciName)
        perm_layer_TvV.(osciName{iOsc}) = NaN([1 nperms]);
    end
    perm_layer_CvV = struct;
    for iOsc = 1:length(osciName)
        perm_layer_CvV.(osciName{iOsc}) = NaN([1 nperms]);
    end
    
    % generate the permuations and take clustermass values
    for iComp = 1:3
        % put the whole group in one container
        if iComp == 1
            contAll = vertcat(KIT,KIC);
            grp1 = grpsizeT;
            grp2 = grpsizeC;
        elseif iComp == 2
            contAll = vertcat(KIT,KIV);
            grp1 = grpsizeT;
            grp2 = grpsizeV;
        elseif iComp == 3
            contAll = vertcat(KIC,KIV);
            grp1 = grpsizeC;
            grp2 = grpsizeV;
        end
        
        for iperm = 1:nperms
            % determine random list order to pull
            order = randperm(length(contAll));
            
            % pull based on random list order
            G1Perm = nan(grp1,size(contAll{1,1},1),size(contAll{1,1},2));
            G2Perm = nan(grp2,size(contAll{1,1},1),size(contAll{1,1},2));
            for ii = 1:grp2
                G2Perm(ii,:,:)=contAll{order(ii)};
            end
            G2Perm = abs(G2Perm);
            
            for ii = grp2+1:grp1+grp2
                G1Perm(ii-grp2,:,:)=contAll{order(ii)};
            end
            G1Perm = abs(G1Perm);
            
            % did the permute.
            
            per1_mean3 = nanmean(G1Perm,1);
            per1_std3 = nanstd(G1Perm,0,1);
            
            per2_mean3 = nanmean(G2Perm,1);
            per2_std3 = nanstd(G2Perm,0,1);
            
            % t test %%%
            per_t = (per1_mean3 - per2_mean3)./...
                sqrt((per1_std3.^2/grpsizeT) + (per2_std3.^2/grpsizeC));
            per_clusters = abs(per_t);
            per_clusters(per_clusters < Tthresh_TvC) = NaN;
            per_clusters(per_clusters >= Tthresh_TvC) = 1;
            
            % check cluster mass for 300 ms from tone onset
            per_clustermass = nansum(nansum(per_clusters));
            
            if iComp == 1
                mass_clustermass_TvC(iperm) = per_clustermass;
                for iOsc = 1:length(osciName)
                    hold_permlayer = per_clusters(:,osciRows{iOsc},:);
                    hold_permlayer = nansum(nansum(hold_permlayer));
                    perm_layer_TvC.(osciName{iOsc})(iperm) = hold_permlayer;
                end
            elseif iComp == 2
                mass_clustermass_TvV(iperm) = per_clustermass;
                for iOsc = 1:length(osciName)
                    hold_permlayer = per_clusters(:,osciRows{iOsc},:);
                    hold_permlayer = nansum(nansum(hold_permlayer));
                    perm_layer_TvV.(osciName{iOsc})(iperm) = hold_permlayer;
                end
            elseif iComp == 3
                mass_clustermass_CvV(iperm) = per_clustermass;
                for iOsc = 1:length(osciName)
                    hold_permlayer = per_clusters(:,osciRows{iOsc},:);
                    hold_permlayer = nansum(nansum(hold_permlayer));
                    perm_layer_CvV.(osciName{iOsc})(iperm) = hold_permlayer;
                end
            end
        end
    end
    
    cd(homedir); cd DATA; cd Spectral; mkdir('Crypt_MagPerm_SP'); cd('Crypt_MagPerm_SP');
    %% Check Significance of full clustermass
    
    % In how many instances is the clustermass of the permutation MORE than
    % the observed clustermass
    sig_mass_TvC = sum(mass_clustermass_TvC>obs_clustermass_TvC,2); %compare how many instances the sum of the permutation is greater than observed
    pVal_TvC = sig_mass_TvC/nperms;
    permMean_TvC = mean(mass_clustermass_TvC);
    permSTD_TvC = std(mass_clustermass_TvC);
    sig_mass_TvC(sig_mass_TvC <= pthresh) = true; %if it's not more than 50, it's siginificant for that frequency
    sig_mass_TvC(sig_mass_TvC > pthresh) = false;
    
    figure('Name',['Spont between Observed cluster vs Permutation ' layer ' ' measurement{iMea}]);
    subplot(131)
    boxplot(mass_clustermass_TvC); hold on;
    title('TvC')
    if sig_mass_TvC == 0
        plot(1,obs_clustermass_TvC,'ro','LineWidth',4)
    else
        plot(1,obs_clustermass_TvC,'go','LineWidth',4)
    end
    legend(['p = ' num2str(pVal_TvC)])
    
    sig_mass_TvV = sum(mass_clustermass_TvV>obs_clustermass_TvV,2); %compare how many instances the sum of the permutation is greater than observed
    pVal_TvV = sig_mass_TvV/nperms;
    permMean_TvV = mean(mass_clustermass_TvV);
    permSTD_TvV = std(mass_clustermass_TvV);
    sig_mass_TvV(sig_mass_TvV <= pthresh) = true; %if it's not more than 50, it's siginificant for that frequency
    sig_mass_TvV(sig_mass_TvV > pthresh) = false;
    
    subplot(132)
    boxplot(mass_clustermass_TvV); hold on;
    title('TvV')
    if sig_mass_TvV == 0
        plot(1,obs_clustermass_TvV,'ro','LineWidth',4)
    else
        plot(1,obs_clustermass_TvV,'go','LineWidth',4)
    end
    legend(['p = ' num2str(pVal_TvV)])
    
    % In how many instances is the clustermass of the permutation MORE than
    % the observed clustermass
    sig_mass_CvV = sum(mass_clustermass_CvV>obs_clustermass_CvV,2); %compare how many instances the sum of the permutation is greater than observed
    pVal_CvV = sig_mass_CvV/nperms;
    permMean_CvV = mean(mass_clustermass_CvV);
    permSTD_CvV = std(mass_clustermass_CvV);
    sig_mass_CvV(sig_mass_CvV <= pthresh) = true; %if it's not more than 50, it's siginificant for that frequency
    sig_mass_CvV(sig_mass_CvV > pthresh) = false;
    
    subplot(133)
    boxplot(mass_clustermass_CvV); hold on;
    title('CvV')
    if sig_mass_CvV == 0
        plot(1,obs_clustermass_CvV,'ro','LineWidth',4)
    else
        plot(1,obs_clustermass_CvV,'go','LineWidth',4)
    end
    legend(['p = ' num2str(pVal_CvV)])
    
    h = gcf;
    h.Renderer = 'Painters';
    set(h, 'PaperType', 'A4');
    set(h, 'PaperOrientation', 'landscape');
    set(h, 'PaperUnits', 'centimeters');
    h_title = ['Spont between Perm Full ' layer ' ' measurement{iMea}];
    sgtitle(h_title)
    savefig(h_title)
    saveas(gcf, [h_title '.png'])
    close(h)
    save([h_title '.mat'],...
        'pVal_TvC','permMean_TvC','permSTD_TvC','pVal_TvV','permMean_TvV',...
        'permSTD_TvV','pVal_CvV','permMean_CvV','permSTD_CvV')
    
    %% Check Significance of layers' clustermass
    
    figure('Name',['Obs vs Perm ' layer ' ' ...
        measurement{iMea} ' ocsillations'],'units','normalized','outerposition',[0 0 1 1])
    pVal_TvC = struct; permMean_TvC = struct; permSTD_TvC = struct;
    pVal_TvV = struct; permMean_TvV = struct; permSTD_TvV = struct;
    pVal_CvV = struct; permMean_CvV = struct; permSTD_CvV = struct;
    
    for iOsc = 1:length(osciName)
        % In how many instances is the clustermass of the permutation MORE than
        % the observed clustermass
        sig_mass_TvC = sum(perm_layer_TvC.(osciName{iOsc})>obs_layer_TvC.(osciName{iOsc}),2); %compare how many instances the sum of the permutation is greater than observed
        pVal_TvC.(osciName{iOsc}) = sig_mass_TvC/nperms;
        permMean_TvC.(osciName{iOsc}) = mean(perm_layer_TvC.(osciName{iOsc}));
        permSTD_TvC.(osciName{iOsc}) = std(perm_layer_TvC.(osciName{iOsc}));
        sig_mass_TvC(sig_mass_TvC <= pthresh) = true; %if it's not more than 50, it's siginificant for that frequency
        sig_mass_TvC(sig_mass_TvC > pthresh) = false;
        
        subplot(3,6,iOsc)
        boxplot(perm_layer_TvC.(osciName{iOsc})); hold on;
        title(['TvC ' osciName{iOsc}])
        if sig_mass_TvC == 0
            plot(1,obs_layer_TvC.(osciName{iOsc}),'ro','LineWidth',4)
        else
            plot(1,obs_layer_TvC.(osciName{iOsc}),'go','LineWidth',4)
        end
        legend(['p = ' num2str(pVal_TvC.(osciName{iOsc}))])
        
        %compare how many instances the sum of the permutation is greater than observed
        sig_mass_TvV = sum(perm_layer_TvV.(osciName{iOsc})>obs_layer_TvV.(osciName{iOsc}),2);
        pVal_TvV.(osciName{iOsc}) = sig_mass_TvV/nperms;
        permMean_TvV.(osciName{iOsc}) = mean(perm_layer_TvV.(osciName{iOsc}));
        permSTD_TvV.(osciName{iOsc}) = std(perm_layer_TvV.(osciName{iOsc}));
        sig_mass_TvV(sig_mass_TvV <= pthresh) = true; %if it's not more than 50, it's siginificant for that frequency
        sig_mass_TvV(sig_mass_TvV > pthresh) = false;
        
        subplot(3,6,iOsc+6)
        boxplot(perm_layer_TvV.(osciName{iOsc})); hold on;
        title(['TvV ' osciName{iOsc}])
        if sig_mass_TvV == 0
            plot(1,obs_layer_TvV.(osciName{iOsc}),'ro','LineWidth',4)
        else
            plot(1,obs_layer_TvV.(osciName{iOsc}),'go','LineWidth',4)
        end
        legend(['p = ' num2str(pVal_TvV.(osciName{iOsc}))])
        
        %compare how many instances the sum of the permutation is greater than observed
        sig_mass_CvV = sum(perm_layer_CvV.(osciName{iOsc})>obs_layer_CvV.(osciName{iOsc}),2);
        pVal_CvV.(osciName{iOsc}) = sig_mass_CvV/nperms;
        permMean_CvV.(osciName{iOsc}) = mean(perm_layer_CvV.(osciName{iOsc}));
        permSTD_CvV.(osciName{iOsc}) = std(perm_layer_CvV.(osciName{iOsc}));
        sig_mass_CvV(sig_mass_CvV <= pthresh) = true; %if it's not more than 50, it's siginificant for that frequency
        sig_mass_CvV(sig_mass_CvV > pthresh) = false;
        
        subplot(3,6,iOsc+12)
        boxplot(perm_layer_CvV.(osciName{iOsc})); hold on;
        title(['CvV ' osciName{iOsc}])
        if sig_mass_CvV == 0
            plot(1,obs_layer_CvV.(osciName{iOsc}),'ro','LineWidth',4)
        else
            plot(1,obs_layer_CvV.(osciName{iOsc}),'go','LineWidth',4)
        end
        legend(['p = ' num2str(pVal_CvV.(osciName{iOsc}))])
        
    end
    
    h = gcf;
    h.Renderer = 'Painters';
    set(h, 'PaperType', 'A4');
    set(h, 'PaperOrientation', 'landscape');
    set(h, 'PaperUnits', 'centimeters');
    h_title = ['Spont between Perm Osci ' layer ' ' measurement{iMea}];
    sgtitle(h_title)
    savefig(h_title)
    saveas(gcf, [h_title '.png'])
    close(h)
    save([h_title '.mat'],...
        'pVal_TvC','permMean_TvC','permSTD_TvC','pVal_TvV','permMean_TvV',...
        'permSTD_TvV','pVal_CvV','permMean_CvV','permSTD_CvV')
end
