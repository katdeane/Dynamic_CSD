function SpectralPerm_Spontwithin(layer,homedir)

% Input:    Layer to analyze and home directory,
%           Loads wavelet transform (WT) mat file from Data/Spectral
% Output:   Figures for means and observed difference of measurements pre 1 vs post 1,
%           pre 2 vs post 2, pre 1 vs pre 2, and pre 1 vs end comparison and
%           figures for observed t values, and clusters in figs/Crypt_MagPerm_SP;
%           boxplot and significance of permutation test in figs/Crypt_MagPerm_SP

% note:     Permutation done within group on spontaneous measurements before
%           and after the 2 laser stimuli and at the very end of recording

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
if ~exist('layer','var')
    layer = 'IV'; % run the granular layer if not specified
end


cd (homedir),cd DATA;

% determine amount of permutations and detection level of significance
nperms = 1000;
pthresh = nperms*(0.05); %if bonferroni wanted, *(0.05/7)

% frequencies can be found in WT_CL.freq{1} to clarify the following
% rows choices; actual intended rows commented
theta = (49:54);        %(4:7);
alpha = (44:48);        %(8:12);
beta_low = (39:43);     %(13:18);
beta_high = (34:38);    %(19:30);
gamma_low = (26:33);    %(31:60);
gamma_high = (19:25);   %(61:100);

osciName    = {'theta' 'alpha' 'beta_low' 'beta_high' 'gamma_low' 'gamma_high'};
osciRows    = {theta alpha beta_low beta_high gamma_low gamma_high};
group       = {'KIC','KIT','KIV'};
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

for iGro = 1:length(group)
    %Pull out conditions and limit them to 1300ms
    preL1 = table2cell(WT_SP(contains(WT_SP.group,group{iGro}) & ...
        startsWith(WT_SP.measurement,measurement{1}),1));
    preL1 =  cellfun(@(x) x(:,1:params.limit),preL1,'UniformOutput',false);
    
    postL1 = table2cell(WT_SP(contains(WT_SP.group,group{iGro}) & ...
        startsWith(WT_SP.measurement,measurement{2}),1));
    postL1 =  cellfun(@(x) x(:,1:params.limit),postL1,'UniformOutput',false);
    
    preL2 = table2cell(WT_SP(contains(WT_SP.group,group{iGro}) & ...
        startsWith(WT_SP.measurement,measurement{3}),1));
    preL2 =  cellfun(@(x) x(:,1:params.limit),preL2,'UniformOutput',false);
    
    postL2 = table2cell(WT_SP(contains(WT_SP.group,group{iGro}) & ...
        startsWith(WT_SP.measurement,measurement{4}),1));
    postL2 =  cellfun(@(x) x(:,1:params.limit),postL2,'UniformOutput',false);
    
    postAll = table2cell(WT_SP(contains(WT_SP.group,group{iGro}) & ...
        startsWith(WT_SP.measurement,measurement{5}),1));
    postAll =  cellfun(@(x) x(:,1:params.limit),postAll,'UniformOutput',false);
    
    %Stack the individual animals' data (animal#x54x600)
    Pre1 = nan(size(preL1,1),size(preL1{1,1},1),size(preL1{1,1},2));
    Pos1 = nan(size(postL1,1),size(postL1{1,1},1),size(postL1{1,1},2));
    Pre2 = nan(size(preL2,1),size(preL2{1,1},1),size(preL2{1,1},2));
    Pos2 = nan(size(postL2,1),size(postL2{1,1},1),size(postL2{1,1},2));
    Post = nan(size(postAll,1),size(postAll{1,1},1),size(postAll{1,1},2));
    
    for ii = 1:length(preL1)
        Pre1(ii,:,:)= preL1{ii};
    end
    Pre1 = abs(Pre1);
    
    for ii = 1:length(postL1)
        Pos1(ii,:,:)= postL1{ii};
    end
    Pos1 = abs(Pos1);
    
    for ii = 1:length(preL2)
        Pre2(ii,:,:)= preL2{ii};
    end
    Pre2 = abs(Pre2);
    
    for ii = 1:length(postL2)
        Pos2(ii,:,:)= postL2{ii};
    end
    Pos2 = abs(Pos2);
    
    for ii = 1:length(postAll)
        Post(ii,:,:)= postAll{ii};
    end
    Post = abs(Post);
    
    grpsize = length(preL1);
      
    clear KIT08CasePre1 KIT08CasepreL1
    if contains(group{iGro},'KIT') && length(preL1) ~= length(postL1) %KIT08 missing Post L 1 measurement
        KIT08CasePre1  = Pre1(1:2,:,:);
        KIT08CasepreL1 = preL1(1:2);
    end
    %% Degrees of Freedom and t Threshold
    
    df = grpsize-1;
    
    if df == 0; continue; end
    
    Tchart = [12.71,4.303,3.182,2.776,2.571,2.447,2.365,2.306,2.262,2.228,...
        2.201,2.179,2.160,2.145,2.131,2.120,2.110,2.101,2.093,2.086];
    if df > 20
        error('Tchart only goes up to df == 20, to expand check this link: http://www.ttable.org/')
    end
    
    Tthresh = Tchart(df);
    
    %% Permutation Step 1 - Observed Differences
    
    obsPre1_mean = nanmean(Pre1,1);
    obsPos1_mean = nanmean(Pos1,1);
    obsPre2_mean = nanmean(Pre2,1);
    obsPos2_mean = nanmean(Pos2,1);
    obsPost_mean = nanmean(Post,1);
    
    % pre 1 vs post 1
    if ~exist('KIT08CasePre1','var') %KIT08 missing Post L 1 measurement
        meanofdiff_Pr1vPo1  = nanmean(Pre1-Pos1,1);
        stdofdiff_Pr1vPo1   = nanstd(Pre1-Pos1,0,1);
        groupsize_1         = length(postL1);
    else
        meanofdiff_Pr1vPo1  = nanmean(KIT08CasePre1-Pos1,1);
        stdofdiff_Pr1vPo1   = nanstd(KIT08CasePre1-Pos1,0,1);
        groupsize_1         = length(postL1);
    end
    % pre 2 vs post 2
    meanofdiff_Pr2vPo2  = nanmean(Pre2-Pos2,1);
    stdofdiff_Pr2vPo2   = nanstd(Pre2-Pos2,0,1);
    % pre 1 vs pre 2
    meanofdiff_Pr1vPr2  = nanmean(Pre1-Pre2,1);
    stdofdiff_Pr1vPr2   = nanstd(Pre1-Pre2,0,1);
    % pre 1 vs end 
    meanofdiff_Pr1vEnd  = nanmean(Pre1-Post,1);
    stdofdiff_Pr1vEnd   = nanstd(Pre1-Post,0,1);
    
    %% Permutation Step 2 - t test
    % find the t values along all data points for each frequency bin
    obs_t_Pr1vPo1 = meanofdiff_Pr1vPo1./(stdofdiff_Pr1vPo1./sqrt(groupsize_1));
    obs_clusters_Pr1vPo1 = abs(obs_t_Pr1vPo1);
    obs_clusters_Pr1vPo1(obs_clusters_Pr1vPo1 < Tthresh) = NaN;
    obs_clusters_Pr1vPo1(obs_clusters_Pr1vPo1 >= Tthresh) = 1;
    % check cluster mass KD: this will need to be improved to include a shorter
    % time window for detection (after each stimulus click)
    obs_clustermass_Pr1vPo1 = nansum(nansum(obs_clusters_Pr1vPo1));
    
    obs_t_Pr2vPo2 = meanofdiff_Pr2vPo2./(stdofdiff_Pr2vPo2./sqrt(grpsize));
    obs_clusters_Pr2vPo2 = abs(obs_t_Pr2vPo2);
    obs_clusters_Pr2vPo2(obs_clusters_Pr2vPo2 < Tthresh) = NaN;
    obs_clusters_Pr2vPo2(obs_clusters_Pr2vPo2 >= Tthresh) = 1;
    obs_clustermass_Pr2vPo2 = nansum(nansum(obs_clusters_Pr2vPo2));
    
    obs_t_Pr1vPr2 = meanofdiff_Pr1vPr2./(stdofdiff_Pr1vPr2./sqrt(grpsize));
    obs_clusters_Pr1vPr2 = abs(obs_t_Pr1vPr2);
    obs_clusters_Pr1vPr2(obs_clusters_Pr1vPr2 < Tthresh) = NaN;
    obs_clusters_Pr1vPr2(obs_clusters_Pr1vPr2 >= Tthresh) = 1;
    obs_clustermass_Pr1vPr2 = nansum(nansum(obs_clusters_Pr1vPr2));
    
    obs_t_Pr1vEnd = meanofdiff_Pr1vEnd./(stdofdiff_Pr1vEnd./sqrt(grpsize));
    obs_clusters_Pr1vEnd = abs(obs_t_Pr1vEnd);
    obs_clusters_Pr1vEnd(obs_clusters_Pr1vEnd < Tthresh) = NaN;
    obs_clusters_Pr1vEnd(obs_clusters_Pr1vEnd >= Tthresh) = 1;
    obs_clustermass_Pr1vEnd = nansum(nansum(obs_clusters_Pr1vEnd));
    
    %% Oscillation frequency bands:
    % % pull out clusters
    
    obs_layer_Pr1vPo1 = struct;
    obs_layer_Pr2vPo2 = struct;
    obs_layer_Pr1vPr2 = struct;
    obs_layer_Pr1vEnd = struct;
    for iOsc = 1:length(osciName)
        obs_layer_Pr1vPo1.(osciName{iOsc}) = obs_clusters_Pr1vPo1(:,osciRows{iOsc},:);
        obs_layer_Pr1vPo1.(osciName{iOsc}) = nansum(nansum(obs_layer_Pr1vPo1.(osciName{iOsc})));
        obs_layer_Pr2vPo2.(osciName{iOsc}) = obs_clusters_Pr2vPo2(:,osciRows{iOsc},:);
        obs_layer_Pr2vPo2.(osciName{iOsc}) = nansum(nansum(obs_layer_Pr2vPo2.(osciName{iOsc})));
        obs_layer_Pr1vPr2.(osciName{iOsc}) = obs_clusters_Pr1vPr2(:,osciRows{iOsc},:);
        obs_layer_Pr1vPr2.(osciName{iOsc}) = nansum(nansum(obs_layer_Pr1vPr2.(osciName{iOsc})));
        obs_layer_Pr1vEnd.(osciName{iOsc}) = obs_clusters_Pr1vEnd(:,osciRows{iOsc},:);
        obs_layer_Pr1vEnd.(osciName{iOsc}) = nansum(nansum(obs_layer_Pr1vEnd.(osciName{iOsc})));
    end
    
    cd(homedir); cd figs; mkdir('Crypt_MagPerm_SP'); cd('Crypt_MagPerm_SP');
    %% dif figs
    [X,Y]=meshgrid(WT_SP.freq{1}(19:54),params.startTime*1000:(params.limit-201));
    figure('Name','Observed Spectral Power and Differences',...
        'units','normalized','outerposition',[0 0 1 1]);
    
    % observed values
    pre1Fig = subplot(251);
    surf(Y,X,squeeze(obsPre1_mean(:,19:54,:))','EdgeColor','None'); view(2);
    set(gca,'YScale','log'); title(measurement{1})
    yticks([0 10 20 30 40 50 60 80 100])
    xlim([-200 1080])
    clim = get(gca,'clim');
    
    pos1Fig = subplot(252);
    surf(Y,X,squeeze(obsPos1_mean(:,19:54,:))','EdgeColor','None'); view(2);
    set(gca,'YScale','log'); title(measurement{2})
    yticks([0 10 20 30 40 50 60 80 100])
    xlim([-200 1080])
    clim = [clim; get(gca,'clim')]; %#ok<*AGROW>
    
    pre2Fig = subplot(253);
    surf(Y,X,squeeze(obsPre2_mean(:,19:54,:))','EdgeColor','None'); view(2);
    set(gca,'YScale','log'); title(measurement{3})
    yticks([0 10 20 30 40 50 60 80 100])
    xlim([-200 1080])
    clim = [clim; get(gca,'clim')];
    
    pos2Fig = subplot(254);
    surf(Y,X,squeeze(obsPos2_mean(:,19:54,:))','EdgeColor','None'); view(2);
    set(gca,'YScale','log'); title(measurement{4})
    yticks([0 10 20 30 40 50 60 80 100])
    xlim([-200 1080])
    clim = [clim; get(gca,'clim')];
    
    postFig = subplot(255);
    surf(Y,X,squeeze(obsPost_mean(:,19:54,:))','EdgeColor','None'); view(2);
    set(gca,'YScale','log'); title(measurement{5})
    yticks([0 10 20 30 40 50 60 80 100])
    xlim([-200 1080])
    clim = [clim; get(gca,'clim')];
    
    % observed difference
    dif1Fig = subplot(256);
    surf(Y,X,squeeze(meanofdiff_Pr1vPo1(:,19:54,:))','EdgeColor','None'); view(2);
    set(gca,'YScale','log'); title('Mean of Diff Pre 1 v Post 1')
    yticks([0 10 20 30 40 50 60 80 100])
    xlim([-200 1080])
    clim = [clim; get(gca,'clim')];
    
    dif2Fig = subplot(257);
    surf(Y,X,squeeze(meanofdiff_Pr2vPo2(:,19:54,:))','EdgeColor','None'); view(2);
    set(gca,'YScale','log'); title('Mean of Diff Pre 2 v Post 2')
    yticks([0 10 20 30 40 50 60 80 100])
    xlim([-200 1080])
    clim = [clim; get(gca,'clim')];
    
    dif3Fig = subplot(258);
    surf(Y,X,squeeze(meanofdiff_Pr1vPr2(:,19:54,:))','EdgeColor','None'); view(2);
    set(gca,'YScale','log'); title('Mean of Diff Pre 1 v Pre 2')
    yticks([0 10 20 30 40 50 60 80 100])
    xlim([-200 1080])
    clim = [clim; get(gca,'clim')];
    
    dif4Fig = subplot(259);
    surf(Y,X,squeeze(meanofdiff_Pr1vEnd(:,19:54,:))','EdgeColor','None'); view(2);
    set(gca,'YScale','log'); title('Mean of Diff Pre 1 v End')
    yticks([0 10 20 30 40 50 60 80 100])
    xlim([-200 1080])
    clim = [clim; get(gca,'clim')];
    
    newC = [min(clim(:)) max(clim(:))];
    
    set(pos1Fig,'Clim',newC);
    set(pre2Fig,'Clim',newC);
    set(pos2Fig,'Clim',newC);
    set(postFig,'Clim',newC);
    set(dif1Fig,'Clim',newC);
    set(dif2Fig,'Clim',newC);
    set(dif3Fig,'Clim',newC);
    set(dif4Fig,'Clim',newC);
    set(pre1Fig,'Clim',newC);colorbar;
    
    h = gcf;
    h.Renderer = 'Painters';
    set(h, 'PaperType', 'A4');
    set(h, 'PaperOrientation', 'landscape');
    set(h, 'PaperUnits', 'centimeters');
    h_title = ['Spont within Spectral Power of ' layer ' ' group{iGro}];
    sgtitle(h_title);
    savefig(h_title)
    saveas(gcf, [h_title '.png'])
    close(h)
    
    %% t fig
    
    [X,Y]=meshgrid(WT_SP.freq{1}(19:54),params.startTime*1000:(params.limit-201));
    figure('Name','Observed t Values BF', ...
        'units','normalized','outerposition',[0 0 1 1]);
    colormap('winter')
    
    subplot(241);
    surf(Y,X,squeeze(obs_t_Pr1vPo1(:,19:54,:))','EdgeColor','None'); view(2);
    set(gca,'YScale','log'); title('t Value Mat pre 1 vs post 1')
    yticks([0 10 20 30 40 50 60 80 100])
    xlim([-200 1080])
    colorbar
    
    subplot(242);
    surf(Y,X,squeeze(obs_t_Pr2vPo2(:,19:54,:))','EdgeColor','None'); view(2);
    set(gca,'YScale','log'); title('t Value Mat pre 2 vs post 2')
    yticks([0 10 20 30 40 50 60 80 100])
    xlim([-200 1080])
    colorbar
    
    subplot(243);
    surf(Y,X,squeeze(obs_t_Pr1vPr2(:,19:54,:))','EdgeColor','None'); view(2);
    set(gca,'YScale','log'); title('t Value Mat pre 1 vs pre 2')
    yticks([0 10 20 30 40 50 60 80 100])
    xlim([-200 1080])
    colorbar
    
    subplot(244);
    surf(Y,X,squeeze(obs_t_Pr1vEnd(:,19:54,:))','EdgeColor','None'); view(2);
    set(gca,'YScale','log'); title('t Value Mat pre 1 vs end')
    yticks([0 10 20 30 40 50 60 80 100])
    xlim([-200 1080])
    colorbar
    
    subplot(245);
    surf(Y,X,squeeze(obs_clusters_Pr1vPo1(:,19:54,:))','EdgeColor','None'); view(2);
    set(gca,'YScale','log'); title('Clusters where p<0.05')
    yticks([0 10 20 30 40 50 60 80 100])
    xlim([-200 1080])
    colorbar
    
    subplot(246);
    surf(Y,X,squeeze(obs_clusters_Pr2vPo2(:,19:54,:))','EdgeColor','None'); view(2);
    set(gca,'YScale','log'); title('Clusters where p<0.05')
    yticks([0 10 20 30 40 50 60 80 100])
    xlim([-200 1080])
    colorbar
    
    subplot(247);
    surf(Y,X,squeeze(obs_clusters_Pr1vPr2(:,19:54,:))','EdgeColor','None'); view(2);
    set(gca,'YScale','log'); title('Clusters where p<0.05')
    yticks([0 10 20 30 40 50 60 80 100])
    xlim([-200 1080])
    colorbar
    
    subplot(248);
    surf(Y,X,squeeze(obs_clusters_Pr1vEnd(:,19:54,:))','EdgeColor','None'); view(2);
    set(gca,'YScale','log'); title('Clusters where p<0.05')
    yticks([0 10 20 30 40 50 60 80 100])
    xlim([-200 1080])
    colorbar
    
    h = gcf;
    h.Renderer = 'Painters';
    set(h, 'PaperType', 'A4');
    set(h, 'PaperOrientation', 'landscape');
    set(h, 'PaperUnits', 'centimeters');
    h_title = ['Spont within Observed t and p ' layer ' ' group{iGro}];
    sgtitle(h_title);
    savefig(h_title)
    saveas(gcf, [h_title '.png'])
    close(h)
    
    
    %% Permutation Step 3 - do the permute
    % preallocate
    mass_clustermass_Pr1vPo1 = NaN([1 nperms]);
    mass_clustermass_Pr2vPo2 = NaN([1 nperms]);
    mass_clustermass_Pr1vPr2 = NaN([1 nperms]);
    mass_clustermass_Pr1vEnd = NaN([1 nperms]);
    
    perm_layer_Pr1vPo1 = struct;
    perm_layer_Pr2vPo2 = struct;
    perm_layer_Pr1vPr2 = struct;
    perm_layer_Pr1vEnd = struct;
    for iOsc = 1:length(osciName)
        perm_layer_Pr1vPo1.(osciName{iOsc}) = NaN([1 nperms]);
        perm_layer_Pr2vPo2.(osciName{iOsc}) = NaN([1 nperms]);
        perm_layer_Pr1vPr2.(osciName{iOsc}) = NaN([1 nperms]);
        perm_layer_Pr1vEnd.(osciName{iOsc}) = NaN([1 nperms]);
    end
    
    for iComp = 1:4 % specify comparison and data going in
        if iComp == 1
            if ~exist('KIT08CasepreL1','var') % KIT08 missing Post L 1
                contAll = vertcat(preL1,postL1);
            else
                contAll = vertcat(KIT08CasepreL1,postL1);
            end
        elseif iComp == 2
            contAll = vertcat(preL1,preL2);
        elseif iComp == 3
            contAll = vertcat(preL1,postL2);
        elseif iComp == 4
            contAll = vertcat(preL1,postAll);
        end
        
        grpsize = size(contAll,1)/2;
        % generate the permuations and take clustermass values
        for iperm = 1:nperms
            % determine random list order to pull
            order = randperm(length(contAll));
            
            % pull based on random list order
            M1Perm = nan(grpsize,size(contAll{1,1},1),size(contAll{1,1},2));
            M2Perm = nan(grpsize,size(contAll{1,1},1),size(contAll{1,1},2));
            for ii = 1:grpsize
                M2Perm(ii,:,:)=contAll{order(ii)};
            end
            M2Perm = abs(M2Perm);
            
            for ii = (grpsize)+1:grpsize+grpsize
                M1Perm(ii-grpsize,:,:)=contAll{order(ii)};
            end
            M1Perm = abs(M1Perm);
            
            % did the permute.
            
            perm_meanofdiff = nanmean(M1Perm-M2Perm,1);
            perm_stdofdiff  = nanstd(M1Perm-M2Perm,0,1);
            
            % t test %%%
            per_t = perm_meanofdiff./(perm_stdofdiff./sqrt(grpsize));
            per_clusters = abs(per_t);
            per_clusters(per_clusters < Tthresh) = NaN;
            per_clusters(per_clusters >= Tthresh) = 1;
            
            % check cluster mass for 300 ms from tone onset
            per_clustermass = nansum(nansum(per_clusters));
            
            if iComp == 1
                mass_clustermass_Pr1vPo1(iperm) = per_clustermass;
            elseif iComp == 2
                mass_clustermass_Pr2vPo2(iperm) = per_clustermass;
            elseif iComp == 3
                mass_clustermass_Pr1vPr2(iperm) = per_clustermass;
            elseif iComp == 4
                mass_clustermass_Pr1vEnd(iperm) = per_clustermass;
            end
            
            for iOsc = 1:length(osciName)
                hold_permlayer = per_clusters(:,osciRows{iOsc},:);
                hold_permlayer = nansum(nansum(hold_permlayer));
                if iComp == 1
                    perm_layer_Pr1vPo1.(osciName{iOsc})(iperm) = hold_permlayer;
                elseif iComp == 2
                    perm_layer_Pr2vPo2.(osciName{iOsc})(iperm) = hold_permlayer;
                elseif iComp == 3
                    perm_layer_Pr1vPr2.(osciName{iOsc})(iperm) = hold_permlayer;
                elseif iComp == 4
                    perm_layer_Pr1vEnd.(osciName{iOsc})(iperm) = hold_permlayer;
                end
            end
        end
    end
    
    cd(homedir); cd DATA; cd Spectral; mkdir('Crypt_MagPerm_SP'); cd('Crypt_MagPerm_SP');
    %% Check Significance of full clustermass
    
    figure('Name',['Spont within Obs cluster vs Perm ' layer ' ' group{iGro}]);
    % In how many instances is the clustermass of the permutation MORE than
    % the observed clustermass
    sig_mass_Pr1vPo1 = sum(mass_clustermass_Pr1vPo1>obs_clustermass_Pr1vPo1,2); %compare how many instances the sum of the permutation is greater than observed
    pVal_Pr1vPo1 = sig_mass_Pr1vPo1/nperms;
    permMean_Pr1vPo1 = mean(mass_clustermass_Pr1vPo1);
    permSTD_Pr1vPo1 = std(mass_clustermass_Pr1vPo1);
    sig_mass_Pr1vPo1(sig_mass_Pr1vPo1 <= pthresh) = true; %if it's not more than 50, it's siginificant for that frequency
    sig_mass_Pr1vPo1(sig_mass_Pr1vPo1 > pthresh) = false;
    
    subplot(141)
    boxplot(mass_clustermass_Pr1vPo1); hold on;
    title([measurement{1} ' v ' measurement{2}])
    if sig_mass_Pr1vPo1 == 0
        plot(1,obs_clustermass_Pr1vPo1,'ro','LineWidth',4)
    else
        plot(1,obs_clustermass_Pr1vPo1,'go','LineWidth',4)
    end
    legend(['p = ' num2str(pVal_Pr1vPo1)])
    
    sig_mass_Pr2vPo2 = sum(mass_clustermass_Pr2vPo2>obs_clustermass_Pr2vPo2,2); %compare how many instances the sum of the permutation is greater than observed
    pVal_Pr2vPo2 = sig_mass_Pr2vPo2/nperms;
    permMean_Pr2vPo2 = mean(mass_clustermass_Pr2vPo2);
    permSTD_Pr2vPo2 = std(mass_clustermass_Pr2vPo2);
    sig_mass_Pr2vPo2(sig_mass_Pr2vPo2 <= pthresh) = true; %if it's not more than 50, it's siginificant for that frequency
    sig_mass_Pr2vPo2(sig_mass_Pr2vPo2 > pthresh) = false;
    
    subplot(142)
    boxplot(mass_clustermass_Pr2vPo2); hold on;
    title([measurement{3} ' v ' measurement{4}])
    if sig_mass_Pr2vPo2 == 0
        plot(1,obs_clustermass_Pr2vPo2,'ro','LineWidth',4)
    else
        plot(1,obs_clustermass_Pr2vPo2,'go','LineWidth',4)
    end
    legend(['p = ' num2str(pVal_Pr2vPo2)])
    
    sig_mass_Pr1vPr2 = sum(mass_clustermass_Pr1vPr2>obs_clustermass_Pr1vPr2,2); %compare how many instances the sum of the permutation is greater than observed
    pVal_Pr1vPr2 = sig_mass_Pr1vPr2/nperms;
    permMean_Pr1vPr2 = mean(mass_clustermass_Pr1vPr2);
    permSTD_Pr1vPr2 = std(mass_clustermass_Pr1vPr2);
    sig_mass_Pr1vPr2(sig_mass_Pr1vPr2 <= pthresh) = true; %if it's not more than 50, it's siginificant for that frequency
    sig_mass_Pr1vPr2(sig_mass_Pr1vPr2 > pthresh) = false;
    
    subplot(143)
    boxplot(mass_clustermass_Pr1vPr2); hold on;
    title([measurement{1} ' v ' measurement{3}])
    if sig_mass_Pr1vPr2 == 0
        plot(1,obs_clustermass_Pr1vPr2,'ro','LineWidth',4)
    else
        plot(1,obs_clustermass_Pr1vPr2,'go','LineWidth',4)
    end
    legend(['p = ' num2str(pVal_Pr1vPr2)])
    
    sig_mass_Pr1vEnd = sum(mass_clustermass_Pr1vEnd>obs_clustermass_Pr1vEnd,2); %compare how many instances the sum of the permutation is greater than observed
    pVal_Pr1vEnd = sig_mass_Pr1vEnd/nperms;
    permMean_Pr1vEnd = mean(mass_clustermass_Pr1vEnd);
    permSTD_Pr1vEnd = std(mass_clustermass_Pr1vEnd);
    sig_mass_Pr1vEnd(sig_mass_Pr1vEnd <= pthresh) = true; %if it's not more than 50, it's siginificant for that frequency
    sig_mass_Pr1vEnd(sig_mass_Pr1vEnd > pthresh) = false;
    
    subplot(144)
    boxplot(mass_clustermass_Pr1vEnd); hold on;
    title([measurement{1} ' v ' measurement{5}])
    if sig_mass_Pr1vEnd == 0
        plot(1,obs_clustermass_Pr1vEnd,'ro','LineWidth',4)
    else
        plot(1,obs_clustermass_Pr1vEnd,'go','LineWidth',4)
    end
    legend(['p = ' num2str(pVal_Pr1vEnd)])
    
    h = gcf;
    h.Renderer = 'Painters';
    set(h, 'PaperType', 'A4');
    set(h, 'PaperOrientation', 'landscape');
    set(h, 'PaperUnits', 'centimeters');
    h_title = ['Spont within Obs cluster vs Perm ' layer ' ' group{iGro}];
    sgtitle(h_title);
    savefig(h_title)
    saveas(gcf, [h_title '.png'])
    close(h)
    save(['Spont within Perm Full ' layer ' ' group{iGro} '.mat'],...
        'pVal_Pr1vPo1','permMean_Pr1vPo1','permSTD_Pr1vPo1','pVal_Pr2vPo2',...
        'permMean_Pr2vPo2','permSTD_Pr2vPo2','pVal_Pr1vPr2','permMean_Pr1vPr2',...
        'permSTD_Pr1vPr2','pVal_Pr1vEnd','permMean_Pr1vEnd','permSTD_Pr1vEnd')
    
    %% Check Significance of layers' clustermass
    
    figure('Name',['Spont within Obs vs Perm ' layer ' ' group{iGro} ' ocsillations']...
        ,'units','normalized','outerposition',[0 0 1 1])
    pVal_Pr1vPo1 = struct; permMean_Pr1vPo1 = struct; permSTD_Pr1vPo1 = struct;
    pVal_Pr2vPo2 = struct; permMean_Pr2vPo2 = struct; permSTD_Pr2vPo2 = struct;
    pVal_Pr1vPr2 = struct; permMean_Pr1vPr2 = struct; permSTD_Pr1vPr2 = struct;
    pVal_Pr1vEnd = struct; permMean_Pr1vEnd = struct; permSTD_Pr1vEnd = struct;
    
    for iOsc = 1:length(osciName)
        % In how many instances is the clustermass of the permutation MORE than
        % the observed clustermass
        sig_mass_Pr1vPo1 = sum(perm_layer_Pr1vPo1.(osciName{iOsc})>obs_layer_Pr1vPo1.(osciName{iOsc}),2); %compare how many instances the sum of the permutation is greater than observed
        pVal_Pr1vPo1.(osciName{iOsc}) = sig_mass_Pr1vPo1/nperms;
        permMean_Pr1vPo1.(osciName{iOsc}) = mean(perm_layer_Pr1vPo1.(osciName{iOsc}));
        permSTD_Pr1vPo1.(osciName{iOsc}) = std(perm_layer_Pr1vPo1.(osciName{iOsc}));
        sig_mass_Pr1vPo1(sig_mass_Pr1vPo1 <= pthresh) = true; %if it's not more than 50, it's siginificant for that frequency
        sig_mass_Pr1vPo1(sig_mass_Pr1vPo1 > pthresh) = false;
        
        subplot(4,6,iOsc)
        boxplot(perm_layer_Pr1vPo1.(osciName{iOsc})); hold on;
        title([measurement{1} 'v' measurement{2} ' ' osciName{iOsc}])
        if sig_mass_Pr1vPo1 == 0
            plot(1,obs_layer_Pr1vPo1.(osciName{iOsc}),'ro','LineWidth',4)
        else
            plot(1,obs_layer_Pr1vPo1.(osciName{iOsc}),'go','LineWidth',4)
        end
        legend(['p = ' num2str(pVal_Pr1vPo1.(osciName{iOsc}))])
        
        sig_mass_Pr2vPo2 = sum(perm_layer_Pr2vPo2.(osciName{iOsc})>obs_layer_Pr2vPo2.(osciName{iOsc}),2); %compare how many instances the sum of the permutation is greater than observed
        pVal_Pr2vPo2.(osciName{iOsc}) = sig_mass_Pr2vPo2/nperms;
        permMean_Pr2vPo2.(osciName{iOsc}) = mean(perm_layer_Pr2vPo2.(osciName{iOsc}));
        permSTD_Pr2vPo2.(osciName{iOsc}) = std(perm_layer_Pr2vPo2.(osciName{iOsc}));
        sig_mass_Pr2vPo2(sig_mass_Pr2vPo2 <= pthresh) = true; %if it's not more than 50, it's siginificant for that frequency
        sig_mass_Pr2vPo2(sig_mass_Pr2vPo2 > pthresh) = false;
        
        subplot(4,6,iOsc+6)
        boxplot(perm_layer_Pr2vPo2.(osciName{iOsc})); hold on;
        title([measurement{1} 'v' measurement{3} ' ' osciName{iOsc}])
        if sig_mass_Pr2vPo2 == 0
            plot(1,obs_layer_Pr2vPo2.(osciName{iOsc}),'ro','LineWidth',4)
        else
            plot(1,obs_layer_Pr2vPo2.(osciName{iOsc}),'go','LineWidth',4)
        end
        legend(['p = ' num2str(pVal_Pr2vPo2.(osciName{iOsc}))])
        
        sig_mass_Pr1vPr2 = sum(perm_layer_Pr1vPr2.(osciName{iOsc})>obs_layer_Pr1vPr2.(osciName{iOsc}),2); %compare how many instances the sum of the permutation is greater than observed
        pVal_Pr1vPr2.(osciName{iOsc}) = sig_mass_Pr1vPr2/nperms;
        permMean_Pr1vPr2.(osciName{iOsc}) = mean(perm_layer_Pr1vPr2.(osciName{iOsc}));
        permSTD_Pr1vPr2.(osciName{iOsc}) = std(perm_layer_Pr1vPr2.(osciName{iOsc}));
        sig_mass_Pr1vPr2(sig_mass_Pr1vPr2 <= pthresh) = true; %if it's not more than 50, it's siginificant for that frequency
        sig_mass_Pr1vPr2(sig_mass_Pr1vPr2 > pthresh) = false;
        
        subplot(4,6,iOsc+12)
        boxplot(perm_layer_Pr1vPr2.(osciName{iOsc})); hold on;
        title([measurement{1} 'v' measurement{4} ' ' osciName{iOsc}])
        if sig_mass_Pr1vPr2 == 0
            plot(1,obs_layer_Pr1vPr2.(osciName{iOsc}),'ro','LineWidth',4)
        else
            plot(1,obs_layer_Pr1vPr2.(osciName{iOsc}),'go','LineWidth',4)
        end
        legend(['p = ' num2str(pVal_Pr1vPr2.(osciName{iOsc}))])
        
        sig_mass_Pr1vEnd = sum(perm_layer_Pr1vEnd.(osciName{iOsc})>obs_layer_Pr1vEnd.(osciName{iOsc}),2); %compare how many instances the sum of the permutation is greater than observed
        pVal_Pr1vEnd.(osciName{iOsc}) = sig_mass_Pr1vEnd/nperms;
        permMean_Pr1vEnd.(osciName{iOsc}) = mean(perm_layer_Pr1vEnd.(osciName{iOsc}));
        permSTD_Pr1vEnd.(osciName{iOsc}) = std(perm_layer_Pr1vEnd.(osciName{iOsc}));
        sig_mass_Pr1vEnd(sig_mass_Pr1vEnd <= pthresh) = true; %if it's not more than 50, it's siginificant for that frequency
        sig_mass_Pr1vEnd(sig_mass_Pr1vEnd > pthresh) = false;
        
        subplot(4,6,iOsc+18)
        boxplot(perm_layer_Pr1vEnd.(osciName{iOsc})); hold on;
        title([measurement{1} 'v' measurement{5} ' ' osciName{iOsc}])
        if sig_mass_Pr1vEnd == 0
            plot(1,obs_layer_Pr1vEnd.(osciName{iOsc}),'ro','LineWidth',4)
        else
            plot(1,obs_layer_Pr1vEnd.(osciName{iOsc}),'go','LineWidth',4)
        end
        legend(['p = ' num2str(pVal_Pr1vEnd.(osciName{iOsc}))])
    end
    
    h = gcf;
    h.Renderer = 'Painters';
    set(h, 'PaperType', 'A4');
    set(h, 'PaperOrientation', 'landscape');
    set(h, 'PaperUnits', 'centimeters');
    h_title = ['Spont within Perm Osci ' layer ' ' group{iGro}];
    sgtitle(h_title)
    savefig(h_title)
    saveas(gcf, [h_title '.png'])
    close(h)
    save(['Spont within Perm Osci ' layer ' ' group{iGro} '.mat'],...
        'pVal_Pr1vPo1','permMean_Pr1vPo1','permSTD_Pr1vPo1','pVal_Pr2vPo2',...
        'permMean_Pr2vPo2','permSTD_Pr2vPo2','pVal_Pr1vPr2','permMean_Pr1vPr2',...
        'permSTD_Pr1vPr2','pVal_Pr1vEnd','permMean_Pr1vEnd','permSTD_Pr1vEnd')
    
end
