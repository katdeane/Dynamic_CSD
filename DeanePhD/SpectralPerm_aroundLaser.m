function SpectralPerm_aroundLaser(layer,homedir,type)

% Input:    Layer to analyze, home directory, type of stimuli (CL or AM).
%           Loads wavelet transform (WT) mat file from Data/Spectral 
% Output:   Figures for means and observed difference of measurements pre vs post 1,
%           pre vs post 2, pre vs post 3, and pre vs post 4 comparison and 
%           figures for observed t values, and clusters in figs/['Crypt_MagPerm_' type]; 
%           boxplot and significance of permutation test in figs/['Crypt_MagPerm_' type]

% note:     Permutation done within group on each measurement after laser
%           to the measurement before the laser

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
if ~exist('type','var')
    type = 'CL'; % run the click measurements
end

cd (homedir),cd DATA;

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
stimfreq    = [2,5,10,20,40];

%% Load in and seperate Data
disp(['loading mat file WT of ' type])
tic
if contains(type,'CL')
    measurement = {'preCL_1','CL_1','CL_2','CL_3','CL_4'};
    load('WT_CL.mat', 'WT_CL')
    WT = WT_CL(startsWith(WT_CL.layer,layer) & endsWith(WT_CL.layer,layer),:);
    clear WT_CL
elseif contains(type,'AM')
    measurement = {'preAM_1','AM_1','AM_2','AM_3','AM_4'};
    load('WT_AM.mat','WT_AM')
    WT = WT_AM(startsWith(WT_AM.layer,layer) & endsWith(WT_AM.layer,layer),:);
    clear WT_AM
else
     error('Input "type" needs to be either "CL" or "AM" for clicks or amplitude modulations respectively')
end
toc

% varNames = unique(wtTable.layer);
params.startTime = -0.2; % seconds
params.limit = 1300;

for iGro = 1:length(group)
    for iStim = 1:length(stimfreq)
        %Pull out conditions and limit them to 1300ms
        preL = table2cell(WT(contains(WT.group,group{iGro}) & ...
            WT.stimulus==stimfreq(iStim) & startsWith(WT.measurement,measurement{1}),1));
        preL =  cellfun(@(x) x(:,1:params.limit),preL,'UniformOutput',false);
        
        postL1 = table2cell(WT(contains(WT.group,group{iGro}) & ...
            WT.stimulus==stimfreq(iStim) & startsWith(WT.measurement,measurement{2}),1));
        postL1 =  cellfun(@(x) x(:,1:params.limit),postL1,'UniformOutput',false);
        
        postL2 = table2cell(WT(contains(WT.group,group{iGro}) & ...
            WT.stimulus==stimfreq(iStim) & startsWith(WT.measurement,measurement{3}),1));
        postL2 =  cellfun(@(x) x(:,1:params.limit),postL2,'UniformOutput',false);
        
        postL3 = table2cell(WT(contains(WT.group,group{iGro}) & ...
            WT.stimulus==stimfreq(iStim) & startsWith(WT.measurement,measurement{4}),1));
        postL3 =  cellfun(@(x) x(:,1:params.limit),postL3,'UniformOutput',false);
        
        postL4 = table2cell(WT(contains(WT.group,group{iGro}) & ...
            WT.stimulus==stimfreq(iStim) & startsWith(WT.measurement,measurement{5}),1));
        postL4 =  cellfun(@(x) x(:,1:params.limit),postL4,'UniformOutput',false);
        
        %Stack the individual animals' data (animal#x54x600)
        Pre = nan(size(preL,1),size(preL{1,1},1),size(preL{1,1},2));
        Po1 = nan(size(postL1,1),size(postL1{1,1},1),size(postL1{1,1},2));
        Po2 = nan(size(postL2,1),size(postL2{1,1},1),size(postL2{1,1},2));
        Po3 = nan(size(postL3,1),size(postL3{1,1},1),size(postL3{1,1},2));
        Po4 = nan(size(postL4,1),size(postL4{1,1},1),size(postL4{1,1},2));
        
        for ii = 1:length(preL)
            Pre(ii,:,:)= preL{ii};
        end
        Pre = abs(Pre);
        
        for ii = 1:length(postL1)
            Po1(ii,:,:)= postL1{ii};
        end
        Po1 = abs(Po1);
        
        for ii = 1:length(postL2)
            Po2(ii,:,:)= postL2{ii};
        end
        Po2 = abs(Po2);
        
        for ii = 1:length(postL3)
            Po3(ii,:,:)= postL3{ii};
        end
        Po3 = abs(Po3);
        
        for ii = 1:length(postL4)
            Po4(ii,:,:)= postL4{ii};
        end
        Po4 = abs(Po4);
        
        grpsize = length(preL);
        
        if length(preL) ~= length(postL1)
            keyboard %this is to catch if measurments are missing
        end
        
        clear KIT04CasePre KIT04CasepreL
        if contains(group{iGro},'KIT') && length(preL) ~= length(postL4) %KIT04 missing AM_4 measurement
            KIT04CasePre  = vertcat(Pre(1:3,:,:),Pre(5:end,:,:));
            KIT04CasepreL = vertcat(preL(1:3),preL(5:end));
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
        
        obsPre_mean = nanmean(Pre,1);
        obsPo1_mean = nanmean(Po1,1);
        obsPo2_mean = nanmean(Po2,1);
        obsPo3_mean = nanmean(Po3,1);
        obsPo4_mean = nanmean(Po4,1);
        
        meanofdiff_v1 = nanmean(Pre-Po1,1);
        stdofdiff_v1  = nanstd(Pre-Po1,0,1);
        meanofdiff_v2 = nanmean(Pre-Po2,1);
        stdofdiff_v2  = nanstd(Pre-Po2,0,1);
        meanofdiff_v3 = nanmean(Pre-Po3,1);
        stdofdiff_v3  = nanstd(Pre-Po3,0,1);
        if ~exist('KIT04CasePre','var') %KIT04 missing AM_4 measurement
            meanofdiff_v4 = nanmean(Pre-Po4,1);
            stdofdiff_v4  = nanstd(Pre-Po4,0,1);
            grpsize_4     = length(postL4);
        else
            meanofdiff_v4 = nanmean(KIT04CasePre-Po4,1);
            stdofdiff_v4  = nanstd(KIT04CasePre-Po4,0,1);
            grpsize_4     = length(postL4);
        end
        
        %% Permutation Step 2 - t test
        % find the t values along all data points for each frequency bin
        obs_t_v1 = meanofdiff_v1./(stdofdiff_v1./sqrt(grpsize));
        obs_clusters_v1 = abs(obs_t_v1);
        obs_clusters_v1(obs_clusters_v1 < Tthresh) = NaN;
        obs_clusters_v1(obs_clusters_v1 >= Tthresh) = 1;
        % check cluster mass KD: this will need to be improved to include a shorter
        % time window for detection (after each stimulus click)
        obs_clustermass_v1 = nansum(nansum(obs_clusters_v1));
        
        obs_t_v2 = meanofdiff_v2./(stdofdiff_v2./sqrt(grpsize));
        obs_clusters_v2 = abs(obs_t_v2);
        obs_clusters_v2(obs_clusters_v2 < Tthresh) = NaN;
        obs_clusters_v2(obs_clusters_v2 >= Tthresh) = 1;
        obs_clustermass_v2 = nansum(nansum(obs_clusters_v2));
        
        obs_t_v3 = meanofdiff_v3./(stdofdiff_v3./sqrt(grpsize));
        obs_clusters_v3 = abs(obs_t_v3);
        obs_clusters_v3(obs_clusters_v3 < Tthresh) = NaN;
        obs_clusters_v3(obs_clusters_v3 >= Tthresh) = 1;
        obs_clustermass_v3 = nansum(nansum(obs_clusters_v3));
        
        obs_t_v4 = meanofdiff_v4./(stdofdiff_v4./sqrt(grpsize_4));
        obs_clusters_v4 = abs(obs_t_v4);
        obs_clusters_v4(obs_clusters_v4 < Tthresh) = NaN;
        obs_clusters_v4(obs_clusters_v4 >= Tthresh) = 1;
        obs_clustermass_v4 = nansum(nansum(obs_clusters_v4));
        
        %% Oscillation frequency bands:
        % % pull out clusters
        
        obs_layer_v1 = struct;
        obs_layer_v2 = struct;
        obs_layer_v3 = struct;
        obs_layer_v4 = struct;
        for iOsc = 1:length(osciName)
            obs_layer_v1.(osciName{iOsc}) = obs_clusters_v1(:,osciRows{iOsc},:);
            obs_layer_v1.(osciName{iOsc}) = nansum(nansum(obs_layer_v1.(osciName{iOsc})));
            obs_layer_v2.(osciName{iOsc}) = obs_clusters_v2(:,osciRows{iOsc},:);
            obs_layer_v2.(osciName{iOsc}) = nansum(nansum(obs_layer_v2.(osciName{iOsc})));
            obs_layer_v3.(osciName{iOsc}) = obs_clusters_v3(:,osciRows{iOsc},:);
            obs_layer_v3.(osciName{iOsc}) = nansum(nansum(obs_layer_v3.(osciName{iOsc})));
            obs_layer_v4.(osciName{iOsc}) = obs_clusters_v4(:,osciRows{iOsc},:);
            obs_layer_v4.(osciName{iOsc}) = nansum(nansum(obs_layer_v4.(osciName{iOsc})));
        end
        
        cd(homedir); cd figs; mkdir(['Crypt_MagPerm_' type]); cd(['Crypt_MagPerm_' type]);
        %% dif figs
        [X,Y]=meshgrid(WT.freq{1}(19:54),params.startTime*1000:(params.limit-201));
        figure('Name','Observed Spectral Power and Differences',...
            'units','normalized','outerposition',[0 0 1 1]);
        
        % observed values
        preFig = subplot(341);
        surf(Y,X,squeeze(obsPre_mean(:,19:54,:))','EdgeColor','None'); view(2);
        set(gca,'YScale','log'); title(measurement{1})
        yticks([0 10 20 30 40 50 60 80 100])
        xlim([-200 1080])
        clim = get(gca,'clim');
        
        Po1Fig = subplot(345);
        surf(Y,X,squeeze(obsPo1_mean(:,19:54,:))','EdgeColor','None'); view(2);
        set(gca,'YScale','log'); title(measurement{2})
        yticks([0 10 20 30 40 50 60 80 100])
        xlim([-200 1080])
        clim = [clim; get(gca,'clim')]; %#ok<*AGROW>
        
        Po2Fig = subplot(346);
        surf(Y,X,squeeze(obsPo2_mean(:,19:54,:))','EdgeColor','None'); view(2);
        set(gca,'YScale','log'); title(measurement{3})
        yticks([0 10 20 30 40 50 60 80 100])
        xlim([-200 1080])
        clim = [clim; get(gca,'clim')];
        
        Po3Fig = subplot(347);
        surf(Y,X,squeeze(obsPo3_mean(:,19:54,:))','EdgeColor','None'); view(2);
        set(gca,'YScale','log'); title(measurement{4})
        yticks([0 10 20 30 40 50 60 80 100])
        xlim([-200 1080])
        clim = [clim; get(gca,'clim')];
        
        Po4Fig = subplot(348);
        surf(Y,X,squeeze(obsPo4_mean(:,19:54,:))','EdgeColor','None'); view(2);
        set(gca,'YScale','log'); title(measurement{5})
        yticks([0 10 20 30 40 50 60 80 100])
        xlim([-200 1080])
        clim = [clim; get(gca,'clim')];
        
        % observed difference
        dif1Fig = subplot(349);
        surf(Y,X,squeeze(meanofdiff_v1(:,19:54,:))','EdgeColor','None'); view(2);
        set(gca,'YScale','log'); title('Mean of Diff Pre v Post 1')
        yticks([0 10 20 30 40 50 60 80 100])
        xlim([-200 1080])
        clim = [clim; get(gca,'clim')];
        
        dif2Fig = subplot(3,4,10);
        surf(Y,X,squeeze(meanofdiff_v2(:,19:54,:))','EdgeColor','None'); view(2);
        set(gca,'YScale','log'); title('Mean of Diff Pre v Post 2')
        yticks([0 10 20 30 40 50 60 80 100])
        xlim([-200 1080])
        clim = [clim; get(gca,'clim')];
        
        dif3Fig = subplot(3,4,11);
        surf(Y,X,squeeze(meanofdiff_v3(:,19:54,:))','EdgeColor','None'); view(2);
        set(gca,'YScale','log'); title('Mean of Diff Pre v Post 3')
        yticks([0 10 20 30 40 50 60 80 100])
        xlim([-200 1080])
        clim = [clim; get(gca,'clim')];
        
        dif4Fig = subplot(3,4,12);
        surf(Y,X,squeeze(meanofdiff_v4(:,19:54,:))','EdgeColor','None'); view(2);
        set(gca,'YScale','log'); title('Mean of Diff Pre v Post 4')
        yticks([0 10 20 30 40 50 60 80 100])
        xlim([-200 1080])
        clim = [clim; get(gca,'clim')];
        
        newC = [min(clim(:)) max(clim(:))];
        
        set(Po1Fig,'Clim',newC);
        set(Po2Fig,'Clim',newC);
        set(Po3Fig,'Clim',newC);
        set(Po4Fig,'Clim',newC);
        set(dif1Fig,'Clim',newC);
        set(dif2Fig,'Clim',newC);
        set(dif3Fig,'Clim',newC);
        set(dif4Fig,'Clim',newC);
        set(preFig,'Clim',newC);colorbar;
        
        h = gcf;
        h.Renderer = 'Painters';
        set(h, 'PaperType', 'A4');
        set(h, 'PaperOrientation', 'landscape');
        set(h, 'PaperUnits', 'centimeters');
        h_title = ['aroundLaser Spectral Power of ' type ' ' layer ' ' num2str(stimfreq(iStim)) ' ' group{iGro}];
        sgtitle(h_title)
        savefig(h_title)
        saveas(gcf, [h_title '.png'])
        close(h)
        
        %% t fig
        
        [X,Y]=meshgrid(WT.freq{1}(19:54),params.startTime*1000:(params.limit-201));
        figure('Name','Observed t Values BF', ...
            'units','normalized','outerposition',[0 0 1 1]);
        colormap('winter')
        
        subplot(241);
        surf(Y,X,squeeze(obs_t_v1(:,19:54,:))','EdgeColor','None'); view(2);
        set(gca,'YScale','log'); title('t Value Mat')
        yticks([0 10 20 30 40 50 60 80 100])
        xlim([-200 1080])
        colorbar
        
        subplot(242);
        surf(Y,X,squeeze(obs_t_v2(:,19:54,:))','EdgeColor','None'); view(2);
        set(gca,'YScale','log'); title('t Value Mat')
        yticks([0 10 20 30 40 50 60 80 100])
        xlim([-200 1080])
        colorbar
        
        subplot(243);
        surf(Y,X,squeeze(obs_t_v3(:,19:54,:))','EdgeColor','None'); view(2);
        set(gca,'YScale','log'); title('t Value Mat')
        yticks([0 10 20 30 40 50 60 80 100])
        xlim([-200 1080])
        colorbar
        
        subplot(244);
        surf(Y,X,squeeze(obs_t_v4(:,19:54,:))','EdgeColor','None'); view(2);
        set(gca,'YScale','log'); title('t Value Mat')
        yticks([0 10 20 30 40 50 60 80 100])
        xlim([-200 1080])
        colorbar
        
        subplot(245);
        surf(Y,X,squeeze(obs_clusters_v1(:,19:54,:))','EdgeColor','None'); view(2);
        set(gca,'YScale','log'); title('Clusters where p<0.05')
        yticks([0 10 20 30 40 50 60 80 100])
        xlim([-200 1080])
        colorbar
        
        subplot(246);
        surf(Y,X,squeeze(obs_clusters_v2(:,19:54,:))','EdgeColor','None'); view(2);
        set(gca,'YScale','log'); title('Clusters where p<0.05')
        yticks([0 10 20 30 40 50 60 80 100])
        xlim([-200 1080])
        colorbar
        
        subplot(247);
        surf(Y,X,squeeze(obs_clusters_v3(:,19:54,:))','EdgeColor','None'); view(2);
        set(gca,'YScale','log'); title('Clusters where p<0.05')
        yticks([0 10 20 30 40 50 60 80 100])
        xlim([-200 1080])
        colorbar
        
        subplot(248);
        surf(Y,X,squeeze(obs_clusters_v4(:,19:54,:))','EdgeColor','None'); view(2);
        set(gca,'YScale','log'); title('Clusters where p<0.05')
        yticks([0 10 20 30 40 50 60 80 100])
        xlim([-200 1080])
        colorbar
        
        h = gcf;
        h.Renderer = 'Painters';
        set(h, 'PaperType', 'A4');
        set(h, 'PaperOrientation', 'landscape');
        set(h, 'PaperUnits', 'centimeters');
        h_title = ['aroundLaser Observed t and p ' type ' ' layer ' ' num2str(stimfreq(iStim)) ' ' group{iGro} ];
        sgtitle(h_title)
        savefig(h_title)
        saveas(gcf, [h_title '.png'])
        close(h)
        

        %% Permutation Step 3 - do the permute
        % preallocate
        mass_clustermass_v1 = NaN([1 nperms]);
        mass_clustermass_v2 = NaN([1 nperms]);
        mass_clustermass_v3 = NaN([1 nperms]);
        mass_clustermass_v4 = NaN([1 nperms]);
        
        perm_layer_v1 = struct;
        perm_layer_v2 = struct;
        perm_layer_v3 = struct;
        perm_layer_v4 = struct;
        for iOsc = 1:length(osciName)
            perm_layer_v1.(osciName{iOsc}) = NaN([1 nperms]);
            perm_layer_v2.(osciName{iOsc}) = NaN([1 nperms]);
            perm_layer_v3.(osciName{iOsc}) = NaN([1 nperms]);
            perm_layer_v4.(osciName{iOsc}) = NaN([1 nperms]);
        end
        
        for iComp = 1:4 % specify comparison and data going in
            if iComp == 1
                contAll = vertcat(preL,postL1);
            elseif iComp == 2
                contAll = vertcat(preL,postL2);
            elseif iComp == 3
                contAll = vertcat(preL,postL3);
            elseif iComp == 4
                if ~exist('KIT04CasepreL','var') % KIT04 missing AM_4 
                    contAll = vertcat(preL,postL4);
                else 
                    contAll = vertcat(KIT04CasepreL,postL4);
                end
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
                    mass_clustermass_v1(iperm) = per_clustermass;
                elseif iComp == 2
                    mass_clustermass_v2(iperm) = per_clustermass;
                elseif iComp == 3
                    mass_clustermass_v3(iperm) = per_clustermass;
                elseif iComp == 4
                    mass_clustermass_v4(iperm) = per_clustermass;
                end
                
                for iOsc = 1:length(osciName)
                    hold_permlayer = per_clusters(:,osciRows{iOsc},:);
                    hold_permlayer = nansum(nansum(hold_permlayer));
                    if iComp == 1
                        perm_layer_v1.(osciName{iOsc})(iperm) = hold_permlayer;
                    elseif iComp == 2
                        perm_layer_v2.(osciName{iOsc})(iperm) = hold_permlayer;
                    elseif iComp == 3
                        perm_layer_v3.(osciName{iOsc})(iperm) = hold_permlayer;
                    elseif iComp == 4
                        perm_layer_v4.(osciName{iOsc})(iperm) = hold_permlayer;
                    end
                end
            end
        end
        
        cd(homedir); cd DATA; cd Spectral; mkdir(['Crypt_MagPerm_' type]); cd(['Crypt_MagPerm_' type]);
        %% Check Significance of full clustermass
        
        figure('Name',['aroundLaser Obs cluster vs Perm ' type ' ' layer ' clickfreq ' num2str(stimfreq(iStim)) ' ' group{iGro}]);
        % In how many instances is the clustermass of the permutation MORE than
        % the observed clustermass
        sig_mass_v1 = sum(mass_clustermass_v1>obs_clustermass_v1,2); %compare how many instances the sum of the permutation is greater than observed
        pVal_v1 = sig_mass_v1/nperms;
        permMean_v1 = mean(mass_clustermass_v1);
        permSTD_v1 = std(mass_clustermass_v1);
        sig_mass_v1(sig_mass_v1 <= pthresh) = true; %if it's not more than 50, it's siginificant for that frequency
        sig_mass_v1(sig_mass_v1 > pthresh) = false;
        
        subplot(141)
        boxplot(mass_clustermass_v1); hold on;
        title([measurement{1} 'v' measurement{2}])
        if sig_mass_v1 == 0
            plot(1,obs_clustermass_v1,'ro','LineWidth',4)
        else
            plot(1,obs_clustermass_v1,'go','LineWidth',4)
        end
        legend(['p = ' num2str(pVal_v1)])
        
        sig_mass_v2 = sum(mass_clustermass_v2>obs_clustermass_v2,2); %compare how many instances the sum of the permutation is greater than observed
        pVal_v2 = sig_mass_v2/nperms;
        permMean_v2 = mean(mass_clustermass_v2);
        permSTD_v2 = std(mass_clustermass_v2);
        sig_mass_v2(sig_mass_v2 <= pthresh) = true; %if it's not more than 50, it's siginificant for that frequency
        sig_mass_v2(sig_mass_v2 > pthresh) = false;
        
        subplot(142)
        boxplot(mass_clustermass_v2); hold on;
        title([measurement{1} 'v' measurement{3}])
        if sig_mass_v2 == 0
            plot(1,obs_clustermass_v2,'ro','LineWidth',4)
        else
            plot(1,obs_clustermass_v2,'go','LineWidth',4)
        end
        legend(['p = ' num2str(pVal_v2)])
        
        sig_mass_v3 = sum(mass_clustermass_v3>obs_clustermass_v3,2); %compare how many instances the sum of the permutation is greater than observed
        pVal_v3 = sig_mass_v3/nperms;
        permMean_v3 = mean(mass_clustermass_v3);
        permSTD_v3 = std(mass_clustermass_v3);
        sig_mass_v3(sig_mass_v3 <= pthresh) = true; %if it's not more than 50, it's siginificant for that frequency
        sig_mass_v3(sig_mass_v3 > pthresh) = false;
        
        subplot(143)
        boxplot(mass_clustermass_v3); hold on;
        title([measurement{1} 'v' measurement{4}])
        if sig_mass_v3 == 0
            plot(1,obs_clustermass_v3,'ro','LineWidth',4)
        else
            plot(1,obs_clustermass_v3,'go','LineWidth',4)
        end
        legend(['p = ' num2str(pVal_v3)])
        
        sig_mass_v4 = sum(mass_clustermass_v4>obs_clustermass_v4,2); %compare how many instances the sum of the permutation is greater than observed
        pVal_v4 = sig_mass_v4/nperms;
        permMean_v4 = mean(mass_clustermass_v4);
        permSTD_v4 = std(mass_clustermass_v4);
        sig_mass_v4(sig_mass_v4 <= pthresh) = true; %if it's not more than 50, it's siginificant for that frequency
        sig_mass_v4(sig_mass_v4 > pthresh) = false;
        
        subplot(144)
        boxplot(mass_clustermass_v4); hold on;
        title([measurement{1} 'v' measurement{5}])
        if sig_mass_v4 == 0
            plot(1,obs_clustermass_v4,'ro','LineWidth',4)
        else
            plot(1,obs_clustermass_v4,'go','LineWidth',4)
        end
        legend(['p = ' num2str(pVal_v4)])
        
        h = gcf;
        h.Renderer = 'Painters';
        set(h, 'PaperType', 'A4');
        set(h, 'PaperOrientation', 'landscape');
        set(h, 'PaperUnits', 'centimeters');
        h_title = ['aroundLaser Perm Full ' type ' ' layer ' ' num2str(stimfreq(iStim)) ' ' group{iGro}];
        sgtitle(h_title)
        savefig(h_title)
        saveas(gcf, [h_title '.png'])
        close(h)
        save([h_title '.mat'],...
            'pVal_v1','permMean_v1','permSTD_v1','pVal_v2','permMean_v2','permSTD_v2',...
            'pVal_v3','permMean_v3','permSTD_v3','pVal_v4','permMean_v4','permSTD_v4')
        
        %% Check Significance of layers' clustermass
        
        figure('Name',['aroundLaser Obs vs Perm ' type ' ' layer ' clickfreq ' num2str(stimfreq(iStim)) ' ' ...
            group{iGro} ' ocsillations'],'units','normalized','outerposition',[0 0 1 1])
        pVal_v1 = struct; permMean_v1 = struct; permSTD_v1 = struct;
        pVal_v2 = struct; permMean_v2 = struct; permSTD_v2 = struct;
        pVal_v3 = struct; permMean_v3 = struct; permSTD_v3 = struct;
        pVal_v4 = struct; permMean_v4 = struct; permSTD_v4 = struct;
        
        for iOsc = 1:length(osciName)
            % In how many instances is the clustermass of the permutation MORE than
            % the observed clustermass
            sig_mass_v1 = sum(perm_layer_v1.(osciName{iOsc})>obs_layer_v1.(osciName{iOsc}),2); %compare how many instances the sum of the permutation is greater than observed
            pVal_v1.(osciName{iOsc}) = sig_mass_v1/nperms;
            permMean_v1.(osciName{iOsc}) = mean(perm_layer_v1.(osciName{iOsc}));
            permSTD_v1.(osciName{iOsc}) = std(perm_layer_v1.(osciName{iOsc}));
            sig_mass_v1(sig_mass_v1 <= pthresh) = true; %if it's not more than 50, it's siginificant for that frequency
            sig_mass_v1(sig_mass_v1 > pthresh) = false;
            
            subplot(4,6,iOsc)
            boxplot(perm_layer_v1.(osciName{iOsc})); hold on;
            title([measurement{1} 'v' measurement{2} ' ' osciName{iOsc}])
            if sig_mass_v1 == 0
                plot(1,obs_layer_v1.(osciName{iOsc}),'ro','LineWidth',4)
            else
                plot(1,obs_layer_v1.(osciName{iOsc}),'go','LineWidth',4)
            end
            legend(['p = ' num2str(pVal_v1.(osciName{iOsc}))])
            
            sig_mass_v2 = sum(perm_layer_v2.(osciName{iOsc})>obs_layer_v2.(osciName{iOsc}),2); %compare how many instances the sum of the permutation is greater than observed
            pVal_v2.(osciName{iOsc}) = sig_mass_v2/nperms;
            permMean_v2.(osciName{iOsc}) = mean(perm_layer_v2.(osciName{iOsc}));
            permSTD_v2.(osciName{iOsc}) = std(perm_layer_v2.(osciName{iOsc}));
            sig_mass_v2(sig_mass_v2 <= pthresh) = true; %if it's not more than 50, it's siginificant for that frequency
            sig_mass_v2(sig_mass_v2 > pthresh) = false;
            
            subplot(4,6,iOsc+6)
            boxplot(perm_layer_v2.(osciName{iOsc})); hold on;
            title([measurement{1} 'v' measurement{3} ' ' osciName{iOsc}])
            if sig_mass_v2 == 0
                plot(1,obs_layer_v2.(osciName{iOsc}),'ro','LineWidth',4)
            else
                plot(1,obs_layer_v2.(osciName{iOsc}),'go','LineWidth',4)
            end
            legend(['p = ' num2str(pVal_v2.(osciName{iOsc}))])
            
            sig_mass_v3 = sum(perm_layer_v3.(osciName{iOsc})>obs_layer_v3.(osciName{iOsc}),2); %compare how many instances the sum of the permutation is greater than observed
            pVal_v3.(osciName{iOsc}) = sig_mass_v3/nperms;
            permMean_v3.(osciName{iOsc}) = mean(perm_layer_v3.(osciName{iOsc}));
            permSTD_v3.(osciName{iOsc}) = std(perm_layer_v3.(osciName{iOsc}));
            sig_mass_v3(sig_mass_v3 <= pthresh) = true; %if it's not more than 50, it's siginificant for that frequency
            sig_mass_v3(sig_mass_v3 > pthresh) = false;
            
            subplot(4,6,iOsc+12)
            boxplot(perm_layer_v3.(osciName{iOsc})); hold on;
            title([measurement{1} 'v' measurement{4} ' ' osciName{iOsc}])
            if sig_mass_v3 == 0
                plot(1,obs_layer_v3.(osciName{iOsc}),'ro','LineWidth',4)
            else
                plot(1,obs_layer_v3.(osciName{iOsc}),'go','LineWidth',4)
            end
            legend(['p = ' num2str(pVal_v3.(osciName{iOsc}))])
            
            sig_mass_v4 = sum(perm_layer_v4.(osciName{iOsc})>obs_layer_v4.(osciName{iOsc}),2); %compare how many instances the sum of the permutation is greater than observed
            pVal_v4.(osciName{iOsc}) = sig_mass_v4/nperms;
            permMean_v4.(osciName{iOsc}) = mean(perm_layer_v4.(osciName{iOsc}));
            permSTD_v4.(osciName{iOsc}) = std(perm_layer_v4.(osciName{iOsc}));
            sig_mass_v4(sig_mass_v4 <= pthresh) = true; %if it's not more than 50, it's siginificant for that frequency
            sig_mass_v4(sig_mass_v4 > pthresh) = false;
            
            subplot(4,6,iOsc+18)
            boxplot(perm_layer_v4.(osciName{iOsc})); hold on;
            title([measurement{1} 'v' measurement{5} ' ' osciName{iOsc}])
            if sig_mass_v4 == 0
                plot(1,obs_layer_v4.(osciName{iOsc}),'ro','LineWidth',4)
            else
                plot(1,obs_layer_v4.(osciName{iOsc}),'go','LineWidth',4)
            end
            legend(['p = ' num2str(pVal_v4.(osciName{iOsc}))])
        end
        
        h = gcf;
        h.Renderer = 'Painters';
        set(h, 'PaperType', 'A4');
        set(h, 'PaperOrientation', 'landscape');
        set(h, 'PaperUnits', 'centimeters');
        h_title = ['aroundLaser Perm Osci ' type ' ' layer ' ' osciName{iOsc} ' ' num2str(stimfreq(iStim)) ' ' group{iGro}];
        savefig(h_title)
        sgtitle(h_title)
        saveas(gcf, [h_title '.png'])
        close(h)
        save([h_title '.mat'],...
            'pVal_v1','permMean_v1','permSTD_v1','pVal_v2','permMean_v2','permSTD_v2',...
            'pVal_v3','permMean_v3','permSTD_v3','pVal_v4','permMean_v4','permSTD_v4')
    end
end
