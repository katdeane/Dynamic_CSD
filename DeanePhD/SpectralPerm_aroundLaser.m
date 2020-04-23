function SpectralPerm_aroundLaser(layer,homedir)

% Input:    Layer to analyze, (possible input: 2,5,10,20,40)
%           Needs WT_CL.mat from /Data/Spectral
% Output:   Figures for means and observed difference of group
%           comparison; figures for observed t values, clusters
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
if ~exist('layer','var')
    layer = 'I_II'; % run the granular layer if not specified
end

cd (homedir),cd DATA;

nperms = 1000;
pthresh = nperms*(0.05); %if bonferroni wanted, *(0.05/7)

% frequencies can be found in wtTable.freq{1} to clarify the following
% rows choices; actual intended rows commented
theta = (49:54);        %(4:7);
alpha = (44:48);        %(8:12);
beta_low = (39:43);     %(13:18);
beta_high = (34:38);    %(19:30);
gamma_low = (26:33);    %(31:60);
gamma_high = (19:25);   %(61:100);

osciName    = {'theta' 'alpha' 'beta_low' 'beta_high' 'gamma_low' 'gamma_high'};
osciRows    = {theta alpha beta_low beta_high gamma_low gamma_high};
measurement = {'preCL_1','CL_1'};
group       = {'KIC','KIT','KIV'};
stimfreq    = [2,5,10,20,40];

%% Load in and seperate Data
load('WT_CL.mat', 'WT_CL')

% varNames = unique(wtTable.layer);
params.startTime = -0.2; % seconds
params.limit = 1300;

% check if layer or full mat needed
if ~strcmp(layer, 'ALL')
    %Pull out layer
    WT_CL = WT_CL(startsWith(WT_CL.layer,layer) & endsWith(WT_CL.layer,layer),:);
end

for iGro = 1:length(group)
    for iStim = 1:length(stimfreq)
        %Pull out conditions and limit them to 1300ms
        Measure1 = table2cell(WT_CL(contains(WT_CL.group,group{iGro}) & ...
            WT_CL.stimulus==stimfreq(iStim) & startsWith(WT_CL.measurement,measurement{1}),1));
        Measure1 =  cellfun(@(x) x(:,1:params.limit),Measure1,'UniformOutput',false);
        
        Measure2 = table2cell(WT_CL(contains(WT_CL.group,group{iGro}) & ...
            WT_CL.stimulus==stimfreq(iStim) & startsWith(WT_CL.measurement,measurement{2}),1));
        Measure2 =  cellfun(@(x) x(:,1:params.limit),Measure2,'UniformOutput',false);
        
        %Stack the individual animals' data (animal#x54x600)
        M1 = nan(size(Measure1,1),size(Measure1{1,1},1),size(Measure1{1,1},2));
        M2 = nan(size(Measure2,1),size(Measure2{1,1},1),size(Measure2{1,1},2));
        
        for ii = 1:length(Measure1)
            M1(ii,:,:)= Measure1{ii};
        end
        M1 = abs(M1);
        
        for ii = 1:length(Measure2)
            M2(ii,:,:)= Measure2{ii};
        end
        M2 = abs(M2);
        
        if ~strcmp(layer, 'ALL')
            grpsize = length(Measure1);
        else
            keyboard; % need to verify if this step is necessary/correct
            grpsize = length(Measure1)/7;
        end
        
        if length(Measure1) ~= length(Measure2)
            keyboard %decide what to do where measurments are missing
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
        
        obs1_mean = nanmean(M1,1);
        obs2_mean = nanmean(M2,1);
        
        meanofdiff = nanmean(M1-M2,1);
        stdofdiff  = nanstd(M1-M2,0,1);
        
        %% Permutation Step 2 - t test
        % find the t values along all data points for each frequency bin
        obs_t = meanofdiff./(stdofdiff./sqrt(grpsize));
        obs_clusters = abs(obs_t);
        obs_clusters(obs_clusters < Tthresh) = NaN;
        obs_clusters(obs_clusters >= Tthresh) = 1;
        
        % check cluster mass KD: this will need to be improved to include a shorter
        % time window for detection (after each stimulus click)
        obs_clustermass = nansum(nansum(obs_clusters));
        
        %% Oscillation frequency bands:
        % % pull out clusters
        
        obs_layer = struct;
        for iOsc = 1:length(osciName)
            obs_layer.(osciName{iOsc}) = obs_clusters(:,osciRows{iOsc},:);
            obs_layer.(osciName{iOsc}) = nansum(nansum(obs_layer.(osciName{iOsc})));
        end
        
        cd(homedir); cd figs; mkdir('Crypt_MagPerm'); cd('Crypt_MagPerm');
        %% dif figs
        [X,Y]=meshgrid(WT_CL.freq{1}(19:54),params.startTime*1000:(params.limit-201));
        figure('Name','Observed Spectral Power and Differences',...
            'units','normalized','outerposition',[0 0 1 1]);
        
        % observed values
        obs1Fig = subplot(231);
        surf(Y,X,squeeze(meanofdiff(:,19:54,:))','EdgeColor','None'); view(2);
        set(gca,'YScale','log'); title(measurement{1})
        yticks([0 10 20 30 40 50 60 80 100])
        clim = get(gca,'clim');
        
        obs2Fig = subplot(232);
        surf(Y,X,squeeze(obs2_mean(:,19:54,:))','EdgeColor','None'); view(2);
        set(gca,'YScale','log'); title(measurement{2})
        yticks([0 10 20 30 40 50 60 80 100])
        clim = [clim; get(gca,'clim')];
        
        % observed difference
        diffFig = subplot(233);
        surf(Y,X,squeeze(meanofdiff(:,19:54,:))','EdgeColor','None'); view(2);
        set(gca,'YScale','log'); title('Mean of Difference')
        yticks([0 10 20 30 40 50 60 80 100 200 300 500])
        clim = [clim; get(gca,'clim')];
        
        newC = [min(clim(:)) max(clim(:))];
        
        set(obs1Fig,'Clim',newC);colorbar;
        set(obs2Fig,'Clim',newC);colorbar;
        set(diffFig,'Clim',newC);colorbar;
        
        subplot(234);
        surf(Y,X,squeeze(obs_t(:,19:54,:))','EdgeColor','None'); view(2);
        set(gca,'YScale','log'); title('t Value Mat')
        yticks([0 10 20 30 40 50 60 80 100])
        colorbar
        
        subplot(235);
        surf(Y,X,squeeze(obs_clusters(:,19:54,:))','EdgeColor','None'); view(2);
        set(gca,'YScale','log'); title('Clusters where p<0.05')
        yticks([0 10 20 30 40 50 60 80 100])
        colorbar
        
        h = gcf;
        h.Renderer = 'Painters';
        set(h, 'PaperType', 'A4');
        set(h, 'PaperOrientation', 'landscape');
        set(h, 'PaperUnits', 'centimeters');
        savefig(['aroundLaser Spectral Power of ' layer ' clickfrq ' num2str(stimfreq(iStim)) ' ' group{iGro}])
        saveas(gcf, ['aroundLaser Spectral Power of ' layer ' clickfrq ' num2str(stimfreq(iStim)) ' ' group{iGro} '.png'])
        close(h)
        
        
        %% Permutation Step 3 - do the permute
        % preallocate
        mass_clustermass = NaN([1 nperms]);
        
        perm_layer = struct;
        for iOsc = 1:length(osciName)
            perm_layer.(osciName{iOsc}) = NaN([1 nperms]);
        end
        
        contAll = vertcat(Measure1,Measure2);
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
            
            perm1_mean = nanmean(M1Perm,1);
            perm2_mean = nanmean(M2Perm,1);
            
            perm_meanofdiff = nanmean(M1Perm-M2Perm,1);
            perm_stdofdiff  = nanstd(M1Perm-M2Perm,0,1);
            
            % t test %%%
            per_t = perm_meanofdiff./(perm_stdofdiff./sqrt(grpsize));
            per_clusters = abs(per_t);
            per_clusters(per_clusters < Tthresh) = NaN;
            per_clusters(per_clusters >= Tthresh) = 1;
            
            % check cluster mass for 300 ms from tone onset
            per_clustermass = nansum(nansum(per_clusters));
            
            mass_clustermass(iperm) = per_clustermass;
            for iOsc = 1:length(osciName)
                hold_permlayer = per_clusters(:,osciRows{iOsc},:);
                hold_permlayer = nansum(nansum(hold_permlayer));
                perm_layer.(osciName{iOsc})(iperm) = hold_permlayer;
            end
            
        end
        
        cd(homedir); cd DATA; cd Spectral; mkdir('Crypt_MagPerm'); cd('Crypt_MagPerm');
        %% Check Significance of full clustermass
        
        % In how many instances is the clustermass of the permutation MORE than
        % the observed clustermass
        sig_mass = sum(mass_clustermass>obs_clustermass,2); %compare how many instances the sum of the permutation is greater than observed
        pVal = sig_mass/nperms;
        permMean = mean(mass_clustermass);
        permSTD = std(mass_clustermass);
        sig_mass(sig_mass <= pthresh) = true; %if it's not more than 50, it's siginificant for that frequency
        sig_mass(sig_mass > pthresh) = false;
        
        figure('Name',['aroundLaser Obs cluster vs Perm ' layer ' clickfreq ' num2str(stimfreq(iStim)) ' ' group{iGro}]);
        boxplot(mass_clustermass); hold on;
        title([measurement{1} 'v' measurement{1}])
        if sig_mass == 0
            plot(1,obs_clustermass,'ro','LineWidth',4)
        else
            plot(1,obs_clustermass,'go','LineWidth',4)
        end
        legend(['p = ' num2str(pVal)])
        
        h = gcf;
        h.Renderer = 'Painters';
        set(h, 'PaperType', 'A4');
        set(h, 'PaperOrientation', 'landscape');
        set(h, 'PaperUnits', 'centimeters');
        savefig(['aroundLaser Perm Full ' layer ' clickfreq ' num2str(stimfreq(iStim)) ' ' group{iGro}])
        saveas(gcf, ['aroundLaser Perm Full ' layer ' clickfreq ' num2str(stimfreq(iStim)) ' ' group{iGro} '.png'])
        close(h)
        save(['aroundLaser Perm Full ' layer ' clickfreq ' num2str(stimfreq(iStim)) ' ' group{iGro} '.mat'],...
            'pVal','permMean','permSTD')
        
        %% Check Significance of layers' clustermass
        
        figure('Name',['aroundLaser Obs vs Perm ' layer ' clickfreq ' num2str(stimfreq(iStim)) ' ' ...
            group{iGro} ' ocsillations'],'units','normalized','outerposition',[0 0 1 1])
        pVal = struct; permMean = struct; permSTD = struct;
        
        for iOsc = 1:length(osciName)
            % In how many instances is the clustermass of the permutation MORE than
            % the observed clustermass
            sig_mass = sum(perm_layer.(osciName{iOsc})>obs_layer.(osciName{iOsc}),2); %compare how many instances the sum of the permutation is greater than observed
            pVal.(osciName{iOsc}) = sig_mass/nperms;
            permMean.(osciName{iOsc}) = mean(perm_layer.(osciName{iOsc}));
            permSTD.(osciName{iOsc}) = std(perm_layer.(osciName{iOsc}));
            sig_mass(sig_mass <= pthresh) = true; %if it's not more than 50, it's siginificant for that frequency
            sig_mass(sig_mass > pthresh) = false;
            
            subplot(1,6,iOsc)
            boxplot(perm_layer.(osciName{iOsc})); hold on;
            title([measurement{1} 'v' measurement{2} ' ' osciName{iOsc}])
            if sig_mass == 0
                plot(1,obs_layer.(osciName{iOsc}),'ro','LineWidth',4)
            else
                plot(1,obs_layer.(osciName{iOsc}),'go','LineWidth',4)
            end
            legend(['p = ' num2str(pVal.(osciName{iOsc}))])
        end
        
        h = gcf;
        h.Renderer = 'Painters';
        set(h, 'PaperType', 'A4');
        set(h, 'PaperOrientation', 'landscape');
        set(h, 'PaperUnits', 'centimeters');
        savefig(['aroundLaser Perm Osci ' layer ' ' osciName{iOsc} ' ' num2str(stimfreq(iStim)) ' ' group{iGro}])
        saveas(gcf, ['aroundLaser Perm Osci ' layer ' ' osciName{iOsc} ' ' num2str(stimfreq(iStim)) ' ' group{iGro} '.png'])
        close(h)
        save(['aroundLaser Perm Osci ' layer ' ' osciName{iOsc} ' ' num2str(stimfreq(iStim)) ' ' group{iGro} '.mat'],...
            'pVal','permMean','permSTD')
    end
end
