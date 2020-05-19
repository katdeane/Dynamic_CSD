function SpectralPerm_PhaseCoWithin(homedir,Meas1,Meas2)

% Input:    home directory, 2 measurements to compare, between groups
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
    elseif exist('D:\Dynamic_CSD','dir') == 7
        cd('D:\Dynamic_CSD')
    end
    
    homedir = pwd;
    addpath(genpath(homedir));
end

cd (homedir),cd DATA;

nperms = 1000;
pthresh = nperms*(0.05);

% frequencies can be found in oscifreq to clarify the following
% rows choices; actual intended rows commented
theta = (49:54);        %(4:7);
alpha = (44:48);        %(8:12);
beta_low = (39:43);     %(13:18);
beta_high = (34:38);    %(19:30);
gamma_low = (26:33);    %(31:60);
gamma_high = (19:25);   %(61:100);

osciName = {'theta' 'alpha' 'beta_low' 'beta_high' 'gamma_low' 'gamma_high'};
osciRows = {theta alpha beta_low beta_high gamma_low gamma_high};
layer    = {'I_II','IV','V','VI'};
stimfrq  = [2,5,10,20,40];

%% Load in the two groups for comparison

params.startTime = -0.2; % seconds
timelimit     = 1300;

for iLay = 1:length(layer)
    for iSti = 1:length(stimfrq)
        % this needs to be loaded in and cut per loop because the parfor
        % later requires a lot of memory
        disp(['Loading ' Meas1(1:end-4) ' and ' Meas2(1:end-4) ' for ' ...
            layer{iLay} ' ' num2str(stimfrq(iSti))]);
        tic
        load(Meas1,'WT_st');
        FullDat = WT_st; clear WT_st
        load(Meas2,'WT_st');
        FullDat = [FullDat;WT_st]; clear WT_st %#ok<AGROW>
        toc

        oscifreq   = FullDat.freq{1,1};
        
        FullDat = table2cell(FullDat(strcmp(FullDat.layer,layer{iLay})...
            & FullDat.stimulus == stimfrq(iSti),1:5));
        
        fullgroup = unique(FullDat(:,2),'stable');
        measure   = unique(FullDat(:,5),'stable');
        grpsize   = length(fullgroup);
        
        % Get phase coherence for groups
        Grp1All = nan(54,timelimit,grpsize);
        Grp2All = nan(54,timelimit,grpsize);
        
        % through each animal
        for iAn = 1: grpsize
            
            % pull out the data for this animal's first measurement
            boohold = (strcmp({FullDat{:,2}},fullgroup{iAn}) ...
                & strcmp({FullDat{:,5}},measure{1})); %#ok<*CCAT1>
            grp1 = {FullDat{boohold,1}}';
            for igrp = 1:length(grp1)
                grp1{igrp} = grp1{igrp}(:,1:timelimit);
            end
            % set up output cells for transformed data
            transGrp1 = cell(size(grp1));
            
            if isempty(grp1)
                continue
            end
            
            % take the phase z/abs(z) for each single trial
            for iTrial = 1:length(grp1)
                curTrial = grp1{iTrial};
                transGrp1{iTrial} = curTrial./abs(curTrial);
            end
            
            % for the subject:
            % get the mean of single trials and the absolute for that mean
            Grp1All(:,:,iAn) = abs(mean(cat(3,transGrp1{:}),3));
            
            
            % pull out the data for this animal's second measurement
            boohold = (strcmp({FullDat{:,2}},fullgroup{iAn}) ...
                & strcmp({FullDat{:,5}},measure{2})); %#ok<*CCAT1>
            grp2 = {FullDat{boohold,1}}';
            for igrp = 1:length(grp2)
                grp2{igrp} = grp2{igrp}(:,1:timelimit);
            end
            % set up output cells for transformed data
            transGrp2 = cell(size(grp2));
            
            if isempty(grp2)
                continue
            end
            
            % take the phase z/abs(z) for each single trial
            for iTrial = 1:length(grp2)
                curTrial = grp2{iTrial};
                transGrp2{iTrial} = curTrial./abs(curTrial);
            end
            
            % get the mean of single trials and the absolute for that mean
            Grp2All(:,:,iAn) = abs(mean(cat(3,transGrp2{:}),3));
            

        end
        
        clear grp1 grp2 transGrp1 transGrp2 curTrial
        
        % for the group:
        % take the mean across subjects
        grp1mean = nanmean(Grp1All,3);
        grp2mean = nanmean(Grp2All,3);
        
        diffmeans = grp2mean - grp1mean;
        
        %% Mann Whitney U Test (ranksum)
        
        % shift dimensions so that subject is first, then y and x axis
        % remove nan rows (animal KIT04 doesn't have layer II)
        shifted1 = shiftdim(Grp1All,2);
        whereNan = find(isnan(shifted1(:,1,1)));
        if ~isempty(whereNan)
            shifted1 = vertcat(shifted1(1:whereNan-1,:,:),shifted1(whereNan+1:end,:,:));
        end
        shifted2 = shiftdim(Grp2All,2);
        whereNan = find(isnan(shifted2(:,1,1)));
        if ~isempty(whereNan)
            shifted2 = vertcat(shifted2(1:whereNan-1,:,:),shifted2(whereNan+1:end,:,:));
        end
        
        stats = mwwtest(shifted1,shifted2);
        obs_Cmass = squeeze(stats.p{2});
        obs_Cmass(obs_Cmass <= 0.05) = true;
        obs_Cmass(obs_Cmass ~= 1) = false;
        
        % effect size of mwwtest is r = abs(z/sqrt(n1+n2)) / 0.1 is small, 0.3 is
        % medium, 0.5 is large
        obs_Esize = abs((squeeze(stats.Z))./sqrt(grpsize+grpsize));
        
        % full matrix
        obs_clustermass = nansum(nansum(obs_Cmass));
        
        % layer specific
        obs_layer = struct;
        
        for iOsc = 1:length(osciName)
            obs_layer.(osciName{iOsc}) = obs_Cmass(osciRows{iOsc},:);
            
            % % sum clusters (twice to get final value)
            for i = 1:2
                obs_layer.(osciName{iOsc}) = nansum(obs_layer.(osciName{iOsc}));
            end
        end
        
        clear stats shifted1 shifted2 Grp1All Grp2All
        %% cluster fig
        %% Fig
        cd(homedir); cd figs; mkdir('Crypt_PhCPerm'); cd('Crypt_PhCPerm');
        
        [X,Y]=meshgrid(oscifreq(19:54),params.startTime*1000:(timelimit-201));
        figure('Name',['Observed Phase Difference' Meas1(1:end-4)...
            ' vs ' Meas2(1:end-4)],'Position',[100 100 1065 700]);
        gr1Fig = subplot(231);
        surf(Y',X',grp1mean(19:54,:),'EdgeColor','None'); view(2);
        set(gca,'YScale','log'); title(Meas1(1:end-4))
        yticks([0 10 20 30 40 50 60 80 100])
        colorbar
        % clim = get(gca,'clim');
        
        gr2Fig = subplot(232);
        surf(Y',X',grp2mean(19:54,:),'EdgeColor','None'); view(2);
        set(gca,'YScale','log'); title(Meas2(1:end-4))
        yticks([4 7 12 18 30 60 100])
        colorbar
        clim = get(gca,'clim');
        % clim = [clim; get(gca,'clim')];
        
        diffFig = subplot(233);
        surf(Y',X',diffmeans(19:54,:),'EdgeColor','None'); view(2);
        set(gca,'YScale','log'); title('Observed Diff')
        yticks([4 7 12 18 30 60 100])
        colorbar
        clim = [clim; get(gca,'clim')];
        
        newC = [min(clim(:)) max(clim(:))];
        
        %figure(awakeFig);
        set(gr2Fig,'Clim',newC);colorbar;
        %figure(ketFig);
        set(gr1Fig,'Clim',newC);colorbar;
        %figure(muscFig);
        set(diffFig,'Clim',newC);colorbar;
        
        subplot(234);
        surf(Y,X,obs_Esize(19:54,:)','EdgeColor','None'); view(2);
        set(gca,'YScale','log'); title('Effect size mat')
        yticks([0 10 20 30 40 50 60 80 100])
        colorbar
        
        subplot(235);
        surf(Y,X,obs_Cmass(19:54,:)','EdgeColor','None'); view(2);
        set(gca,'YScale','log'); title('Clusters where p>0.05')
        yticks([0 10 20 30 40 50 60 80 100])
        
        h = gcf;
        h.Renderer = 'Painters';
        set(h, 'PaperType', 'A4');
        set(h, 'PaperOrientation', 'landscape');
        set(h, 'PaperUnits', 'centimeters');
        htitle = ['Observed Phase Difference ' Meas1(1:end-4)...
            ' vs ' Meas2(1:end-4) ' layer ' layer{iLay} ' clickfrq ' num2str(stimfrq(iSti))];
        sgtitle(htitle)
        savefig(htitle)
        saveas(h, [htitle '.png'])
        close(h)

        clear obs_Cmass obs_Esize diffmeans grp1mean grp2mean X Y ...
            gr1Fig gr2Fig diffFig newC boohold
        %% Permutation Step 3 - do the permute
        
        % set up some containers
        mass_clustermass = NaN([1 nperms]);
        
        Theta  = NaN([1 nperms]);
        Alpha  = NaN([1 nperms]);
        BetaL  = NaN([1 nperms]);
        BetaH  = NaN([1 nperms]);
        GammaL = NaN([1 nperms]);
        GammaH = NaN([1 nperms]);
        
        disp('Permuting...')
        tic
        % echo loop to find observed values
        parfor iperm = 1:nperms
            
            % randomize the order of measurement selection
            order = randi(2,1,grpsize); % for the second "group", take the opposite
            
            
            % Get phase coherence for groups
            permbox1 = nan(54,1300,grpsize);
            permbox2 = nan(54,1300,grpsize);
                        
            % through each animal
            for iAn = 1: grpsize
                
                % pull out the data for the random animal
                boohold = (strcmp({FullDat{:,2}},fullgroup{iAn}) ...
                    & strcmp({FullDat{:,5}},measure{order(iAn)})); %#ok<*CCAT1>
                perm1 = {FullDat{boohold,1}}';
                for igrp = 1:length(perm1)
                    perm1{igrp} = perm1{igrp}(:,1:1300);
                end
                % set up output cells for transformed data
                transPerm1 = cell(size(perm1));
                
                if isempty(perm1)
                    continue
                end
                
                % take the phase z/abs(z) for each single trial
                for iTrial = 1:length(perm1)
                    curTrial = perm1{iTrial};
                    transPerm1{iTrial} = curTrial./abs(curTrial);
                end
                
                % for the subject:
                % get the mean of single trials and the absolute for that mean
                curAn = abs(mean(cat(3,transPerm1{:}),3));
                permbox1(:,:,iAn) = curAn;
                
                % get the other measurement for this animal.
                if order(iAn) == 1
                    OtherMeas = 2
                elseif order(iAn) == 2
                    OtherMeas = 1
                end
                
                boohold = (strcmp({FullDat{:,2}},fullgroup{iAn}) ...
                    & strcmp({FullDat{:,5}},measure{OtherMeas}));
                perm2 = {FullDat{boohold,1}}';
                for igrp = 1:length(perm2)
                    perm2{igrp} = perm2{igrp}(:,1:1300);
                end
                % set up output cells for transformed data
                transPerm2 = cell(size(perm2));
                
                if isempty(perm2)
                    continue
                end
                
                % take the phase z/abs(z) for each single trial
                for iTrial = 1:length(perm2)
                    curTrial = perm2{iTrial};
                    transPerm2{iTrial} = curTrial./abs(curTrial);
                end
                
                % get the mean of single trials and the absolute for that mean
                curAn = abs(mean(cat(3,transPerm2{:}),3));
                permbox2(:,:,iAn) = curAn;
                

            end
            
            % Mann Whitney U Test (ranksum) %%
            
            % shift dimensions so that subject is first, then y and x axis
            % remove nan rows (animal KIT04 doesn't have layer II)
            shifted1 = shiftdim(permbox1,2);
            whereNan = find(isnan(shifted1(:,1,1)));
            if ~isempty(whereNan)
                shifted1 = vertcat(shifted1(1:whereNan-1,:,:),shifted1(whereNan+1:end,:,:));
            end
            shifted2 = shiftdim(permbox2,2);
            whereNan = find(isnan(shifted2(:,1,1)));
            if ~isempty(whereNan)
                shifted2 = vertcat(shifted2(1:whereNan-1,:,:),shifted2(whereNan+1:end,:,:));
            end
            
            stats = mwwtest(shifted1,shifted2);
            perm_Cmass = squeeze(stats.p{2});
            perm_Cmass(perm_Cmass <= 0.05) = true;
            perm_Cmass(perm_Cmass ~= 1) = false;
            
            % check cluster mass for 300 ms from tone onset
            per_clustermass = nansum(nansum(perm_Cmass));
            mass_clustermass(iperm) = per_clustermass;
            
            % for layer specific: %%%
            % % pull out clusters
            Theta(iperm)  = nansum(nansum(perm_Cmass(osciRows{1},:)));
            Alpha(iperm)  = nansum(nansum(perm_Cmass(osciRows{2},:)));
            BetaL(iperm)  = nansum(nansum(perm_Cmass(osciRows{3},:)));
            BetaH(iperm)  = nansum(nansum(perm_Cmass(osciRows{4},:)));
            GammaL(iperm) = nansum(nansum(perm_Cmass(osciRows{5},:)));
            GammaH(iperm) = nansum(nansum(perm_Cmass(osciRows{6},:)));
            
        end
        toc
        
        perm_layer.theta      = Theta;
        perm_layer.alpha      = Alpha;
        perm_layer.beta_low   = BetaL;
        perm_layer.beta_high  = BetaH;
        perm_layer.gamma_low  = GammaL;
        perm_layer.gamma_high = GammaH;
        
        cd(homedir); cd DATA; cd Spectral; mkdir('Crypt_PhCPerm'); cd('Crypt_PhCPerm');
        %% Check Significance of full clustermass
        
        figure('Name',['Obs cluster vs Perm Phase Co of ' Meas1 ' vs ' ...
            Meas2 ' ' layer{iLay} ' clickfreq ' num2str(stimfrq(iSti))]);
        % In how many instances is the clustermass of the permutation MORE than
        % the observed clustermass
        sig_mass = sum(mass_clustermass>obs_clustermass,2); %compare how many instances the sum of the permutation is greater than observed
        pVal = sig_mass/nperms;
        permMean = mean(mass_clustermass);
        permSTD = std(mass_clustermass);
        sig_mass(sig_mass <= pthresh) = true; %if it's not more than 50, it's siginificant for that frequency
        sig_mass(sig_mass > pthresh) = false;

        boxplot(mass_clustermass); hold on;
        title([Meas1(1:end-4) ' v ' Meas2(1:end-4) ' of ' layer{iLay} ' ' num2str(stimfrq(iSti))])
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
        htitle = ['Obs cluster vs Perm Phase Co of ' Meas1(1:end-4) ' vs ' ...
            Meas2(1:end-4) ' ' layer{iLay} ' clickfreq ' num2str(stimfrq(iSti))];
        savefig(htitle)
        saveas(h, [htitle '.png'])
        close(h)
        save(['Permutation' Meas1(1:end-4) ' vs ' Meas2(1:end-4) ' ' ...
            layer{iLay} ' clickfreq ' num2str(stimfrq(iSti)) ' Full.mat'],...
            'pVal','permMean','permSTD')
        
        %% Check Significance of layers' clustermass
        
        figure('Name',['Obs cluster vs Perm Phase Co of ' Meas1 ' vs ' ...
            Meas2 ' ' layer{iLay} ' clickfreq ' num2str(stimfrq(iSti)) ...
            ' ocsillations'],'units','normalized','outerposition',[0 0 1 1])
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
            title(osciName{iOsc})
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
        htitle = ['Obs cluster vs Perm Phase Co of ' Meas1(1:end-4) ' vs ' ...
            Meas2(1:end-4) ' ' layer{iLay} ' clickfreq ' num2str(stimfrq(iSti)) ' oscillations'];
        sgtitle(htitle)
        savefig(htitle)
        saveas(h, [htitle '.png'])
        close(h)
        save(['Permutation' Meas1(1:end-4) ' vs ' Meas2(1:end-4) ' ' ...
            layer{iLay} ' clickfreq ' num2str(stimfrq(iSti)) ' Oscillations.mat'],...
            'pVal','permMean','permSTD')
        
        clear FullDat
    end %stim frequency
end %layer