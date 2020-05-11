function SpectralPerm_PhaseCoherence(homedir,Meas1,Meas2)

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
layer    = {'II','IV','V','VI'};
stimfrq  = [2,5,10,20,40];

%% Load in the two groups for comparison
disp(['Loading ' Meas1(1:end-4) ' and ' Meas2(1:end-4)]); 
tic
load(Meas1,'WT_st');
M1Dat = WT_st; clear WT_st
load(Meas2,'WT_st');
M2Dat = WT_st; clear WT_st
toc

grp1Name  = unique(M1Dat.animal,'stable');
grp2Name  = unique(M2Dat.animal,'stable');
grpsize1  = length(grp1Name);
grpsize2  = length(grp2Name);
fullgroup = vertcat(grp1Name,grp2Name);

params.startTime = -0.2; % seconds
params.limit     = 1300;

% Get phase coherence for groups
Grp1All = nan(54,params.limit,grpsize1);
Grp2All = nan(54,params.limit,grpsize2);

for iLay = 1:length(layer)
    for iSti = 1:length(stimfrq)
        % through each animal
        for iAn = 1: grpsize1 + grpsize2
            
            if iAn <= grpsize1
                % pull out the data for this group and this animal
                grp1 = table2cell(M1Dat(contains(M1Dat.animal,fullgroup{iAn})...
                    & contains(M1Dat.layer,layer{iLay}) ...
                    & M1Dat.stimulus == stimfrq(iSti),1));
                grp1 = cellfun(@(x) x(:,1:params.limit),grp1,'UniformOutput',false);
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
                
            else
                grp2 = table2cell(M2Dat(contains(M2Dat.animal,fullgroup{iAn})...
                    & contains(M2Dat.layer,layer{iLay}) ...
                    & M2Dat.stimulus == stimfrq(iSti),1));
                grp2 =  cellfun(@(x) x(:,1:params.limit),grp2,'UniformOutput',false);
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
                Grp2All(:,:,iAn-grpsize1) = abs(mean(cat(3,transGrp2{:}),3));
                
            end
        end
        
        % for the group:
        % take the mean across subjects
        grp1mean = nanmean(Grp1All,3);
        grp2mean = nanmean(Grp2All,3);
        
        diffmeans = grp2mean - grp1mean;
        
        %% Mann Whitney U Test (ranksum)
        
        shifted1 = shiftdim(Grp1All,1);
        shifted2 = shiftdim(Grp2All,1);
        obs_Esize = zeros(size(grp1mean));
        obs_Cmass = zeros(size(grp1mean));
        
        % this test works on vectors. so we're doing it pointwise. Awyis.
        for row = 1:size(Grp1All,1)
            for col = 1:size(Grp1All,2)
                
                % take out the point of interest
                point1 = shifted1(col,:,row);
                point1 = point1(~isnan(point1));
                point2 = shifted2(col,:,row);
                point2 = point2(~isnan(point2));
                
                % two-tailed test (stats.p(2)) like with the ttest2 used in mag plots, the
                % samples are large enough to use method 'normal approximation'
                stats = mwwtest(point1,point2);
                if stats.p(2) <= 0.05
                    sigtest = 1;
                else
                    sigtest = 0;
                end
                
                % effect size of mwwtest is r = abs(z/sqrt(n1+n2)) / 0.1 is small, 0.3 is
                % medium, 0.5 is large
                Esize = abs(stats.Z/sqrt(grpsize1+grpsize2));
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
        
        for iOsc = 1:length(osciName)
            obs_layer.(osciName{iOsc}) = obs_Cmass(osciRows{iOsc},:);
            
            % % sum clusters (twice to get final value)
            for i = 1:2
                obs_layer.(osciName{iOsc}) = nansum(obs_layer.(osciName{iOsc}));
            end
        end
        
        
        %% cluster fig
        %% Fig
        cd(homedir); cd figs; mkdir('Spectral_PhCPerm'); cd('Spectral_PhCPerm');
        
        [X,Y]=meshgrid(M1Dat.freq{1}(19:54),params.startTime*1000:(params.limit-201));
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

        
        %% Permutation Step 3 - do the permute
        
        % set up some containers
        mass_clustermass = NaN([1 nperms]);
        disp('Permuting')
        tic
        perm_layer = struct;
        for iOsc = 1:length(osciName)
            perm_layer.(osciName{iOsc}) = NaN([1 nperms]);
        end
        
        % echo loop to find observed values
        for iperm = 1:nperms
            
            % randomize the order of animal selection
            order = randperm(grpsize1+grpsize2);
            FullDat = [M1Dat;M2Dat];
            
            % Get phase coherence for groups
            permbox1 = nan(54,params.limit,grpsize1);
            permbox2 = nan(54,params.limit,grpsize2);
                        
            % through each animal
            for iAn = 1: grpsize1 + grpsize2
                
                if iAn <= grpsize1
                    % pull out the data for the random animal
                    perm1 = table2cell(FullDat(contains(FullDat.animal,fullgroup{order(iAn)})...
                        & contains(FullDat.layer,layer{iLay}) ...
                        & FullDat.stimulus == stimfrq(iSti),1));
                    perm1 = cellfun(@(x) x(:,1:params.limit),perm1,'UniformOutput',false);
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
                    
                else
                    perm2 = table2cell(FullDat(contains(FullDat.animal,fullgroup{order(iAn)})...
                        & contains(FullDat.layer,layer{iLay}) ...
                        & FullDat.stimulus == stimfrq(iSti),1));
                    perm2 = cellfun(@(x) x(:,1:params.limit),perm2,'UniformOutput',false);
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
                    permbox2(:,:,iAn-grpsize1) = curAn;
                    
                end
            end
            
            % Mann Whitney U Test (ranksum) %%
            
            shifted1 = shiftdim(permbox1,1);
            shifted2 = shiftdim(permbox2,1);
            perm_Cmass = zeros(size(grp1mean));
            
            % pointwise
            for row = 1:size(permbox1,1)
                for col = 1:size(permbox1,2)
                    
                    % take out the point of interest
                    point1 = shifted1(col,:,row);
                    point2 = shifted2(col,:,row);
                    point1 = point1(~isnan(point1));
                    point2 = point2(~isnan(point2));
                    
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
            
            for iOsc = 1:length(osciName)
                hold_permlayer = perm_Cmass(osciRows{iOsc},:);
                
                % % sum clusters (twice to get final value)
                hold_permlayer = nansum(nansum(hold_permlayer));
                perm_layer.(osciName{iOsc})(iperm) = hold_permlayer;
            end
            
        end
        toc
        
        cd(homedir); cd DATA; cd Spectral; mkdir('Spectral_PhCPerm'); cd('Spectral_PhCPerm');
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
    end %stim frequency
end %layer