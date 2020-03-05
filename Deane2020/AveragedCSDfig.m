%% Averaged CSD

% The purpose of this script is to provide an averaged CSD for visual
% representation of each analysis group. 

% ttest2 with unequal variance provides H = 1 in each place the null hypothesis can be
% rejected with confidence of 0.05

%Input:     is DATA; specifically (not automatically) named per Kat's MT groups
%Output:    is in figure folder AvgCSDs; figures only for representation of
%           characteristic profile
%           Welch's t test averaged line graph with significance 

clear
%% standard operations
warning('OFF');
dbstop if error

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
cd (home),cd DATA;

%load in data and call up the list of animal names
load('AnesthetizedPre_Data.mat')
Anesthetized = Data; clear Data; ANnames = fieldnames(Anesthetized);
load('Awake10dB_Data.mat')
Awake = Data; clear Data; AWnames = fieldnames(Awake);
load('Muscimol_Data.mat')
Muscimol = Data; clear Data; Mnames = fieldnames(Muscimol);
load('ANChronic_Data.mat')
AnChronic = Data; clear Data; AnCnames = fieldnames(AnChronic);

% Ticks = {'min_three',  'min_two', 'min_one', 'BF', 'plus_one', 'plus_two', 'plus_three'};
Ticknames = {'BF - 3','BF - 2','BF - 1','Best Frequency','BF + 1', 'BF + 2', 'BF + 3'};
ANLFPcont = {}; AWLFPcont = {}; MLFPcont = {}; AnCLFPcont = {};
t_thresh = 1.86;

%% line up CSD's with +-3 around the BF for each group
for iA = 1:length(ANnames) %for each animal
    BF = find((Anesthetized.(ANnames{iA}).GS_BF) == (Anesthetized.(ANnames{iA}).Frqz')); %get the BF Position
    
    for itick = -3:length(Ticknames)-4 %through each freqency around charted BF
        if BF+itick > 0 && BF+itick <= length(Anesthetized.(ANnames{iA}).CSD) 
            ANLFPcont{iA,itick+4} = Anesthetized.(ANnames{iA}).CSD{BF+itick}(1:28,1:600); %had to cut to smallest form of each matrix (28x600)
        else
            ANLFPcont{iA,itick+4} = NaN(28,600); %in case the frequency around the BF doesn't exist
        end
    end
end

for iA = 1:length(AWnames) %for each animal
        BF = find((Awake.(AWnames{iA}).GS_BF) == (Awake.(AWnames{iA}).Frqz'))'; %get the BF Position
        
        for itick = -3:length(Ticknames)-4 %through each freqency around charted BF
            if BF+itick > 0 && BF+itick <= length(Awake.(AWnames{iA}).CSD)
                AWLFPcont{iA,itick+4} = Awake.(AWnames{iA}).CSD{BF+itick}(1:28,1:600);
            else
                AWLFPcont{iA,itick+4} = NaN(28,600);
            end
        end
end

for iA = 1:length(Mnames) %for each animal
    BF = find((Muscimol.(Mnames{iA}).GS_BF) == (Muscimol.(Mnames{iA}).Frqz')); %get the BF Position
    
    for itick = -3:length(Ticknames)-4 %through each freqency around charted BF
        if BF+itick > 0 && BF+itick <= length(Muscimol.(Mnames{iA}).CSD)
            MLFPcont{iA,itick+4} = Muscimol.(Mnames{iA}).CSD{BF+itick}(1:28,1:600);
        else
            MLFPcont{iA,itick+4} = NaN(28,600);
        end
    end
end

for iA = 1:length(AnCnames) %for each animal
    BF = find((AnChronic.(AnCnames{iA}).GS_BF) == (AnChronic.(AnCnames{iA}).Frqz')); %get the BF Position
    
    for itick = -3:length(Ticknames)-4 %through each freqency around charted BF
        if BF+itick > 0 && BF+itick <= length(AnChronic.(AnCnames{iA}).CSD)
            AnCLFPcont{iA,itick+4} = AnChronic.(AnCnames{iA}).CSD{BF+itick}(1:28,1:600);
        else
            AnCLFPcont{iA,itick+4} = NaN(28,600);
        end
    end
end

% check containers for size and shape verification (should be number of
% animals by number of defined ticks with each cell having 28x600)


%% Average CSDs
ANmeans = {};
AWmeans = {};
Mmeans = {};
AnCmeans = {};
difmeans = {};

for itick = 1:length(Ticknames) %creates 1x7 CSDs with 4 being the BF
    
    ancat = cat(3,ANLFPcont{1,itick},ANLFPcont{2,itick},ANLFPcont{3,itick},ANLFPcont{4,itick}...
        ,ANLFPcont{5,itick},ANLFPcont{6,itick},ANLFPcont{7,itick},ANLFPcont{8,itick},ANLFPcont{9,itick}...
        ,ANLFPcont{10,itick},ANLFPcont{11,itick});
    
    ANmeans{itick} = nanmean(ancat,3);
    
    awcat = cat(3,AWLFPcont{1,itick},AWLFPcont{2,itick},AWLFPcont{3,itick},AWLFPcont{4,itick}...
        ,AWLFPcont{5,itick},AWLFPcont{6,itick},AWLFPcont{7,itick},AWLFPcont{8,itick},AWLFPcont{9,itick}); %9 animals
    
    AWmeans{itick} = nanmean(awcat,3);
    
    mcat = cat(3,MLFPcont{1,itick},MLFPcont{2,itick},MLFPcont{3,itick},MLFPcont{4,itick}...
        ,MLFPcont{5,itick},MLFPcont{6,itick},MLFPcont{7,itick},MLFPcont{8,itick},MLFPcont{9,itick}...
        ,MLFPcont{10,itick},MLFPcont{11,itick}); %11 animals
    
    Mmeans{itick} = nanmean(mcat,3);
    
    anccat = cat(3,AnCLFPcont{1,itick},AnCLFPcont{2,itick},AnCLFPcont{3,itick},AnCLFPcont{4,itick}...
        ,AnCLFPcont{5,itick},AnCLFPcont{6,itick},AnCLFPcont{7,itick},AnCLFPcont{8,itick},AnCLFPcont{9,itick}); %9 animals
    
    AnCmeans{itick} = nanmean(anccat,3);
    
    % means of differences per animal for t test
    difcat = cat(3,difmat{1,itick},difmat{2,itick},difmat{3,itick},difmat{4,itick}...
        ,difmat{5,itick},difmat{6,itick},difmat{7,itick},difmat{8,itick},difmat{9,itick}); %9 animals
    
    difmeans{itick} = nanmean(difcat,3);
    difstd{itick} = nanstd(difcat,0,3); % currently using 0 for W but don't know what it should be KD
    
end

% check means for size and shape verification (should be 1 by number of
% defined ticks with each cell still having 28x600)


%% Produce CSD figures


figure('Name','Anesthetized Average CSD')
for itick = 1:length(Ticknames)
    subplot(2,round(length(Ticknames)/2),itick)                  
    imagesc(ANmeans{1,itick}(:,200:500))
    caxis([-0.0005 0.0005])
    colormap('jet')
    title(Ticknames{itick})
end

h = gcf;
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(h,'Anesthetized Average CSD','compact')
close (h)

figure('Name','Awake Average CSD')
for itick = 1:length(Ticknames)
    subplot(2,round(length(Ticknames)/2),itick)                  
    imagesc(AWmeans{1,itick}(:,200:500))
    caxis([-0.0005 0.0005])
    colormap('jet')
    title(Ticknames{itick})
end

h = gcf;
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(h,'Awake Average CSD','compact')
close (h)

figure('Name','Muscimol Average CSD')
for itick = 1:length(Ticknames)
    subplot(2,round(length(Ticknames)/2),itick)                  
    imagesc(Mmeans{1,itick}(:,200:500))
    caxis([-0.0005 0.0005])
    colormap('jet')
    title(Ticknames{itick})
end

h = gcf;
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(h,'Muscimol Average CSD','compact')
close (h)

figure('Name','Anesthetized Chronic Average CSD')
for itick = 1:length(Ticknames)
    subplot(2,round(length(Ticknames)/2),itick)                  
    imagesc(AnCmeans{1,itick}(:,200:500))
    caxis([-0.0005 0.0005])
    colormap('jet')
    title(Ticknames{itick})
end

h = gcf;
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(h,'Anesthetized Chronic Average CSD','compact')
close (h)


%% Stats of matrices; Best Frequency; ttest2

% Best Frequency full matrix
[H,P,~,STATS] = mod_ttest2(AnCmeans{1,4},AWmeans{1,4},'vartype','unequal'); %help ttest2 to understand function
matHfullBF = H;matPfullBF = P; matSfullBF = STATS.tstat;
% getting the average lines from the mats
zeroline = ones(1,301); zeroline(zeroline == 1) = 0;
fullmeanBF = mean(matSfullBF(:,200:500)); fullPBF = mean(matPfullBF(:,200:500));
LIVmeanBF = mean(matSfullBF(8:11,200:500)); IVPBF = mean(matPfullBF(8:11,200:500));
LVbmeanBF = mean(matSfullBF(19:22,200:500)); VPBF = mean(matPfullBF(19:22,200:500));
LVIameanBF = mean(matSfullBF(23:26,200:500)); VIPBF = mean(matPfullBF(23:26,200:500));
extrafull = fullmeanBF; extrafull(fullPBF > 0.05) = NaN;
extraIV = LIVmeanBF; extraIV(IVPBF > 0.05) = NaN;
extraV = LVbmeanBF; extraV(VPBF > 0.05) = NaN;
extraVI = LVIameanBF; extraVI(VIPBF > 0.05) = NaN;

% Best Frequency full column stats
[H,P,~,STATS] = ttest2(AnCmeans{1,4},AWmeans{1,4},'vartype','unequal');
colHfullBF = H;colPfullBF = P; colSfullBF = STATS.tstat;
% Best Frequency full row stats and pulling the layer rows out
[H,P,~,STATS] = ttest2(AnCmeans{1,4},AWmeans{1,4},'vartype','unequal','dim',2);
rowHfullBF = H; rowPfullBF = P; rowSfullBF = STATS.tstat;
ivrows = rowHfullBF(8:11)'; 
vbrows = rowHfullBF(19:22)';
viarows = rowHfullBF(23:26)';

% Best Frequency L IV column stats
[H,P,~,STATS] = ttest2(AnCmeans{1,4}(8:11,1:600),AWmeans{1,4}(8:11,1:600),'vartype','unequal');
colHIVBF = H;colPIVBF = P; colSIVBF = STATS.tstat;
% Best Frequency L Vb column stats
[H,P,~,STATS] = ttest2(AnCmeans{1,4}(19:22,1:600),AWmeans{1,4}(19:22,1:600),'vartype','unequal');
colHVBF = H;colPVBF = P; colSVBF = STATS.tstat;
% Best Frequency L VIa column stats
[H,P,~,STATS] = ttest2(AnCmeans{1,4}(23:26,1:600),AWmeans{1,4}(23:26,1:600),'vartype','unequal');
colHVIBF = H;colPVIBF = P; colSVIBF = STATS.tstat;


%plotting it all (Best Frequency edition)
figure; 
subplot(1,3,1) %t-value matrix
imagesc(matSfullBF(:,200:500));
colormap(jet);
caxis([-10 10]); %not all numbers in range but majority are
colorbar('SouthOutside')
title('t mat (BF)');

subplot(1,3,2) %p-value matrix
imagesc(matPfullBF(:,200:500));
caxis([0 .1]); %makes p=0.05 as if 0 on colorscale (green)
colorbar('SouthOutside')
title('p mat (BF)');

subplot(1,3,3) %average lines directly from the t-matrix
plot(fullmeanBF); hold on; plot(extrafull,'LineWidth',2);
plot(LIVmeanBF); plot(extraIV,'LineWidth',2);
plot(LVbmeanBF); plot(extraV,'LineWidth',2);
plot(LVIameanBF); plot(extraVI,'LineWidth',2);
plot(zeroline);  
ylim([-10,15])
title('layer-wise relative to csd')
legend('full','full p<0.05','LIV','IV p<0.05','LVb','V p<0.05','LVIa','VI p<0.05')

h = gcf;
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(h,'Stats Layout for BF','compact')
close (h)

figure; %average lines from seperate analyses
plot(colSfullBF(:,200:500)); hold on; plot(colSIVBF(:,200:500)); plot(colSVBF(:,200:500)); plot(colSVIBF(:,200:500));
plot(zeroline);  
title('layer-wise seperately analyzed BF')
colHfullBF(colHfullBF == 0) = NaN; colHfullBF(colHfullBF == 1) = 10; 
colHIVBF(colHIVBF == 0) = NaN; colHIVBF(colHIVBF == 1) = 11;
colHVBF(colHVBF == 0) = NaN; colHVBF(colHVBF == 1) = 12;
colHVIBF(colHVIBF == 0) = NaN; colHVIBF(colHVIBF == 1) = 13;
plot(1:301, colHfullBF(:,200:500), '*'); %Full H = 1
plot(1:301, colHIVBF(:,200:500), '*'); %Full H = 1
plot(1:301, colHVBF(:,200:500), '*'); %Full H = 1
plot(1:301, colHVIBF(:,200:500), '*'); %Full H = 1
legend('full','LIV','LVb','LVIa','0','full p<0.05','IV p<0.05','V p<0.05','VI p<0.05')

h = gcf;
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(h,'layer specific analysis BF','compact')
close (h)


%% Stats of matrices; Off Best Frequency; ttest2

% Best Frequency -2 full matrix
[H,P,~,STATS] = mod_ttest2(AnCmeans{1,2},AWmeans{1,2},'vartype','unequal'); %help ttest2 to understand function
matHfulloff = H;matPfulloff = P; matSfulloff = STATS.tstat;
% getting the average lines from the mats
zeroline = ones(1,301); zeroline(zeroline == 1) = 0;
fullmeanoff = mean(matSfulloff(:,200:500)); fullPoff = mean(matPfulloff(:,200:500));
LIVmeanoff = mean(matSfulloff(8:11,200:500)); IVPoff = mean(matPfulloff(8:11,200:500));
LVbmeanoff = mean(matSfulloff(19:22,200:500)); VPoff = mean(matPfulloff(19:22,200:500));
LVIameanoff = mean(matSfulloff(23:26,200:500)); VIPoff = mean(matPfulloff(23:26,200:500));
extraofffull = fullmeanoff; extraofffull(fullPoff > 0.05) = NaN;
extraoffIV = LIVmeanoff; extraoffIV(IVPoff > 0.05) = NaN;
extraoffV = LVbmeanoff; extraoffV(VPoff > 0.05) = NaN;
extraoffVI = LVIameanoff; extraoffVI(VIPoff > 0.05) = NaN;

% Best Frequency -2 full column stats
[H,P,~,STATS] = ttest2(AnCmeans{1,2},AWmeans{1,2},'vartype','unequal');
colHfulloff = H;colPfulloff = P; colSfulloff = STATS.tstat;
% Best Frequency -2 full row stats and pulling the layer rows out
[H,P,~,STATS] = ttest2(AnCmeans{1,2},AWmeans{1,2},'vartype','unequal','dim',2);
rowHfulloff = H; rowPfulloff = P; rowSfulloff = STATS.tstat;
ivrowsoff = rowHfulloff(8:11)'; 
vbrowsoff = rowHfulloff(19:22)';
viarowsoff = rowHfulloff(23:26)';

% Best Frequency -2 L IV column stats
[H,P,~,STATS] = ttest2(AnCmeans{1,2}(8:11,1:600),AWmeans{1,2}(8:11,1:600),'vartype','unequal');
colHIVoff = H;colPIVoff = P; colSIVoff = STATS.tstat;
% Best Frequency -2 L Vb column stats
[H,P,~,STATS] = ttest2(AnCmeans{1,2}(19:22,1:600),AWmeans{1,2}(19:22,1:600),'vartype','unequal');
colHVoff = H;colPVoff = P; colSVoff = STATS.tstat;
% Best Frequency -2 L VIa column stats
[H,P,~,STATS] = ttest2(AnCmeans{1,2}(23:26,1:600),AWmeans{1,2}(23:26,1:600),'vartype','unequal');
colHVIoff = H;colPVIoff = P; colSVIoff = STATS.tstat;


%plotting it all (Best Frequency -2 edition)
figure; 
subplot(1,3,1) %t-value matrix
imagesc(matSfulloff(:,200:500));
colormap(jet);
caxis([-10 10]); %not all numbers in range but majority are
colorbar('SouthOutside')
title('t mat (off BF)');

subplot(1,3,2) %p-value matrix
imagesc(matPfulloff(:,200:500));
caxis([0 .1]); %makes p=0.05 as if 0 on colorscale (green)
colorbar('SouthOutside')
title('p mat (off BF)');

subplot(1,3,3) %average lines directly from the t-matrix
plot(fullmeanoff); hold on; plot(extraofffull,'LineWidth',2);
plot(LIVmeanoff); plot(extraoffIV,'LineWidth',2);
plot(LVbmeanoff); plot(extraoffV,'LineWidth',2);
plot(LVIameanoff); plot(extraoffVI,'LineWidth',2);
plot(zeroline);  
ylim([-10,15])
title('layer-wise relative to csd')
legend('full','full p<0.05','LIV','IV p<0.05','LVb','V p<0.05','LVIa','VI p<0.05')

h = gcf;
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(h,'Stats Layout for off BF','compact')
close (h)

figure; %average lines from seperate analyses
plot(colSfulloff(:,200:500)); hold on; plot(colSIVoff(:,200:500)); plot(colSVoff(:,200:500)); plot(colSVIoff(:,200:500));
plot(zeroline);  
title('layer-wise seperately analyzed')
colHfulloff(colHfulloff == 0) = NaN; colHfulloff(colHfulloff == 1) = 10; 
colHIVoff(colHIVoff == 0) = NaN; colHIVoff(colHIVoff == 1) = 11;
colHVoff(colHVoff == 0) = NaN; colHVoff(colHVoff == 1) = 12;
colHVIoff(colHVIoff == 0) = NaN; colHVIoff(colHVIoff == 1) = 13;
plot(1:301, colHfulloff(:,200:500), '*'); 
plot(1:301, colHIVoff(:,200:500), '*'); 
plot(1:301, colHVoff(:,200:500), '*'); 
plot(1:301, colHVIoff(:,200:500), '*'); 
legend('full','LIV','LVb','LVIa','0','full p<0.05','IV p<0.05','V p<0.05','VI p<0.05')

h = gcf;
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(h,'layer specific analysis off BF','compact')
close (h)



filename = 'stats_output_ALL';
save(filename,...
    'matHfullBF','matPfullBF','matSfullBF',...
    'matHfulloff','matPfulloff','matSfulloff',...
    'colHfullBF','colPfullBF','colSfullBF',...
    'colHfulloff','colPfulloff','colSfulloff',...
    'rowHfullBF','rowPfullBF','rowSfullBF',...
    'rowHfulloff','rowPfulloff','rowSfulloff',...
    'colHIVBF','colPIVBF','colSIVBF',...
    'colHIVoff','colPIVoff','colSIVoff',...
    'colHVBF','colPVBF','colSVBF',...
    'colHVoff','colPVoff','colSVoff',...
    'colHVIBF','colPVIBF','colSVIBF',...
    'colHVIoff','colPVIoff','colSVIoff'...
    )
