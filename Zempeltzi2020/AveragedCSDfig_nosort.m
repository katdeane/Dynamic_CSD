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
load('Awake10dB_Data.mat')
Awake = Data; clear Data; AWnames = fieldnames(Awake);

Ticks = {'250',  '500', '1000', '2000', '4000', '8000', '16000'};
% Ticknames = {'BF - 3','BF - 2','BF - 1','Best Frequency','BF + 1', 'BF + 2', 'BF + 3'};
AWLFPcont = {}; 
t_thresh = 1.86;

%% line up CSD's with +-3 around the BF for each group

for iA = 1:length(AWnames) %for each animal
        for itick = 1:length(Ticks) %through each freqency around charted BF
            AWLFPcont{iA,itick} = Awake.(AWnames{iA}).CSD{itick}(1:28,1:600);
        end
end


% check containers for size and shape verification (should be number of
% animals by number of defined ticks with each cell having 28x600)


%% Average CSDs
AWmeans = {};

for itick = 1:length(Ticks) %creates 1x7 CSDs with 4 being the BF
     
    awcat = cat(3,AWLFPcont{1,itick},AWLFPcont{2,itick},AWLFPcont{3,itick},AWLFPcont{4,itick}...
        ,AWLFPcont{5,itick},AWLFPcont{6,itick},AWLFPcont{7,itick},AWLFPcont{8,itick},AWLFPcont{9,itick}); %9 animals
    
    AWmeans{itick} = nanmean(awcat,3);
    
end

% check means for size and shape verification (should be 1 by number of
% defined ticks with each cell still having 28x600)


%% Produce CSD figures


figure('Name','Awake Average CSD')
for itick = 1:length(Ticks)
    subplot(2,round(length(Ticks)/2),itick)                  
    imagesc(AWmeans{1,itick}(:,200:500))
    caxis([-0.0005 0.0005])
    colormap('jet')
    title(Ticks{itick})
end

h = gcf;
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(h,'Awake Average CSD','compact')
close (h)



