%% README
%This script takes the Data.mat from the DATA folder and pulls out the bandwidths
%(based on fwhm) and the tuning curve (based on 1 std) from the groups in
%that folder and outputs one bar graph for each.

%currently works best with single condition groups

cd('D:\MyCode\Dynamic_CSD_Analysis');
warning('OFF');
dbstop if error

home = pwd;
addpath(genpath(home));

cd DATA;
input = dir('*.mat');
entries = length(input);

%% Run
bandmean = [];
bandsem = [];
tuningmean = [];
tuningsem = [];
ticks = {};

for i1 = 1:entries
    load (input(i1).name) %loads data into workspace
    ticks = [ticks input(i1).name(1:end-9)];
    names = fieldnames(Data); %list of animal names
    NumAnimals = length(names);
    
    tuningwidth = [];
    bandwidth = [];
    
    for iB = 1:NumAnimals
        bw = Data.(names{iB}).Bandwidth;
        bandwidth = [bandwidth bw];
        tw = Data.(names{iB}).Tuningwidth;
        tuningwidth = [tuningwidth tw];
    end
    
    bandmean = [bandmean nanmean(bandwidth)];
    bandsem = [bandsem nanstd(bandwidth,1)/sqrt(sum(~isnan(bandmean)))];
    
    tuningmean = [tuningmean nanmean(tuningwidth)];
    tuningsem = [tuningsem nanstd(tuningwidth,1)/sqrt(sum(~isnan(tuningmean)))];
    
end

cd(home); cd figs; %opens figure folder to store future images
mkdir('Group Bandwidth'); %adds folder to directory

h=figure;

title('Mean Bandwidth with SEM','FontSize',15,'FontWeight','bold');
bar(bandmean); hold on;
errorbar(1:entries,bandmean, bandsem, '.');
set(gca,'XTick',1:entries); set(gca,'XTickLabel',ticks,'FontSize',8);
set(gca,'XTickLabelRotation',45);
ylabel('Bandwidth','FontSize',15,'FontWeight','bold');

cd([home '\figs\' 'Group Bandwidth']);
savefig(h,[ticks{1} '_' ticks{2} ' Mean Bandwidth with SEM']); close all

h=figure;

title('Mean Tuningwidth with SEM','FontSize',15,'FontWeight','bold');
bar(tuningmean); hold on;
errorbar(1:entries,tuningmean, tuningsem, '.');
set(gca,'XTick',1:entries); set(gca,'XTickLabel',ticks,'FontSize',8);
set(gca,'XTickLabelRotation',45);
ylabel('Tuningwidth','FontSize',15,'FontWeight','bold');

savefig(h,[ticks{1} '_' ticks{2} ' Mean Tuningwidth with SEM']); close all
