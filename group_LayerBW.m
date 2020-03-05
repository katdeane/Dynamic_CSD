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
entrynames = {'AnChronic','Ketamine','Awake','Muscimol'};
layers = {'I_IIE','IVE','VbE','VIaE'};

cd([home '\figs\' 'Group Bandwidth']);


for i1 = 1:entries
    load (input(i1).name) %loads data into workspace
    names = fieldnames(Data);
    sinkbandwidth = nan(length(names),length(layers));
    
    for isub = 1:length(names)
        for ilayer = 1:length(layers)
            bandwidth = nansum(~isnan(horzcat(Data.(names{isub}).SinkRMS.(layers{ilayer}))));
            sinkbandwidth(isub, ilayer) = bandwidth;
        end
    end

    h = figure;
    title([entrynames{i1} ' Bandwidth of Sink Response'],'FontSize',15,'FontWeight','bold');
    boxplot(sinkbandwidth)
    ylabel('Bandwidth [Octaves]','FontSize',15,'FontWeight','bold');
    savefig(h,[entrynames{i1} ' Bandwidth of Sink Response']); 
end

close all
