%% BF Dist
% This code needs some tuning - however it currently loads in whatever
% scripts are in the group folder to create a gramm plot of the
% distributions of BFs within each group and a gramm plot of all groups overlaid. 

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

%% Run
input = dir('*.mat');
entries = length(input);
cd(home),cd figs;
mkdir('BF Distribution'); cd('BF Distribution');

fullBFlist = []; fullnamelist = []; fullgrouplist = [];
freq = {'125' '250' '500' '1000' '2000' '4000' '8000' '16000'};

for i1 = 1:entries    
    load(input(i1).name); %load Data.mat into workspace
    group = (input(i1).name(1:2));
    
    names = fieldnames(Data); %list of animal names
    NumAnimals = length(names); %groupsize
    namelist = 1:NumAnimals;
    grouplist = cell(1,NumAnimals); %this code can be used to make other gramm plots more automatic
    grouplist(:) = {group};
    
    BFlist = nan(1,NumAnimals);
    for iA = 1:NumAnimals
        BFlist(iA) = Data.(names{iA}).GS_BF;
    end
    
    clear g
    g = gramm('x',log2(BFlist),'y',namelist); %log2 logarithmic doubling (for log scale/octaves)
    g.stat_bin('geom','overlaid_bar','fill','transparent');
    g.set_names('x','Frequency');
    g.set_title([(group) ' BF distribution']);
    g.draw();
    g.export('file_name',[(group) ' BF distribution'], 'file_type','pdf');
    g.export('file_name',[(group) ' BF distribution'], 'file_type','png');
    close all;
    
    fullBFlist = [fullBFlist BFlist];
    fullnamelist = [fullnamelist namelist];
    fullgrouplist = [fullgrouplist grouplist]; 
end

clear g 
g = gramm('x',log2(fullBFlist),'y',fullnamelist,'color',fullgrouplist);
g.stat_bin('geom','overlaid_bar','fill','transparent');
g.set_names('x','Frequency');
g.set_title('BF distribution');
g.axe_property('XTickLabel',freq);
g.set_color_options('map','matlab');
g.draw();
g.export('file_name','BF distribution', 'file_type','pdf');
g.export('file_name','BF distribution', 'file_type','png');
close all;