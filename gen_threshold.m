clear all
%% README 

% This script generates an average threshold over multiple files specific
% in the Condition. i.e. if 3 Pre's are present and the Condition checks
% these, the baseline threshold taken into the sink_dura code will be based
% on the averaged baseline data of all 3 pre measurements. 

% Note: this take the layers specified so there may be mismatch if you
% change the layers and don't re-generate thresholds. 

% Note: If this is not run, the sink_dura code will generate baseline
% threshold on each CSD it analyzes individually - this simply gives it a
% more stable baseline because of the averaging but it is not necessary

%% standard operations
warning('OFF');
dbstop if error %stops execution of program if and allows user to examine the workspace
cd('D:\MyCode\Dynamic_CSD_Analysis');
home = pwd; %pwd: print working directory, gives full current directory title of "home"
addpath(genpath(home)); %adds all subdirectories of 'home' to path for access

%% Start
cd (home),cd groups;

input = dir('*.m'); %dir *.m lists all program files in the current directory (stores 5 fields of data but struc array contains #x1)
entries = length(input); %number of input entries

%how many groups?
for i1 = 1:entries
    run(input(i1).name);
    
    Indexer = struct; %creating empty structure to fill
    Condition = {'Pre'}; %this is the condition within which we want to determine the average threshold level
    
    %how many animals?
    for i2 = 1:length(animals)
        cd (home); cd raw;    %we are switching folders to raw to load the .mat's
        
        clear name; clear Baseline;
        std_BL = [];
        mean_BL = [];
        Baseline = struct;
        name = animals{i2}; %assigns the i2th entry in the list in the current group to variable 'name'
        measurementlist = Cond.condition{i2};
        disp(['Analyzing animal: ' name])
        
        %how many post measurements in that animal?
        for i3 = 1:length(measurementlist)
            
            %for each measurement of each animal in each group...
            measurement = measurementlist{i3};
            
            load ([name '_' measurement]);clear avgFP;   %loads info from mat file (Data, DS, Header, P, Spik_list, SWEEP)
            
            BL = Header.t_pre*P.Fs_AD(1); %BL-baseline %t_pre is the time before the tone %Fs_AD - Sampling frequency of channels (they are all the same so we use first value)
            
            cd (home),cd subfunc;
            
            [~, AvgCSD, ~,~,~,~,~,~] =...
                SingleTrialCSD_reduced(SWEEP, str2num(channels{i2}),BL);
            
            %to create a full matrix by horizontally stacking the BL's of each AvgCSD
            global_BL = cell2mat(cellfun(@(x) x(:,1:BL),AvgCSD, 'UniformOutput',0));
            STD = nanstd(global_BL(:)); %the standard deviation of the baseline
            AVG = nanmean(global_BL(:)); %the mean of the baseline
            
            std_BL = [std_BL STD];
            mean_BL = [mean_BL AVG];
            
        end
        threshstd = mean(std_BL);
        threshmean = mean(mean_BL);
        
        Baseline.threshmean = threshmean;
        Baseline.threshstd = threshstd;
        
        cd (home); cd groups;
        save([name '_Baseline'],'Baseline');
    end
end
