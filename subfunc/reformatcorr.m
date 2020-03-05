%% README

%   This script is designed to format data into .csv table based on
%   specific Sinks and Parameters. Input is from the DATA folder

%   Note: It is not currently needed for it's original specific purchase so
%   feel free to edit it as needed for output. If we need a long term
%   script doing this then we can create a permanent template in the master
%   which can be then altered in branches for specific projects

clear;
warning('OFF');
dbstop if error;

cd('D:\MyCode\Dynamic_CSD_Analysis');
home = pwd;
addpath(genpath(home));

cd DATA;
input = dir('*.mat');
entries = length(input);

%% Run
for i1 = 1:entries
    clear RMS_table tempRMS_table Duration_table PeakAmplitude_table
    disp(['Making tables for: ' (input(i1).name(1:end-9))])
    tic
    load (input(i1).name) %loads data into workspace
    
    cd tables; %opens figure folder to store future images
    mkdir(['Group_' input(i1).name(1:end-9)]); %adds folder to directory
    cd(['Group_' input(i1).name(1:end-9)])
    
    for animal = fieldnames(Data)'; %cycles through all animals in structure 1 by 1
        
        clear L1_2 L4 L5a L5b L6 next
        tones = Data.(animal{1}).Frqz';
        if size(tones,1)==1
            tones = tones';
        end
        current_animal = [animal animal animal animal animal animal animal animal]'; %for table
        current_animal = current_animal(1:length(tones));
        
        %% SinkRMS
        %All EARLY layers
        L1_2 = vertcat(Data.(animal{1}).SinkRMS.I_IIE);
        L4 = vertcat(Data.(animal{1}).SinkRMS.IVE); %seperate layer 4
        L5a = vertcat(Data.(animal{1}).SinkRMS.VaE);
        L5b = vertcat(Data.(animal{1}).SinkRMS.VbE);
        L6 = vertcat(Data.(animal{1}).SinkRMS.VIE);
        
        %make or add to the table
        if ~exist('RMS_table','var')
            RMS_table = table(current_animal,L1_2,L4,L5a,L5b,L6,tones);
        else
            next = table(current_animal,L1_2,L4,L5a,L5b,L6,tones);
            RMS_table = [RMS_table;next];
        end
    end
    writetable(RMS_table,'RMS_table.csv')
    
    for animal = fieldnames(Data)';
        %% tempSinkRMS
        clear L1_2 L4 L5a L5b L6 next 
        tones = Data.(animal{1}).Frqz';
        if size(tones,1)==1
            tones = tones';
        end
        current_animal = [animal animal animal animal animal animal animal animal]'; %for table
        current_animal = current_animal(1:length(tones));

        %All EARLY layers
        L1_2 = vertcat(Data.(animal{1}).tempSinkRMS.I_IIE)';
        L4 = vertcat(Data.(animal{1}).tempSinkRMS.IVE)'; %seperate layer 4
        L5a = vertcat(Data.(animal{1}).tempSinkRMS.VaE)';
        L5b = vertcat(Data.(animal{1}).tempSinkRMS.VbE)';
        L6 = vertcat(Data.(animal{1}).tempSinkRMS.VIE)';
        
        %make or add to the table
        if ~exist('tempRMS_table','var')
            tempRMS_table = table(current_animal,L1_2,L4,L5a,L5b,L6,tones);
        else
            next = table(current_animal,L1_2,L4,L5a,L5b,L6,tones);
            tempRMS_table = [tempRMS_table;next];
        end
    end
    writetable(tempRMS_table,'tempRMS_table.csv')
    
    for animal = fieldnames(Data)';
        %% SinkDur
        clear L1_2 L4 L5a L5b L6 next
        tones = Data.(animal{1}).Frqz';
        if size(tones,1)==1
            tones = tones';
        end
        current_animal = [animal animal animal animal animal animal animal animal]'; %for table
        current_animal = current_animal(1:length(tones));
        
        %All EARLY layers
        L1_2 = vertcat(Data.(animal{1}).SinkDur.I_IIE)';
        L4 = vertcat(Data.(animal{1}).SinkDur.IVE)'; %seperate layer 4
        L5a = vertcat(Data.(animal{1}).SinkDur.VaE)';
        L5b = vertcat(Data.(animal{1}).SinkDur.VbE)';
        L6 = vertcat(Data.(animal{1}).SinkDur.VIE)';
        
        %make or add to the table
        if ~exist('Duration_table','var')
            Duration_table = table(current_animal,L1_2,L4,L5a,L5b,L6,tones);
        else
            next = table(current_animal,L1_2,L4,L5a,L5b,L6,tones);
            Duration_table = [Duration_table;next];
        end
    end
    writetable(Duration_table,'Duration_table.csv')
    
    for animal = fieldnames(Data)';
        %% PeakAmplitude
        clear L1_2 L4 L5a L5b L6 next
        tones = Data.(animal{1}).Frqz';
        if size(tones,1)==1
            tones = tones';
        end
        current_animal = [animal animal animal animal animal animal animal animal]'; %for table
        current_animal = current_animal(1:length(tones));
        
        %All EARLY layers
        L1_2 = vertcat(Data.(animal{1}).PeakAmplitude.I_IIE);
        L4 = vertcat(Data.(animal{1}).PeakAmplitude.IVE); %seperate layer 4
        L5a = vertcat(Data.(animal{1}).PeakAmplitude.VaE);
        L5b = vertcat(Data.(animal{1}).PeakAmplitude.VbE);
        L6 = vertcat(Data.(animal{1}).PeakAmplitude.VIE);
        
        %make or add to the table
        if ~exist('PeakAmplitude_table','var')
            PeakAmplitude_table = table(current_animal,L1_2,L4,L5a,L5b,L6,tones);
        else
            next = table(current_animal,L1_2,L4,L5a,L5b,L6,tones);
            PeakAmplitude_table = [PeakAmplitude_table;next];
        end
    end
    writetable(PeakAmplitude_table,'PeakAmplitude_table.csv')

    cd(home);cd DATA;
end