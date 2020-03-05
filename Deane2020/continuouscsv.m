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
% 

Ticknames = {'BF - 3','BF - 2','BF - 1','Best Frequency','BF + 1', 'BF + 2', 'BF + 3'};
ANLFPcont = {}; AWLFPcont = {}; MLFPcont = {}; AnCLFPcont = {};
%% line up CSD's with +-3 around the BF for each group

for iA = 1:length(ANnames) %for each animal
    BF = find((Anesthetized.(ANnames{iA}).GS_BF) == (Anesthetized.(ANnames{iA}).Frqz')); %get the BF Position
    
    for itick = -3:length(Ticknames)-4 %through each freqency around charted BF
        if BF+itick > 0 && BF+itick <= length(Anesthetized.(ANnames{iA}).singletrialCSD) 
            ANLFPcont{iA,itick+4} = Anesthetized.(ANnames{iA}).singletrialCSD{BF+itick}; %had to cut to smallest form of each matrix (28x600)
        else
            ANLFPcont{iA,itick+4} = NaN(28,600,40); %in case the frequency around the BF doesn't exist
        end
    end
end

for iA = 1:length(AWnames) %for each animal
        BF = find((Awake.(AWnames{iA}).GS_BF) == (Awake.(AWnames{iA}).Frqz'))'; %get the BF Position
        
        for itick = -3:length(Ticknames)-4 %through each freqency around charted BF
            if BF+itick > 0 && BF+itick <= length(Awake.(AWnames{iA}).singletrialCSD)
                AWLFPcont{iA,itick+4} = Awake.(AWnames{iA}).singletrialCSD{BF+itick};
            else
                AWLFPcont{iA,itick+4} = NaN(28,600,50);
            end
        end
end

for iA = 1:length(Mnames) %for each animal
    BF = find((Muscimol.(Mnames{iA}).GS_BF) == (Muscimol.(Mnames{iA}).Frqz')); %get the BF Position
    
    for itick = -3:length(Ticknames)-4 %through each freqency around charted BF
        if BF+itick > 0 && BF+itick <= length(Muscimol.(Mnames{iA}).singletrialCSD)
            MLFPcont{iA,itick+4} = Muscimol.(Mnames{iA}).singletrialCSD{BF+itick};
        else
            MLFPcont{iA,itick+4} = NaN(28,600,50);
        end
    end
end

for iA = 1:length(AnCnames) %for each animal
    BF = find((AnChronic.(AnCnames{iA}).GS_BF) == (AnChronic.(AnCnames{iA}).Frqz')); %get the BF Position
    
    for itick = -3:length(Ticknames)-4 %through each freqency around charted BF
        if BF+itick > 0 && BF+itick <= length(AnChronic.(AnCnames{iA}).singletrialCSD)
            AnCLFPcont{iA,itick+4} = AnChronic.(AnCnames{iA}).singletrialCSD{BF+itick};
        else
            AnCLFPcont{iA,itick+4} = NaN(28,600,50);
        end
    end
end

% check containers for size and shape verification (should be number of
% animals by number of defined ticks with each cell having 28x600)

%% Get the BF out

% for this particular animal for now
GXL03BF = [];
for i = 1:size(ANLFPcont{1,4},3)
    % one single trial csd at a time
    sheet = ANLFPcont{1,4}(:,:,i);
    if i == 1
        GXL03BF = sheet;
        Single = sheet;
    else
        GXL03BF = horzcat(GXL03BF, sheet);   
    end
end

filename1 = 'GXL03_BF.csv';

csvwrite(filename1, GXL03BF)

filename2 = 'SingleCSD.csv';

csvwrite(filename2, Single)









