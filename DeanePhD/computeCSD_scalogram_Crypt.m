function computeCSD_scalogram_Crypt(homedir)

% Input:    Dynamic_CSD\DATA -> *DATA.mat; manually called through function
%           runCwtCsd_mice.m and Dynamic_CSD\groups -> *.m (* = group name)
% Output:   Runs average trial CWT analysis using the Wavelet toolbox. 
%           tables 'WT_AM.mat' 'WT_CL.mat' and 'WT_Tuning.mat' with all 
%           data -> Dynamic_CSD\DATA\Spectral
%           Runs single trial wavelet analysis and saves out tables per
%           group and measurement into -> Dynamic_CSD\DATA\Spectral

%% standard operations
warning('OFF');
dbstop if error

% Change directory to your working folder
if ~exist('homedir','var')
    if exist('D:\MyCode\Dynamic_CSD','dir') == 7
        cd('D:\MyCode\Dynamic_CSD');
    elseif exist('C:\Users\kedea\Documents\Work Stuff\Dynamic_CSD','dir') == 7
        cd('C:\Users\kedea\Documents\Work Stuff\Dynamic_CSD')
    elseif exist('D:\Dynamic_CSD','dir') == 7
        cd('D:\Dynamic_CSD')
    end
    
    homedir = pwd;
    addpath(genpath(homedir));
end

cd(homedir)
%% INIT PARAMETERS

params.sampleRate = 1000; % Hz
params.startTime = -0.2; % seconds
params.timeLimits = [-0.2 0.399]; % seconds
params.frequencyLimits = [5 params.sampleRate/2]; % Hz
params.voicesPerOctave = 8;
params.timeBandWidth = 54;
params.layers = {'I_II','IV','V','VI'};
params.rel2BFlist = [0 -2];

cd(homedir); cd DATA; mkdir('Spectral'); cd('Spectral');
%% CWT analysis Average

% disp('Running average WTs for KIC group')
% tic
% [KIC_Tune,KIC_CL,KIC_AM,KIC_SP] = runCwtCsd_CryptAV('KIC',params,homedir);
% toc
% disp('... for KIT group')
% tic
% [KIT_Tune,KIT_CL,KIT_AM,KIT_SP] = runCwtCsd_CryptAV('KIT',params,homedir);
% toc
% disp('... for KIV group')
% tic
% [KIV_Tune,KIV_CL,KIV_AM,KIV_SP] = runCwtCsd_CryptAV('KIV',params,homedir);
% toc
% % Reorganize data into Gramm-compatible structure
% 
% organize 3 tables all with each group
% WT_Tuning = [struct2table(KIC_Tune); struct2table(KIT_Tune); struct2table(KIV_Tune)];
% WT_CL = [struct2table(KIC_CL); struct2table(KIT_CL); struct2table(KIV_CL)];
% WT_AM = [struct2table(KIC_AM); struct2table(KIT_AM); struct2table(KIV_AM)];
% WT_SP = [struct2table(KIC_SP); struct2table(KIT_SP); struct2table(KIV_SP)];
% 
% save('WT_Tuning.mat','WT_Tuning')
% save('WT_CL.mat','WT_CL')
% save('WT_AM.mat','WT_AM')
% save('WT_SP.mat','WT_SP')
% clear WT_Tuning WT_CL WT_AM WT_SP KIC_Tune KIC_CL KIC_AM KIC_SP KIT_Tune ...
%     KIT_CL KIT_AM KIT_SP KIV_Tune KIV_CL KIV_AM KIV_SP

%% CWT analysis Single Trial

% these need to be sorted by specific measurement and group because they are
% otherwise to large. Permutation comparisons will always have to call in
% the two appropiate measurements/groups for comparison.

group = {'KIC','KIT','KIV'};
measurements = {'preCL_1','CL_1','CL_2','CL_3','CL_4','preAM_1','AM_1',...
    'AM_2','AM_3','AM_4','spPre1_1','spPost1_1','spPre2_1','spPost2_1','spEnd_1'}; 
% 'Pre_1','Pre_2','Pre_3','Pre_4'
% 'preAMtono_1','preAMtono_2','preAMtono_3','preAMtono_4'
% 'preCLtono_1','preCLtono_2','preCLtono_3','preCLtono_4'
% 'CLtono_1','AMtono_1'

for imeas = 1:length(measurements)
    for iGro = 1:length(group)
        
        disp(['Running single trial WTs of ' measurements{imeas} ...
            ' for ' group{iGro} ' group'])
        
        tic
        [datOut] = runCwtCsd_CryptST(group{iGro},params,homedir,...
            measurements{imeas});        
        
        WT_st = struct2table(datOut);
        save([group{iGro} '_' measurements{imeas} '.mat'],'WT_st')
        clear WT_st datOut
        toc
    end
end

