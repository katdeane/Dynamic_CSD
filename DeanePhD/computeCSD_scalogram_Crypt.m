function computeCSD_scalogram_Crypt(homedir)

% Input:    Dynamic_CSD\DATA -> *DATA.mat; manually called through function
% Output:   Runs CWT analysis using the Wavelet toolbox. tables 'WT_AM.mat' 
%           'WT_CL.mat' and 'WT_Tuning.mat' with all data -> Dynamic_CSD\DATA\Spectral

%% standard operations
warning('OFF');
dbstop if error

% Change directory to your working folder
if ~exist('homedir','var')
    if exist('D:\MyCode\Dynamic_CSD','dir') == 7
        cd('D:\MyCode\Dynamic_CSD');
    elseif exist('C:\Users\kedea\Documents\Work Stuff\Dynamic_CSD','dir') == 7
        cd('C:\Users\kedea\Documents\Work Stuff\Dynamic_CSD')
    end
    
    homedir = pwd;
    addpath(genpath(homedir));
end

cd(homedir)
%% INIT PARAMETERS

%close all
params.sampleRate = 1000; % Hz
params.startTime = -0.2; % seconds
params.timeLimits = [-0.2 0.399]; % seconds
params.frequencyLimits = [5 params.sampleRate/2]; % Hz
params.voicesPerOctave = 8;
params.timeBandWidth = 54;
params.layers = {'I_II','IV','V','VI'};
params.rel2BFlist = [0 -2];

%% CWT analysis
tic
[KIC_Tune,KIC_CL,KIC_AM,KIC_SP] = runCwtCsd_Crypt('KIC',params,homedir);
toc
tic
[KIT_Tune,KIT_CL,KIT_AM,KIT_SP] = runCwtCsd_Crypt('KIT',params,homedir);
toc
tic
[KIV_Tune,KIV_CL,KIV_AM,KIV_SP] = runCwtCsd_Crypt('KIV',params,homedir);
toc
%% Reorganize data into Gramm-compatible structure

cd(homedir); cd DATA; mkdir('Spectral'); cd('Spectral');

% organize 3 tables all with each group
WT_Tuning = [struct2table(KIC_Tune); struct2table(KIT_Tune); struct2table(KIV_Tune)];
WT_CL = [struct2table(KIC_CL); struct2table(KIT_CL); struct2table(KIV_CL)];
WT_AM = [struct2table(KIC_AM); struct2table(KIT_AM); struct2table(KIV_AM)];
WT_SP = [struct2table(KIC_SP); struct2table(KIT_SP); struct2table(KIV_SP)];

save('WT_Tuning.mat','WT_Tuning')
save('WT_CL.mat','WT_CL')
save('WT_AM.mat','WT_AM')
save('WT_SP.mat','WT_SP')
