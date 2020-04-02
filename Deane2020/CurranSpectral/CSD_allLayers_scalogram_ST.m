function CSD_allLayers_scalogram_ST(homedir)
% Input:    Dynamic_CSD\DATA -> *DATA.mat; manually called
% Output:   Runs CWT analysis using the Wavelet toolbox. figures of
%           animal-wise scalograms -> Dynamic_CSD\figs\Spectral_AngPlot
%           and table 'STscalograms_*.mat' with data -> Dynamic_CSD\DATA\Spectral
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
params.sampleRate = 1000; % Hz
params.startTime = -0.2; % seconds
params.timeLimits = [-0.2 0.399]; % seconds
params.frequencyLimits = [5 params.sampleRate/2]; % Hz
params.voicesPerOctave = 8;
params.timeBandWidth = 54;

layers   = {'I_IIE','IVE','VbE','VIaE'};
rel2BF   = [0 -2];
savelist = {'STscalograms_IIBF','STscalograms_IIoffBF','STscalograms_IVBF',...
    'STscalograms_IVoffBF','STscalograms_VbBF','STscalograms_VboffBF',...
    'STscalograms_VIaBF','STscalograms_VIaoffBF'};

count = 1;
for i2BF = 1:length(rel2BF)
    for iLay = 1:length(layers)
        % Change these for specific runs % MAKE THIS A LOOP
        params.layers = {layers{iLay}};
        params.rel2BFlist = rel2BF(i2BF);  % BF = 0; offBF = -2
               
        %% Run CWT
        
        [awake] = runCwtCsd('Awake10dB',params,homedir,'single');
        [anest] = runCwtCsd('AnesthetizedPre',params,homedir,'single');
        [musc]  = runCwtCsd('Muscimol',params,homedir,'single');
        
        %% Reorganize data into Gramm-compatible structure
        cd(homedir); cd DATA; mkdir('Spectral'); cd('Spectral');
        savename = savelist{count};
        wtTable = [struct2table(awake); struct2table(anest); struct2table(musc)];
        save(savename,'wtTable')
        count = count + 1;
    end
end