%% PIPELINE DEANE 2020 \m/(>.<)\m/

% By: Katrina Deane; katrina.deane@lin-magdeburg.de

% This code is intended only for use with the data specifically for
% Deane et al. 2020 (Anesthetized, Awake10dB, and Muscimol).
% Any other data run through here will require manual edits and seperate scripts.

% Please ensure that any external files are not in main folders groups or
% DATA. Files generated from the output folder are called manually; all
% other input is dynamic and will attempt to run whatever is inside the
% main folder.

% If everything is properly called and output, this entire pipeline
% generates all statistics and plots used for Deane et al 2020. Please cite
% us if you use these scripts :)

%% Start
clear; clc;

% add your directory to the list OR comment this step and simply start in
% the correct directory 
if exist('D:\MyCode\Dynamic_CSD','dir') == 7
    cd('D:\MyCode\Dynamic_CSD');
elseif exist('C:\Users\kedea\Documents\Work Stuff\Dynamic_CSD','dir') == 7
    cd('C:\Users\kedea\Documents\Work Stuff\Dynamic_CSD')
end

homedir = pwd;
addpath(genpath(homedir));

% uncomment on first run: 
% mkdir('groups') % -> for input
% mkdir('raw') % -> for input
% mkdir('figs') % -> for output
% mkdir('DATA') % -> for output

%% Group creation

%Input:     sink_dura_single.m and several other edited functions
%Output:    Figures of all single animals in "Single..." folder 
%           DATA.mat files in DATA folder
Dynamic_CSD_Inf_single(homedir)

%% Group Sorting

%Input:     D:\MyCode\Dynamic_CSD_Analysis\DATA -> *DATA.mat; (bin, zscore,
%           mirror)
%Output:    Figures of groups in "Group..." folder 
%           .mat files in DATA/Output folder
%           AVG data
GroupAnalysis_fnc_singlecondition(1,0,0,homedir); 

%% Permutation analysis of CSD profiles

cd(homedir);
diary permutation_readout 

%Input:     Dynamic_CSD\DATA\; specifically (not automatically) called
%Output:    is placed in figs - PermuteCSDs; observed difference matrix,
%           observed t matrix, observed cluster matrix in .fig / the
%           observed vs permustation curves across frequency bins in .pdf
%           and .png (grammplots), ALSO the stats output calculated in .mat
PermutationTestCSD(homedir)

diary off
%% Tuning curves

cd(homedir);
diary sinktuning_readout_omitnan

% Input:    Dynamic_CSD\DATA\Output -> *.mat 
% IMPORTANT:Command window needs to be copy/pasted to carry info from
%           post-hoc tests. This is done with the diary function here
% Output:   *.mat files in DATA\Stats\
Stats_Loop_5Freq(homedir) %St with +-2 & cohen's d

diary off

% Input:     Dynamic_CSD\DATA\Output -> *.mat 
% Output:    gramm plots for all tuning curves in "Group grammplots ST"
grammplots_ST(homedir)

% USE (i.e) grammcombo_4sinks.m to visualize multiple tuning curves. Go
% there directly to control which sinks and in which paramter you would
% like to output
open('grammcombo_4sinks.m')

%% AvRec

%Input:     D:\MyCode\Dynamic_CSD_Analysis\DATA -> *DATA.mat; 
%Output:    a carry over .mat and pictures in D:\MyCode\Dynamic_CSD_Analysis\DATA\avrec_compare
snglAVREC_BF(homedir)

cd(homedir);
diary avrec_readout_newcutoff 

%STATS:
% Input:    Dynamic_CSD\DATA -> *DATA.mat 
% IMPORTANT:Command window needs to be copy/pasted to carry info from
%           post-hoc tests. This is done with the diary function here
% Output:   Single frequency bin and tuning curve statistics (p and cohen's
%           d values) into \DATA\Stats*
avrec_relres_stats(homedir)

diary off

% folder:   Dynamic_CSD\DATA\avrec_compare
% Input:    stats from single trial avrec and relres 
% Output:   BrownScyth statistical comparison of variances; box plots and
%           stats saved as pictures
BrownScyth(homedir)

%% Spectral Analysis

% Input:    Dynamic_CSD\DATA -> *DATA.mat; manually called
% Output:   Runs CWT analysis using the Wavelet toolbox. figures of
%           animal-wise scalograms -> Dynamic_CSD\figs\Spectral_MagPlot
%           and table 'scalograms.mat' with all data -> Dynamic_CSD\DATA\Spectral
CSD_allLayers_scalogram(homedir)

% Input:    Dynamic_CSD\DATA -> *DATA.mat; manually called
% Output:   Runs CWT analysis using the Wavelet toolbox. figures of
%           animal-wise scalograms -> Dynamic_CSD\figs\Spectral_AngPlot
%           and table 'STscalograms_*.mat' with data -> Dynamic_CSD\DATA\Spectral
CSD_allLayers_scalogram_ST(homedir)

 %% Run Permutation Magnitude
rel2BFin = [0 -2];
layerin = {'I_IIE','IVE','VbE','VIaE';};
for rel2BF = 1:length(rel2BFin)
    for layer = 1:length(layerin)

    % Input:    Layer to analyze, relative to BF
    %           Dynamic_CSD\DATA\scalograms.mat
    % Output:   Figures for means and observed difference of awake/ketamine
    %           comparison; figures for observed t values, clusters, ttest line
    %           output; boxplot and significance of permutation test ->
    %           Dynamic_CSD\figs\Spectral_MagPerm and
    %           Dynamic_CSD\DATA\Spectral\MagPerm
    PermutationTestScalogram(layerin{layer},rel2BFin(rel2BF),homedir)
        

    % Input:    Layer to analyze, relative to BF
    %           Dynamic_CSD\DATA\scalograms.mat
    % Output:   Figures for cohen's D of observed difference of awake/ketamine
    %           comparison;  -> Dynamic_CSD\figs\Spectral_MagPerm
    cohensDScalogram(layerin{layer},rel2BFin(rel2BF),homedir)

    end
end

%% Run Permutation Phase Coherence
% Input:    matrix to load, layer and rel2BF as a name
%           Dynamic_CSD\DATA\STscalograms_*.mat
% Output:   Figures for means and observed difference of awake/ketamine
%           comparison; figures for observed t values, clusters, ttest line
%           output; boxplot and significance of permutation test -> Pictures 
%           folder -> phase permutations

% please note that this code takes an enourmously long time to run
PermutationTest_PhaseCoh('STscalograms_IIBF.mat',' I_II BF',homedir)
PermutationTest_PhaseCoh('STscalograms_IVBF.mat',' IV BF',homedir)
PermutationTest_PhaseCoh('STscalograms_VbBF.mat',' Vb BF',homedir)
PermutationTest_PhaseCoh('STscalograms_VIaBF.mat',' VIa BF',homedir)

PermutationTest_PhaseCoh('STscalograms_IIoffBF.mat',' I_II off BF',homedir)
PermutationTest_PhaseCoh('STscalograms_IVoffBF.mat',' IV off BF',homedir)
PermutationTest_PhaseCoh('STscalograms_VboffBF.mat',' Vb off BF',homedir)
PermutationTest_PhaseCoh('STscalograms_VIaoffBF.mat',' VIa off BF',homedir)

%% End
