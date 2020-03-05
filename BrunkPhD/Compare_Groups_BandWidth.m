clear,dbstop if error,warning ('off','all');
Home = pwd; FigFolder = [Home '\Figs_Group_comparison'];
cd ('K:\CSD_dynamic_analysis\DATA')
addpath('K:\CSD_dynamic_analysis\subfunc')

Opto = load('Input_ALL_OPTO_7post_Data.mat');
High = load('Input_HighP_7post_Data.mat');
Control = load('Input_Control_7post_Data.mat');
YFP = load('Input_YFP_7post_Data.mat');

Groups = {'Control' 'YFP' 'Opto' 'High'};

para = {'SinkRMS','tempSinkRMS','SinkPeakAmp','Full_RMS_AVREC','Early_RMS_AVREC',...
    'Late_RMS_AVREC','Full_RMS_RELRES','Early_RMS_RELRES','Late_RMS_RELRES',...
    'Full_RMS_ABSRES','Early_RMS_ABSRES','Late_RMS_ABSRES'};

Cond = [3, 4, 10];
Cond2 = {'Pre','Combi','Post'};


sqrt(mean(Data(1).GOT13.CSD{1}(:,1:200).^2))