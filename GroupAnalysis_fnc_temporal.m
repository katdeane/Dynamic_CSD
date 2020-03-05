function X = GroupAnalysis_fnc_temporal(bin, zscore, mirror)
% This function asks if it should bin (1 or 0), take zscore (1 or 0), or
% mirror (1 or 0)

% This group analysis generates temporal changes, shifts away from the BF,
% and %MB what else is this producing?

warning('OFF');
dbstop if error
nThresh = 0.25; % Threshold for animals

% Change directory to your working folder
try cd('D:\MyCode\Dynamic_CSD_Analysis');
catch
    cd('D:\Dynamic_CSD_Analysis');
end
home = pwd; %pwd: print working directory, gives full current directory title of "home"
addpath(genpath(home)); %adds all subdirectories of 'home' to path for access

cd DATA

input = dir('*.mat');
entries = length(input);

para = { 'SinkRMS', 'tempSinkRMS', 'init_Peak_BF_tune','SinkPeakAmp',...
    'SinkPeakLate','IntCurFlow', 'SinkDur',...
    'Sinkonset', 'Sinkoffset','Full_RMS_AVREC', 'Early_RMS_AVREC',...
    'Late_RMS_AVREC', 'Full_RMS_RELRES', 'Early_RMS_RELRES',...
    'Late_RMS_RELRES','Full_RMS_ABSRES','Early_RMS_ABSRES','Late_RMS_ABSRES'};%,'BW_full','BW_early','BW_late'

pass_on ={'Full_Single_RMS_AVREC','Early_Single_RMS_AVREC',...
    'Late_Single_RMS_AVREC','Full_Single_RMS_RELRES',...
    'Early_Single_RMS_RELRES','Late_Single_RMS_RELRES',...
    'Full_Single_RMS_ABSRES','Early_Single_RMS_ABSRES',...
    'Late_Single_RMS_ABSRES',};% all the single entries

% SINKs = {'VbE', 'IVE', 'VIE', 'VaE', 'I_IIE','VbL', 'IVL', 'VIL', 'VaL', 'I_IIL'};
SINKs = {'VbE','IVE','VIaE','VIbE','VaE','I_IIE','VbL','IVL','VIaL','VIbL','VaL','I_IIL'};
% SINKs = {'I_IIL','IVE','VaE','VbE', 'VIE','InfE', 'VIL'}; 
for i1 = 1:entries

    load (input(i1).name)
    CurAn = (input(i1).name);
    
    names = fieldnames(Data);
    NumAnimals = length(names);
    
    BF_Pos = [];
    DimDat = length(Data);
    
    %to determine where the BF should be placed (i.e. the length of the longest frequency list)
    for iNum = 1:NumAnimals        
        for idim = 1:DimDat
            if ~isempty(Data(idim).(names{iNum}))
                if isempty(BF_Pos)
                    BF_Pos = length(Data(idim).(names{iNum}).Frqz);
                elseif length(Data(idim).(names{iNum}).Frqz) > BF_Pos
                    BF_Pos = length(Data(idim).(names{iNum}).Frqz);
                end
            end
        end        
    end
    
    NumFreq = (BF_Pos*2)-1; %so that BF_Pos is center

    % Granular BF based
    [Data_GSBased,~] = Groupsorting(Data,names,para,DimDat,SINKs,NumFreq,'GS',mirror,nThresh,bin);

    % Selftuning based
    [Data_STBased,ST_Shift] = Groupsorting(Data,names,para,DimDat,SINKs,NumFreq,'ST',mirror,nThresh,bin);


     if strcmp(Data(3).(names{1}).Condition,'Pre_3')
                PRE = 3;               
     else
                PRE = 1;   
     end
     
    if zscore == 1 
    [Zscrd_DATA_GS,~] = ParaZScore(Data_GSBased, para, SINKs,PRE,0,'bin');
    [Zscrd_DATA_ST,~] = ParaZScore(Data_STBased, para, SINKs,PRE,0,'bin');
    Data_GSBased_BU = Data_GSBased; Data_GSBased = Zscrd_DATA_GS;
    Data_STBased_BU = Data_STBased; Data_STBased = Zscrd_DATA_ST;
    end
    
    cd (home), cd figs;
    mkdir(['Group_' input(i1).name(1:end-2)]);
    
    %% Plot ST shift single animals
    cd (home), cd subfunc;
    SParam = {'SinkRMS' 'tempSinkRMS'};    
    plotshift_single(Data,names,SINKs,ST_Shift,NumAnimals,home,CurAn,SParam)

    %% Plot ST_Shift group
    cd (home), cd subfunc;
    plotshift_group(Data,names,SINKs,ST_Shift,home,CurAn)

    %% Plot GS tuning of single animals
    cd (home), cd subfunc;
    plotGS_single(Data,Data_GSBased,names,SINKs,para,zscore,nThresh,bin,mirror,home,CurAn,BF_Pos)

    %% Plot ST tuning single animal
    cd (home), cd subfunc;
    plotST_single(Data,Data_STBased,names,SINKs,para,zscore,nThresh,bin,mirror,home,CurAn,BF_Pos)
    
    %% Plot Group Tuning Curves ST based
    cd (home), cd subfunc;
    plotST_group(Data,Data_STBased,names,SINKs,para,zscore,nThresh,bin,mirror,home,CurAn)
    
    %% Group TUning Curves GS based  
    cd (home), cd subfunc;
    plotGS_group(Data,Data_GSBased,names,SINKs,para,zscore,nThresh,bin,mirror,home,CurAn)

    %% Pass on without any sorting
    PASS = struct;   
    for i3= 1:length(pass_on)
        PASS.(pass_on{i3})= [];
    end
    for i3= 1:length(pass_on)
        
        for i4 = 1:size(Data,2)
            
            clear curPar
            try
                curPar= nan(BF_Pos,size(Data(i4).(names{1}).(pass_on{i3}),2),size(Data,2));
            catch
            end
            for i5 = 1:length(names)
                try
                    curPar(1:size(Data(i4).(names{i5}).(pass_on{i3}),1),...
                        1:size(Data(i4).(names{1}).(pass_on{i3}),2),i5)=...
                        Data(i4).(names{i5}).(pass_on{i3})(:,1:size(Data(i4).(names{1}).(pass_on{i3}),2));
                catch
                    try
                        curPar(1:size(Data(i4).(names{i5}).(pass_on{i3}),1),...
                            1:size(Data(i4).(names{i5}).(pass_on{i3}),2),i5)=...
                            Data(i4).(names{i5}).(pass_on{i3})(:,1:size(Data(i4).(names{i5}).(pass_on{i3}),2));
                    end
                end
            end
            try PASS(i4).(pass_on{i3})= curPar;
            end
        end
        
    end
 
 %% Save and Exit
 cd(home); cd DATA; clear Data;
 Data.GS_based = Data_GSBased;
 Data.ST_based = Data_STBased;
 Data.Shifts = ST_Shift;
 Data.pass_on =PASS;
 Data.names = names;
 Data.BF_Pos = BF_Pos;
 Data.Sinks = SINKs;
 
 cd Output
 
 save(['Output_' input(i1).name(1:end-4) '_Threshold_' num2str(nThresh*100) '_Zscore_' num2str(zscore) '_binned_' num2str(bin) '_mirror_' num2str(mirror)],'Data')
 clear Data_GSBased Data Data_STBased ST_Shift
 cd (home); cd Data;
 
 cd (home);
end

