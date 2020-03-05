%% Dynamic CSD for sinks I_II through VIb
%   This script takes input from the groups and raw folders. It calculates 
%   and stores CSD information in Data struct which is
%   saved in the DATA folder. 
% 
%   6a and 6b are detected here and there is no single trial data
%   IMPORTANT: DO NOT change sink list here. If you need another set of
%   sinks then create a SEPARATE and UNIQUELY NAMED script. 
%   Note: we are keeping scripts that do not run single trial stuff for now
%   because it is faster without. We will eventually remove these if we
%   permenantly move to using single-trial data across the board
% 
%   CHANGE if needed: add your working directory to the try/catch; change
%   condition to run 

clear;
%% Choose Condition

% Condition = {'Pre' 'post' 'att' 'FM' 'between' 'muscimol' 'musatt' 'musFM'};
Condition = {'condition'};
% Condition = {'Pre' 'post'}; %MUST HAVE CAPITAL P FOR STRING MATCH IN NEXT CODE

%% standard operations
warning('OFF');
dbstop if error %stops execution of program if and allows user to examine the workspace

% Change directory to your working folder
if exist('D:\MyCode\Dynamic_CSD_Analysis','dir') == 7
    cd('D:\MyCode\Dynamic_CSD_Analysis');
elseif exist('D:\Dynamic_CSD_Analysis','dir') == 7
    cd('D:\Dynamic_CSD_Analysis');
elseif exist('C:\Users\kedea\Documents\Dynamic_CSD_Analysis','dir') == 7
    cd('C:\Users\kedea\Documents\Dynamic_CSD_Analysis')
end

home = pwd; %pwd: print working directory, gives full current directory title of "home"
addpath(genpath(home)); %adds all subdirectories of 'home' to path for access
cd (home),cd groups; 

%% Load in 
input = dir('*.m'); %dir *.m lists all program files in the current directory (stores 5 fields of data but struc array contains #x1)
entries = length(input); %number of input entries

    
for i1 = 1:entries
    run(input(i1).name);
    
    cd (home); cd figs; %opens figure folder to store future images
    mkdir(['Single_' input(i1).name(1:end-2)]); %adds folder to directory - name without '.m' (-2 from the end)
    
    Data = struct; %creating empty structure to fill
    Indexer = struct; %creating empty structure to fill
    

    %% Condition and Indexer
    for iC = 1:length(Condition) %input from script (length of condition: 'Pre' 'Combo' 'Post1'=3) %i2 is condition we are in
        for iA = 1:length(animals) %input from groups (in workspace) %i3 is animal we are on
            counter = 0;
            counter = counter + length(Cond.(Condition{iC}){iA}); %Cond has all conditions recorded, condition(current cond type selection for current animal selection)
            % counter = the number of entries in current condition for the current animal
            if iA == 1 %for the first animal
                Indexer(1).(Condition{iC}) = counter; %the indexer structure is the size of the current counter
            elseif Indexer.(Condition{iC}) < counter %if the next animal has a structure larger than the current counter, it will make the indexer that size
                Indexer(1).(Condition{iC}) = counter;
            end
            %An indexer is now made for the i2th condition that is the size of the largest amongst the animals
        end
        % Indexer now contains all conditions at their greatest sizes
    end
    
    count = 1; %new variable 'count'
    for iC = 1:length(Condition) %i2 is still which condition we are in and it's reset to 1
        if iC == 1 %first condition type
            Indexer(2).(Condition{iC}) =1; %adding a second field to the condition i2
        else
            Indexer(2).(Condition{iC}) = count;
        end
        count = count + Indexer(1).(Condition{iC}); %current count (e.g. 1)+size of condition i2 in index(e.g. 3) = (e.g. 4)
    end
    %now the indexer has 2 columns beyond the cond names, indicating the condition size in the first and the starting count in the second
    
    %%
    for iA = 1:length(animals) 
        cd (home); cd raw;    
        name = animals{iA}; %assigns the i2th entry in the list in the current group to variable 'name'
        
        for iC = 1:length(Condition) 
            clear SWEEP;
            for i4 = 1:length(Cond.(Condition{iC}){iA}) %i4 goes through the number of entries in condition i3 for animal i4
                if i4 == 1 
                    CondIDX = Indexer(2).(Condition{iC}); %new count variable calling the start count number for the first condition entry (which equals 1)
                else 
                    CondIDX = Indexer(2).(Condition{iC})+i4-1;
                end
                
                measurement = Cond.(Condition{iC}){iA}{i4}; 
                %e.g. 003 is the first entry of the first animal in the first condition
                
                if ~isempty(measurement) % ~ is the logical argument 'not' so this says if the measurement variable is not empty than true
                    disp(['Analyzing animal: ' name '_' measurement])
                    tic
                    
                    try 
                        load ([name '_' measurement]);clear avgFP;   %loads info from mat file (Data, DS, Header, P, Spik_list, SWEEP)
                    catch
                        fprintf('the name or measurement does not exist/n')
                    end
                    
                    cd (home),cd groups;
                    try load([name '_Baseline']); %this Baseline is determined by gen_threshold.m from multiple recordings of one animal
                    catch
                        Baseline = []; %if a baseline wasn't taken, create an empty variable
                    end
                    
                    % all of the above is to indicate which animal and
                    % condition is being analyzed
                    if exist('SWEEP','var') %Acute recordings
                        clear DATA frqz
                        BL = Header.t_pre*P.Fs_AD(1); %BL-baseline %t_pre is the time before the tone %Fs_AD - Sampling frequency of channels (they are all the same so we use first value)
                        tone = Header.t_sig(1)*P.Fs_AD(1); %t_sig is duration of stimulus * sampling rate = 200
                        frqz = Header.stimlist(:,1); %stimlist contains tone frequencies in all rows (:), first column (:,1)
                        Fs = P.Fs_AD(1); %sampling rate
                        frqz(find(frqz == inf))=[];
                        frqz(find(frqz == 0))=[];
                    else %Chronic recordings
                        cd (home); cd subfunc;
                        %To apply channel interpolation
                        if ~isempty(DATA.ART)
                            [DATA.LFP]=chaninterp(DATA.LFP, 'linextra', DATA.ART.chan,...
                                DATA.ELECPOS);
                        end
                        %To select attenuation specific recordings
                        if ~isempty(dB_lev{iA})
                            if strcmp(dB_lev{iA},'10dB')
                                [DATA]=attenuation10(DATA);
                            elseif strcmp(dB_lev{iA},'30dB')
                                [DATA]=attenuation30(DATA);
                            end
                        end
                        %To create SWEEP
                        [BL, tone, Fs, SWEEP, frqz]=reformat(DATA);
                    end
                    
                    %% CSD full
                    cd (home),cd subfunc;

                    [AvgFP, AvgCSD,AvgRecCSD, SingleTrialAvgRecCSD, AvgRelResCSD,...
                        SingleTrialAvgRelResCSD, AvgAbsResCSD, AvgLayerRelRes] =...
                        SingleTrialCSD_reduced(SWEEP, str2num(channels{iA}),BL);

                    
                    
                    %delete empty columns to have the correct amount of stimuli
                    %present (needed for attenuation 30)
                    SingleTrialAvgRecCSD = SingleTrialAvgRecCSD(~cellfun('isempty',SingleTrialAvgRecCSD'));
                    SingleTrialAvgRelResCSD = SingleTrialAvgRelResCSD(~cellfun('isempty',SingleTrialAvgRelResCSD'));
                    AvgCSD = AvgCSD(~cellfun('isempty', AvgCSD'));
                    AvgFP = AvgFP(~cellfun('isempty', AvgFP'));
                    AvgRecCSD = AvgRecCSD(~cellfun('isempty', AvgRecCSD')); AvgCSD=AvgCSD(1:length(frqz)); 
                    AvgRelResCSD = AvgRelResCSD(~cellfun('isempty', AvgRelResCSD'));
                    
                    % Sink durations
                    L.I_IIE = str2num(Layer.I_IIE{iA}); L.I_IIL = str2num(Layer.I_IIL{iA}); 
                    L.IVE = str2num(Layer.IVE{iA}); L.IVL = str2num(Layer.IVL{iA});
                    L.VaE = str2num(Layer.VaE{iA}); L.VaL = str2num(Layer.VaL{iA}); 
                    L.VbE = str2num(Layer.VbE{iA}); L.VbL = str2num(Layer.VbL{iA});
                    L.VIaE = str2num(Layer.VIaE{iA}); L.VIaL = str2num(Layer.VIaL{iA});
                    L.VIbE = str2num(Layer.VIbE{iA}); L.VIbL = str2num(Layer.VIbL{iA});
                    
                    curChan = str2num(channels{iA});
                    
                    %Generate Sink Boxes
                    [DUR,RMS,PAMP,PLAT,INT] = sink_dura_group(L,AvgCSD,BL,tone,SWEEP,curChan,Baseline);
                    
                    %Calculate CSD
                    %columns of starting and ending values minus the baseline
                    DUR2.I_IIE = reshape(round([DUR(:).I_IIE]),[],length([DUR(:).I_IIE])/2)'-BL;
                    DUR2.I_IIL = reshape(round([DUR(:).I_IIL]),[],length([DUR(:).I_IIL])/2)'-BL;
                    DUR2.IVE = reshape(round([DUR(:).IVE]),[],length([DUR(:).IVE])/2)'-BL;
                    DUR2.IVL = reshape(round([DUR(:).IVL]),[],length([DUR(:).IVL])/2)'-BL;
                    DUR2.VaE = reshape(round([DUR(:).VaE]),[],length([DUR(:).VaE])/2)'-BL;
                    DUR2.VaL = reshape(round([DUR(:).VaL]),[],length([DUR(:).VaL])/2)'-BL;
                    DUR2.VbE = reshape(round([DUR(:).VbE]),[],length([DUR(:).VbE])/2)'-BL;
                    DUR2.VbL = reshape(round([DUR(:).VbL]),[],length([DUR(:).VbL])/2)'-BL;
                    DUR2.VIaE = reshape(round([DUR(:).VIaE]),[],length([DUR(:).VIaE])/2)'-BL;
                    DUR2.VIaL = reshape(round([DUR(:).VIaL]),[],length([DUR(:).VIaL])/2)'-BL;
                    DUR2.VIbE = reshape(round([DUR(:).VIbE]),[],length([DUR(:).VIbE])/2)'-BL;
                    DUR2.VIbL = reshape(round([DUR(:).VIbL]),[],length([DUR(:).VIbL])/2)'-BL;
                    
                    %rows of values for onsets and offsets
                    onset.I_IIE = DUR2.I_IIE(:,1)';   offset.I_IIE = DUR2.I_IIE(:,2)';
                    onset.I_IIL = DUR2.I_IIL(:,1)';   offset.I_IIL = DUR2.I_IIL(:,2)';
                    onset.IVE = DUR2.IVE(:,1)';       offset.IVE = DUR2.IVE(:,2)';
                    onset.IVL = DUR2.IVL(:,1)';       offset.IVL = DUR2.IVL(:,2)';
                    onset.VaE = DUR2.VaE(:,1)';       offset.VaE = DUR2.VaE(:,2)';
                    onset.VaL = DUR2.VaL(:,1)';       offset.VaL = DUR2.VaL(:,2)';
                    onset.VbE = DUR2.VbE(:,1)';       offset.VbE = DUR2.VbE(:,2)';
                    onset.VbL = DUR2.VbL(:,1)';       offset.VbL = DUR2.VbL(:,2)';
                    onset.VIaE = DUR2.VIaE(:,1)';     offset.VIaE = DUR2.VIaE(:,2)';
                    onset.VIaL = DUR2.VIaL(:,1)';     offset.VIaL = DUR2.VIaL(:,2)';
                    onset.VIbE = DUR2.VIbE(:,1)';     offset.VIbE = DUR2.VIbE(:,2)';
                    onset.VIbL = DUR2.VIbL(:,1)';     offset.VIbL = DUR2.VIbL(:,2)';
                    
                    %subracts onset value from offset value to recreate duration (total time of sink)
                    DUR2.I_IIE = DUR2.I_IIE(:,2)' - DUR2.I_IIE(:,1)';
                    DUR2.I_IIL = DUR2.I_IIL(:,2)' - DUR2.I_IIL(:,1)';
                    DUR2.IVE = DUR2.IVE(:,2)' - DUR2.IVE(:,1)';
                    DUR2.IVL = DUR2.IVL(:,2)' - DUR2.IVL(:,1)';
                    DUR2.VaE = DUR2.VaE(:,2)' - DUR2.VaE(:,1)';
                    DUR2.VaL = DUR2.VaL(:,2)' - DUR2.VaL(:,1)';
                    DUR2.VbE = DUR2.VbE(:,2)' - DUR2.VbE(:,1)';
                    DUR2.VbL = DUR2.VbL(:,2)' - DUR2.VbL(:,1)';
                    DUR2.VIaE = DUR2.VIaE(:,2)' - DUR2.VIaE(:,1)';
                    DUR2.VIaL = DUR2.VIaL(:,2)' - DUR2.VIaL(:,1)';
                    DUR2.VIbE = DUR2.VIbE(:,2)' - DUR2.VIbE(:,1)';
                    DUR2.VIbL = DUR2.VIbL(:,2)' - DUR2.VIbL(:,1)';
                    
                    threshold_std = 2; 
                    threshold_dur = 0.005;
                    Latency_HighCutoff = 0.2;
                    Latency_LowCutoff = 0.015;
                    state = 1:(length(AvgRecCSD{1})-BL); %KD this idicates the length after the stimulus 
                    

                    %calculate full RMS of AVREC and RELRES (1-max Duration)
                    [full_BW,~, ~, full_RMS_AvgRecCSD,~, ~,~, ~] =...
                        ExtractCSDBasedPar(1,AvgRecCSD,BL,...
                        Fs, 1, 0, threshold_std, threshold_dur,...
                        Latency_HighCutoff, Latency_LowCutoff,state);
                    
                    [~,~, ~, full_RMS_RelResCSD,~, ~, ~] =...
                        ExtractCSDBasedPar(1,AvgRelResCSD,BL,...
                        Fs, 0, 0, threshold_std, threshold_dur,...
                        Latency_HighCutoff, Latency_LowCutoff,state);
                    
                    [~,~, ~, full_RMS_AbsResCSD,~, ~, ~] =...
                        ExtractCSDBasedPar(1,AvgAbsResCSD,BL,...
                        Fs, 0, 0, threshold_std, threshold_dur,...
                        Latency_HighCutoff, Latency_LowCutoff,state);
                    
                    %calculate early RMS of AVREC and RELRES (1-50 ms)
                    state = 1:50;         
                    [e_BW,~, ~, e_RMS_AvgRecCSD,~, ~, ~, ~] =...
                        ExtractCSDBasedPar(1,AvgRecCSD,BL,...
                        Fs, 1, 0, threshold_std, threshold_dur,...
                        Latency_HighCutoff, Latency_LowCutoff,state);
                    
                    [~,~, ~, e_RMS_RelResCSD, ~, ~, ~] =...
                        ExtractCSDBasedPar(1,AvgRelResCSD,BL,...
                        Fs, 0, 0, threshold_std, threshold_dur,...
                        Latency_HighCutoff, Latency_LowCutoff,state);
                    
                    [~,~, ~, e_RMS_AbsResCSD, ~, ~, ~] =...
                        ExtractCSDBasedPar(1,AvgAbsResCSD,BL,...
                        Fs, 0, 0, threshold_std, threshold_dur,...
                        Latency_HighCutoff, Latency_LowCutoff,state);
                    
                    %calculate late RMS of AVREC and RELRES (81-300 ms)
                    state = 81:300;
                    [l_BW,~, ~, l_RMS_AvgRecCSD,~, ~, ~, ~] =...
                        ExtractCSDBasedPar(1,AvgRecCSD,BL,...
                        Fs, 1, 0, threshold_std, threshold_dur,...
                        Latency_HighCutoff, Latency_LowCutoff,state);
                    
                    [~,~, ~, l_RMS_RelResCSD, ~, ~, ~] =...
                        ExtractCSDBasedPar(1,AvgRelResCSD,BL,...
                        Fs, 0, 0, threshold_std, threshold_dur,...
                        Latency_HighCutoff, Latency_LowCutoff,state);
                    
                    [~,~, ~, l_RMS_AbsResCSD, ~, ~, ~] =...
                        ExtractCSDBasedPar(1,AvgAbsResCSD,BL,...
                        Fs, 0, 0, threshold_std, threshold_dur,...
                        Latency_HighCutoff, Latency_LowCutoff,state);

                    toc

                    %rename data for plots
                    if ~isempty(dB_lev{iA})
                        name = [animals{iA} '_' dB_lev{iA}];
                    end
                    
                    %% BANDWIDTH and TUNINGWIDTH
                    rmscurve = [];
                    
                    for iB = 1:size(AvgCSD,2)-1
                        rawCSD = (nanmean(AvgCSD{iB}(:,200:300)))*-1;
                        rootms = rms(rawCSD);
                        rmscurve = [rmscurve rootms];
                    end
                    
                    halfMax = (min(rmscurve) + max(rmscurve)) / 2;
                    bandwidth = length(find(rmscurve >= halfMax));
                    
                    firststd = std(rmscurve);
                    tuningwidth = length(find(rmscurve >= firststd));
                                                          
                    %% Plots and Figures %%
                    figure('Name',[name ' ' measurement ': ' Condition{iC}])
                    tic
                    disp('Plotting CSD with sink detections')
                    for i5 = 1:length(AvgCSD)
                        subplot(2,round(length(AvgCSD)/2),i5)
                        imagesc(AvgCSD{i5})
                        if i5 == 1 %full title for first image
                            title ([name ' ' measurement ': ' Condition{iC} ' ' num2str(i5) ' ' num2str(frqz(i5)) ' Hz'])
                        else %abbreviated version for images to follow on same figure
                            title ([num2str(frqz(i5)) ' Hz'])
                        end
                        
                        colormap (jet)
                        caxis([-0.0005 0.0005])

                        
                        hold on
                        % Layer I_IIE
                        y =[(max(L.I_IIE)+0.5),(max(L.I_IIE)+0.5),(min(L.I_IIE)-0.5),(min(L.I_IIE)-0.5),(max(L.I_IIE)+0.5)];
                        if isempty(y); y = [NaN NaN NaN NaN NaN]; end
                        x = [DUR(i5).I_IIE(1), DUR(i5).I_IIE(2),DUR(i5).I_IIE(2),DUR(i5).I_IIE(1),DUR(i5).I_IIE(1)];
                        plot(x,y,'black','LineWidth',2)
                        
                        % Layer I_IIL
                        y =[(max(L.I_IIL)+0.5),(max(L.I_IIL)+0.5),(min(L.I_IIL)-0.5),(min(L.I_IIL)-0.5),(max(L.I_IIL)+0.5)];
                        if isempty(y); y = [NaN NaN NaN NaN NaN]; end
                        x = [DUR(i5).I_IIL(1), DUR(i5).I_IIL(2),DUR(i5).I_IIL(2),DUR(i5).I_IIL(1),DUR(i5).I_IIL(1)];
                        plot(x,y,'white','LineWidth',2)
                        
                        % Layer IVE
                        y =[(max(L.IVE)+0.5),(max(L.IVE)+0.5),(min(L.IVE)-0.5),(min(L.IVE)-0.5),(max(L.IVE)+0.5)];
                        x = [DUR(i5).IVE(1), DUR(i5).IVE(2),DUR(i5).IVE(2),DUR(i5).IVE(1),DUR(i5).IVE(1)];
                        plot(x,y,'black','LineWidth',2)
                        
                        % Layer IVL
                        y =[(max(L.IVL)+0.5),(max(L.IVL)+0.5),(min(L.IVL)-0.5),(min(L.IVL)-0.5),(max(L.IVL)+0.5)];
                        x = [DUR(i5).IVL(1), DUR(i5).IVL(2),DUR(i5).IVL(2),DUR(i5).IVL(1),DUR(i5).IVL(1)];
                        plot(x,y,'white','LineWidth',2)
                        
                        % Layer VaE
                        y =[(max(L.VaE)+0.5),(max(L.VaE)+0.5),(min(L.VaE)-0.5),(min(L.VaE)-0.5),(max(L.VaE)+0.5)];
                        x = [DUR(i5).VaE(1), DUR(i5).VaE(2),DUR(i5).VaE(2),DUR(i5).VaE(1),DUR(i5).VaE(1)];
                        plot(x,y,'black','LineWidth',2)
                        
                        % Layer VaL
                        y =[(max(L.VaL)+0.5),(max(L.VaL)+0.5),(min(L.VaL)-0.5),(min(L.VaL)-0.5),(max(L.VaL)+0.5)];
                        x = [DUR(i5).VaL(1), DUR(i5).VaL(2),DUR(i5).VaL(2),DUR(i5).VaL(1),DUR(i5).VaL(1)];
                        plot(x,y,'white','LineWidth',2)
                        
                        % Layer VbE
                        y =[(max(L.VbE)+0.5),(max(L.VbE)+0.5),(min(L.VbE)-0.5),(min(L.VbE)-0.5),(max(L.VbE)+0.5)];
                        x = [DUR(i5).VbE(1), DUR(i5).VbE(2),DUR(i5).VbE(2),DUR(i5).VbE(1),DUR(i5).VbE(1)];
                        plot(x,y,'black','LineWidth',2)
                        
                        % Layer VbL
                        y =[(max(L.VbL)+0.5),(max(L.VbL)+0.5),(min(L.VbL)-0.5),(min(L.VbL)-0.5),(max(L.VbL)+0.5)];
                        x = [DUR(i5).VbL(1), DUR(i5).VbL(2),DUR(i5).VbL(2),DUR(i5).VbL(1),DUR(i5).VbL(1)];
                        plot(x,y,'white','LineWidth',2)
                        
                        % Layer VIaE
                        y =[(max(L.VIaE)+0.5),(max(L.VIaE)+0.5),(min(L.VIaE)-0.5),(min(L.VIaE)-0.5),(max(L.VIaE)+0.5)];
                        x = [DUR(i5).VIaE(1), DUR(i5).VIaE(2),DUR(i5).VIaE(2),DUR(i5).VIaE(1),DUR(i5).VIaE(1)];
                        plot(x,y,'black','LineWidth',2)
                        
                        % Layer VIaL
                        y =[(max(L.VIaL)+0.5),(max(L.VIaL)+0.5),(min(L.VIaL)-0.5),(min(L.VIaL)-0.5),(max(L.VIaL)+0.5)];
                        x = [DUR(i5).VIaL(1), DUR(i5).VIaL(2),DUR(i5).VIaL(2),DUR(i5).VIaL(1),DUR(i5).VIaL(1)];
                        plot(x,y,'white','LineWidth',2)
                        
                        % Layer VIbE
                        y =[(max(L.VIbE)+0.5),(max(L.VIbE)+0.5),(min(L.VIbE)-0.5),(min(L.VIbE)-0.5),(max(L.VIbE)+0.5)];
                        x = [DUR(i5).VIbE(1), DUR(i5).VIbE(2),DUR(i5).VIbE(2),DUR(i5).VIbE(1),DUR(i5).VIbE(1)];
                        plot(x,y,'black','LineWidth',2)
                        
                        % Layer VIbL
                        y =[(max(L.VIbL)+0.5),(max(L.VIbL)+0.5),(min(L.VIbL)-0.5),(min(L.VIbL)-0.5),(max(L.VIbL)+0.5)];
                        x = [DUR(i5).VIbL(1), DUR(i5).VIbL(2),DUR(i5).VIbL(2),DUR(i5).VIbL(1),DUR(i5).VIbL(1)];
                        plot(x,y,'white','LineWidth',2)
                        
                        hold off
                    end
                    
                    cd([home '\figs\' 'Single_' input(i1).name(1:end-2)])
                    h = gcf; %returns the handle of the current figure
                    set(h, 'PaperType', 'A4');
                    set(h, 'PaperOrientation', 'landscape');
                    set(h, 'PaperUnits', 'centimeters');
                    savefig(h,[name '_' measurement '_CSDs' ],'compact')
                    close (h)
                    
                    Order = {'I_IIE','IVE','VaE','VbE','VIaE','VIaE',...
                             'I_IIL','IVL','VaL','VbL','VIaL','VIaL'};
                    h = figure('Name',[name ' ' measurement ': ' Condition{iC} ' ' num2str(i4)]);
                    for i5 = 1:(length(Order)/2)
                        subplot(2,(round((length(Order)/2)/2)),i5)
                        time = [DUR(:).(Order{i5})];
                        time = reshape(time,[],length(time)/2);
                        time = time-BL;
                        hold on
                        plot (time(1,:),'b--o','LineWidth',2);
                        plot (time(2,:),'r--o','LineWidth',2);
                        hold off
                        title(Order{i5})
                        ylim([0 (length(AvgCSD{1})-BL)])
                        xlim([0 (length(AvgCSD)+1)])
                        ylabel ('Time in ms')
                        xlabel ('Index of Stimuli')
                        if i5 == 1
                            legend('Sink onset','Sink offset','Location','best')
                        end
                    end
                    set(h, 'PaperType', 'A4');
                    set(h, 'PaperOrientation', 'landscape');
                    set(h, 'PaperUnits', 'centimeters');
                    savefig(h,[name '_' measurement '_Sink_onset+offset_Early' ],'compact')
                    close (h)
                    
                    h = figure('Name',[name ' ' measurement ': ' Condition{iC} ' ' num2str(i4)]);
                    for i5 = (length(Order)/2)+1:(length(Order)/2)*2
                        subplot(2,(round((length(Order)/2)/2)),(i5-(length(Order)/2)))
                        time = [DUR(:).(Order{i5})];
                        time = reshape(time,[],length(time)/2);
                        time = time-BL;
                        hold on
                        plot (time(1,:),'b--o','LineWidth',2);
                        plot (time(2,:),'r--o','LineWidth',2);
                        hold off
                        title(Order{i5})
                        ylim([0 (length(AvgCSD{1})-BL)])
                        xlim([0 (length(AvgCSD)+1)])
                        ylabel ('Time in ms')
                        xlabel ('Index of Stimuli')
                        if i5 == 1
                            legend('Sink onset','Sink offset','Location','best')
                        end
                    end
                    set(h, 'PaperType', 'A4');
                    set(h, 'PaperOrientation', 'landscape');
                    set(h, 'PaperUnits', 'centimeters');
                    savefig(h,[name '_' measurement '_Sink_onset+offset_Late' ],'compact')
                    close (h)
                    
                    BF = find([RMS.IVE] == max([RMS.IVE])); 
                    BF_frqz = frqz(BF);
                    
                    if isnan(BF), BF = 1; end
                    if isempty(BF), BF = 1; end
                    
                    SINK = {'I_IIE','IVE','VaE','VbE','VIaE','VIbE',...
                             'I_IIL','IVL','VaL','VbL','VIaL','VIbL'};
                    h = figure('Name',[name ' ' measurement ': ' Condition{iC} ' ' num2str(i4)]);
                    for i5 = 1:(length(Order)/2)
                        subplot(3,(length(Order)/2),i5)
                        durtime = DUR2.(SINK{i5});
                        plot(durtime,'b--o','LineWidth',2)
                        hold on
                        plot (BF,durtime(BF),'r*','LineWidth',2)
                        hold off
                        xlim([0 (length(AvgCSD)+1)])
                        title (['Sink duration Layer ' SINK{i5} ' ms'])
                        
                        subplot(3,(length(Order)/2),i5+(length(Order)/2))
                        plot([RMS.(SINK{i5})],'b--o','LineWidth',2)
                        hold on
                        plot (BF,RMS(BF).(SINK{i5}),'r*','LineWidth',2)
                        hold off
                        xlim([0 (length(AvgCSD)+1)])
                        title (['RMS Layer ' SINK{i5} ' V/mm?'])
                        
                        subplot(3,(length(Order)/2),i5+((length(Order)/2)*2))
                        relRMS = [RMS.(SINK{i5})]./durtime;
                        plot(relRMS,'b--o','LineWidth',2)
                        hold on
                        plot (BF,relRMS(BF),'r*','LineWidth',2)
                        hold off
                        xlim([0 (length(AvgCSD)+1)])
                        title (['RMS / ms L ' SINK{i5} ' mV/mm? ms'])
                        RMSTIME.(SINK{i5}) = relRMS;
                        
                    end
                    
                    set(h, 'PaperType', 'A4');
                    set(h, 'PaperOrientation', 'landscape');
                    set(h, 'PaperUnits', 'centimeters');
                    savefig(h,[name '_' measurement '_Duration_RMS_Ratio_Early' ],'compact')
                    close (h)
                    
                    h = figure('Name',[name ' ' measurement ': ' Condition{iC} ' ' num2str(i4)]);
                    for i5 = (length(Order)/2)+1:(length(Order)/2)*2
                        subplot(3,(length(Order)/2),i5-(length(Order)/2))
                        durtime = DUR2.(SINK{i5});
                        plot(durtime,'b--o','LineWidth',2)
                        hold on
                        plot (BF,durtime(BF),'r*','LineWidth',2)
                        hold off
                        xlim([0 (length(AvgCSD)+1)])
                        title (['Sink duration Layer ' SINK{i5} ' ms'])
                        
                        subplot(3,(length(Order)/2),i5)
                        plot([RMS.(SINK{i5})],'b--o','LineWidth',2)
                        hold on
                        plot (BF,RMS(BF).(SINK{i5}),'r*','LineWidth',2)
                        hold off
                        xlim([0 (length(AvgCSD)+1)])
                        title (['RMS Layer ' SINK{i5} ' V/mm?'])
                        
                        subplot(3,(length(Order)/2),i5+(length(Order)/2))
                        relRMS = [RMS.(SINK{i5})]./durtime;
                        plot(relRMS,'b--o','LineWidth',2)
                        hold on
                        plot (BF,relRMS(BF),'r*','LineWidth',2)
                        hold off
                        xlim([0 (length(AvgCSD)+1)])
                        title (['RMS / ms L ' SINK{i5} ' mV/mm? ms'])
                        RMSTIME.(SINK{i5}) = relRMS;
                        
                    end
                    
                    set(h, 'PaperType', 'A4');
                    set(h, 'PaperOrientation', 'landscape');
                    set(h, 'PaperUnits', 'centimeters');
                    savefig(h,[name '_' measurement '_Duration_RMS_Ratio_Late' ],'compact')
                    close (h)
                    
                    % BF normalized
                    h = figure('Name',[name ' ' measurement ': ' Condition{iC} ' ' num2str(i4)]);
                    for i5 = 1:(length(Order)/2)
                        subplot(3,(length(Order)/2),i5)
                        durtime = DUR2.(SINK{i5});
                        plot(durtime./durtime(BF),'LineWidth',2)
                        hold on
                        plot (BF,durtime(BF)./durtime(BF),'r*','LineWidth',2)
                        hold off
                        xlim([0 (length(AvgCSD)+1)])
                        title (['Sink duration Layer ' SINK{i5} ' ms'])
                        
                        subplot(3,(length(Order)/2),i5+(length(Order)/2))
                        RS =[RMS.(SINK{i5})];
                        plot(RS./RS(BF),'LineWidth',2)
                        hold on
                        plot (BF,RMS(BF).(SINK{i5})/RS(BF),'r*','LineWidth',2)
                        hold off
                        xlim([0 (length(AvgCSD)+1)])
                        title (['RMS Layer ' SINK{i5} ' mV/mm?'])
                        
                        subplot(3,(length(Order)/2),i5+(length(Order)/2))
                        relRMS = [RMS.(SINK{i5})].*durtime;
                        plot(relRMS./relRMS(BF),'LineWidth',2)
                        hold on
                        plot (BF,relRMS(BF)/relRMS(BF),'r*','LineWidth',2)
                        hold off
                        xlim([0 (length(AvgCSD)+1)])
                        title (['RMS * ms L ' SINK{i5} ' mV/mm? ms'])
                        
                    end
                    
                    set(h, 'PaperType', 'A4');
                    set(h, 'PaperOrientation', 'landscape');
                    set(h, 'PaperUnits', 'centimeters');
                    savefig(h,[name '_' measurement '_Duration_RMS_Ratio_BF_normalized_Early' ],'compact')
                    close (h)
                    
                    h = figure('Name',[name ' ' measurement ': ' Condition{iC} ' ' num2str(i4)]);
                    for i5 = (length(Order)/2)+1:(length(Order)/2)*2
                        subplot(3,(length(Order)/2),i5-(length(Order)/2))
                        durtime = DUR2.(SINK{i5});
                        plot(durtime./durtime(BF),'LineWidth',2)
                        hold on
                        plot (BF,durtime(BF)./durtime(BF),'r*','LineWidth',2)
                        hold off
                        xlim([0 (length(AvgCSD)+1)])
                        title (['Sink duration Layer ' SINK{i5} ' ms'])
                        
                        subplot(3,(length(Order)/2),i5)
                        RS =[RMS.(SINK{i5})];
                        plot(RS./RS(BF),'LineWidth',2)
                        hold on
                        plot (BF,RMS(BF).(SINK{i5})/RS(BF),'r*','LineWidth',2)
                        hold off
                        xlim([0 (length(AvgCSD)+1)])
                        title (['RMS Layer ' SINK{i5} ' mV/mm?'])
                        
                        subplot(3,(length(Order)/2),i5+(length(Order)/2))
                        relRMS = [RMS.(SINK{i5})].*durtime;
                        plot(relRMS./relRMS(BF),'LineWidth',2)
                        hold on
                        plot (BF,relRMS(BF)/relRMS(BF),'r*','LineWidth',2)
                        hold off
                        xlim([0 (length(AvgCSD)+1)])
                        title (['RMS * ms L ' SINK{i5} ' mV/mm? ms'])
                        
                    end
                    
                    set(h, 'PaperType', 'A4');
                    set(h, 'PaperOrientation', 'landscape');
                    set(h, 'PaperUnits', 'centimeters');
                    savefig(h,[name '_' measurement '_Duration_RMS_Ratio_BF_normalized_Late' ],'compact')
                    close (h)
                    
                    
                    h = figure('Name',[name ' ' measurement ': ' Condition{iC} ' ' num2str(i4)]);
                    for i5 = 1: length(SINK)
                        subplot(3,length(SINK),i5)
                        durtime = DUR2.(SINK{i5});
                        plot(durtime,'b--o','LineWidth',2)
                        hold on
                        plot (BF,durtime(BF),'r*','LineWidth',2)
                        hold off
                        xlim([0 (length(AvgCSD)+1)])
                        title (['Sink duration Layer ' SINK{i5} ' ms'])
                        
                        subplot(3,length(SINK),i5+length(SINK))
                        plot([INT.(SINK{i5})],'b--o','LineWidth',2)
                        hold on
                        plot (BF,INT(BF).(SINK{i5}),'r*','LineWidth',2)
                        hold off
                        xlim([0 (length(AvgCSD)+1)])
                        title (['Current flow Layer ' SINK{i5} ' mV*ms/mm?'])
                        
                        subplot(3,length(SINK),i5+2*length(SINK))
                        relRMS = [INT.(SINK{i5})]./durtime;
                        plot(relRMS,'b--o','LineWidth',2)
                        hold on
                        plot (BF,relRMS(BF),'r*','LineWidth',2)
                        hold off
                        xlim([0 (length(AvgCSD)+1)])
                        title (['mean CF L ' SINK{i5} ' mV/mm?'])
                        RMSTIME.(SINK{i5}) = relRMS;
                        
                    end
                    
                    set(h, 'PaperType', 'A4');
                    set(h, 'PaperOrientation', 'landscape');
                    set(h, 'PaperUnits', 'centimeters');
                    savefig(h,[name '_' measurement '_Current_flow' ],'compact')
                    close (h)
                    
                    
                    all = [AvgRecCSD{:}];
                    all = all(BL+15:BL+50,:);
                    initPeakTune = max(all);
                    
                    %% Save and Quit
                    Data(CondIDX).(name).measurement =[name '_' measurement];
                    Data(CondIDX).(name).Condition = [Condition{iC} '_' num2str(i4)];
                    Data(CondIDX).(name).BL = BL;
                    Data(CondIDX).(name).StimDur = tone;
                    Data(CondIDX).(name).Frqz = frqz';
                    Data(CondIDX).(name).GS_BF = frqz(BF);
                    Data(CondIDX).(name).Bandwidth = bandwidth;
                    Data(CondIDX).(name).Tuningwidth = tuningwidth;
                    Data(CondIDX).(name).init_Peak_BF = frqz(find(initPeakTune == max(initPeakTune(1:size(frqz,1)))));
                    Data(CondIDX).(name).init_Peak_BF_tune = initPeakTune;
                    Data(CondIDX).(name).SinkPeakAmp =PAMP;
                    Data(CondIDX).(name).SinkPeakLate =PLAT;
                    Data(CondIDX).(name).IntCurFlow = INT;
                    Data(CondIDX).(name).SinkDur = DUR2;
                    Data(CondIDX).(name).Sinkonset = onset;
                    Data(CondIDX).(name).Sinkoffset = offset;
                    Data(CondIDX).(name).SinkRMS = RMS;
                    Data(CondIDX).(name).tempSinkRMS = RMSTIME;
                    Data(CondIDX).(name).LFP = AvgFP;
                    Data(CondIDX).(name).CSD = AvgCSD;
                    Data(CondIDX).(name).BW_full = full_BW;
                    Data(CondIDX).(name).BW_early = e_BW;
                    Data(CondIDX).(name).BW_late = l_BW;
                    Data(CondIDX).(name).LayerRelRes =AvgLayerRelRes;
                    Data(CondIDX).(name).AVREC_raw = AvgRecCSD;
                    Data(CondIDX).(name).SingleTrial_AVREC_raw = SingleTrialAvgRecCSD;
                    Data(CondIDX).(name).SingleTrial_RELRES_raw = SingleTrialAvgRelResCSD;
                    Data(CondIDX).(name).RELRES_raw = AvgRelResCSD;
                    Data(CondIDX).(name).Full_RMS_AVREC = full_RMS_AvgRecCSD;
                    Data(CondIDX).(name).Early_RMS_AVREC = e_RMS_AvgRecCSD;
                    Data(CondIDX).(name).Late_RMS_AVREC = l_RMS_AvgRecCSD;
                    Data(CondIDX).(name).Early_RMS_RELRES = e_RMS_RelResCSD*100;
                    Data(CondIDX).(name).Late_RMS_RELRES = l_RMS_RelResCSD*100;
                    Data(CondIDX).(name).Full_RMS_ABSRES= full_RMS_AbsResCSD;
                    Data(CondIDX).(name).Early_RMS_ABSRES = e_RMS_AbsResCSD;
                    Data(CondIDX).(name).Late_RMS_ABSRES = l_RMS_AbsResCSD;                   

                    
                    figure('Name',[name ' ' measurement ': ' Condition{iC} ' ' num2str(i4)]);
                    plot(1:length(frqz),[Data(CondIDX).(name).SinkRMS.I_IIE],'LineWidth',2),... %'Color','black'
                        hold on,...
                        plot(1:length(frqz),[Data(CondIDX).(name).SinkRMS.I_IIL],'LineWidth',2),...
                        plot(1:length(frqz),[Data(CondIDX).(name).SinkRMS.IVE],'LineWidth',2),...
                        plot(1:length(frqz),[Data(CondIDX).(name).SinkRMS.IVL],'LineWidth',2),...
                        plot(1:length(frqz),[Data(CondIDX).(name).SinkRMS.VaE],'LineWidth',2),...
                        plot(1:length(frqz),[Data(CondIDX).(name).SinkRMS.VaL],'LineWidth',2),...
                        plot(1:length(frqz),[Data(CondIDX).(name).SinkRMS.VbE],'LineWidth',2),...
                        plot(1:length(frqz),[Data(CondIDX).(name).SinkRMS.VbL],'LineWidth',2),...
                        plot(1:length(frqz),[Data(CondIDX).(name).SinkRMS.VIaE],'LineWidth',2),...
                        plot(1:length(frqz),[Data(CondIDX).(name).SinkRMS.VIaL],'LineWidth',2),...
                        plot(1:length(frqz),[Data(CondIDX).(name).SinkRMS.VIbE],'LineWidth',2),...
                        plot(1:length(frqz),[Data(CondIDX).(name).SinkRMS.VIbL],'LineWidth',2),...
                        legend('I_IIE','I_IIL','IVE','IVL','VaE','VaL',...
                               'VbE','VbL','VIaE','VIaL','VIbE','VIbL');
                    hold off
                    h = gcf;
                    set(h, 'PaperType', 'A4');
                    set(h, 'PaperOrientation', 'landscape');
                    set(h, 'PaperUnits', 'centimeters');
                    savefig(h,[name '_' measurement '_RMS Sink tuning' ],'compact')
                    close (h)
                    
                    figure('Name',[name ' ' measurement ': ' Condition{iC} ' ' num2str(i4)]);
                    plot(1:length(frqz),[Data(CondIDX).(name).tempSinkRMS.I_IIE],'LineWidth',2),...%'Color','black'
                        hold on,...
                        plot(1:length(frqz),[Data(CondIDX).(name).tempSinkRMS.I_IIL],'LineWidth',2),...
                        plot(1:length(frqz),[Data(CondIDX).(name).tempSinkRMS.IVE],'LineWidth',2),...
                        plot(1:length(frqz),[Data(CondIDX).(name).tempSinkRMS.IVL],'LineWidth',2),...
                        plot(1:length(frqz),[Data(CondIDX).(name).tempSinkRMS.VaE],'LineWidth',2),...
                        plot(1:length(frqz),[Data(CondIDX).(name).tempSinkRMS.VaL],'LineWidth',2),...
                        plot(1:length(frqz),[Data(CondIDX).(name).tempSinkRMS.VbE],'LineWidth',2),...
                        plot(1:length(frqz),[Data(CondIDX).(name).tempSinkRMS.VbL],'LineWidth',2),...
                        plot(1:length(frqz),[Data(CondIDX).(name).tempSinkRMS.VIaE],'LineWidth',2),...
                        plot(1:length(frqz),[Data(CondIDX).(name).tempSinkRMS.VIaL],'LineWidth',2),...
                        plot(1:length(frqz),[Data(CondIDX).(name).tempSinkRMS.VIbE],'LineWidth',2),...
                        plot(1:length(frqz),[Data(CondIDX).(name).tempSinkRMS.VIbL],'LineWidth',2),...
                        legend('I_IIE','I_IIL','IVE','IVL','VaE','VaL',...
                               'VbE','VbL','VIaE','VIaL','VIbE','VIbL');
                    hold off
                    h = gcf;
                    set(h, 'PaperType', 'A4');
                    set(h, 'PaperOrientation', 'landscape');
                    set(h, 'PaperUnits', 'centimeters');
                    savefig(h,[name '_' measurement '_temporal RMS Sink tuning' ],'compact')
                    close (h)
                    toc
                end
            end
        end
    end
    

    cd ([home '/DATA'])
    save([input(i1).name(1:end-2) '_Data'],'Data');
    clear Data
end
