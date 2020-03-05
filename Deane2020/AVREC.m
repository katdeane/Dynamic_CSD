%% README

%  This script takes input from the groups and raw folders. It generates
%  CSD profiles and then averages the AVREC across all measurements within
%  the group per animal. Output is saved in pics and DATA folders. This is
%  very basic and can be adapted.

%  This is a modification of the AVRECFM code

clear
%% Inputs
SO = 1; %sinks only: positive values are set as nan during calculations

%% standard operations
warning('OFF');
dbstop if error %stops execution of program if and allows user to examine the workspace
cd('D:\MyCode\Dynamic_CSD_Analysis');
home = pwd; %pwd: print working directory, gives full current directory title of "home"
addpath(genpath(home)); %adds all subdirectories of 'home' to path for access

%% code
cd (home),cd groups; %cd: change directory

input = dir('*.m'); %dir *.m lists all program files in the current directory (stores 5 fields of data but struc array contains #x1)
entries = length(input); %number of input entries

for i1 = 1:entries
    run(input(i1).name);
    
    cd (home); cd figs; %opens figure folder to store future images
    mkdir(['Single_' input(i1).name(1:end-2)]); %adds folder to directory - name without '.m' (-2 from the end)
    
    avgrec = struct; %creating empty structure to fill
    Indexer = struct; %creating empty structure to fill
%     Condition = {'pre' 'post' 'att' 'FM' 'between' 'muscimol' 'musatt' 'musFM'};
%     Condition = {'naiveawake' 'naivepre' 'muscimol' };
    Condition = {'condition'};
%     Condition = {'Pre' 'post'};
    
    for i2 = 1:length(Condition) %input from script (length of condition: 'Pre' 'Combo' 'Post1'=3) %i2 is condition we are in
        
        for i3 = 1:length(animals) %input from groups (in workspace) %i3 is animal we are on
            counter = 0;
            counter = counter + length(Cond.(Condition{i2}){i3}); %Cond has all conditions recorded, condition(current cond type selection for current animal selection)
            % counter = the number of entries in current condition for the current animal
            if i3 == 1 %for the first animal
                Indexer(1).(Condition{i2}) = counter; %the indexer structure is the size of the current counter
            elseif Indexer.(Condition{i2}) < counter %if the next animal has a structure larger than the current counter, it will make the indexer that size
                Indexer(1).(Condition{i2}) = counter;
            end
            %An indexer is now made for the i2th condition that is the size of the largest amongst the animals
        end
        % Indexer now contains all conditions at their greatest sizes
    end
    
    count = 1; %new variable 'count'
    for i2 = 1:length(Condition) %i2 is still which condition we are in and it's reset to 1
        
        if i2 == 1 %first condition type
            Indexer(2).(Condition{i2}) =1; %adding a second field to the condition i2
        else
            Indexer(2).(Condition{i2}) = count;
        end
        count = count + Indexer(1).(Condition{i2}); %current count (e.g. 1)+size of condition i2 in index(e.g. 3) = (e.g. 4)
    end
    %now the indexer has 2 columns beyond the cond names, indicating the condition size in the first and the starting count in the second
    
    for i2 = 1:length(animals) %i2 is now for animals!
        cd (home); cd raw;    %we are switching folders to raw to load the .mat's
        
        for i3 = 1:length(Condition) %i3 is now condition!
            clear name; clear SWEEP
            name = animals{i2}; %assigns the i2th entry in the list in the current group to variable 'name'
            for i4 = 1:length(Cond.(Condition{i3}){i2}) %i4 goes through the number of entries in condition i3 for animal i4
                %now we are in group i1, animal i2, condition i3, entry i4
                if i4 == 1 %the first entry
                    CondIDX = Indexer(2).(Condition{i3}); %new count variable calling the start count number for the first condition entry (which equals 1)
                else %entries following
                    CondIDX = Indexer(2).(Condition{i3})+i4-1;
                end
                
                measurement = Cond.(Condition{i3}){i2}{i4}; %new variable= the i3 condition from the Cond struct and further the i2 animal and its i4 entry
                %e.g. 003 is the first entry of the first animal in the first condition
                
                if ~isempty(measurement) % ~ is the logical argument 'not' so this says if the measurement variable is not empty than true
                    disp(['Analyzing animal: ' name '_' measurement])
                    tic
                    
                    try %allows the code to continue running if statements in try/catch do not succeed
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
                    if exist('SWEEP','var') %Micha's recordings
                        clear DATA frqz
                        BL = Header.t_pre*P.Fs_AD(1); %BL-baseline %t_pre is the time before the tone %Fs_AD - Sampling frequency of channels (they are all the same so we use first value)
                        tone = Header.t_sig(1)*P.Fs_AD(1); %t_sig is duration of stimulus * sampling rate = 200
                        frqz = Header.stimlist(:,1); %stimlist contains tone frequencies in all rows (:), first column (:,1)
                        Fs = P.Fs_AD(1); %sampling rate
                    else %Marina's recordings
                        cd (home); cd subfunc;
                        %To apply channel interpolation
                        if ~isempty(DATA.ART)
                            [DATA.LFP]=chaninterp(DATA.LFP, 'linextra', DATA.ART.chan,...
                                DATA.ELECPOS);
                        end
                        %To select attenuation specific recordings
                        if ~isempty(dB_lev{i2})
                            if strcmp(dB_lev{i2},'10dB')
                                [DATA]=attenuation10(DATA);
                            elseif strcmp(dB_lev{i2},'30dB')
                                [DATA]=attenuation30(DATA);
                            end
                        end
                        %To create SWEEP
                        [BL, tone, Fs, SWEEP, frqz]=reformat(DATA);
                    end
                    
                    cd (home),cd subfunc;
                    %let's do this
                    [~, AvgCSD, AvgRecCSD, ~, ~, ~, ~, ~] =...
                        SingleTrialCSD_reduced(SWEEP, str2num(channels{i2}),BL);
                    
                    %delete empty columns to have the correct amount of stimuli
                    %present (needed for attenuation 30)
                    AvgCSD = AvgCSD(~cellfun('isempty', AvgCSD'));
                    AvgRecCSD = AvgRecCSD(~cellfun('isempty', AvgRecCSD'));
                    
                    % Sink durations
                    L.I_IIE = str2num(Layer.I_IIE{i2}); L.I_IIL = str2num(Layer.I_IIL{i2}); 
                    L.IVE = str2num(Layer.IVE{i2}); L.IVL = str2num(Layer.IVL{i2});
                    L.VaE = str2num(Layer.VaE{i2}); L.VaL = str2num(Layer.VaL{i2}); 
                    L.VbE = str2num(Layer.VbE{i2}); L.VbL = str2num(Layer.VbL{i2});
                    L.VIaE = str2num(Layer.VIaE{i2}); L.VIaL = str2num(Layer.VIaL{i2});
                    L.VIbE = str2num(Layer.VIbE{i2}); L.VIbL = str2num(Layer.VIbL{i2});
                    
                    curChan = str2num(channels{i2});
                    
                    %Generate Sink Boxes
                    [DUR,RMS,PAMP,PLAT] = sink_dura_group(L,AvgCSD,BL,tone,SWEEP,curChan,Baseline);
                    
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
                    state = 1:(length(AvgRecCSD{1})-BL); 
                    
                    %calculate full RMS of AVREC  (1-max Duration)
                    [~, ~, full_RMS_AvgRecCSD,~, ~,~, ~] =...
                        ExtractCSDBasedPar(1,AvgRecCSD,BL,...
                        Fs, 1, 0, threshold_std, threshold_dur,...
                        Latency_HighCutoff, Latency_LowCutoff,state);
                    
                    
                    %calculate early RMS of AVREC (1-50 ms)
                    state = 1:50;
                    [~, ~, e_RMS_AvgRecCSD,~, ~, ~, ~] =...
                        ExtractCSDBasedPar(1,AvgRecCSD,BL,...
                        Fs, 1, 0, threshold_std, threshold_dur,...
                        Latency_HighCutoff, Latency_LowCutoff,state);
                    
                    
                    %calculate late RMS of AVREC  (81-300 ms)
                    state = 51:200;
                    [~, ~, l_RMS_AvgRecCSD,~, ~, ~, ~] =...
                        ExtractCSDBasedPar(1,AvgRecCSD,BL,...
                        Fs, 1, 0, threshold_std, threshold_dur,...
                        Latency_HighCutoff, Latency_LowCutoff,state);
                    
                    toc

                    %rename data for plots
                    if ~isempty(dB_lev{i2})
                        name = [animals{i2} '_' dB_lev{i2}];
                    end
                   
                    %% Plots and Figures %%
                    figure('Name',[name ' ' measurement ': ' Condition{i3}])
                    tic
                    disp('Plotting CSD with sink detections')
                    for i5 = 1:length(AvgCSD)
                        subplot(2,round(length(AvgCSD)/2),i5)
                        imagesc(AvgCSD{i5})
                        if i5 == 1 %full title for first image
                            title ([name ' ' measurement ': ' Condition{i3} ' ' num2str(i5) ' ' num2str(frqz(i5)) ' Hz'])
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
                    
                    % the avgrec mat will be saved. First through seventh
                    % go in order of lowest to highest frequency
                    avgrec.(name).first = AvgRecCSD{1,1}';
                    avgrec.(name).second = AvgRecCSD{1,2}';
                    avgrec.(name).third = AvgRecCSD{1,3}';
                    avgrec.(name).fourth = AvgRecCSD{1,4}';
                    avgrec.(name).fifth = AvgRecCSD{1,5}';
                    avgrec.(name).sixth = AvgRecCSD{1,6}';
                    avgrec.(name).seventh = AvgRecCSD{1,7}';
%                     avgrec.(name).pause = AvgRecCSD{1,8}';
%                     avgrec.(name).click = AvgRecCSD{1,9}';
                    
                end
            end
        end
    end
    
    %% plotting avgrec
    FNames = fields(avgrec);
    Legend = vertcat(FNames, 'Average');
    
    %up to BF
    h = figure('Name','AvgRec - Frqz 1');
    average = [];
    for i1 = 1:length(animals)
        line = avgrec.(FNames{i1}).first;
        line = line(1:600); %GXL have 600ms and GKD have 800!
        plot(line);
        average = [average line'];
        hold on
    end
    
    average = mean(average,2)';
    plot(average,'LineWidth',2);
    legend(Legend)
    
    set(h, 'PaperType', 'A4');
    set(h, 'PaperOrientation', 'landscape');
    set(h, 'PaperUnits', 'centimeters');
    savefig(h, 'AvgRec - Frqz 1','compact')
    close (h)
    
    %up full sweep
    h = figure('Name','AvgRec - Frqz 2');
    average = [];
    for i1 = 1:length(animals)
        line = avgrec.(FNames{i1}).second;
        line = line(1:600); %GXL have 600ms and GKD have 800!
        plot(line);
        average = [average line'];
        hold on
    end
    
    average = mean(average,2)';
    plot(average,'LineWidth',2);
    legend(Legend)
    
    set(h, 'PaperType', 'A4');
    set(h, 'PaperOrientation', 'landscape');
    set(h, 'PaperUnits', 'centimeters');
    savefig(h, 'AvgRec - Frqz 2','compact')
    close (h)
    
    %BF up
    h = figure('Name','AvgRec - Frqz 3');
    average = [];
    for i1 = 1:length(animals)
        line = avgrec.(FNames{i1}).third;
        line = line(1:600); %GXL have 600ms and GKD have 800!
        plot(line);
        average = [average line'];
        hold on
    end
    
    average = mean(average,2)';
    plot(average,'LineWidth',2);
    legend(Legend)
    
    set(h, 'PaperType', 'A4');
    set(h, 'PaperOrientation', 'landscape');
    set(h, 'PaperUnits', 'centimeters');
    savefig(h, 'AvgRec - Frqz 3','compact')
    close (h)
    
    %down to BF
    h = figure('Name','AvgRec - Frqz 4');
    average = [];
    for i1 = 1:length(animals)
        line = avgrec.(FNames{i1}).fourth;
        line = line(1:600); %GXL have 600ms and GKD have 800!
        plot(line);
        average = [average line'];
        hold on
    end
    
    average = mean(average,2)';
    plot(average,'LineWidth',2);
    legend(Legend)
    
    set(h, 'PaperType', 'A4');
    set(h, 'PaperOrientation', 'landscape');
    set(h, 'PaperUnits', 'centimeters');
    savefig(h, 'AvgRec - Frqz 4','compact')
    close (h)
    
    %Down full sweep
    h = figure('Name','AvgRec - Frqz 5');
    average = [];
    for i1 = 1:length(animals)
        line = avgrec.(FNames{i1}).fifth;
        line = line(1:600); %GXL have 600ms and GKD have 800!
        plot(line);
        average = [average line'];
        hold on
    end
    
    average = mean(average,2)';
    plot(average,'LineWidth',2);
    legend(Legend)
    
    set(h, 'PaperType', 'A4');
    set(h, 'PaperOrientation', 'landscape');
    set(h, 'PaperUnits', 'centimeters');
    savefig(h, 'AvgRec - Frqz 5','compact')
    close (h)
    
    %BF down
    h = figure('Name','AvgRec - Frqz 6');
    average = [];
    for i1 = 1:length(animals)
        line = avgrec.(FNames{i1}).sixth;
        line = line(1:600); %GXL have 600ms and GKD have 800!
        plot(line);
        average = [average line'];
        hold on
    end
    
    average = mean(average,2)';
    plot(average,'LineWidth',2);
    legend(Legend)
    
    set(h, 'PaperType', 'A4');
    set(h, 'PaperOrientation', 'landscape');
    set(h, 'PaperUnits', 'centimeters');
    savefig(h, 'AvgRec - Frqz 6','compact')
    close (h)
    
    %BF down
    h = figure('Name','AvgRec - Frqz 7');
    average = [];
    for i1 = 1:length(animals)
        line = avgrec.(FNames{i1}).seventh;
        line = line(1:600); %GXL have 600ms and GKD have 800!
        plot(line);
        average = [average line'];
        hold on
    end
    
    average = mean(average,2)';
    plot(average,'LineWidth',2);
    legend(Legend)
    
    set(h, 'PaperType', 'A4');
    set(h, 'PaperOrientation', 'landscape');
    set(h, 'PaperUnits', 'centimeters');
    savefig(h, 'AvgRec - Frqz 7','compact')
    close (h)
    
%     %BF down
%     h = figure('Name','AvgRec - Pause');
%     average = [];
%     for i1 = 1:length(animals);
%         line = avgrec.(FNames{i1}).pause;
%         line = line(1:600); %GXL have 600ms and GKD have 800!
%         plot(line);
%         average = [average line'];
%         hold on
%     end
%     
%     average = mean(average,2)';
%     plot(average,'LineWidth',2);
%     legend(Legend)
%     
%     set(h, 'PaperType', 'A4');
%     set(h, 'PaperOrientation', 'landscape');
%     set(h, 'PaperUnits', 'centimeters');
%     savefig(h, 'AvgRec - Pause','compact')
%     close (h)
    
%     %BF down
%     h = figure('Name','AvgRec - Click');
%     average = [];
%     for i1 = 1:length(animals);
%         line = avgrec.(FNames{i1}).click;
%         plot(line);
%         average = [average line'];
%         hold on
%     end
%     
%     average = mean(average,2)';
%     plot(average,'LineWidth',2);
%     legend(Legend)
%     
%     set(h, 'PaperType', 'A4');
%     set(h, 'PaperOrientation', 'landscape');
%     set(h, 'PaperUnits', 'centimeters');
%     savefig(h, 'AvgRec - Click','compact')
%     close (h)
    
    
    cd ([home '/DATA'])
    save([input(1).name(1:end-2) '_avgrec'],'avgrec');
    clear avgrec
end
