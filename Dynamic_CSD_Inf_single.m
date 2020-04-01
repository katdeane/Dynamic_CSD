function Dynamic_CSD_Inf_single(homedir)
%% Dynamic CSD for sinks I_II through Inf; incl. single

%   This script takes input from the groups and raw folders. It calculates 
%   and stores CSD information in Data struct which is
%   saved in the DATA folder.
% 
%   6a and 6b and Inf are detected here and there is single trial data
%   IMPORTANT: DO NOT change sink list here. If you need another set of
%   sinks then create a SEPARATE and UNIQUELY NAMED script.
%   Note: Currently the graphs being generated are limited to select sinks
%   although all sinks are being stored in the data structure
% 
%   CHANGE if needed: add your working directory to the try/catch; change
%   condition to run

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

cd (homedir),cd groups;

%% Load in
input = dir('*.m');
entries = length(input);


for i1 = 1:entries    
    
    run(input(i1).name);
    
    %% Choose Condition    
%     Condition = {'TM1' 'TM2' 'TM3'};
    Condition = {'condition'};
    %If Pre condition, must have capital P to match string in group code 
    
    %% Condition and Indexer   
    Data = struct;
    cd (homedir); cd figs;
    mkdir(['Single_Spike_' input(i1).name(1:end-2)]);
    
    Indexer = imakeIndexer(Condition,animals,Cond);
    %%
    
    for iA = 1:length(animals)
        cd (homedir); cd raw;
        name = animals{iA};
        
        for iC = 1:length(Condition)
            for i4 = 1:length(Cond.(Condition{iC}){iA})
                if i4 == 1
                    CondIDX = Indexer(2).(Condition{iC});
                else
                    CondIDX = Indexer(2).(Condition{iC})+i4-1;
                end
                tic
                
                measurement = Cond.(Condition{iC}){iA}{i4};
                if ~isempty(measurement)
                    disp(['Analyzing animal: ' name '_' measurement])
                    clear SWEEP
                    try
                        load ([name '_' measurement]);clear avgFP;
                    catch
                        fprintf('the name or measurement does not exist/n')
                    end
                    
                    cd (homedir),cd groups;
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
                        frqz(find(frqz == inf))=[]; % takes click out of analysis
                        frqz(find(frqz == 0))=[]; % takes pause out of analysis
                    else %Chronic recordings
                        cd (homedir); cd subfunc;
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
                    cd (homedir),cd subfunc;
                    
                    %note: channel order is given twice because the whole
                    %CSD needs to be checked now (and not any one layer)
                    [SingleLayer_AVREC,AvgLayer_AVREC,AvgFP, SingleTrialFP, AvgCSD,...
                        SingleTrialCSD, AvgRecCSD, SingleTrialAvgRecCSD,...
                        SingleTrialRelResCSD, AvgRelResCSD,AvgAbsResCSD,...
                        SingleTrialAbsResCSD, SingleLayerRelRes, AvgLayerRelRes] =...
                        SingleTrialCSD_full(SWEEP, str2num(channels{iA}),1:length(str2num(channels{iA})),BL);
                                       
                    %delete empty columns to have the correct amount of stimuli
                    %present (needed for attenuation 30)
                    SingleLayer_AVREC = SingleLayer_AVREC(~cellfun('isempty',SingleLayer_AVREC'));
                    SingleTrialFP = SingleTrialFP(~cellfun('isempty',SingleTrialFP'));
                    SingleTrialCSD = SingleTrialCSD(~cellfun('isempty',SingleTrialCSD'));
                    SingleTrialRelResCSD = SingleTrialRelResCSD(~cellfun('isempty',SingleTrialRelResCSD'));
                    AvgAbsResCSD = AvgAbsResCSD(~cellfun('isempty',AvgAbsResCSD'));
                    SingleTrialAbsResCSD = SingleTrialAbsResCSD(~cellfun('isempty',SingleTrialAbsResCSD'));
                    SingleLayerRelRes = SingleLayerRelRes(~cellfun('isempty',SingleLayerRelRes'));
                    AvgLayerRelRes = AvgLayerRelRes(~cellfun('isempty',AvgLayerRelRes'));
                    SingleTrialAvgRecCSD = SingleTrialAvgRecCSD(~cellfun('isempty',SingleTrialAvgRecCSD'));
                    AvgRelResCSD = AvgRelResCSD(~cellfun('isempty', AvgRelResCSD'));
                    AvgCSD = AvgCSD(~cellfun('isempty', AvgCSD')); AvgCSD=AvgCSD(1:length(frqz));
                    AvgFP = AvgFP(~cellfun('isempty', AvgFP'));
                    AvgRecCSD = AvgRecCSD(~cellfun('isempty', AvgRecCSD'));
                    
                    
                    % Sink durations
                    L.I_IIE = str2num(Layer.I_IIE{iA}); L.I_IIL = str2num(Layer.I_IIL{iA}); 
                    L.IVE = str2num(Layer.IVE{iA}); L.IVL = str2num(Layer.IVL{iA});
                    L.VaE = str2num(Layer.VaE{iA}); L.VaL = str2num(Layer.VaL{iA}); 
                    L.VbE = str2num(Layer.VbE{iA}); L.VbL = str2num(Layer.VbL{iA});
                    L.VIaE = str2num(Layer.VIaE{iA}); L.VIaL = str2num(Layer.VIaL{iA});
                    L.VIbE = str2num(Layer.VIbE{iA}); L.VIbL = str2num(Layer.VIbL{iA});
                    L.InfE = str2num(Layer.InfE{iA}); L.InfL = str2num(Layer.InfL{iA});
                    
                    curChan = str2num(channels{iA});
                    
                    %Generate Sink Boxes
                    [DUR,RMS,SINGLE_RMS,SINT,PAMP,SINGLE_PAMP,PLAT,SINGLE_PLAT,INT] =...
                        sink_dura_single(L,AvgCSD,SingleTrialCSD,BL,Baseline);
                    
                    toc
                    
                    % Calculate CSD %                    
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
                    DUR2.InfE = reshape(round([DUR(:).InfE]),[],length([DUR(:).InfE])/2)'-BL;
                    DUR2.InfL = reshape(round([DUR(:).InfL]),[],length([DUR(:).InfL])/2)'-BL;
                    
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
                    onset.InfE = DUR2.InfE(:,1)';     offset.InfE = DUR2.InfE(:,2)';
                    onset.InfL = DUR2.InfL(:,1)';     offset.InfL = DUR2.InfL(:,2)';
                    
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
                    DUR2.InfE = DUR2.InfE(:,2)' - DUR2.InfE(:,1)'; 
                    DUR2.InfL = DUR2.InfL(:,2)' - DUR2.InfL(:,1)'; 
                                       
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
                    
                    %% Plots and Figures
                    figure('Name',[name ' ' measurement ': ' Condition{iC}])
                    tic
                    disp('Plotting CSD with sink detections')
                    for i5 = 1:length(AvgCSD)
                        subplot(2,round(length(AvgCSD)/2),i5)
                        imagesc(AvgCSD{i5})
                        if i5 == 1
                            title ([name ' ' measurement ': ' Condition{iC} ' ' num2str(i5) ' ' num2str(frqz(i5)) ' Hz'])
                        else
                            title ([num2str(frqz(i5)) ' Hz'])
                        end
                        
                        colormap (jet)
                        
                        
                        caxis([-0.0005 0.0005])
                        
                        
                        
                        hold on
                        % Layer I_II
                        y =[(max(L.I_IIE)+0.5),(max(L.I_IIE)+0.5),(min(L.I_IIE)-0.5),(min(L.I_IIE)-0.5),(max(L.I_IIE)+0.5)];
                        if isempty(y); y = [NaN NaN NaN NaN NaN]; end %in case the upper layer is not there
                        x = [DUR(i5).I_IIE(1), DUR(i5).I_IIE(2),DUR(i5).I_IIE(2),DUR(i5).I_IIE(1),DUR(i5).I_IIE(1)];
                        plot(x,y,'black','LineWidth',2)
                        
                        y =[(max(L.I_IIL)+0.5),(max(L.I_IIL)+0.5),(min(L.I_IIL)-0.5),(min(L.I_IIL)-0.5),(max(L.I_IIL)+0.5)];
                        if isempty(y); y = [NaN NaN NaN NaN NaN]; end
                        x = [DUR(i5).I_IIL(1), DUR(i5).I_IIL(2),DUR(i5).I_IIL(2),DUR(i5).I_IIL(1),DUR(i5).I_IIL(1)];
                        plot(x,y,'white','LineWidth',2)
                        
                        % Layer IV
                        y =[(max(L.IVE)+0.5),(max(L.IVE)+0.5),(min(L.IVE)-0.5),(min(L.IVE)-0.5),(max(L.IVE)+0.5)];
                        x = [DUR(i5).IVE(1), DUR(i5).IVE(2),DUR(i5).IVE(2),DUR(i5).IVE(1),DUR(i5).IVE(1)];
                        plot(x,y,'black','LineWidth',2)
                        
                        y =[(max(L.IVL)+0.5),(max(L.IVL)+0.5),(min(L.IVL)-0.5),(min(L.IVL)-0.5),(max(L.IVL)+0.5)];
                        x = [DUR(i5).IVL(1), DUR(i5).IVL(2),DUR(i5).IVL(2),DUR(i5).IVL(1),DUR(i5).IVL(1)];
                        plot(x,y,'white','LineWidth',2)
                        
                        % Layer Va
                        y =[(max(L.VaE)+0.5),(max(L.VaE)+0.5),(min(L.VaE)-0.5),(min(L.VaE)-0.5),(max(L.VaE)+0.5)];
                        x = [DUR(i5).VaE(1), DUR(i5).VaE(2),DUR(i5).VaE(2),DUR(i5).VaE(1),DUR(i5).VaE(1)];
                        plot(x,y,'black','LineWidth',2)
                        
                        y =[(max(L.VaL)+0.5),(max(L.VaL)+0.5),(min(L.VaL)-0.5),(min(L.VaL)-0.5),(max(L.VaL)+0.5)];
                        x = [DUR(i5).VaL(1), DUR(i5).VaL(2),DUR(i5).VaL(2),DUR(i5).VaL(1),DUR(i5).VaL(1)];
                        plot(x,y,'white','LineWidth',2)
                        
                        % Layer Vb
                        y =[(max(L.VbE)+0.5),(max(L.VbE)+0.5),(min(L.VbE)-0.5),(min(L.VbE)-0.5),(max(L.VbE)+0.5)];
                        x = [DUR(i5).VbE(1), DUR(i5).VbE(2),DUR(i5).VbE(2),DUR(i5).VbE(1),DUR(i5).VbE(1)];
                        plot(x,y,'black','LineWidth',2)
                        
                        y =[(max(L.VbL)+0.5),(max(L.VbL)+0.5),(min(L.VbL)-0.5),(min(L.VbL)-0.5),(max(L.VbL)+0.5)];
                        x = [DUR(i5).VbL(1), DUR(i5).VbL(2),DUR(i5).VbL(2),DUR(i5).VbL(1),DUR(i5).VbL(1)];
                        plot(x,y,'white','LineWidth',2)
                        
                        % Layer VIa
                        y =[(max(L.VIaE)+0.5),(max(L.VIaE)+0.5),(min(L.VIaE)-0.5),(min(L.VIaE)-0.5),(max(L.VIaE)+0.5)];
                        x = [DUR(i5).VIaE(1), DUR(i5).VIaE(2),DUR(i5).VIaE(2),DUR(i5).VIaE(1),DUR(i5).VIaE(1)];
                        plot(x,y,'black','LineWidth',2)
                        
                        y =[(max(L.VIaL)+0.5),(max(L.VIaL)+0.5),(min(L.VIaL)-0.5),(min(L.VIaL)-0.5),(max(L.VIaL)+0.5)];
                        x = [DUR(i5).VIaL(1), DUR(i5).VIaL(2),DUR(i5).VIaL(2),DUR(i5).VIaL(1),DUR(i5).VIaL(1)];
                        plot(x,y,'white','LineWidth',2)
                        
                        % Layer VIb
                        y =[(max(L.VIbE)+0.5),(max(L.VIbE)+0.5),(min(L.VIbE)-0.5),(min(L.VIbE)-0.5),(max(L.VIbE)+0.5)];
                        x = [DUR(i5).VIbE(1), DUR(i5).VIbE(2),DUR(i5).VIbE(2),DUR(i5).VIbE(1),DUR(i5).VIbE(1)];
                        plot(x,y,'black','LineWidth',2)
                        
                        y =[(max(L.VIbL)+0.5),(max(L.VIbL)+0.5),(min(L.VIbL)-0.5),(min(L.VIbL)-0.5),(max(L.VIbL)+0.5)];
                        x = [DUR(i5).VIbL(1), DUR(i5).VIbL(2),DUR(i5).VIbL(2),DUR(i5).VIbL(1),DUR(i5).VIbL(1)];
                        plot(x,y,'white','LineWidth',2)
                        
                        % Layer Infragranular
                        y =[(max(L.InfE)+0.5),(max(L.InfE)+0.5),(min(L.InfE)-0.5),(min(L.InfE)-0.5),(max(L.InfE)+0.5)];
                        x = [DUR(i5).InfE(1), DUR(i5).InfE(2),DUR(i5).InfE(2),DUR(i5).InfE(1),DUR(i5).InfE(1)];
                        plot(x,y,'black','LineWidth',2)
                        
                        y =[(max(L.InfL)+0.5),(max(L.InfL)+0.5),(min(L.InfL)-0.5),(min(L.InfL)-0.5),(max(L.InfL)+0.5)];
                        x = [DUR(i5).InfL(1), DUR(i5).InfL(2),DUR(i5).InfL(2),DUR(i5).InfL(1),DUR(i5).InfL(1)];
                        plot(x,y,'white','LineWidth',2)
                        
                        hold off
                        
                    end
                    toc
                    
                    cd([homedir '\figs\']); mkdir(['Single_' input(i1).name(1:end-2)]);
                    cd(['Single_' input(i1).name(1:end-2)])
                    h = gcf;
                    set(h, 'PaperType', 'A4');
                    set(h, 'PaperOrientation', 'landscape');
                    set(h, 'PaperUnits', 'centimeters');
                    savefig(h,[name '_' measurement '_CSDs' ],'compact')
                    close (h)

                    Order = {'I_IIL','IVE','VaE','VbE','VIaE','InfE','VIaL'};
                    h = figure('Name',[name ' ' measurement ': ' Condition{iC} ' ' num2str(i4)]);
                    for i5 = 1:length(Order)
                        subplot(1,7,i5)
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
                    savefig(h,[name '_' measurement '_Sink_onset+offset' ],'compact')
                    close (h)

                    BF = find([PAMP.IVE] == max([PAMP.IVE]));
                    BF_frqz = frqz(BF);
                    
                    if isnan(BF), BF = 1; end
                    if isempty(BF), BF = 1; end
                    
                    SINK = {'VbE', 'IVE', 'VIaE', 'VIbE', 'VaE', 'I_IIE','InfE',...
                        'VbL', 'IVL', 'VIaL', 'VIbL', 'VaL', 'I_IIL','InfL'};
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
                        plot([RMS.(SINK{i5})],'b--o','LineWidth',2)
                        hold on
                        plot (BF,RMS(BF).(SINK{i5}),'r*','LineWidth',2)
                        hold off
                        xlim([0 (length(AvgCSD)+1)])
                        title (['RMS Layer ' SINK{i5} ' mV/mm?'])
                        
                        subplot(3,length(SINK),i5+2*length(SINK))
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
                    savefig(h,[name '_' measurement '_Duration_RMS_Ratio' ],'compact')
                    close (h)
                    
                    % Integral current flow V*ms/mm?
                    h = figure('Name',[name ' ' measurement ': ' Condition{iC} ' ' num2str(i4)]);
                    for i5 = 1: length(SINK)
                        subplot(3,length(SINK),i5)
                        durtime = DUR2.(SINK{i5});
                        plot(durtime./durtime(BF),'LineWidth',2)
                        hold on
                        plot (BF,durtime(BF)./durtime(BF),'r*','LineWidth',2)
                        hold off
                        xlim([0 (length(AvgCSD)+1)])
                        title (['Sink duration Layer ' SINK{i5} ' ms'])
                        
                        subplot(3,length(SINK),i5+length(SINK))
                        RS =[RMS.(SINK{i5})];
                        plot(RS./RS(BF),'LineWidth',2)
                        hold on
                        plot (BF,RMS(BF).(SINK{i5})/RS(BF),'r*','LineWidth',2)
                        hold off
                        xlim([0 (length(AvgCSD)+1)])
                        title (['RMS Layer ' SINK{i5} ' mV/mm?'])
                        
                        subplot(3,length(SINK),i5+2*length(SINK))
                        relRMS = [RMS.(SINK{i5})]./durtime;
                        plot(relRMS./relRMS(BF),'LineWidth',2)
                        hold on
                        plot (BF,relRMS(BF)/relRMS(BF),'r*','LineWidth',2)
                        hold off
                        xlim([0 (length(AvgCSD)+1)])
                        title (['RMS / ms L ' SINK{i5} ' mV/mm? ms'])
                        
                    end
                    
                    set(h, 'PaperType', 'A4');
                    set(h, 'PaperOrientation', 'landscape');
                    set(h, 'PaperUnits', 'centimeters');
                    savefig(h,[name '_' measurement 'Duration_RMS_Ratio_BF_normalized' ],'compact')
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
                    Data(CondIDX).(name).SingleSinkPeakAmp =SINGLE_PAMP;
                    Data(CondIDX).(name).SinkPeakLate =PLAT;
                    Data(CondIDX).(name).SingleSinkPeakLat =SINGLE_PLAT;
                    Data(CondIDX).(name).IntCurFlow = INT;
                    Data(CondIDX).(name).SinkDur = DUR2;
                    Data(CondIDX).(name).Sinkonset = onset;
                    Data(CondIDX).(name).Sinkoffset = offset;
                    Data(CondIDX).(name).SinkRMS = RMS;
                    Data(CondIDX).(name).SinkINT = SINT;
                    Data(CondIDX).(name).SingleSinkRMS = SINGLE_RMS;
                    Data(CondIDX).(name).tempSinkRMS = RMSTIME;
                    Data(CondIDX).(name).singletrialLFP = SingleTrialFP;
                    Data(CondIDX).(name).LFP = AvgFP;
                    Data(CondIDX).(name).singletrialCSD =SingleTrialCSD;
                    Data(CondIDX).(name).CSD = AvgCSD;
                    Data(CondIDX).(name).LayerRelRes =AvgLayerRelRes;
                    Data(CondIDX).(name).AVREC_raw = AvgRecCSD;
                    Data(CondIDX).(name).SingleTrial_AVREC_raw = SingleTrialAvgRecCSD;
                    Data(CondIDX).(name).SingleTrial_RelRes_raw = SingleTrialRelResCSD;
                    Data(CondIDX).(name).SingleTrial_AbsRes_raw = SingleTrialAbsResCSD;
                    Data(CondIDX).(name).RELRES_raw = AvgRelResCSD;
                    Data(CondIDX).(name).ABSRES_raw = AvgAbsResCSD;
                   
                    figure('Name',[name ' ' measurement ': ' Condition{iC} ' ' num2str(i4)]);
                    plot(1:length(frqz),[Data(CondIDX).(name).SinkRMS.I_IIL],'LineWidth',2,'Color','black'),...
                        hold on,...
                        plot(1:length(frqz),[Data(CondIDX).(name).SinkRMS.IVE],'LineWidth',2,'Color','red'),...
                        plot(1:length(frqz),[Data(CondIDX).(name).SinkRMS.VaE],'LineWidth',2,'Color','green'),...
                        plot(1:length(frqz),[Data(CondIDX).(name).SinkRMS.VbE],'LineWidth',2,'Color','blue'),...
                        plot(1:length(frqz),[Data(CondIDX).(name).SinkRMS.VIaE],'LineWidth',2,'Color','magenta'),...
                        plot(1:length(frqz),[Data(CondIDX).(name).SinkRMS.InfE],'LineWidth',2,'Color','yellow'),...
                        plot(1:length(frqz),[Data(CondIDX).(name).SinkRMS.VIaL],'LineWidth',2,'Color','cyan'),...
                        legend('I/IIL', 'IVE', 'VaE', 'VbE', 'VIaE','IG', 'VIaL')
                    hold off
                    h = gcf;
                    set(gca,'XTickLabel',Data(CondIDX).(name).Frqz,'FontSize',12);
                    set(h, 'PaperType', 'A4');
                    set(h, 'PaperOrientation', 'landscape');
                    set(h, 'PaperUnits', 'centimeters');
                    savefig(h,[name '_' measurement '_RMS Sink tuning' ],'compact')
                    close (h)
                    
                    figure('Name',[name ' ' measurement ': ' Condition{iC} ' ' num2str(i4)]);
                    plot(1:length(frqz),[Data(CondIDX).(name).tempSinkRMS.I_IIL],'LineWidth',2,'Color','black'),...
                        hold on,...
                        plot(1:length(frqz),[Data(CondIDX).(name).tempSinkRMS.IVE],'LineWidth',2,'Color','red'),...
                        plot(1:length(frqz),[Data(CondIDX).(name).tempSinkRMS.VaE],'LineWidth',2,'Color','green'),...
                        plot(1:length(frqz),[Data(CondIDX).(name).tempSinkRMS.VbE],'LineWidth',2,'Color','blue'),...
                        plot(1:length(frqz),[Data(CondIDX).(name).tempSinkRMS.VIaE],'LineWidth',2,'Color','magenta'),...
                        plot(1:length(frqz),[Data(CondIDX).(name).tempSinkRMS.InfE],'LineWidth',2,'Color','yellow'),...
                        plot(1:length(frqz),[Data(CondIDX).(name).tempSinkRMS.VIaL],'LineWidth',2,'Color','cyan'),...
                        legend('I/IIL', 'IVE', 'VaE', 'VbE', 'VIaE', 'GS','VIaL')
                    hold off
                    h = gcf;
                    set(gca,'XTickLabel',Data(CondIDX).(name).Frqz,'FontSize',12);
                    set(h, 'PaperType', 'A4');
                    set(h, 'PaperOrientation', 'landscape');
                    set(h, 'PaperUnits', 'centimeters');
                    savefig(h,[name '_' measurement '_temporal RMS Sink tuning' ],'compact')
                    close (h)
                end
            end
        end
    end
    cd ([homedir '/DATA'])
    save([input(i1).name(1:end-2) '_Data'],'Data');
    clear Data
end
cd(homedir)
toc