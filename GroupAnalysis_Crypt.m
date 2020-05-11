function GroupAnalysis_Crypt(homedir)
%Input:     D:\MyCode\Dynamic_CSD_Analysis\DATA -> *DATA.mat; (bin,mirror)
%Output:    Figures of groups in "Group..." folder 
%           .mat files in DATA/output folder
%           AVG data

% Data out:     Tuning struct contains sorted tuning of all tonotopies per
%               per layer per parameter (for FIRST sink in layer if it
%               falls between 0:65 ms, pause&click not included); Click
%               struct containes sorted click response per stimulus
%               frequency per layer per parameter; Normalized Clicks
%               normalize all sinks within one animal/one layer/one click 
%               frequency to the first ~isnan in the pre-laser condition -
%               if no detected sinks in pre-laser condition, animal layer
%               stim type (2hz, 5hz, etc.) is nanned. 

% Figures out:  SinkRMS of tonotopy tuning curves through the recording
%               session; Boxplots of consecutive click responses per
%               measurement and layer; Boxplots of consecutive click
%               responses normalized to the first detected sink of the
%               pre-CL

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
cd(homedir);

cd DATA
nThresh = 0.25; 
input = dir('*.mat');
entries = length(input);

para = { 'SinkRMS','SinkPeakAmp','SinkPeakLate','SinkDur','Sinkonset','Sinkoffset'};
layer = {'I_II','IV','V','VI'};
clfrqz = {'twoHz','fiveHz','tenHz','twentyHz','fortyHz'};

%% make list of animals for each group:

Groups = {'KIC','KIT','KIV'};
% KIC is knock in control, no transduction 
% KIT is knock in transduced, TREATED group
% KIV is knock in transduced, viral control
KIClist = {}; KITlist = {}; KIVlist = {};
for iAn = 1:entries
    if strcmp(input(iAn).name(1:3),'KIC')
        KIClist = horzcat(KIClist, input(iAn).name);
    elseif strcmp(input(iAn).name(1:3),'KIT')
        KITlist = horzcat(KITlist, input(iAn).name);
    elseif strcmp(input(iAn).name(1:3),'KIV')
        KIVlist = horzcat(KIVlist, input(iAn).name);
    end
end

%% Pulling and Sorting

for iG = 1:length(Groups)
    
    clear Tuning
    clear Clicks
    clear ClickNorm
    
    % Determine which list of animals to use
    if strcmp(Groups{iG},'KIC')
        Anlist = KIClist;
    elseif strcmp(Groups{iG},'KIT')
        Anlist = KITlist;
    elseif strcmp(Groups{iG},'KIV')
        Anlist = KIVlist;
    end
    
    % get list of conditions for this group
    cd(homedir); cd groups
    run([Groups{iG} '.m'])
    clear channels animals Layer
    conditions = fieldnames(Cond);
    
    % classical tonotopy tuning curves:
    % Self-tuning so that each layer has its own BF from sinkRMS
    NumFreq = 15;
    BF_Pos = 8;
    Tuning = struct;
    
    for iGsub = 1:length(Anlist)
        
        cd(homedir); cd DATA
        load (Anlist{iGsub})
        CurAn = (Anlist{iGsub}(1:5));
        
        for imeas = 1:length(Data)
            condname = Data(imeas).Condition;
            
            if isempty(Data(imeas).Condition)
                continue
            end
            
            if contains(Data(imeas).Condition,'Pre_') || contains(Data(imeas).Condition,'tono')
                %% classical tonotopy tuning curves:
                % Self-tuning so that each layer has its own BF from sinkRMS
                
                for ipar = 1:length(para)
                    
                    for ilay = 1:length(layer)
                        
                        % create a container to collect the sorted features
                        % based on the layer's sink rms BF:
                        Tuning.(layer{ilay}).(condname).(para{ipar})(iGsub,:) = nan(1,NumFreq);
                        
                        % Layer BF
                        if ilay == 1
                            rmsBF = Data(imeas).BF_II;
                        elseif ilay == 2
                            rmsBF = Data(imeas).BF_IV;
                        elseif ilay == 3
                            rmsBF = Data(imeas).BF_V;
                        elseif ilay == 4
                            rmsBF = Data(imeas).BF_VI;
                        end
                        
                        if isempty(rmsBF)
                            % if this layer doesn't have a BF, use BF_IV
                            BF = find(Data(imeas).Frqz == Data(imeas).BF_IV);
                        else
                            BF = find(Data(imeas).Frqz == rmsBF);
                        end
                        
                        for istim = 1:size(Data(imeas).(para{ipar}),2)
                            % don't add a pause or click to the tuning
                            if Data(imeas).Frqz(istim) == 0 || isinf(Data(imeas).Frqz(istim))
                                continue
                            end
                            % check that the first sink is between 0:65 ms
                            if 0 > Data(imeas).Sinkonset(istim).(layer{ilay})(1) < 65
                                % if it is, take this value for the
                                % parameter.layer; Center the BF on BF_Pos
                                Tuning.(layer{ilay}).(condname).(para{ipar})(iGsub,BF_Pos-BF+istim) = ...
                                    Data(imeas).(para{ipar})(istim).(layer{ilay})(1);
                            end
                        end %stimulus
                    end %layer
                end %parameter
                
            elseif contains(Data(imeas).Condition,'CL_')
                %% Clicks
                
                for ipar = 1:length(para)
                    
                    for ilay = 1:length(layer)
                        
                        for istim = 1:size(Data(imeas).(para{ipar}),2)
                            
                            % 2Hz
                            if istim == 1
                                % create a container to collect the
                                % consecutive sinks after clicks 
                                % NOTE: currently supragranular sinks are not being taken as they are ~80 ms
                                Clicks.(layer{ilay}).(condname).(clfrqz{istim}).(para{ipar})(iGsub,:)...
                                    = consec_sinks(Data(imeas).(para{ipar})(istim).(layer{ilay}),...
                                    Data(imeas).Sinkonset(istim).(layer{ilay}),...
                                    2, 1000, 1); % 2 sinks, 1000 dur, 1 ms after click start detection
                                
%                           % 5Hz
                            elseif istim == 2
                                
                                Clicks.(layer{ilay}).(condname).(clfrqz{istim}).(para{ipar})(iGsub,:)...
                                    = consec_sinks(Data(imeas).(para{ipar})(istim).(layer{ilay}),...
                                    Data(imeas).Sinkonset(istim).(layer{ilay}),...
                                    5, 1000, 1); % 5 sinks, 1000 dur, 1 ms after click start detection
                                
                            % 10Hz
                            elseif istim == 3
                                
                                Clicks.(layer{ilay}).(condname).(clfrqz{istim}).(para{ipar})(iGsub,:)...
                                    = consec_sinks(Data(imeas).(para{ipar})(istim).(layer{ilay}),...
                                    Data(imeas).Sinkonset(istim).(layer{ilay}),...
                                    10, 1000, 1); % 10 sinks, 1000 dur, 1 ms after click start detection
                                
                            % 20Hz
                            elseif istim == 4
                                
                                Clicks.(layer{ilay}).(condname).(clfrqz{istim}).(para{ipar})(iGsub,:)...
                                    = consec_sinks(Data(imeas).(para{ipar})(istim).(layer{ilay}),...
                                    Data(imeas).Sinkonset(istim).(layer{ilay}),...
                                    20, 1000, 1); % 20 sinks, 1000 dur, 1 ms after click start detection
                                
                            % 40Hz
                            elseif istim == 5
                                
                               Clicks.(layer{ilay}).(condname).(clfrqz{istim}).(para{ipar})(iGsub,:)...
                                    = consec_sinks(Data(imeas).(para{ipar})(istim).(layer{ilay}),...
                                    Data(imeas).Sinkonset(istim).(layer{ilay}),...
                                    40, 1000, 1); % 40 sinks, 1000 dur, 1 ms after click start detection
                            
                            end %which stim
                            
                            if contains(para{ipar},'RMS') || contains(para{ipar},'PeakAmp')
                                
                                preclicks = Clicks.(layer{ilay}).preCL_1.(clfrqz{istim}).(para{ipar})(iGsub,:);
                                nonanpre = preclicks(~isnan(preclicks));
                                curclicks = Clicks.(layer{ilay}).(condname).(clfrqz{istim}).(para{ipar})(iGsub,:);
                                
                                if ~isempty(nonanpre)
                                    ClickNorm.(layer{ilay}).(condname).(clfrqz{istim}).(para{ipar})(iGsub,:) = ...
                                        curclicks ./ nonanpre(:,1);
                                else
                                    ClickNorm.(layer{ilay}).(condname).(clfrqz{istim}).(para{ipar})(iGsub,:) = ...
                                        curclicks ./ preclicks(:,1);
                                end
                            end
                        end %stimulus
                    end %layer
                end % parameter
                
            elseif contains(Data(imeas).Condition,'AM_')
                %% Amplitude Modulation
                
                for ipar = 1:length(para)
                    
                    for ilay = 1:length(layer)
                        
                        for istim = 1:size(Data(imeas).(para{ipar}),2)
                            
                            % 2Hz
                            if istim == 1
                                % create a container to collect the
                                % consecutive sinks after clicks 
                                % NOTE: currently supragranular sinks are not being taken as they are ~80 ms
                                AMs.(layer{ilay}).(condname).(clfrqz{istim}).(para{ipar})(iGsub,:)...
                                    = consec_sinks(Data(imeas).(para{ipar})(istim).(layer{ilay}),...
                                    Data(imeas).Sinkonset(istim).(layer{ilay}),...
                                    2, 1000, 1); % 2 sinks, 1000 dur, 1 ms after click start detection
                                
%                           % 5Hz
                            elseif istim == 2
                                
                                AMs.(layer{ilay}).(condname).(clfrqz{istim}).(para{ipar})(iGsub,:)...
                                    = consec_sinks(Data(imeas).(para{ipar})(istim).(layer{ilay}),...
                                    Data(imeas).Sinkonset(istim).(layer{ilay}),...
                                    5, 1000, 1); % 5 sinks, 1000 dur, 1 ms after click start detection
                                
                            % 10Hz
                            elseif istim == 3
                                
                                AMs.(layer{ilay}).(condname).(clfrqz{istim}).(para{ipar})(iGsub,:)...
                                    = consec_sinks(Data(imeas).(para{ipar})(istim).(layer{ilay}),...
                                    Data(imeas).Sinkonset(istim).(layer{ilay}),...
                                    10, 1000, 1); % 10 sinks, 1000 dur, 1 ms after click start detection
                                
                            % 20Hz
                            elseif istim == 4
                                
                                AMs.(layer{ilay}).(condname).(clfrqz{istim}).(para{ipar})(iGsub,:)...
                                    = consec_sinks(Data(imeas).(para{ipar})(istim).(layer{ilay}),...
                                    Data(imeas).Sinkonset(istim).(layer{ilay}),...
                                    20, 1000, 1); % 20 sinks, 1000 dur, 1 ms after click start detection
                                
                            % 40Hz
                            elseif istim == 5
                                
                               AMs.(layer{ilay}).(condname).(clfrqz{istim}).(para{ipar})(iGsub,:)...
                                    = consec_sinks(Data(imeas).(para{ipar})(istim).(layer{ilay}),...
                                    Data(imeas).Sinkonset(istim).(layer{ilay}),...
                                    40, 1000, 1); % 40 sinks, 1000 dur, 1 ms after click start detection
                            
                            end %which stim
                            
                            if contains(para{ipar},'RMS') || contains(para{ipar},'PeakAmp')
                                
                                preclicks = AMs.(layer{ilay}).preAM_1.(clfrqz{istim}).(para{ipar})(iGsub,:);
                                nonanpre = preclicks(~isnan(preclicks));
                                curclicks = AMs.(layer{ilay}).(condname).(clfrqz{istim}).(para{ipar})(iGsub,:);
                                
                                if ~isempty(nonanpre)
                                    AMNorm.(layer{ilay}).(condname).(clfrqz{istim}).(para{ipar})(iGsub,:) = ...
                                        curclicks ./ nonanpre(:,1);
                                else
                                    AMNorm.(layer{ilay}).(condname).(clfrqz{istim}).(para{ipar})(iGsub,:) = ...
                                        curclicks ./ preclicks(:,1);
                                end
                            end
                        end %stimulus
                    end %layer
                end % parameter
            end %which meas in which struct
        end %measurement
    end % animal
    
    %% Some Plots
    
    cd(homedir);cd figs
    mkdir(['Group_' (Groups{iG})]); cd(['Group_' (Groups{iG})])
    bfposticks = {'-7' '-6' '-5' '-4' '-3' '-2' '-1' 'BF' '+1' '+2' '+3' '+4' '+5' '+6'};
    take_tono = {'Pre_4','preAMtono_4','preCLtono_4','CLtono_1','AMtono_1'};
    
    % Layer-wise Sink RMS of tonotopies through recording day
    y_min = nan(1,length(take_tono));
    y_max = nan(1,length(take_tono));
     
    h = figure('Name','Tonos_SinkRMS','Position',[500 500 1000 700]); 
    for itono = 1:length(take_tono)
        
        subplot(2,3,itono) %not currently able to take more than 6 subplots
        IImean = nanmean(Tuning.I_II.(take_tono{itono}).SinkRMS);
        IIstd = nanstd(Tuning.I_II.(take_tono{itono}).SinkRMS);
        IVmean = nanmean(Tuning.IV.(take_tono{itono}).SinkRMS);
        IVstd = nanstd(Tuning.IV.(take_tono{itono}).SinkRMS);
        Vmean = nanmean(Tuning.V.(take_tono{itono}).SinkRMS);
        Vstd = nanstd(Tuning.V.(take_tono{itono}).SinkRMS);
        VImean = nanmean(Tuning.VI.(take_tono{itono}).SinkRMS);
        VIstd = nanstd(Tuning.VI.(take_tono{itono}).SinkRMS);
                
        errorbar(IImean,IIstd,'LineWidth',2); hold on
        errorbar(IVmean,IVstd,'LineWidth',2);
        errorbar(Vmean,Vstd,'LineWidth',2);
        errorbar(VImean,VIstd,'LineWidth',2);
        set(gca,'XTick',1:1:15); set(gca,'XTickLabel',bfposticks,'FontSize',8);
        title([(take_tono{itono}) ' RMS'],'FontSize',10,'FontWeight','bold');
        
        if itono == 1
            legend(layer)
        end
        y_min(itono) = min(ylim);
        y_max(itono) = max(ylim);
        AX(itono) = subplot(2,3,itono);
        
    end
    
    for itono = 1:length(take_tono)
        set(AX(itono),'YLim',[min(y_min) max(y_max)])
    end
    
    savefig(h,[(Groups{iG}) '_Tonos_SinkRMS']); close all
    clear AX y_min y_max
    
    % Click boxplots
    take_click = {'preCL_1','CL_1','CL_2','CL_3','CL_4'}; %pre = before laser
    parastr = { 'SinkRMS','SinkPeakAmp'};
    for istim = 1:length(clfrqz)
        for ipar = 1:length(parastr)
            % 1 fig per frequency; rows are layer, columns are conditions
            h = figure('Name',['Clicks_' (parastr{ipar}) '_' (clfrqz{istim})],'Position',[100 300 2300 1000]);
            pos = 0;
            
            for ilay = 1:length(layer)
                for icond = 1:length(take_click)
                    
                    % boxplot per condition per layer
                    pos = pos+1;
                    subplot(length(layer),length(take_click),pos)
                    AX(pos) = subplot(length(layer),length(take_click),pos);
                    
                    boxplot(Clicks.(layer{ilay}).(take_click{icond}).(clfrqz{istim}).(parastr{ipar}))
                    hold on
                    if pos <= length(take_click)
                        title(take_click{icond})
                    end
                    if pos == 1 || pos == 6 || pos == 11 || pos == 16
                        ylabel(layer{ilay},'FontSize',12,'FontWeight','bold','Color','b')
                    end
                                     
                end %layer
            end %condition/measurement
            
            % in order to standardize the y axis. Throw out the upper and
            % lower outliers and then set each subplot to the y min and max
            if contains(parastr{ipar},'RMS')
                for pos = 1:length(take_click)*length(layer)
                    set(AX(pos),'YLim',[-0.002 0.007])
                end
            elseif contains(parastr{ipar},'PeakAmp')
                for pos = 1:length(take_click)*length(layer)
                    set(AX(pos),'YLim',[0 0.01])
                end
            end
            savefig(h,[(Groups{iG}) '_Clicks_' (parastr{ipar}) '_' (clfrqz{istim})],'compact'); 
        end %parameter
    end %click stimulus
    
    % Normalized Click boxplots
    take_click = {'preCL_1','CL_1','CL_2','CL_3','CL_4'}; %pre = before laser
    
    for istim = 1:length(clfrqz)
        for ipar = 1:length(parastr)
            % 1 fig per frequency; rows are layer, columns are conditions
            h = figure('Name',['NormClicks_' (parastr{ipar}) '_' (clfrqz{istim})],'Position',[100 300 2300 1000]);
            pos = 0;
            
            for ilay = 1:length(layer)
                for icond = 1:length(take_click)
                    
                    % boxplot per condition per layer
                    pos = pos+1;
                    subplot(length(layer),length(take_click),pos)
                    AX(pos) = subplot(length(layer),length(take_click),pos);
                    
                    boxplot(ClickNorm.(layer{ilay}).(take_click{icond}).(clfrqz{istim}).(parastr{ipar}))
                    hold on
                    if pos <= length(take_click)
                        title(take_click{icond})
                    end
                    if pos == 1 || pos == 6 || pos == 11 || pos == 16
                        ylabel(layer{ilay},'FontSize',12,'FontWeight','bold','Color','b')
                    end
                                     
                end %layer
            end %condition/measurement
            
            % in order to standardize the y axis. Throw out the upper and
            % lower outliers and then set each subplot to the y min and max
            if contains(parastr{ipar},'RMS')
                for pos = 1:length(take_click)*length(layer)
                    set(AX(pos),'YLim',[-1 4])
                end
            elseif contains(parastr{ipar},'PeakAmp')
                for pos = 1:length(take_click)*length(layer)
                    set(AX(pos),'YLim',[-1 4])
                end
            end
            savefig(h,[(Groups{iG}) '_NormClicks_' (parastr{ipar}) '_' (clfrqz{istim})],'compact'); 
        end %parameter
    end %click stimulus
    
    %% Save it!
    close all
    cd(homedir); cd DATA; cd Output;
    save([(Groups{iG}) '_Tuning_Avg'],'Tuning');
    save([(Groups{iG}) '_ClickSinks'],'Clicks','ClickNorm')
    save([(Groups{iG}) '_AMSinks'],'AMs','AMNorm')
    cd(homedir);
    
end % group

