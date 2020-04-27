function spectralplots_scatter(homedir,type)

% Input:        WT_CL.mat generated through computeCSD_scalogram_mice.m
% Output:       Plots showing power the oscillatory frequencies sorted by
%               group but showing individual peaks and then by osci freq 
%               showing individual animal peaks -> figs/Group_Spectral_Plots

% this script could potentially be used to take out single animal features
% for statistics but cannot be used for single trial statistics

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
if ~exist('type','var')
    type = 'CL'; % run the click measurements
end

disp(['loading mat file WT of ' type])
tic
if contains(type,'CL')
    load('WT_CL.mat','WT_CL')
    WT = WT_CL;
    clear WT_CL
elseif contains(type,'AM')
    load('WT_AM.mat','WT_AM')
    WT = WT_AM;
    clear WT_AM
else
    error('Input "type" needs to be either "CL" or "AM" for clicks or amplitude modulations respectively')
end
toc


% frequencies can be found in wtTable.freq{1} to clarify the following
% rows choices; actual intended rows commented
theta = (49:54);        %(4:7);
alpha = (44:48);        %(8:12);
beta_low = (39:43);     %(13:18);
beta_high = (34:38);    %(19:30);
gamma_low = (26:33);    %(31:60);
gamma_high = (19:25);   %(61:100);

osciName  = {'theta (4:7)' 'alpha (8:12)' 'low beta (13:18)' ...
    'high beta (19:30)' 'low gamma (31:60)' 'high gamma (61:100)'};
osciRows  = {theta alpha beta_low beta_high gamma_low gamma_high};
% color picker: http://doc.instantreality.org/tools/color_calculator/
osciColor = {[0.007, 0.196, 1],[0.854, 0.039, 1],[1, 0.196, 0.007],[1, 0.819, 0.039],[0.058, 0.788, 0.031],[0.039, 0.968, 1]};

Meas  = unique(WT.measurement,'stable');
Stim  = unique(WT.stimulus,'stable');
Layer = unique(WT.layer,'stable');
Group = unique(WT.group,'stable');
Animal= unique(WT.animal,'stable');

cd(homedir); cd figs; mkdir(['Group_Spectral_Plots_' type]); cd(['Group_Spectral_Plots_' type]);

%% plots with power; oscillatory bands in one window
for iGro = 1:length(Group)
    % loop through layers
    for iLay = 1:length(Layer)
        % loop through stimulation frequencies
        for iSti = 1:length(Stim)
            % loop through measurement type
            
            h = figure('units','normalized','outerposition',[0 0 1 1]); 
            hold on
            hAx=gobjects(length(Meas),1); 
            ylm = nan(length(Meas),2);
            
            for iMea = 1:length(Meas)
                
                %pull out which group, layer, stim, and measurement
                curWT = WT(startsWith(WT.group,Group{iGro}) & ...
                    startsWith(WT.measurement,Meas{iMea}) & ...
                    WT.stimulus == Stim(iSti) & ...
                    startsWith(WT.layer,Layer{iLay}),:);
                
                hAx(iMea) = subplot(floor(length(Meas)/2),ceil(length(Meas)/2),iMea);
                hold on
                Y = nan(1,size(curWT,1)*Stim(iSti));
                X = nan(1,size(curWT,1)*Stim(iSti));
                icount = 1:Stim(iSti):length(Y)+41;
                for iOsc = 1:length(osciRows)
                    for iWT = 1:size(curWT,1)
                    % loop through oscillation frequencies
                    curOsci = curWT.scalogram{iWT,1}(osciRows{iOsc},:);
                    % get power
                    curOsci = abs(curOsci).^2;
                    
                    [Y(icount(iWT):icount(iWT+1)-1),...
                        X(icount(iWT):icount(iWT+1)-1)] = ...
                        consec_peaks(curOsci, Stim(iSti), 1200, 200);
                                        
                    end %which scalogram
                    plot(X,Y,'*','color',osciColor{iOsc});
                   
                end %oscillation
                
                clicktime = 200;
                for icli = 1:Stim(iSti) 
                    xline(clicktime,'--r');
                    clicktime = clicktime + 1000/Stim(iSti);
                end
                legend(osciName)
                title(Meas{iMea})
                ylm(iMea,:) = hAx(iMea).YLim;
                
            end %measurement
            
            ylim(hAx,[min(ylm(:,1)) max(ylm(:,2))]);
            xlim(hAx,[0 1300]);
            
            h_title = ['Spectral Power of animals ' type ' in group ' Group{iGro} ' layer ' ...
                Layer{iLay} ' for click freq ' num2str(Stim(iSti))];
            sgtitle(h_title)
            saveas(h,h_title)
            print(h,'-dpng',h_title)
            print(h,'-dpdf',h_title)
            close(h)
        end %stimulus
    end %layer
end %groups


%% plots with power; groups in one window

osciName  = {'Theta' 'Alpha' 'Low beta' 'High beta' ...
    'Low gamma' 'High gamma'};
osciColor = {[0.007, 0.196, 1],[0.254, 0.098, 0.901],[0.486, 0.054, 0.866],...
    [0.854, 0.039, 1],[1, 0.039, 0.784],[1, 0.039, 0.372],[1, 0.039, 0.062],...
    [1, 0.196, 0.007],[1, 0.658, 0.007],[1, 0.819, 0.039],[0.815, 1, 0.039],...
    [0.415, 1, 0.039],[0.058, 0.788, 0.031],[0.031, 0.788, 0.294],[0.031, 0.788, 0.533],...
    [0.031, 0.788, 0.772],[0.039, 0.968, 1]};

for iOsc = 1:length(osciRows)
    % loop through layers
    for iLay = 1:length(Layer)
        % loop through stimulation frequencies
        for iSti = 1:length(Stim)
            
            h = figure('units','normalized','outerposition',[0 0 1 1]); 
            hold on
            hAx=gobjects(length(Meas),1); 
            ylm = nan(length(Meas),2);
            
            
            % loop through measurement type
            for iMea = 1:length(Meas)
                icount = 1;
                
                %pull out which layer, stim, and measurement
                curWT = WT(startsWith(WT.measurement,Meas{iMea}) & ...
                    WT.stimulus == Stim(iSti) & ...
                    startsWith(WT.layer,Layer{iLay}) &...
                    endsWith(WT.layer,Layer{iLay}),:);
                
                hAx(iMea) = subplot(floor(length(Meas)/2),ceil(length(Meas)/2),iMea);
                hold on

                % loop through groups to make curves
                for iGro = 1:length(Group)
                    
                    % pull out which group
                    groupWT = curWT(startsWith(curWT.group,Group{iGro}),:);
                    
                    for iWT = 1:size(groupWT,1)
                        
                        curOsci = groupWT.scalogram{iWT,1}(osciRows{iOsc},:);
                        % get power
                        curOsci = abs(curOsci).^2;
                        
                        [Y,X] = consec_peaks(curOsci, Stim(iSti), 1200, 200);
                                               
                        plot(X,Y,'*','color',osciColor{icount});
                        icount = icount+1;
                        
                    end %which scalogram
                end %group
                
                clicktime = 200;
                for icli = 1:Stim(iSti) 
                    xline(clicktime,'--r');
                    clicktime = clicktime + 1000/Stim(iSti);
                end
                legend(Animal)
                title(Meas{iMea})
                ylm(iMea,:) = hAx(iMea).YLim;
                
            end %measurement
            
            ylim(hAx,[min(ylm(:,1)) max(ylm(:,2))]);
            xlim(hAx,[0 1300]);
            
            h_title = [osciName{iOsc} ' peak power of animals ' type ' for layer ' ...
                Layer{iLay} ' click freq ' num2str(Stim(iSti))];
            sgtitle(h_title)
            saveas(h,h_title)
            print(h,'-dpng',h_title)
%             print(h,'-dpdf',h_title)
            close(h)
        end %stimulus
    end %layer
end %groups

