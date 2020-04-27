function spectralplots(homedir,type)

% Input:        WT_CL.mat generated through computeCSD_scalogram_mice.m
% Output:       Plots showing power the oscillatory frequencies sorted by
%               group and then by osci freq -> figs/Group_Spectral_Plots

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

% frequencies can be found in WT_CL.freq{1} to clarify the following
% rows choices; actual intended rows commented
theta       = (49:54);   %(4:7);
alpha       = (44:48);   %(8:12);
beta_low    = (39:43);   %(13:18);
beta_high   = (34:38);   %(19:30);
gamma_low   = (26:33);   %(31:60);
gamma_high  = (19:25);   %(61:100);

osciName  = {'theta (4:7)' 'alpha (8:12)' 'low beta (13:18)' ...
    'high beta (19:30)' 'low gamma (31:60)' 'high gamma (61:100)'};
osciRows  = {theta alpha beta_low beta_high gamma_low gamma_high};
osciColor = {[0.007, 0.196, 1],[0.854, 0.039, 1],[1, 0.196, 0.007],...
    [1, 0.819, 0.039],[0.058, 0.788, 0.031],[0.039, 0.968, 1]};
Meas    = unique(WT.measurement,'stable'); Meas = Meas(1:5);
Stim    = unique(WT.stimulus,'stable');
Layer   = unique(WT.layer,'stable');
Group   = unique(WT.group,'stable');

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
                    startsWith(WT.layer,Layer{iLay}) &...
                    endsWith(WT.layer,Layer{iLay}),:);
                
                hAx(iMea) = subplot(2,ceil(length(Meas)/2),iMea);
                holdPower = nan(size(curWT,1),length(curWT.scalogram{1,1}));
                for iOsc = 1:length(osciRows)
                    for iWT = 1:size(curWT,1)
                    % loop through oscillation frequencies
                    curOsci = curWT.scalogram{iWT,1}(osciRows{iOsc},:);
                    % get power
                    curOsci = abs(curOsci).^2;
                    mean_power = mean(curOsci);
                    
                    holdPower(iWT,:) = mean_power;
                    
                    end %which scalogram
                    
                    y = mean(holdPower,1);
                    x = 1:length(y);
                    err = std(holdPower,1,1);
                    
                    p = shadedErrorBar(x,y,err,'lineprops',{'color',osciColor{iOsc}});
                    % turn off legend for shade and upper/lower edges
                    set(get(get(p.patch,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                    set(get(get(p.edge(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                    set(get(get(p.edge(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                end %oscillation
                
                legend(osciName)
                title(Meas{iMea})
                ylm(iMea,:) = hAx(iMea).YLim;
                
            end %measurement
            
            ylim(hAx,[min(ylm(:,1)) max(ylm(:,2))]);
            
            h_title = ['Spectral Power of group ' type ' ' Group{iGro} ' layer ' ...
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
osciRows  = {theta alpha beta_low beta_high gamma_low gamma_high};
osciColor = {[0.007, 0.196, 0.949],[0.949, 0.196, 0.007],[0.058, 0.788, 0.031]};

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
                
                %pull out which layer, stim, and measurement
                curWT = WT(startsWith(WT.measurement,Meas{iMea}) & ...
                    WT.stimulus == Stim(iSti) & ...
                    startsWith(WT.layer,Layer{iLay}) &...
                    endsWith(WT.layer,Layer{iLay}),:);
                
                hAx(iMea) = subplot(floor(length(Meas)/2),ceil(length(Meas)/2),iMea);
                % loop through groups to make curves
                for iGro = 1:length(Group)
                    
                    % pull out which group
                    groupWT = curWT(startsWith(curWT.group,Group{iGro}),:);
                    holdPower = nan(size(groupWT,1),length(groupWT.scalogram{1,1}));
                    
                    for iWT = 1:size(groupWT,1)
                        
                        curOsci = groupWT.scalogram{iWT,1}(osciRows{iOsc},:);
                        % get power
                        curOsci = abs(curOsci).^2;
                        mean_power = mean(curOsci);
                        
                        holdPower(iWT,:) = mean_power;
                        
                    end %which scalogram
                    
                    y = mean(holdPower,1);
                    x = 1:length(y);
                    err = std(holdPower,1,1);
                    
                    p = shadedErrorBar(x,y,err,'lineprops',{'color',osciColor{iGro}});
                    % turn off legend for shade and upper/lower edges
                    set(get(get(p.patch,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                    set(get(get(p.edge(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                    set(get(get(p.edge(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                end %oscillation
                
                legend(Group)
                title(Meas{iMea})
                ylm(iMea,:) = hAx(iMea).YLim;
                
            end %measurement
            
            ylim(hAx,[min(ylm(:,1)) max(ylm(:,2))]);
            
            h_title = [osciName{iOsc} ' power of groups ' type ' for layer ' ...
                Layer{iLay} ' click freq ' num2str(Stim(iSti))];
            sgtitle(h_title)
            saveas(h,h_title)
            print(h,'-dpng',h_title)
            print(h,'-dpdf',h_title)
            close(h)
        end %stimulus
    end %layer
end %groups
