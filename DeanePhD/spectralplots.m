clear

% work in progress;
% Input:        WT_CL.mat generated through computeCSD_scalogram_mice.m
% Output:       Plots showing power the oscillatory frequencies sorted by
%               group and then by osci freq -> figs/Group_Spectral_Plots

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

disp('loading mat file WT_CL')
tic
load('WT_CL.mat','WT_CL')
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
osciColor = {[0.007, 0.196, 1],[0.854, 0.039, 1],[1, 0.196, 0.007],[1, 0.819, 0.039],[0.058, 0.788, 0.031],[0.039, 0.968, 1]};
Meas  = unique(WT_CL.measurement,'stable');
Stim  = unique(WT_CL.stimulus,'stable');
Layer = unique(WT_CL.layer,'stable');
Group = unique(WT_CL.group,'stable');

cd(homedir); cd figs; mkdir Group_Spectral_Plots; cd Group_Spectral_Plots;

%% plots with power; oscillatory bands in one window
for iGro = 1:length(Group)
    % loop through layers
    for iLay = 1:length(Layer)
        % loop through stimulation frequencies
        for iSti = 1:length(Stim)
            % loop through measurement type
            
            h = figure('units','normalized','outerposition',[0 0 1 1]); 
%             title(['Spectral Power of group ' Group{iGro} ' layer ' ...
%                 Layer{iLay} ' for click freq ' num2str(Stim(iSti))])
            hold on
            
            for iMea = 1:length(Meas)
                
                %pull out which group, layer, stim, and measurement
                curWT = WT_CL(startsWith(WT_CL.group,Group{iGro}) & ...
                    startsWith(WT_CL.measurement,Meas{iMea}) & ...
                    WT_CL.stimulus == Stim(iSti) & ...
                    startsWith(WT_CL.layer,Layer{iLay}),:);
                
                subplot(floor(length(Meas)/2),ceil(length(Meas)/2),iMea)
                holdPower = nan(size(curWT,1),length(curWT.scalogram{1,1}));
                for iOsc = 1:length(osciRows)
                    for iSca = 1:size(curWT,1)
                    % loop through oscillation frequencies
                    curOsci = curWT.scalogram{iSca,1}(osciRows{iOsc},:);
                    % get power
                    curOsci = abs(curOsci).^2;
                    mean_power = mean(curOsci);
                    
                    holdPower(iSca,:) = mean_power;
                    
                    % possible points to pull out:
                    % should take the peak after stim selection (consec
                    % peak function may work here)
                    % peak_power = max(max(curOsci));
                    % [row,peaklat] = find(curOsci == peak_power);
                    
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
                
            end %measurement
            h_title = ['Spectral Power of group ' Group{iGro} ' layer ' ...
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

osciName  = {'Low beta' 'High beta' ...
    'Low gamma' 'High gamma'};
osciRows  = {beta_low beta_high gamma_low gamma_high};
osciColor = {[0.007, 0.196, 0.949],[0.949, 0.196, 0.007],[0.058, 0.788, 0.031]};

for iOsc = 1:length(osciRows)
    % loop through layers
    for iLay = 1:length(Layer)
        % loop through stimulation frequencies
        for iSti = 1:length(Stim)
            
            h = figure('units','normalized','outerposition',[0 0 1 1]); 
            hold on
            
            % loop through measurement type
            for iMea = 1:length(Meas)
                
                %pull out which layer, stim, and measurement
                curWT = WT_CL(startsWith(WT_CL.measurement,Meas{iMea}) & ...
                    WT_CL.stimulus == Stim(iSti) & ...
                    startsWith(WT_CL.layer,Layer{iLay}),:);
                
                subplot(floor(length(Meas)/2),ceil(length(Meas)/2),iMea)
                % loop through groups to make curves
                for iGro = 1:length(Group)
                    
                    % pull out which group
                    groupWT = curWT(startsWith(curWT.group,Group{iGro}),:);
                    holdPower = nan(size(groupWT,1),length(groupWT.scalogram{1,1}));
                    
                    for iSca = 1:size(groupWT,1)
                        
                        curOsci = groupWT.scalogram{iSca,1}(osciRows{iOsc},:);
                        % get power
                        curOsci = abs(curOsci).^2;
                        mean_power = mean(curOsci);
                        
                        holdPower(iSca,:) = mean_power;
                        
                        % possible points to pull out:
                        % should take the peak after stim selection (consec
                        % peak function may work here)
                        % peak_power = max(max(curOsci));
                        % [row,peaklat] = find(curOsci == peak_power);
                        
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
                
            end %measurement
            h_title = [osciName{iOsc} ' power of groups for layer ' ...
                Layer{iLay} ' click freq ' num2str(Stim(iSti))];
            sgtitle(h_title)
            saveas(h,h_title)
            print(h,'-dpng',h_title)
            print(h,'-dpdf',h_title)
            close(h)
        end %stimulus
    end %layer
end %groups

    
    
    
    
    
    
    
    
    
    
    
    
    

