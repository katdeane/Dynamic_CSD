%% README

%  This script takes the data mats as input and produces gramm and matlab
%  plots for the RELRES as output.

% Note: Currently only running with one condition/one measurement - needs to be updated for
% more

%% Run
cd('C:\Users\Katrina\Documents\CortXplorers\MyCode\Dynamic CSD');
warning('OFF');
dbstop if error

home = pwd;
addpath(genpath(home));

cd DATA;
input = dir('*.mat');
entries = length(input);
plotdata = struct;
plotdata.BF = []; plotdata.min_one = []; plotdata.min_two = []; plotdata.min_three = [];
plotdata.plus_one = []; plotdata.plus_two = []; plotdata.plus_three = [];
plotdata.time = [];
plotdata.group = []; 

for i1 = 1:entries
    %% Preliminary functions
    disp(['Analyzing Group: ' (input(i1).name(1:end-9))])
    tic
    load (input(i1).name) %loads data into workspace
    
    cd(home); cd figs; %opens figure folder to store future images
    mkdir(['Group_' input(i1).name(1:end-12)]); cd(['Group_' input(i1).name(1:end-9)]);

    names = fieldnames(Data); %list of animal names
    avrec = struct;
    
    time = (1:1:600);
    groupname = (input(i1).name(1:end-12));
    if strcmp ('Awake10dB',groupname)
        groupname = 'Awake'; ANlist = repmat({groupname}, 1, 9);
    elseif strcmp('AnesthetizedPre',groupname)
        groupname = 'Anesthetized'; AWlist = repmat({groupname}, 1, 9);
    else
        Mlist = repmat({groupname}, 1, 9);
    end
    grouplist = repmat({groupname}, 1, 600);
   
   for i2 = 1:length(names);
       animal = names{i2};
       %BF
       BF = find((Data.(names{i2}).GS_BF) == (Data.(names{i2}).Frqz'));
       avrec.(animal).BF =  Data.(animal).RELRES_raw{1,BF}';
       
       plotdata.BF = horzcat(plotdata.BF, avrec.(animal).BF(1:600));
       plotdata.time = horzcat(plotdata.time, time);
       plotdata.group = horzcat(plotdata.group, grouplist);
       %-1
       try 
           avrec.(animal).min_one = Data.(animal).RELRES_raw{1,(BF-1)}';  
       catch
           avrec.(animal).min_one = NaN(1,length(avrec.(animal).BF));
       end
       plotdata.min_one = horzcat(plotdata.min_one, avrec.(animal).min_one(1:600));
       %-2
       try 
           avrec.(animal).min_two = Data.(animal).RELRES_raw{1,(BF-2)}';  
       catch
           avrec.(animal).min_two = NaN(1,length(avrec.(animal).BF));
       end
       plotdata.min_two = horzcat(plotdata.min_two, avrec.(animal).min_two(1:600));
       %-3
       try 
           avrec.(animal).min_three = Data.(animal).RELRES_raw{1,(BF-3)}';  
       catch
           avrec.(animal).min_three = NaN(1,length(avrec.(animal).BF));
       end
       plotdata.min_three = horzcat(plotdata.min_three, avrec.(animal).min_three(1:600));
       %+1
       try 
           avrec.(animal).plus_one = Data.(animal).RELRES_raw{1,(BF+1)}';  
       catch
           avrec.(animal).plus_one = NaN(1,length(avrec.(animal).BF));
       end
       plotdata.plus_one = horzcat(plotdata.plus_one, avrec.(animal).plus_one(1:600));
       %+2
       try 
           avrec.(animal).plus_two = Data.(animal).RELRES_raw{1,(BF+2)}';  
       catch
           avrec.(animal).plus_two = NaN(1,length(avrec.(animal).BF));
       end
       plotdata.plus_two = horzcat(plotdata.plus_two, avrec.(animal).plus_two(1:600));
       %+3
       try 
           avrec.(animal).plus_three = Data.(animal).RELRES_raw{1,(BF+3)}';  
       catch
           avrec.(animal).plus_three = NaN(1,length(avrec.(animal).BF));
       end
       plotdata.plus_three = horzcat(plotdata.plus_three, avrec.(animal).plus_three(1:600));
   end
   
   %plots
   %BF
    h = figure('Name','Relres - BF');
    average = [];
    for i3 = 1:length(names);
        line = avrec.(names{i3}).BF;
        line = line(1:600); %GXL have 600ms and GKD have 800!
        plot(line);
        average = [average line'];
        hold on
    end
    
    average = nanmean(average,2)';
    plot(average,'k','LineWidth',2);
    legend(names, 'Average')
    
    set(h, 'PaperType', 'A4');
    set(h, 'PaperOrientation', 'landscape');
    set(h, 'PaperUnits', 'centimeters');
    savefig(h, 'Relres - BF','compact')
    close (h)
    
    %BF -1
    h = figure('Name','Relres - BF -1');
    average = [];
    for i3 = 1:length(names);
        line = avrec.(names{i3}).min_one;
        line = line(1:600); %GXL have 600ms and GKD have 800!
        plot(line);
        average = [average line'];
        hold on
    end
    
    average = nanmean(average,2)';
    plot(average,'k','LineWidth',2);
    legend(names, 'Average')
    
    set(h, 'PaperType', 'A4');
    set(h, 'PaperOrientation', 'landscape');
    set(h, 'PaperUnits', 'centimeters');
    savefig(h, 'Relres - BF -1','compact')
    close (h)
    
    %BF -2
    h = figure('Name','Relres - BF -2');
    average = [];
    for i3 = 1:length(names);
        line = avrec.(names{i3}).min_two;
        line = line(1:600); %GXL have 600ms and GKD have 800!
        plot(line);
        average = [average line'];
        hold on
    end
    
    average = nanmean(average,2)';
    plot(average,'k','LineWidth',2);
    legend(names, 'Average')
    
    set(h, 'PaperType', 'A4');
    set(h, 'PaperOrientation', 'landscape');
    set(h, 'PaperUnits', 'centimeters');
    savefig(h, 'Relres - BF -2','compact')
    close (h)
    
     %BF -3
    h = figure('Name','Relres - BF -3');
    average = [];
    for i3 = 1:length(names);
        line = avrec.(names{i3}).min_three;
        line = line(1:600); %GXL have 600ms and GKD have 800!
        plot(line);
        average = [average line'];
        hold on
    end
    
    average = nanmean(average,2)';
    plot(average,'k','LineWidth',2);
    legend(names, 'Average')
    
    set(h, 'PaperType', 'A4');
    set(h, 'PaperOrientation', 'landscape');
    set(h, 'PaperUnits', 'centimeters');
    savefig(h, 'Relres - BF -3','compact')
    close (h)
   
    %BF +1
    h = figure('Name','Relres - BF +1');
    average = [];
    for i3 = 1:length(names);
        line = avrec.(names{i3}).plus_one;
        line = line(1:600); %GXL have 600ms and GKD have 800!
        plot(line);
        average = [average line'];
        hold on
    end
    
    average = nanmean(average,2)';
    plot(average,'k','LineWidth',2);
    legend(names, 'Average')
    
    set(h, 'PaperType', 'A4');
    set(h, 'PaperOrientation', 'landscape');
    set(h, 'PaperUnits', 'centimeters');
    savefig(h, 'Relres - BF +1','compact')
    close (h)
    
    %BF +2
    h = figure('Name','Relres - BF +2');
    average = [];
    for i3 = 1:length(names);
        line = avrec.(names{i3}).plus_two;
        line = line(1:600); %GXL have 600ms and GKD have 800!
        plot(line);
        average = [average line'];
        hold on
    end
    
    average = nanmean(average,2)';
    plot(average,'k','LineWidth',2);
    legend(names, 'Average')
    
    set(h, 'PaperType', 'A4');
    set(h, 'PaperOrientation', 'landscape');
    set(h, 'PaperUnits', 'centimeters');
    savefig(h, 'Relres - BF +2','compact')
    close (h)
    
     %BF +3
    h = figure('Name','Relres - BF +3');
    average = [];
    for i3 = 1:length(names);
        line = avrec.(names{i3}).plus_three;
        line = line(1:600); %GXL have 600ms and GKD have 800!
        plot(line);
        average = [average line'];
        hold on
    end
    
    average = nanmean(average,2)';
    plot(average,'k','LineWidth',2);
    legend(names, 'Average')
    
    set(h, 'PaperType', 'A4');
    set(h, 'PaperOrientation', 'landscape');
    set(h, 'PaperUnits', 'centimeters');
    savefig(h, 'Relres - BF +3','compact')
    close (h)
    close all;
end

Ticks = {'BF', 'min_one', 'min_two', 'min_three', 'plus_one', 'plus_two', 'plus_three'};
grammticks = {'d-BF','c-1','b-2','a-3','e+1','f+2','g+3'};
Ticknames = {'Best Frequency', 'BF - 1','BF - 2','BF - 3','BF + 1', 'BF + 2', 'BF + 3'};
%% Gramm Plots
cd(home); cd figs; %opens figure folder to store future images
mkdir('Group Plots Relres'); cd('Group Plots Relres');

for iticks = 1:length(Ticks)
    clear g
    g=gramm('x',plotdata.time,'y',plotdata.(Ticks{iticks}),'color',plotdata.group);
    g.stat_summary('type','sem','geom','area'); %mean and sem shown
    g.set_layout_options('Position',[0 0 0.7 0.7],...
            'legend_pos',[0.71 0.66 0.2 0.2],... %We detach the legend from the plot and move it to the top right
            'margin_height',[0.1 0.1],...
            'margin_width',[0.1 0.1],...
            'redraw',false);
    g.set_names('x','Time (ms)','y','%','color','Group');
    g.axe_property('xlim',[0 600],'ylim',[-0.3 0.2])
    g.set_color_options('map','matlab');
    g.set_title((Ticknames{iticks}));
    g.draw();
    g.export('file_name',(Ticks{iticks}), 'file_type','png');
    g.export('file_name',(Ticks{iticks}), 'file_type','pdf');
    close all;
end

%% This level of statistics holds one data point for each animal rather than the 50 

%% Stats
levels = [2,1];
cd(home); cd DATA
mkdir('Stats Relres'); cd('Stats Relres')

%set up the variables
plotdata.peakamp = []; plotdata.peaklat = []; plotdata.peakrms = []; plotdata.totalrms = []; plotdata.frqz = []; plotdata.smallgroup = [];
ANsize = sum(strcmp('Anesthetized',plotdata.group))/600; %change to before rerun: ANsize = sum(strcmp('AnesthetizedPre',plotdata.group))/600;
AWsize = sum(strcmp('Awake',plotdata.group))/600;
Msize = sum(strcmp('Muscimol',plotdata.group))/600;
winstart = (200:600:size(plotdata.group,2)); winend = (300:600:size(plotdata.group,2));
totalstart = (1:600:size(plotdata.group,2)); totalend = (600:600:size(plotdata.group,2));
sixhun = (0:600:size(plotdata.group,2));

for iticks = 1:length(Ticks)
    %all the variables
    ANpeakamp = []; ANpeaklat = []; ANtotalvalue = []; ANpeakvalue = [];
    AWpeakamp = []; AWpeaklat = []; AWtotalvalue = []; AWpeakvalue = [];
    Mpeakamp = []; Mpeaklat = []; Mtotalvalue = []; Mpeakvalue = [];

    disp(['***********************STATS FOR ' (Ticks{iticks}) '***********************'])
    for istack = 1:ANsize
        ANpeakamp = [ANpeakamp nanmin(plotdata.(Ticks{iticks})(winstart(istack):winend(istack)))];
        if ~isnan(ANpeakamp(istack))
            temppeaklat = (find(plotdata.(Ticks{iticks})(winstart(istack):winend(istack)) == ANpeakamp(istack)));
            ANpeaklat = [ANpeaklat temppeaklat(1)];
        else 
            ANpeaklat = [ANpeaklat NaN];
        end
        ANtotalvalue = [ANtotalvalue nanmean(abs(plotdata.(Ticks{iticks})(totalstart(istack):totalend(istack))))];
        ANpeakvalue = [ANpeakvalue nanmean(abs(plotdata.(Ticks{iticks})(winstart(istack):winend(istack))))];
    end
    
    for istack = ANsize+1:ANsize+AWsize
        AWpeakamp = [AWpeakamp nanmin(plotdata.(Ticks{iticks})(winstart(istack):winend(istack)))];
        if ~isnan(AWpeakamp(istack-(ANsize)))
            temppeaklat = (find(plotdata.(Ticks{iticks})(winstart(istack):winend(istack)) == AWpeakamp(istack-(ANsize))));
            AWpeaklat = [AWpeaklat temppeaklat(1)]; %there may be two values and I want just the first
        else
            AWpeaklat = [AWpeaklat NaN];
        end
        AWtotalvalue = [AWtotalvalue nanmean(abs(plotdata.(Ticks{iticks})(totalstart(istack):totalend(istack))))];
        AWpeakvalue = [AWpeakvalue nanmean(abs(plotdata.(Ticks{iticks})(winstart(istack):winend(istack))))];
    end
    
    for istack = ANsize+AWsize+1:ANsize+AWsize+Msize
        Mpeakamp = [Mpeakamp nanmin(plotdata.(Ticks{iticks})(winstart(istack):winend(istack)))];
        if ~isnan(Mpeakamp(istack-(ANsize+AWsize)))
            temppeaklat = (find(plotdata.(Ticks{iticks})(winstart(istack):winend(istack)) == Mpeakamp(istack-(ANsize+AWsize))));
            Mpeaklat = [Mpeaklat temppeaklat];
        else 
            Mpeaklat = [Mpeaklat NaN];
        end
        Mtotalvalue = [Mtotalvalue nanmean(abs(plotdata.(Ticks{iticks})(totalstart(istack):totalend(istack))))];
        Mpeakvalue = [Mpeakvalue nanmean(abs(plotdata.(Ticks{iticks})(winstart(istack):winend(istack))))];
    end
    
    %need 9 values to be equal to awake group
    ANpeakamp = [ANpeakamp NaN];
    ANpeaklat = [ANpeaklat NaN];
    ANtotalvalue = [ANtotalvalue NaN];
    ANpeakvalue = [ANpeakvalue NaN];
    Mpeakamp = [Mpeakamp NaN];
    Mpeaklat = [Mpeaklat NaN];
    Mtotalvalue = [Mtotalvalue NaN];
    Mpeakvalue = [Mpeakvalue NaN];
    
    % setting up for tuning curves
    plotdata.peakamp = horzcat(plotdata.peakamp, ANpeakamp, AWpeakamp, Mpeakamp);
    plotdata.peaklat = horzcat(plotdata.peaklat, ANpeaklat, AWpeaklat, Mpeaklat);
    plotdata.peakrms = horzcat(plotdata.peakrms, ANpeakvalue, AWpeakvalue, Mpeakvalue);
    plotdata.totalrms = horzcat(plotdata.totalrms, ANtotalvalue, AWtotalvalue, Mtotalvalue);
    plotdata.smallgroup = horzcat(plotdata.smallgroup, ANlist, AWlist, Mlist);
    plotdata.frqz = horzcat(plotdata.frqz, repmat(grammticks(iticks), 1, 27));
    
    disp('********Peak Amplitude********')
    peakamp_ANvsAW = horzcat(ANpeakamp', AWpeakamp');
    peakamp_ANvsM = horzcat(ANpeakamp', Mpeakamp');
    var_peakamp = {'Groups','Peak Amplitude'};
    
    disp('**Anesthetized vs Awake**')
    AnesthetizedvsAwake = teg_repeated_measures_ANOVA(peakamp_ANvsAW, levels, var_peakamp);
    disp('**Anesthetized vs Muscimol**')
    AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(peakamp_ANvsM, levels, var_peakamp);
    
    savefile = [(Ticks{iticks}) ' PeakAmp RelresStats.mat'];
    save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol')
    
    disp('********Peak Latitude********')
    peaklat_ANvsAW = horzcat(ANpeaklat', AWpeaklat');
    peaklat_ANvsM = horzcat(ANpeaklat', Mpeaklat');
    var_peaklat = {'Groups','Peak Latitude'};
    
    disp('**Anesthetized vs Awake**')
    AnesthetizedvsAwake = teg_repeated_measures_ANOVA(peaklat_ANvsAW, levels, var_peaklat);
    disp('**Anesthetized vs Muscimol**')
    AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(peaklat_ANvsM, levels, var_peaklat);
    
    savefile = [(Ticks{iticks}) ' PeakLat RelresStats.mat'];
    save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol')
    
    disp('********Total Value********')
    totalvalue_ANvsAW = horzcat(ANtotalvalue', AWtotalvalue');
    totalvalue_ANvsM = horzcat(ANtotalvalue', Mtotalvalue');
    var_totalvalue = {'Groups','Total'};
    
    disp('**Anesthetized vs Awake**')
    AnesthetizedvsAwake = teg_repeated_measures_ANOVA(totalvalue_ANvsAW, levels, var_totalvalue);
    disp('**Anesthetized vs Muscimol**')
    AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(totalvalue_ANvsM, levels, var_totalvalue);
    
    savefile = [(Ticks{iticks}) ' Total RelresStats.mat'];
    save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol')
    
    disp('********Peak Value********')
    peakvalue_ANvsAW = horzcat(ANpeakvalue', AWpeakvalue');
    peakvalue_ANvsM = horzcat(ANpeakvalue', Mpeakvalue');
    var_peakvalue = {'Groups','Total'};
    
    disp('**Anesthetized vs Awake**')
    AnesthetizedvsAwake = teg_repeated_measures_ANOVA(peakvalue_ANvsAW, levels, var_peakvalue);
    disp('**Anesthetized vs Muscimol**')
    AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(peakvalue_ANvsM, levels, var_peakvalue);
    
    savefile = [(Ticks{iticks}) ' PeakValue RelresStats.mat'];
    save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol')
end

%% PEAK AMP Tuning curve
feature = {'peakamp','peaklat', 'peakrms','totalrms'};
Name = {'Peak Amplitude', 'Peak Latency', 'RMS between 200 and 300 ms', 'RMS of total Relres'};
cd(home); cd figs; mkdir('Group Relres Tuning'); cd('Group Relres Tuning')

for i1 = 1:length(feature);
clear g
figure();

g(1,1)=gramm('x',plotdata.frqz,'y',plotdata.(feature{i1}), 'color', plotdata.smallgroup); %,'y',cars.Acceleration,'color',cars.Cylinders,'subset',cars.Cylinders~=3 & cars.Cylinders~=5
g(1,1).stat_summary('type','sem','geom','errorbar'); %mean and sem shown
g(1,1).stat_summary('type','sem','geom','point'); %mean and sem shown
g(1,1).stat_summary('type','sem','geom','area'); %mean and sem shown
g(1,1).set_text_options('base_size',12) 
if strcmp(feature{i1},'peaklat')
    g(1,1).set_names('x','Tone','y','ms','color','Group');
    g(1,1).axe_property('YLim',[10 90]);
elseif strcmp(feature{i1},'peakamp')
    g(1,1).set_names('x','Tone','y','%','color','Group');
    g(1,1).axe_property('YLim',[-0.4 0.1]);
elseif strcmp(feature{i1},'peakrms')
    g(1,1).set_names('x','Tone','y','%','color','Group');
    g(1,1).axe_property('YLim',[0.02 0.12]);
else
    g(1,1).set_names('x','Tone','y','%','color','Group');
    g(1,1).axe_property('YLim',[0.02 0.1]);
end
g(1,1).set_color_options('map','matlab');
g(1,1).axe_property('XTickLabel',{'-3', '-2', '-1', 'BF', '+1', '+2', '+3'})

g.set_title(Name{i1});
g.draw();
g.export('file_name',['RelRes ' Name{i1}], 'file_type','pdf');
g.export('file_name',['RelRes ' Name{i1}], 'file_type','png');
close all;
end