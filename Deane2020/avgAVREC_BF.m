%% README

%  This script takes the data mats as input and produces gramm and matlab
%  plots for the averaged AVREC as output.

% Note: Currently only running with one condition/one measurement - needs to be updated for
% more

%% Run
cd('D:\MyCode\Dynamic_CSD_Analysis');
warning('OFF');
dbstop if error

home = pwd;
addpath(genpath(home));
Group = {'Anesthetized','Awake','Muscimol'};

cd DATA;
input = dir('*.mat');
entries = length(input);
plotdata = struct;
%for full graph
plotdata.BF = []; plotdata.min_one = []; plotdata.min_two = []; plotdata.min_three = [];
plotdata.plus_one = []; plotdata.plus_two = []; plotdata.plus_three = [];
plotdata.time = [];
plotdata.group = [];
%for tuning curves
grpsz = 11;
plotdata.frqz = [];
plotdata.tgroup = [];
plotdata.peakamp = [];plotdata.peaklate = [];plotdata.smallmean = [];plotdata.totalmean = [];

for i1 = 1:entries
    %% Preliminary functions
    disp(['Analyzing Group: ' (input(i1).name(1:end-9))])
    tic
    load (input(i1).name) %loads data into workspace
    
    cd(home); cd figs; %opens figure folder to store future images
    mkdir(['AvgAVREC_' input(i1).name(1:end-9)]); cd(['AvgAVREC_' input(i1).name(1:end-9)]);

    names = fieldnames(Data); %list of animal names
    avrec = struct;
    
    time = (1:1:600);
    groupname = (input(i1).name(1:end-12));
    if strcmp ('Awake10dB',groupname)
        groupname = 'Awake';
    elseif strcmp('AnesthetizedPre',groupname)
        groupname = 'Anesthetized';
    end
    grouplist = repmat({groupname}, 1, 600);
   
   for i2 = 1:length(names)
       animal = names{i2};
       %BF
       BF = find((Data.(names{i2}).GS_BF) == (Data.(names{i2}).Frqz'));
       avrec.(animal).BF =  Data.(animal).AVREC_raw{1,BF}';
       
       plotdata.BF = horzcat(plotdata.BF, avrec.(animal).BF(1:600));
       plotdata.time = horzcat(plotdata.time, time);
       plotdata.group = horzcat(plotdata.group, grouplist);
       %-1
       try 
           avrec.(animal).min_one = Data.(animal).AVREC_raw{1,(BF-1)}';  
       catch
           avrec.(animal).min_one = NaN(1,length(avrec.(animal).BF));
       end
       plotdata.min_one = horzcat(plotdata.min_one, avrec.(animal).min_one(1:600));
       %-2
       try 
           avrec.(animal).min_two = Data.(animal).AVREC_raw{1,(BF-2)}';  
       catch
           avrec.(animal).min_two = NaN(1,length(avrec.(animal).BF));
       end
       plotdata.min_two = horzcat(plotdata.min_two, avrec.(animal).min_two(1:600));
       %-3
       try 
           avrec.(animal).min_three = Data.(animal).AVREC_raw{1,(BF-3)}';  
       catch
           avrec.(animal).min_three = NaN(1,length(avrec.(animal).BF));
       end
       plotdata.min_three = horzcat(plotdata.min_three, avrec.(animal).min_three(1:600));
       %+1
       try 
           avrec.(animal).plus_one = Data.(animal).AVREC_raw{1,(BF+1)}';  
       catch
           avrec.(animal).plus_one = NaN(1,length(avrec.(animal).BF));
       end
       plotdata.plus_one = horzcat(plotdata.plus_one, avrec.(animal).plus_one(1:600));
       %+2
       try 
           avrec.(animal).plus_two = Data.(animal).AVREC_raw{1,(BF+2)}';  
       catch
           avrec.(animal).plus_two = NaN(1,length(avrec.(animal).BF));
       end
       plotdata.plus_two = horzcat(plotdata.plus_two, avrec.(animal).plus_two(1:600));
       %+3
       try 
           avrec.(animal).plus_three = Data.(animal).AVREC_raw{1,(BF+3)}';  
       catch
           avrec.(animal).plus_three = NaN(1,length(avrec.(animal).BF));
       end
       plotdata.plus_three = horzcat(plotdata.plus_three, avrec.(animal).plus_three(1:600));
   end
   
   %plots
   %BF
    h = figure('Name','AvgRec - BF');
    average = [];
    for i3 = 1:length(names)
        line = avrec.(names{i3}).BF;
        line = line(1:600); %GXL have 600ms and GKD have 800!
        plot(line);
        average = [average line'];
        hold on
    end
    
    average = nanmean(average,2)';
    plot(average,'k','LineWidth',2);
    allnames = vertcat(names,{'Average'});
    legend(allnames)
    
    set(h, 'PaperType', 'A4');
    set(h, 'PaperOrientation', 'landscape');
    set(h, 'PaperUnits', 'centimeters');
    savefig(h, 'AvgRec - BF','compact')
    close (h)
    
    %BF -1
    h = figure('Name','AvgRec - BF -1');
    average = [];
    for i3 = 1:length(names)
        line = avrec.(names{i3}).min_one;
        line = line(1:600); %GXL have 600ms and GKD have 800!
        plot(line);
        average = [average line'];
        hold on
    end
    
    average = nanmean(average,2)';
    plot(average,'k','LineWidth',2);
    legend(allnames)
    
    set(h, 'PaperType', 'A4');
    set(h, 'PaperOrientation', 'landscape');
    set(h, 'PaperUnits', 'centimeters');
    savefig(h, 'AvgRec - BF -1','compact')
    close (h)
    
    %BF -2
    h = figure('Name','AvgRec - BF -2');
    average = [];
    for i3 = 1:length(names)
        line = avrec.(names{i3}).min_two;
        line = line(1:600); %GXL have 600ms and GKD have 800!
        plot(line);
        average = [average line'];
        hold on
    end
    
    average = nanmean(average,2)';
    plot(average,'k','LineWidth',2);
    legend(allnames)
    
    set(h, 'PaperType', 'A4');
    set(h, 'PaperOrientation', 'landscape');
    set(h, 'PaperUnits', 'centimeters');
    savefig(h, 'AvgRec - BF -2','compact')
    close (h)
    
     %BF -3
    h = figure('Name','AvgRec - BF -3');
    average = [];
    for i3 = 1:length(names)
        line = avrec.(names{i3}).min_three;
        line = line(1:600); %GXL have 600ms and GKD have 800!
        plot(line);
        average = [average line'];
        hold on
    end
    
    average = nanmean(average,2)';
    plot(average,'k','LineWidth',2);
    legend(allnames)
    
    set(h, 'PaperType', 'A4');
    set(h, 'PaperOrientation', 'landscape');
    set(h, 'PaperUnits', 'centimeters');
    savefig(h, 'AvgRec - BF -3','compact')
    close (h)
   
    %BF +1
    h = figure('Name','AvgRec - BF +1');
    average = [];
    for i3 = 1:length(names)
        line = avrec.(names{i3}).plus_one;
        line = line(1:600); %GXL have 600ms and GKD have 800!
        plot(line);
        average = [average line'];
        hold on
    end
    
    average = nanmean(average,2)';
    plot(average,'k','LineWidth',2);
    legend(allnames)
    
    set(h, 'PaperType', 'A4');
    set(h, 'PaperOrientation', 'landscape');
    set(h, 'PaperUnits', 'centimeters');
    savefig(h, 'AvgRec - BF +1','compact')
    close (h)
    
    %BF +2
    h = figure('Name','AvgRec - BF +2');
    average = [];
    for i3 = 1:length(names)
        line = avrec.(names{i3}).plus_two;
        line = line(1:600); %GXL have 600ms and GKD have 800!
        plot(line);
        average = [average line'];
        hold on
    end
    
    average = nanmean(average,2)';
    plot(average,'k','LineWidth',2);
    legend(allnames)
    
    set(h, 'PaperType', 'A4');
    set(h, 'PaperOrientation', 'landscape');
    set(h, 'PaperUnits', 'centimeters');
    savefig(h, 'AvgRec - BF +2','compact')
    close (h)
    
     %BF +3
    h = figure('Name','AvgRec - BF +3');
    average = [];
    for i3 = 1:length(names)
        line = avrec.(names{i3}).plus_three;
        line = line(1:600); %GXL have 600ms and GKD have 800!
        plot(line);
        average = [average line'];
        hold on
    end
    
    average = nanmean(average,2)';
    plot(average,'k','LineWidth',2);
    legend(allnames)
    
    set(h, 'PaperType', 'A4');
    set(h, 'PaperOrientation', 'landscape');
    set(h, 'PaperUnits', 'centimeters');
    savefig(h, 'AvgRec - BF +3','compact')
    close (h)
    close all;
end

Ticks = {'min_three', 'min_two','min_one' ,'BF',  'plus_one', 'plus_two', 'plus_three'};
Ticknames = {'BF - 3','BF - 2', 'BF - 1','Best Frequency','BF + 1', 'BF + 2', 'BF + 3'};
%% Gramm Plots
cd(home); cd figs; %opens figure folder to store future images
mkdir('Gramm Plots AvgAvrec'); cd('Gramm Plots AvgAvrec');

for iticks = 1:length(Ticks)
    clear g
    g=gramm('x',plotdata.time,'y',plotdata.(Ticks{iticks}),'color',plotdata.group);
    g.stat_summary('type','sem','geom','area'); %mean and sem shown
    g.set_layout_options('Position',[0 0 0.7 0.7],...
            'legend_pos',[0.71 0.66 0.2 0.2],... %We detach the legend from the plot and move it to the top right
            'margin_height',[0.1 0.1],...
            'margin_width',[0.1 0.1],...
            'redraw',false);
    g.set_names('x','Time (ms)','y','µV','color','Group');
    g.axe_property('ylim',[0 0.002],'xlim',[0 600])
    g.set_color_options('map','matlab');
    g.set_title((Ticknames{iticks}));
    g.draw();
    g.export('file_name',(Ticks{iticks}), 'file_type','png');
    g.export('file_name',(Ticks{iticks}), 'file_type','pdf');
    close all;
end

%% This level of statistics holds one data point for each animal rather than the 50 

grammTicks = {'a-3', 'b-2','c-1' ,'d BF',  'e+1', 'f+2', 'g+3'};
%% Stats
levels = [2,1];
cd(home); cd DATA
mkdir('Stats avgAvRec'); cd('Stats avgAvRec')

%set up the variables
ANsize = sum(strcmp('Anesthetized',plotdata.group))/600; %change to before rerun: ANsize = sum(strcmp('AnesthetizedPre',plotdata.group))/600;
AWsize = sum(strcmp('Awake1',plotdata.group))/600;
Msize = sum(strcmp('Musci',plotdata.group))/600;
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
        [latency, ampl, meanOfdata] = iGetPeakData_avrec(plotdata.(Ticks{iticks})(winstart(istack):winend(istack)));
        if isnan(meanOfdata)
            ANpeakamp = [ANpeakamp NaN];
            ANpeaklat = [ANpeaklat NaN];
        else
            ANpeakamp = [ANpeakamp ampl];
            ANpeaklat = [ANpeaklat latency];
        end
        ANtotalvalue = [ANtotalvalue nanmean(plotdata.(Ticks{iticks})(totalstart(istack):totalend(istack)))];
        ANpeakvalue = [ANpeakvalue meanOfdata];
    end
    
    for istack = ANsize+1:ANsize+AWsize
        [latency, ampl, meanOfdata] = iGetPeakData_avrec(plotdata.(Ticks{iticks})(winstart(istack):winend(istack)));
        if isnan(meanOfdata)
            AWpeakamp = [AWpeakamp NaN];
            AWpeaklat = [AWpeaklat NaN];
        else
            AWpeakamp = [AWpeakamp ampl];
            AWpeaklat = [AWpeaklat latency];
        end
        AWtotalvalue = [AWtotalvalue nanmean(plotdata.(Ticks{iticks})(totalstart(istack):totalend(istack)))];
        AWpeakvalue = [AWpeakvalue meanOfdata];
    end
    
    for istack = ANsize+AWsize+1:ANsize+AWsize+Msize
        [latency, ampl, meanOfdata] = iGetPeakData_avrec(plotdata.(Ticks{iticks})(winstart(istack):winend(istack)));
        if isnan(meanOfdata)
            Mpeakamp = [Mpeakamp NaN];
            Mpeaklat = [Mpeaklat NaN];
        else
            Mpeakamp = [Mpeakamp ampl];
            Mpeaklat = [Mpeaklat latency];
        end
        Mtotalvalue = [Mtotalvalue nanmean(plotdata.(Ticks{iticks})(totalstart(istack):totalend(istack)))];
        Mpeakvalue = [Mpeakvalue meanOfdata];
    end
    
    %apply a cutoff threshold
    if sum(isnan(ANpeakamp)) > ceil(ANsize/4) % 25% at least
        ANpeakamp = NaN(1,ANsize);
        ANpeaklate = NaN(1,ANsize);
    end
    if sum(isnan(AWpeakamp)) > ceil(AWsize/4) % 25% at least
        AWpeakamp = NaN(1,AWsize);
        AWpeaklate = NaN(1,AWsize);
    end
    if sum(isnan(Mpeakamp)) > ceil(Msize/4) % 25% at least
        Mpeakamp = NaN(1,Msize);
        Mpeaklate = NaN(1,Msize);
    end
    %CHANGE and USE this to even out groups if there is a different number of animals
%     ANpeakamp = [ANpeakamp NaN];
%     ANpeaklat = [ANpeaklat NaN];
%     ANtotalvalue = [ANtotalvalue NaN];
%     ANpeakvalue = [ANpeakvalue NaN];
    AWpeakamp = [AWpeakamp NaN NaN];
    AWpeaklat = [AWpeaklat NaN NaN];
    AWtotalvalue = [AWtotalvalue NaN NaN];
    AWpeakvalue = [AWpeakvalue NaN NaN];
%     Mpeakamp = [Mpeakamp NaN];
%     Mpeaklat = [Mpeaklat NaN];
%     Mtotalvalue = [Mtotalvalue NaN];
%     Mpeakvalue = [Mpeakvalue NaN];

    %for tuning graph
    plotdata.peakamp = horzcat(plotdata.peakamp, ANpeakamp, AWpeakamp, Mpeakamp);
    plotdata.peaklate = horzcat(plotdata.peaklate, ANpeaklat, AWpeaklat, Mpeaklat);
    plotdata.smallmean = horzcat(plotdata.smallmean, ANpeakvalue, AWpeakvalue, Mpeakvalue);
    plotdata.totalmean = horzcat(plotdata.totalmean, ANtotalvalue, AWtotalvalue, Mtotalvalue);
    plotdata.tgroup = horzcat(plotdata.tgroup, repmat({Group{1}},1,grpsz), repmat({Group{2}},1,grpsz), repmat({Group{3}},1,grpsz)); %manual insertion of group size - make automatic
    plotdata.frqz = horzcat(plotdata.frqz, repmat({grammTicks{iticks}}, 1, grpsz), repmat({grammTicks{iticks}}, 1, grpsz), repmat({grammTicks{iticks}}, 1, grpsz));
    
    disp('********Peak Amplitude********')
    peakamp_ANvsAW = horzcat(ANpeakamp', AWpeakamp');
    peakamp_ANvsM = horzcat(ANpeakamp', Mpeakamp');
    var_peakamp = {'Groups','Peak Amplitude'};
    
    disp('**Anesthetized vs Awake**')
    AnesthetizedvsAwake = teg_repeated_measures_ANOVA(peakamp_ANvsAW, levels, var_peakamp);
    disp('**Anesthetized vs Muscimol**')
    AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(peakamp_ANvsM, levels, var_peakamp);
    
    savefile = [(Ticks{iticks}) ' PeakAmp AvrecStats.mat'];
    save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol')
    
    disp('********Peak Latitude********')
    peaklat_ANvsAW = horzcat(ANpeaklat', AWpeaklat');
    peaklat_ANvsM = horzcat(ANpeaklat', Mpeaklat');
    var_peaklat = {'Groups','Peak Latitude'};
    
    disp('**Anesthetized vs Awake**')
    AnesthetizedvsAwake = teg_repeated_measures_ANOVA(peaklat_ANvsAW, levels, var_peaklat);
    disp('**Anesthetized vs Muscimol**')
    AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(peaklat_ANvsM, levels, var_peaklat);
    
    savefile = [(Ticks{iticks}) ' PeakLat AvrecStats.mat'];
    save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol')
    
    disp('********Total Value********')
    totalvalue_ANvsAW = horzcat(ANtotalvalue', AWtotalvalue');
    totalvalue_ANvsM = horzcat(ANtotalvalue', Mtotalvalue');
    var_totalvalue = {'Groups','Total'};
    
    disp('**Anesthetized vs Awake**')
    AnesthetizedvsAwake = teg_repeated_measures_ANOVA(totalvalue_ANvsAW, levels, var_totalvalue);
    disp('**Anesthetized vs Muscimol**')
    AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(totalvalue_ANvsM, levels, var_totalvalue);
    
    savefile = [(Ticks{iticks}) ' Total AvrecStats.mat'];
    save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol')
    
    disp('********Peak Value********')
    peakvalue_ANvsAW = horzcat(ANpeakvalue', AWpeakvalue');
    peakvalue_ANvsM = horzcat(ANpeakvalue', Mpeakvalue');
    var_peakvalue = {'Groups','Total'};
    
    disp('**Anesthetized vs Awake**')
    AnesthetizedvsAwake = teg_repeated_measures_ANOVA(peakvalue_ANvsAW, levels, var_peakvalue);
    disp('**Anesthetized vs Muscimol**')
    AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(peakvalue_ANvsM, levels, var_peakvalue);
    
    savefile = [(Ticks{iticks}) ' PeakValue AvrecStats.mat'];
    save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol')
end

%% Tuning curve
feature = {'peakamp','peaklate', 'smallmean','totalmean'};
Name = {'Peak Amplitude', 'Peak Latency', 'RMS between 200 and 300 ms', 'RMS of total Avrec'};
cd(home); cd figs; mkdir('AvgAvRec Tuning'); cd('AvgAvRec Tuning')

for i1 = 1:length(feature)
clear g
figure();

g(1,1)=gramm('x',plotdata.frqz,'y',plotdata.(feature{i1}), 'color', plotdata.tgroup); %,'y',cars.Acceleration,'color',cars.Cylinders,'subset',cars.Cylinders~=3 & cars.Cylinders~=5
g(1,1).stat_summary('type','sem','geom','errorbar'); %mean and sem shown
g(1,1).stat_summary('type','sem','geom','point'); %mean and sem shown
g(1,1).stat_summary('type','sem','geom','area'); %mean and sem shown
g(1,1).set_text_options('base_size',12) 
if strcmp(feature{i1},'peaklate')
    g(1,1).set_names('x','Tone','y','ms','color','Group');
%     g(1,1).axe_property('YLim',[230 270]);
elseif strcmp(feature{i1},'peakamp')
    g(1,1).set_names('x','Tone','y','mV*mm²','color','Group');
%     g(1,1).axe_property('YLim',[0 0.0025]);
else
    g(1,1).set_names('x','Tone','y','mV/mm²','color','Group');
    g(1,1).axe_property('YLim',[0.0002 0.001]);
end
g(1,1).set_color_options('map','matlab');
g(1,1).axe_property('XTickLabel',{'-3', '-2', '-1', 'BF', '+1', '+2', '+3'})

g.set_title(Name{i1});
g.draw();
g.export('file_name', Name{i1}, 'file_type','pdf');
g.export('file_name',Name{i1}, 'file_type','png');
close all;

% box plots, single feature
clear g
figure();

g=gramm('x',plotdata.tgroup,'y',plotdata.(feature{i1}), 'color', plotdata.tgroup);
g.facet_grid(plotdata.frqz,[]);
g.stat_boxplot('notch',true);
if strcmp(feature{i1},'peaklate')
    g(1,1).set_names('x','Group','y','Latency (ms)','color','Group');
elseif strcmp(feature{i1},'peakamp')
    g(1,1).set_names('x','Group','y','Amplituded (mV/mm²)','color','Group');
else
    g(1,1).set_names('x','Tone','y','RMS (mV/mm²)','color','Group');

end
g(1,1).set_color_options('map','matlab');

g.set_title([Name{i1} 'boxplot']);
g.draw();
g.export('file_name',[Name{i1} 'boxplot'], 'file_type','pdf');
g.export('file_name',[Name{i1} 'boxplot'], 'file_type','png');
close all;
end

clear g
figure();

g=gramm('x',plotdata.peaklate,'y',plotdata.peakamp, 'color', plotdata.tgroup);
g.facet_grid(plotdata.frqz,[]);
g.geom_point('alpha',0.5);
g.axe_property('YLim',[0 0.015]);
g.set_color_options('map','matlab');

g.set_title('Peak Latency against Amp');
g.draw();
g.export('file_name','Peak Latency against Amp', 'file_type','pdf');
g.export('file_name','Peak Latency against Amp', 'file_type','png');
close all;

save('AvrecPlotData_group.mat','plotdata')