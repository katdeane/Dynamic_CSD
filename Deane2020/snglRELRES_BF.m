%% README

%  This script is SPECIFIC for running the RelRes stats for Katrina's
%  Master's thesis. It may be used as a template but please make a copy
%  elsewhere for modicifications

clear 
cd('D:\MyCode\Dynamic_CSD_Analysis');
warning('OFF');
dbstop if error

home = pwd;
addpath(genpath(home));

%for teg_repeated measures_ANOVA
levels = [2,1];
rowofnans = NaN(1,50); %IMPORTANT: if animal group sizes are equal this isn't needed. Currently these are commented in the code

Ticks = {'a -3','b -2', 'c -1','d BF', 'e +1', 'f +2', 'g +3'};
% Ticknames = {'Best Frequency', 'BF - 1','BF - 2','BF - 3','BF + 1', 'BF + 2', 'BF + 3'};
Ticknames = {'BF - 3','BF - 2','BF - 1','Best Frequency','BF + 1', 'BF + 2', 'BF + 3'};
% forsave = {'BF', 'BF-1','BF-2','BF-3','BF+1', 'BF+2', 'BF+3'};
forsave = {'BF-3', 'BF-2','BF-1','BF','BF+1', 'BF+2', 'BF+3'};
Group = {'Anesthetized','Awake','Muscimol','Anesthetized Chronic'};
% FromBF = [0,-1,-2,-3,+1,+2,+3];
FromBF = [-3,-2,-1,0,+1,+2,+3];

%% Load in the appropriate files
cd DATA;
load('AnesthetizedPre_Data.mat')
Anesthetized = Data; clear Data;
load('Awake10dB_Data.mat')
Awake = Data; clear Data;
load('Muscimol_Data.mat')
Muscimol = Data; clear Data;
% load('ANChronic_Data.mat')
% AnChronic = Data; clear Data;

%this is from the avrec version of this code
cd avrec_compare;
load('AvrecPlotData_single.mat')
AVRECdata = plottune; clear plottune
cd(home);cd DATA

%% Plots
input = dir('*.mat'); %this is currently redundant, I'm using the loop from avgRelRes_BF so I brought over the load-in code for simplicity
entries = length(input); %^
plotdata = struct;
%for full plots
plotdata.BF = []; plotdata.min_one = []; plotdata.min_two = []; plotdata.min_three = [];
plotdata.plus_one = []; plotdata.plus_two = []; plotdata.plus_three = [];
plotdata.time = [];
plotdata.fgroup = [];
%for tuning curves
plottune = struct;
grpsz = 550;
plottune.frqz = [];
plottune.tgroup = [];
plottune.peakamp = [];plottune.peaklate = [];plottune.smallmean = [];plottune.totalmean = [];

for i1 = 1:entries
    % Preliminary functions
    disp(['Analyzing Group: ' (input(i1).name(1:end-9))])
    tic
    load (input(i1).name) %loads data into workspace
    
    cd(home); cd figs; %opens figure folder to store future images
    mkdir(['SnglRelres_' input(i1).name(1:end-9)]); cd(['SnglRelres_' input(i1).name(1:end-9)]);

    names = fieldnames(Data); %list of animal names
    avrec = struct;
    
    time = (1:1:600);
    groupname = (input(i1).name(1:end-9));
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
       relres.(animal).BF =  nanmean(Data.(animal).SingleTrial_RelRes_raw{1,BF},3)';
       
       plotdata.BF = horzcat(plotdata.BF, relres.(animal).BF(1:600));
       plotdata.time = horzcat(plotdata.time, time);
       plotdata.fgroup = horzcat(plotdata.fgroup, grouplist);
       %-1
       try 
           relres.(animal).min_one = nanmean(Data.(animal).SingleTrial_RelRes_raw{1,BF-1},3)';  
       catch
           relres.(animal).min_one = NaN(1,length(relres.(animal).BF));
       end
       plotdata.min_one = horzcat(plotdata.min_one, relres.(animal).min_one(1:600));
       %-2
       try 
           relres.(animal).min_two = nanmean(Data.(animal).SingleTrial_RelRes_raw{1,BF-2},3)';  
       catch
           relres.(animal).min_two = NaN(1,length(relres.(animal).BF));
       end
       plotdata.min_two = horzcat(plotdata.min_two, relres.(animal).min_two(1:600));
       %-3
       try 
           relres.(animal).min_three = nanmean(Data.(animal).SingleTrial_RelRes_raw{1,BF-3},3)';  
       catch
           relres.(animal).min_three = NaN(1,length(relres.(animal).BF));
       end
       plotdata.min_three = horzcat(plotdata.min_three, relres.(animal).min_three(1:600));
       %+1
       try 
           relres.(animal).plus_one = nanmean(Data.(animal).SingleTrial_RelRes_raw{1,BF+1},3)';  
       catch
           relres.(animal).plus_one = NaN(1,length(relres.(animal).BF));
       end
       plotdata.plus_one = horzcat(plotdata.plus_one, relres.(animal).plus_one(1:600));
       %+2
       try 
           relres.(animal).plus_two = nanmean(Data.(animal).SingleTrial_RelRes_raw{1,BF+2},3)';  
       catch
           relres.(animal).plus_two = NaN(1,length(relres.(animal).BF));
       end
       plotdata.plus_two = horzcat(plotdata.plus_two, relres.(animal).plus_two(1:600));
       %+3
       try 
           relres.(animal).plus_three = nanmean(Data.(animal).SingleTrial_RelRes_raw{1,BF+3},3)';  
       catch
           relres.(animal).plus_three = NaN(1,length(relres.(animal).BF));
       end
       plotdata.plus_three = horzcat(plotdata.plus_three, relres.(animal).plus_three(1:600));
   end
   
%    %plots
%    %BF
%     h = figure('Name','sAvgRec - BF');
%     average = [];
%     for i3 = 1:length(names)
%         line = relres.(names{i3}).BF;
%         line = line(1:600); %GXL have 600ms and GKD have 800!
%         plot(line);
%         average = [average line'];
%         hold on
%     end
%     
%     average = nanmean(average,2)';
%     plot(average,'k','LineWidth',2);
%     allnames = vertcat(names,{'Average'});
%     legend(allnames)
%     
%     set(h, 'PaperType', 'A4');
%     set(h, 'PaperOrientation', 'landscape');
%     set(h, 'PaperUnits', 'centimeters');
%     savefig(h, 'sAvgRec - BF','compact')
%     close (h)
%     
%     %BF -1
%     h = figure('Name','sAvgRec - BF -1');
%     average = [];
%     for i3 = 1:length(names)
%         line = relres.(names{i3}).min_one;
%         line = line(1:600); %GXL have 600ms and GKD have 800!
%         plot(line);
%         average = [average line'];
%         hold on
%     end
%     
%     average = nanmean(average,2)';
%     plot(average,'k','LineWidth',2);
%     legend(allnames)
%     
%     set(h, 'PaperType', 'A4');
%     set(h, 'PaperOrientation', 'landscape');
%     set(h, 'PaperUnits', 'centimeters');
%     savefig(h, 'sAvgRec - BF -1','compact')
%     close (h)
%     
%     %BF -2
%     h = figure('Name','sAvgRec - BF -2');
%     average = [];
%     for i3 = 1:length(names)
%         line = relres.(names{i3}).min_two;
%         line = line(1:600); %GXL have 600ms and GKD have 800!
%         plot(line);
%         average = [average line'];
%         hold on
%     end
%     
%     average = nanmean(average,2)';
%     plot(average,'k','LineWidth',2);
%     legend(allnames)
%     
%     set(h, 'PaperType', 'A4');
%     set(h, 'PaperOrientation', 'landscape');
%     set(h, 'PaperUnits', 'centimeters');
%     savefig(h, 'sAvgRec - BF -2','compact')
%     close (h)
%     
%      %BF -3
%     h = figure('Name','sAvgRec - BF -3');
%     average = [];
%     for i3 = 1:length(names)
%         line = relres.(names{i3}).min_three;
%         line = line(1:600); %GXL have 600ms and GKD have 800!
%         plot(line);
%         average = [average line'];
%         hold on
%     end
%     
%     average = nanmean(average,2)';
%     plot(average,'k','LineWidth',2);
%     legend(allnames)
%     
%     set(h, 'PaperType', 'A4');
%     set(h, 'PaperOrientation', 'landscape');
%     set(h, 'PaperUnits', 'centimeters');
%     savefig(h, 'sAvgRec - BF -3','compact')
%     close (h)
%    
%     %BF +1
%     h = figure('Name','sAvgRec - BF +1');
%     average = [];
%     for i3 = 1:length(names)
%         line = relres.(names{i3}).plus_one;
%         line = line(1:600); %GXL have 600ms and GKD have 800!
%         plot(line);
%         average = [average line'];
%         hold on
%     end
%     
%     average = nanmean(average,2)';
%     plot(average,'k','LineWidth',2);
%     legend(allnames)
%     
%     set(h, 'PaperType', 'A4');
%     set(h, 'PaperOrientation', 'landscape');
%     set(h, 'PaperUnits', 'centimeters');
%     savefig(h, 'sAvgRec - BF +1','compact')
%     close (h)
%     
%     %BF +2
%     h = figure('Name','sAvgRec - BF +2');
%     average = [];
%     for i3 = 1:length(names)
%         line = relres.(names{i3}).plus_two;
%         line = line(1:600); %GXL have 600ms and GKD have 800!
%         plot(line);
%         average = [average line'];
%         hold on
%     end
%     
%     average = nanmean(average,2)';
%     plot(average,'k','LineWidth',2);
%     legend(allnames)
%     
%     set(h, 'PaperType', 'A4');
%     set(h, 'PaperOrientation', 'landscape');
%     set(h, 'PaperUnits', 'centimeters');
%     savefig(h, 'sAvgRec - BF +2','compact')
%     close (h)
%     
%      %BF +3
%     h = figure('Name','sAvgRec - BF +3');
%     average = [];
%     for i3 = 1:length(names)
%         line = relres.(names{i3}).plus_three;
%         line = line(1:600); %GXL have 600ms and GKD have 800!
%         plot(line);
%         average = [average line'];
%         hold on
%     end
%     
%     average = nanmean(average,2)';
%     plot(average,'k','LineWidth',2);
%     legend(allnames)
%     
%     set(h, 'PaperType', 'A4');
%     set(h, 'PaperOrientation', 'landscape');
%     set(h, 'PaperUnits', 'centimeters');
%     savefig(h, 'sAvgRec - BF +3','compact')
%     close (h)
%     close all;
end

Ticks = {'BF', 'min_one', 'min_two', 'min_three', 'plus_one', 'plus_two', 'plus_three'};
Ticknames = {'Best Frequency', 'BF - 1','BF - 2','BF - 3','BF + 1', 'BF + 2', 'BF + 3'};
%% Gramm Plots
cd(home); cd figs; %opens figure folder to store future images
mkdir('Gramm Plots SnglRelres'); cd('Gramm Plots SnglRelres');

for iticks = 1:length(Ticks)
    clear g
    g=gramm('x',plotdata.time,'y',plotdata.(Ticks{iticks}),'color',plotdata.fgroup);
    g.stat_summary('type','sem','geom','area'); %mean and sem shown
    g.set_layout_options('Position',[0 0 0.7 0.7],...
            'legend_pos',[0.71 0.66 0.2 0.2],... %We detach the legend from the plot and move it to the top right
            'margin_height',[0.1 0.1],...
            'margin_width',[0.1 0.1],...
            'redraw',false);
    g.set_names('x','Time (ms)','y','µV','color','Group');
    g.axe_property('ylim',[-0.3 0.2],'xlim',[0 600])
    g.set_color_options('map','matlab');
    g.set_title((Ticknames{iticks}));
    g.draw();
    g.export('file_name',['S' (Ticks{iticks})], 'file_type','png');
    g.export('file_name',['S' (Ticks{iticks})], 'file_type','pdf');
    close all;
end

mkdir('Stats RelRes single'); cd('Stats RelRes single')
%% Run Stats

newTicks = {'a -3','b -2', 'c -1','d BF', 'e +1', 'f +2', 'g +3'};
newTicknames = {'BF - 3', 'BF - 2','BF - 1','Best Frequency','BF + 1', 'BF + 2', 'BF + 3'};
avrec_count = (1:2200:15401);
for ifreq = 1:length(Ticknames)
    disp(['***********************STATS FOR ' (newTicknames{ifreq}) '***********************'])
    ANpeakamp = []; ANpeaklat = []; ANtotalvalue = []; ANpeakvalue = [];
    ANnames = fieldnames(Anesthetized);
    for i1 = 1:length(ANnames)
        animal = ANnames{i1};
        %Find the best frequency
        BF = find((Anesthetized.(ANnames{i1}).GS_BF) == (Anesthetized.(ANnames{i1}).Frqz'));
        
        try %cut out the section of the matrix necessary (BF, BF-1, etc)
            relreslist =  Anesthetized.(animal).SingleTrial_RelRes_raw{1,(BF+FromBF(ifreq))};
        catch %produce NAN if there isn't an entry here
            relreslist = NaN(1,length(relreslist));
        end
        peakampme = []; peaklat = []; totalvalue = []; peakvalue = [];
        
        for i2 = 1:50
            if isnan(relreslist)
                peaklat = [peaklat NaN];
                peakampme = [peakampme NaN];
                totalvalue = [totalvalue NaN];
                peakvalue = [peakvalue NaN];
            else
                [latency, ampl, meanOfdata] = iGetPeakData_relres(relreslist(:,:,i2),[200 300]);
                [~, ~, meanOffulldata] = iGetPeakData_relres(relreslist(:,:,i2));
                peaklat = [peaklat latency];
                peakampme = [peakampme ampl];
                peakvalue = [peakvalue meanOfdata];
                totalvalue = [totalvalue meanOffulldata];
            end
        end
        
        %apply a cutoff threshold
        if sum(isnan(peakampme)) > length(peakampme)/4 % 25% at least
            peakampme = NaN(1,50);
            peaklat = NaN(1,50);
        end
        
        ANpeakamp = [ANpeakamp peakampme];
        ANpeaklat = [ANpeaklat peaklat];
        ANtotalvalue = [ANtotalvalue totalvalue];
        ANpeakvalue = [ANpeakvalue peakvalue];
    end
    
%     ANpeakamp = [ANpeakamp rowofnans];
%     ANpeaklat = [ANpeaklat rowofnans];
%     ANtotalvalue = [ANtotalvalue rowofnans];
%     ANpeakvalue = [ANpeakvalue rowofnans];
    
    
    AWpeakamp = []; AWpeaklat = []; AWtotalvalue = []; AWpeakvalue = [];
    AWnames = fieldnames(Awake);
    for i1 = 1:length(AWnames)
        animal = AWnames{i1};
        BF = find((Awake.(AWnames{i1}).GS_BF) == (Awake.(AWnames{i1}).Frqz'));
        try
            relreslist =  Awake.(animal).SingleTrial_RelRes_raw{:,(BF+FromBF(ifreq))};
        catch
            relreslist = NaN(1,length(relreslist));
        end
        peakampme = []; peaklat = []; totalvalue = []; peakvalue = [];
        for i2 = 1:50
            if isnan(relreslist)
                peaklat = [peaklat NaN];
                peakampme = [peakampme NaN];
                totalvalue = [totalvalue NaN];
                peakvalue = [peakvalue NaN];
            else
                try
                    [latency, ampl, meanOfdata] = iGetPeakData_relres(relreslist(:,:,i2),[200 300]);
                    [~, ~, meanOffulldata] = iGetPeakData_relres(relreslist(:,:,i2));
                    peaklat = [peaklat latency];
                    peakampme = [peakampme ampl];
                    peakvalue = [peakvalue meanOfdata];
                    totalvalue = [totalvalue meanOffulldata];
                catch
                    peaklat = [peaklat NaN];
                    peakampme = [peakampme NaN];
                    totalvalue = [totalvalue NaN];
                    peakvalue = [peakvalue NaN];
                end
            end
        end
        
        %apply a cutoff threshold
        if sum(isnan(peakampme)) > length(peakampme)/4 % 25% at least
            peakampme = NaN(1,50);
            peaklat = NaN(1,50);
        end
        
        AWpeakamp = [AWpeakamp peakampme];
        AWpeaklat = [AWpeaklat peaklat];
        AWtotalvalue = [AWtotalvalue totalvalue];
        AWpeakvalue = [AWpeakvalue peakvalue];
    end
    
    AWpeakamp = [AWpeakamp rowofnans rowofnans];
    AWpeaklat = [AWpeaklat rowofnans rowofnans];
    AWtotalvalue = [AWtotalvalue rowofnans rowofnans];
    AWpeakvalue = [AWpeakvalue rowofnans rowofnans];
    
    Mpeakamp = []; Mpeaklat = []; Mtotalvalue = []; Mpeakvalue = [];
    Mnames = fieldnames(Muscimol);
    for i1 = 1:length(Mnames)
        animal = Mnames{i1};
        BF = find((Muscimol.(Mnames{i1}).GS_BF) == (Muscimol.(Mnames{i1}).Frqz'));
        try
            relreslist =  Muscimol.(animal).SingleTrial_RelRes_raw{:,(BF+FromBF(ifreq))};
        catch
            relreslist = NaN(1,length(relreslist));
        end
        peakampme = []; peaklat = []; totalvalue = []; peakvalue = [];
        for i2 = 1:50
            if isnan(relreslist)
                peaklat = [peaklat NaN];
                peakampme = [peakampme NaN];
                totalvalue = [totalvalue NaN];
                peakvalue = [peakvalue NaN];
            else
                [latency, ampl, meanOfdata] = iGetPeakData_relres(relreslist(:,:,i2),[200 300]);
                [~, ~, meanOffulldata] = iGetPeakData_relres(relreslist(:,:,i2));
                peaklat = [peaklat latency];
                peakampme = [peakampme ampl];
                peakvalue = [peakvalue meanOfdata];
                totalvalue = [totalvalue meanOffulldata];
            end
        end
        
        %apply a cutoff threshold
        if sum(isnan(peakampme)) > length(peakampme)/4 % 25% at least
            peakampme = NaN(1,50);
            peaklat = NaN(1,50);
        end
        
        Mpeakamp = [Mpeakamp peakampme];
        Mpeaklat = [Mpeaklat peaklat];
        Mtotalvalue = [Mtotalvalue totalvalue];
        Mpeakvalue = [Mpeakvalue peakvalue];
    end
    
    
    
    %for later tuning curves
    plottune.peakamp = horzcat(plottune.peakamp, ANpeakamp, AWpeakamp, Mpeakamp);
    plottune.peaklate = horzcat(plottune.peaklate, ANpeaklat, AWpeaklat, Mpeaklat);
    plottune.smallmean = horzcat(plottune.smallmean, ANpeakvalue, AWpeakvalue, Mpeakvalue);
    plottune.totalmean = horzcat(plottune.totalmean, ANtotalvalue, AWtotalvalue, Mtotalvalue);
    plottune.tgroup = horzcat(plottune.tgroup, repmat({Group{1}},1,grpsz), repmat({Group{2}},1,grpsz), repmat({Group{3}},1,grpsz));
    plottune.frqz = horzcat(plottune.frqz, repmat({newTicks{ifreq}}, 1, grpsz), repmat({newTicks{ifreq}}, 1, grpsz), repmat({newTicks{ifreq}}, 1, grpsz));
%     
%     % Peak Amp
%     disp('********Peak Amplitude********')
%     peakamp_ANvsAW = horzcat(ANpeakamp', AWpeakamp');
% %     peakamp_AnCvsAW = horzcat(AnCpeakamp', AWpeakamp');
%     peakamp_ANvsM = horzcat(ANpeakamp', Mpeakamp');
%     var_peakamp = {'Groups','Peak Amplitude'};
%     
%     disp('**Anesthetized vs Awake**')
%     AnesthetizedvsAwake = teg_repeated_measures_ANOVA(peakamp_ANvsAW, levels, var_peakamp);
%     disp('**Anesthetized vs Awake COHEN D**')
%     ANvsAWcohenD = iMakeCohensD(ANpeakamp, AWpeakamp)
% %     disp('**Anesthetized Chronic vs Awake**')
% %     AnChronicvsAwake = teg_repeated_measures_ANOVA(peakamp_AnCvsAW, levels, var_peakamp);
% %     disp('**Anesthetized Chronic vs Awake COHEN D**')
% %     AnCvsAWcohenD = iMakeCohensD(AnCpeakamp, AWpeakamp)
%     disp('**Anesthetized vs Muscimol**')
%     AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(peakamp_ANvsM, levels, var_peakamp);
%     disp('**Anesthetized vs Muscimol COHEN D**')
%     ANvsMcohenD = iMakeCohensD(ANpeakamp, Mpeakamp)
%     
%     savefile = [forsave{ifreq} 'PeakAmp RelresStats.mat'];
%     save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol','ANvsAWcohenD','ANvsMcohenD')
%     
%     % Peak Lat
%     disp('********Peak Latency********')
%     peaklat_ANvsAW = horzcat(ANpeaklat', AWpeaklat');
% %     peaklat_AnCvsAW = horzcat(AnCpeaklat', AWpeaklat');
%     peaklat_ANvsM = horzcat(ANpeaklat', Mpeaklat');
%     var_peaklat = {'Groups','Peak Latency'};
%     
%     disp('**Anesthetized vs Awake**')
%     AnesthetizedvsAwake = teg_repeated_measures_ANOVA(peaklat_ANvsAW, levels, var_peaklat);
%     disp('**Anesthetized vs Awake COHEN D**')
%     ANvsAWcohenD = iMakeCohensD(ANpeaklat, AWpeaklat)
% %     disp('**Anesthetized Chronic vs Awake**')
% %     AnChronicvsAwake = teg_repeated_measures_ANOVA(peaklat_AnCvsAW, levels, var_peaklat);
% %     disp('**Anesthetized Chronic vs Awake COHEN D**')
% %     AnCvsAWcohenD = iMakeCohensD(AnCpeaklat, AWpeaklat)
%     disp('**Anesthetized vs Muscimol**')
%     AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(peaklat_ANvsM, levels, var_peaklat);
%     disp('**Anesthetized vs Muscimol COHEN D**')
%     ANvsMcohenD = iMakeCohensD(ANpeaklat, Mpeaklat)
%     
%     savefile = [forsave{ifreq} ' Peaklat RelresStats.mat'];
%     save(savefile, 'AnesthetizedvsAwake',  'AnesthetizedvsMuscimol','ANvsAWcohenD','ANvsMcohenD')
%     
%     % Total Value
%     disp('********Total Value********')
%     totalvalue_ANvsAW = horzcat(ANtotalvalue', AWtotalvalue');
% %     totalvalue_AnCvsAW = horzcat(AnCtotalvalue', AWtotalvalue');
%     totalvalue_ANvsM = horzcat(ANtotalvalue', Mtotalvalue');
%     var_totalvalue = {'Groups','Total'};
%     
%     disp('**Anesthetized vs Awake**')
%     AnesthetizedvsAwake = teg_repeated_measures_ANOVA(totalvalue_ANvsAW, levels, var_totalvalue);
%     disp('**Anesthetized vs Awake COHEN D**')
%     ANvsAWcohenD = iMakeCohensD(ANtotalvalue, AWtotalvalue)
% %     disp('**Anesthetized Chronic vs Awake**')
% %     AnChronicvsAwake = teg_repeated_measures_ANOVA(totalvalue_AnCvsAW, levels, var_totalvalue);
% %     disp('**Anesthetized Chronic vs Awake COHEN D**')
% %     AnCvsAWcohenD = iMakeCohensD(AnCtotalvalue, AWtotalvalue)
%     disp('**Anesthetized vs Muscimol**')
%     AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(totalvalue_ANvsM, levels, var_totalvalue);
%     disp('**Anesthetized vs Muscimol COHEN D**')
%     ANvsMcohenD = iMakeCohensD(ANtotalvalue, Mtotalvalue)
%     
%     savefile = [forsave{ifreq} ' Total RelresStats.mat'];
%     save(savefile, 'AnesthetizedvsAwake',  'AnesthetizedvsMuscimol','ANvsAWcohenD','ANvsMcohenD')
%     
%     % Peak Value
%     disp('********Peak Value********')
%     peakvalue_ANvsAW = horzcat(ANpeakvalue', AWpeakvalue');
% %     peakvalue_AnCvsAW = horzcat(AnCpeakvalue', AWpeakvalue');
%     peakvalue_ANvsM = horzcat(ANpeakvalue', Mpeakvalue');
%     var_peakvalue = {'Groups','Total'};
%     
%     disp('**Anesthetized vs Awake**')
%     AnesthetizedvsAwake = teg_repeated_measures_ANOVA(peakvalue_ANvsAW, levels, var_peakvalue);
%     disp('**Anesthetized vs Awake COHEN D**')
%     ANvsAWcohenD = iMakeCohensD(ANpeakvalue, AWpeakvalue)
% %     disp('**Anesthetized Chronic vs Awake**')
% %     AnChronicvsAwake = teg_repeated_measures_ANOVA(peakvalue_AnCvsAW, levels, var_peakvalue);
% %     disp('**Anesthetized Chronic vs Awake COHEN D**')
% %     AnCvsAWcohenD = iMakeCohensD(AnCpeakvalue, AWpeakvalue)
%     disp('**Anesthetized vs Muscimol**')
%     AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(peakvalue_ANvsM, levels, var_peakvalue);
%     disp('**Anesthetized vs Muscimol COHEN D**')
%     ANvsMcohenD = iMakeCohensD(ANpeakvalue, Mpeakvalue)
%     
%     savefile = [forsave{ifreq} ' PeakValue RelresStats.mat'];
%     save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol','ANvsAWcohenD','ANvsMcohenD')
end


plottune.peakamp(isnan(AVRECdata.peakamp)) = NaN;
plottune.peaklate(isnan(AVRECdata.peaklate)) = NaN;
plottune.smallmean(isnan(AVRECdata.smallmean)) = NaN;
%% Tuning curve
feature = {'peakamp','peaklate', 'smallmean','totalmean'};
Name = {'Peak Amplitude', 'Peak Latency', 'RMS between 200 and 300 ms', 'RMS of total Relres'};
cd(home); cd figs; mkdir('SnglRelRes Tuning'); cd('SnglRelRes Tuning')

for i1 = 1:length(feature)
% tuning curves
clear g
figure();

g(1,1)=gramm('x',plottune.frqz,'y',plottune.(feature{i1}), 'color', plottune.tgroup); %,'y',cars.Acceleration,'color',cars.Cylinders,'subset',cars.Cylinders~=3 & cars.Cylinders~=5
g(1,1).stat_summary('type','sem','geom','errorbar'); %mean and sem shown
g(1,1).stat_summary('type','sem','geom','point'); %mean and sem shown
g(1,1).stat_summary('type','sem','geom','area'); %mean and sem shown
g(1,1).set_text_options('base_size',12) 
if strcmp(feature{i1},'peaklate')
    g(1,1).set_names('x','Tone','y','ms','color','Group');
    g(1,1).axe_property('YLim',[220 270]);
elseif strcmp(feature{i1},'peakamp')
    g(1,1).set_names('x','Tone','y','mV*mm²','color','Group');
    g(1,1).axe_property('YLim',[-0.55 -0.1]);
else
    g(1,1).set_names('x','Tone','y','ratio','color','Group');
    g(1,1).axe_property('YLim',[0 0.4]);
end
g(1,1).set_color_options('map','matlab');
g(1,1).axe_property('XTickLabel',{'-3', '-2', '-1', 'BF', '+1', '+2', '+3'})

g.set_title(Name{i1});
g.draw();
g.export('file_name',Name{i1}, 'file_type','pdf');
g.export('file_name',Name{i1}, 'file_type','png');
close all;

% box plots, single feature
clear g
figure();

g=gramm('x',plottune.tgroup,'y',plottune.(feature{i1}), 'color', plottune.tgroup);
g.facet_grid(plottune.frqz,[]);
g.stat_boxplot('notch',true);
if strcmp(feature{i1},'peaklate')
    g(1,1).set_names('x','Group','y','Latency (ms)','color','Group');
elseif strcmp(feature{i1},'peakamp')
    g(1,1).set_names('x','Group','y','Percentage (mV/mm²)','color','Group');
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

%all frequencies amp vs latency
clear g
figure();

g=gramm('x',plottune.peaklate,'y',plottune.peakamp, 'color', plottune.tgroup);
g.facet_grid(plottune.frqz,[]);
g.set_point_options('base_size',2);
g.geom_point();
g.set_color_options('map','matlab');

g.set_title('Peak Latency against Amp');
g.draw();
g.export('file_name','Peak Latency against Amp', 'file_type','pdf');
g.export('file_name','Peak Latency against Amp', 'file_type','png');
close all;



BFlogical = strcmp('d BF',plottune.frqz);
min2logical = strcmp('b -2',plottune.frqz);

BF_avrec_late = AVRECdata.peaklate(BFlogical);
BF_avrec_amp = AVRECdata.peakamp(BFlogical);
BF_avrec_group = AVRECdata.tgroup(BFlogical);

BF_relres_late = plottune.peaklate(BFlogical);
BF_relres_amp = plottune.peakamp(BFlogical);
BF_relres_group = plottune.tgroup(BFlogical);

m2_avrec_late = AVRECdata.peaklate(min2logical);
m2_avrec_amp = AVRECdata.peakamp(min2logical);
m2_avrec_group = AVRECdata.tgroup(min2logical);

m2_relres_late = plottune.peaklate(min2logical);
m2_relres_amp = plottune.peakamp(min2logical);
m2_relres_group = plottune.tgroup(min2logical);



%% Relres peak latency against peak amplitute
clear g
figure();

g=gramm('x',BF_relres_late,'y',BF_relres_amp, 'color', BF_relres_group);
g.geom_point();
g.axe_property('YLim',[-1 0.5],'XLim',[200 300]);
g.set_color_options('map','matlab');

g.set_names('x','Peak Latency','y','Peak Amplitude','color','Group');
g.set_title('RelRes Peak Latency against Amp BF');
g.draw();
g.export('file_name','Peak Latency against Amp BF', 'file_type','pdf');
g.export('file_name','Peak Latency against Amp BF', 'file_type','png');
close all;

clear g
figure();

g=gramm('x',m2_relres_late,'y',m2_relres_amp, 'color', m2_avrec_group);
g.geom_point();
g.set_color_options('map','matlab');

g.set_names('x','Peak Latency','y','Peak Amplitude','color','Group');
g.set_title('RelRes Peak Latency against Amp BF-2');
g.draw();
g.export('file_name','Peak Latency against Amp BF-2', 'file_type','pdf');
g.export('file_name','Peak Latency against Amp BF-2', 'file_type','png');
close all;

%% avrec peak latency against relres peak latency
clear g
figure();

g=gramm('x',BF_avrec_late,'y',BF_relres_late, 'color', BF_avrec_group);
g.geom_point(); 
g.geom_abline();
% g.set_point_options('base_size',4);
g.stat_ellipse('type','95percentile','geom','area','patch_opts',{'FaceAlpha',0.1,'LineWidth',2});
g.axe_property('YLim',[200 300],'XLim',[200 300]);
g.set_color_options('map','matlab');
g.set_names('x','Avrec Peak Latency','y','Relres Peak Latency','color','Group');

g.set_title('Avrec against Relres Peak Latency BF');
g.draw();
g.export('file_name','Avrec against Relres Peak Latency BF', 'file_type','pdf');
g.export('file_name','Avrec against Relres Peak Latency BF', 'file_type','png');
close all;

%heatmaps for funsies
clear g
figure();
g=gramm('x',BF_avrec_late,'y',BF_relres_late);
g.facet_grid([],BF_avrec_group);
g.geom_abline();
g.stat_bin2d('nbins',[20 20],'geom','image');
g.axe_property('YLim',[200 300],'XLim',[200 300]);
g.set_names('x','Avrec Peak Latency','y','Relres Peak Latency','color','Count');

g.set_title('Heatmap Avrec against Relres Peak Latency BF');
g.draw();
g.export('file_name','Heatmap Avrec against Relres Peak Latency BF', 'file_type','pdf');
g.export('file_name','Heatmap Avrec against Relres Peak Latency BF', 'file_type','png');
close all;

% BF -2
clear g
figure();

g=gramm('x',m2_avrec_late,'y',m2_relres_late, 'color', m2_avrec_group);
g.geom_point(); 
g.geom_abline();
% g.set_point_options('base_size',4);
g.stat_ellipse('type','95percentile','geom','area','patch_opts',{'FaceAlpha',0.1,'LineWidth',2});
g.axe_property('YLim',[200 300],'XLim',[200 300]);
g.set_color_options('map','matlab');
g.set_names('x','Avrec Peak Latency','y','Relres Peak Latency','color','Group');

g.set_title('Avrec against Relres Peak Latency BF-2');
g.draw();
g.export('file_name','Avrec against Relres Peak Latency BF-2', 'file_type','pdf');
g.export('file_name','Avrec against Relres Peak Latency BF-2', 'file_type','png');
close all;

%heatmaps for funsies
Ketamine_avlate = m2_avrec_late(strcmp('Anesthetized',m2_avrec_group));
Ketamine_rellate = m2_relres_late(strcmp('Anesthetized',m2_relres_group));

Awake_avlate = m2_avrec_late(strcmp('Awake',m2_avrec_group));
Awake_rellate = m2_relres_late(strcmp('Awake',m2_relres_group));

clear g
figure();
g(1,1)=gramm('x',Ketamine_avlate,'y',Ketamine_rellate);
% g.facet_grid([],m2_avrec_group(strcmp('Anesthetized',m2_avrec_group)));
g(1,1).stat_bin2d('nbins',[20 20],'geom','image');
g.geom_abline();
g(1,1).axe_property('YLim',[200 300],'XLim',[200 300]);
g(1,1).set_names('x','Avrec Peak Latency','y','Relres Peak Latency','color','Count');

g(1,2)=gramm('x',Awake_avlate,'y',Awake_rellate);
% g.facet_grid([],m2_avrec_group(strcmp('Anesthetized',m2_avrec_group)));
g(1,2).stat_bin2d('nbins',[20 20],'geom','image');
g.geom_abline();
g(1,2).axe_property('YLim',[200 300],'XLim',[200 300]);
g(1,2).set_names('x','Avrec Peak Latency','y','Relres Peak Latency','color','Count');

g.set_title('Heatmap Avrec against Relres Peak Latency BF-2');
g.draw();
g.export('file_name','Heatmap Avrec against Relres Peak Latency BF-2', 'file_type','pdf');
g.export('file_name','Heatmap Avrec against Relres Peak Latency BF-2', 'file_type','png');
close all;



%% avrec peak amplitude against relres peak amplitude
clear g
figure();

g=gramm('x',BF_avrec_amp,'y',BF_relres_amp, 'color', BF_avrec_group);
g.geom_point(); 
g.stat_ellipse('type','95percentile','geom','area','patch_opts',{'FaceAlpha',0.1,'LineWidth',2});
g.set_color_options('map','matlab');
g.set_names('x','Avrec Peak Amplitude','y','Relres Peak Amplitude','color','Group');

g.set_title('Avrec against Relres Peak Amplitude BF');
g.draw();
g.export('file_name','Avrec against Relres Peak Amplitude BF', 'file_type','pdf');
g.export('file_name','Avrec against Relres Peak Amplitude BF', 'file_type','png');
close all;

%heatmaps for funsies
clear g
figure();
g=gramm('x',BF_avrec_amp,'y',BF_relres_amp);
g.facet_grid([],BF_avrec_group);
g.stat_bin2d('nbins',[20 20],'geom','image');
g.set_names('x','Avrec Peak Amplitude','y','Relres Peak Amplitude','color','Count');

g.set_title('Heatmap Avrec against Relres Peak Amplitude BF');
g.draw();
g.export('file_name','Heatmap Avrec against Relres Peak Amplitude BF', 'file_type','pdf');
g.export('file_name','Heatmap Avrec against Relres Peak Amplitude BF', 'file_type','png');
close all;

% BF -2
clear g
figure();

g=gramm('x',m2_avrec_amp,'y',m2_relres_amp, 'color', m2_avrec_group);
g.geom_point(); 
g.stat_ellipse('type','95percentile','geom','area','patch_opts',{'FaceAlpha',0.1,'LineWidth',2});
% g.axe_property('YLim',[-0.06 0.02],'XLim',[-0.002 0.006]);
g.set_color_options('map','matlab');
g.set_names('x','Avrec Peak Amplitude','y','Relres Peak Amplitude','color','Group');

g.set_title('Avrec against Relres Peak Amplitude BF-2');
g.draw();
g.export('file_name','Avrec against Relres Peak Amplitude BF-2', 'file_type','pdf');
g.export('file_name','Avrec against Relres Peak Amplitude BF-2', 'file_type','png');
close all;


Ketamine_avamp = m2_avrec_amp(strcmp('Anesthetized',m2_avrec_group));
Ketamine_relamp = m2_relres_amp(strcmp('Anesthetized',m2_relres_group));

Awake_avamp = m2_avrec_amp(strcmp('Awake',m2_avrec_group));
Awake_relamp = m2_relres_amp(strcmp('Awake',m2_relres_group));

%heatmaps for funsies
clear g
figure();
g(1,1)=gramm('x',Ketamine_avamp,'y',Ketamine_relamp);
g(1,1).stat_bin2d('nbins',[20 20],'geom','image');
% g(1,1).axe_property('YLim',[-0.07 0.02],'XLim',[-0.002 0.007]);
g(1,1).set_names('x','Avrec Peak Amplitude','y','Relres Peak Amplitude','color','Count');
g(1,1).set_title('Ketamine')

g(1,2)=gramm('x',Awake_avamp,'y',Awake_relamp);
g(1,2).stat_bin2d('nbins',[20 20],'geom','image');
% g(1,2).axe_property('YLim',[-0.07 0.02],'XLim',[-0.002 0.007]);
g(1,2).set_names('x','Avrec Peak Amplitude','y','Relres Peak Amplitude','color','Count');
g(1,2).set_title('Awake')

g.set_title('Heatmap Avrec against Relres Peak Amplitude BF-2');
g.draw();
g.export('file_name','Heatmap Avrec against Relres Peak Amplitude BF-2', 'file_type','pdf');
g.export('file_name','Heatmap Avrec against Relres Peak Amplitude BF-2', 'file_type','png');
close all;

%% Save the data out
cd(home);cd DATA; cd avrec_compare;
save('RelresPlotData_single.mat','plottune')
