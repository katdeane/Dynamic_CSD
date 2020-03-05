%% README

%  This script is SPECIFIC for running the AVREC figures for Katrina's
%  master's thesis and ketamine publicatoin

clear 
cd('D:\MyCode\Dynamic_CSD_Analysis');
warning('OFF');
dbstop if error

home = pwd;
addpath(genpath(home));

%for teg_repeated measures_ANOVA
rowofnans = NaN(1,50); %IMPORTANT: if animal group sizes are equal this isn't needed. Currently these are commented in the code

% forsave = {'BF-2','BF-1','BF','BF+1', 'BF+2'};
Group = {'Anesthetized','Awake','Muscimol','ANChronic'};
FromBF = [-2,-1,0,+1,+2];

%% Load in the appropriate files
cd DATA;
load('AnesthetizedPre_Data.mat')
Anesthetized = Data; clear Data;
load('ANChronic_Data.mat')
ANChronic = Data; clear Data;
load('Awake10dB_Data.mat')
Awake = Data; clear Data;
load('Muscimol_Data.mat')
Muscimol = Data; clear Data;

GroupNames = {'Anesthetized','ANChronic','Awake','Muscimol'};


%% Plots
% The first plots are of the individual animals' AVRECs and a black line
% averageing the group. For each frequency, using Matlab plotting.
% The second plots are of the groups' AVRECs for each frequecny using
% grammplots

%for full plots
plotdata = struct;
plotdata.BF = []; plotdata.min_one = []; plotdata.min_two = []; plotdata.min_three = [];
plotdata.plus_one = []; plotdata.plus_two = []; plotdata.plus_three = [];
plotdata.time = [];
plotdata.fgroup = [];
%for tuning curves
plottune = struct;
grpsz = 550;
plottune.frqz = [];
plottune.tgroup = []; plottune.animal = [];
plottune.peakamp = [];plottune.peaklate = [];plottune.smallmean = [];plottune.totalmean = [];
plottune.peakamp_norm = [];plottune.smallmean_norm = [];

for i1 = 1:length(GroupNames)

    disp(['Analyzing Group: ' GroupNames{i1}])
    tic
    
    clear Data
    if strcmp('Anesthetized',GroupNames{i1})
        Data = Anesthetized;
    elseif strcmp('Awake',GroupNames{i1})
        Data = Awake;
    elseif strcmp('ANChronic',GroupNames{i1})
        Data = ANChronic;
    elseif strcmp('Muscimol',GroupNames{i1})
        Data = Muscimol;
    end
    
    names = fieldnames(Data); %list of animal names
    avrec = struct;
    
    time = (1:1:600);
    groupname = GroupNames{i1};
    grouplist = repmat({groupname}, 1, 600);
   
   for i2 = 1:length(names)
       animal = names{i2};
       %BF
       BF = find((Data.(names{i2}).GS_BF) == (Data.(names{i2}).Frqz'));
       avrec.(animal).BF =  nanmean(Data.(animal).SingleTrial_AVREC_raw{1,BF},3)';
       
       plotdata.BF = horzcat(plotdata.BF, avrec.(animal).BF(1:600));
       plotdata.time = horzcat(plotdata.time, time);
       plotdata.fgroup = horzcat(plotdata.fgroup, grouplist);
       %-1
       try 
           avrec.(animal).min_one = nanmean(Data.(animal).SingleTrial_AVREC_raw{1,BF-1},3)';  
       catch
           avrec.(animal).min_one = NaN(1,length(avrec.(animal).BF));
       end
       plotdata.min_one = horzcat(plotdata.min_one, avrec.(animal).min_one(1:600));
       %-2
       try 
           avrec.(animal).min_two = nanmean(Data.(animal).SingleTrial_AVREC_raw{1,BF-2},3)';  
       catch
           avrec.(animal).min_two = NaN(1,length(avrec.(animal).BF));
       end
       plotdata.min_two = horzcat(plotdata.min_two, avrec.(animal).min_two(1:600));
       %-3
       try 
           avrec.(animal).min_three = nanmean(Data.(animal).SingleTrial_AVREC_raw{1,BF-3},3)';  
       catch
           avrec.(animal).min_three = NaN(1,length(avrec.(animal).BF));
       end
       plotdata.min_three = horzcat(plotdata.min_three, avrec.(animal).min_three(1:600));
       %+1
       try 
           avrec.(animal).plus_one = nanmean(Data.(animal).SingleTrial_AVREC_raw{1,BF+1},3)';  
       catch
           avrec.(animal).plus_one = NaN(1,length(avrec.(animal).BF));
       end
       plotdata.plus_one = horzcat(plotdata.plus_one, avrec.(animal).plus_one(1:600));
       %+2
       try 
           avrec.(animal).plus_two = nanmean(Data.(animal).SingleTrial_AVREC_raw{1,BF+2},3)';  
       catch
           avrec.(animal).plus_two = NaN(1,length(avrec.(animal).BF));
       end
       plotdata.plus_two = horzcat(plotdata.plus_two, avrec.(animal).plus_two(1:600));
       %+3
       try 
           avrec.(animal).plus_three = nanmean(Data.(animal).SingleTrial_AVREC_raw{1,BF+3},3)';  
       catch
           avrec.(animal).plus_three = NaN(1,length(avrec.(animal).BF));
       end
       plotdata.plus_three = horzcat(plotdata.plus_three, avrec.(animal).plus_three(1:600));
   end
   
   % MATLAB PLOTS currently silence to avoid unnecessary replotting
   cd(home); cd figs; %opens figure folder to store future images
   mkdir(['SnglAvrec_' GroupNames{i1}]); cd(['SnglAvrec_' GroupNames{i1}]);
   PlotAnimalAvrec(avrec, names)
    
end

% Gramm Plots currently silenced to avoid unnecessary replotting
Ticks = {'BF', 'min_one', 'min_two', 'min_three', 'plus_one', 'plus_two', 'plus_three'};
Ticknames = {'Best Frequency', 'BF - 1','BF - 2','BF - 3','BF + 1', 'BF + 2', 'BF + 3'};
cd(home); cd figs; %opens figure folder to store future images
mkdir('Gramm Plots SnglAvrec'); cd('Gramm Plots SnglAvrec');

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
    g.axe_property('ylim',[0 0.002],'xlim',[0 600])
    g.set_color_options('map','matlab');
    g.set_title((Ticknames{iticks}));
    g.draw();
    g.export('file_name',['S' (Ticks{iticks})], 'file_type','png');
    g.export('file_name',['S' (Ticks{iticks})], 'file_type','pdf');
    close all;
end

%% Tuning Plots

newTicks = {'b -2', 'c -1','d BF', 'e +1', 'f +2'};
newTicknames = {'BF - 2','BF - 1','Best Frequency','BF + 1', 'BF + 2'};

for ifreq = 1:length(newTicknames)
    disp(['***********************STATS FOR ' (newTicknames{ifreq}) '***********************'])
    
    % Anesthetized Acute
    Kpeakamp = []; Kpeaklat = []; Ktotalvalue = []; Kpeakvalue = [];
    Kpeakamp_norm = []; Kpeakvalue_norm = [];
    Knames = fieldnames(Anesthetized);
    for i1 = 1:length(Knames)
        animal = Knames{i1};
        %Find the best frequency
        BF = find((Anesthetized.(Knames{i1}).GS_BF) == (Anesthetized.(Knames{i1}).Frqz'));
        avreclist_BF =  Anesthetized.(animal).SingleTrial_AVREC_raw{1,BF};
        
        try %cut out the section of the matrix necessary (BF, BF-1, etc)
            avreclist =  Anesthetized.(animal).SingleTrial_AVREC_raw{1,(BF+FromBF(ifreq))};
        catch %produce NAN if there isn't an entry here
            avreclist = NaN(1,length(avreclist));
        end
        
        peakampme = []; peaklat = []; totalvalue = []; peakvalue = [];
        for i2 = 1:50
            if isnan(avreclist)
                peaklat = [peaklat NaN];
                peakampme = [peakampme NaN];
                totalvalue = [totalvalue NaN];
                peakvalue = [peakvalue NaN];
            else
                [latency, ampl, meanOfdata] = iGetPeakData_avrec(avreclist(:,:,i2),[200 300]);
                [~, ~, meanOffulldata] = iGetPeakData_avrec(avreclist(:,:,i2));
                peaklat = [peaklat latency];
                peakampme = [peakampme ampl];
                peakvalue = [peakvalue meanOfdata];
                totalvalue = [totalvalue meanOffulldata];
            end
        end
        
        peakampBF = []; peakvalueBF = [];
        for i2 = 1:50
            if isnan(avreclist_BF)
                peakampBF = [peakampBF NaN];
                peakvalueBF = [peakvalueBF NaN];
            else
                [latency, ampl, meanOfdata] = iGetPeakData_avrec(avreclist_BF(:,:,i2),[200 300]);
                peakampBF = [peakampBF ampl];
                peakvalueBF = [peakvalueBF meanOfdata];
            end
        end
        
        peakamp_norm = peakampme./peakampBF;
        peakvalue_norm = peakvalue./peakvalueBF;
        
        %apply a cutoff threshold
        if sum(isnan(peakampme)) > length(peakampme)*.75 % 25% 
            peakampme = NaN(1,50);
            peaklat = NaN(1,50);
            peakvalue = NaN(1,50);
        end
        
        Kpeakamp = [Kpeakamp peakampme];
        Kpeaklat = [Kpeaklat peaklat];
        Ktotalvalue = [Ktotalvalue totalvalue];
        Kpeakvalue = [Kpeakvalue peakvalue];
        
        Kpeakamp_norm = [Kpeakamp_norm peakamp_norm];
        Kpeakvalue_norm = [Kpeakvalue_norm peakvalue_norm];
    end
    
%     Kpeakamp = [Kpeakamp rowofnans];
%     Kpeaklat = [Kpeaklat rowofnans];
%     Ktotalvalue = [Ktotalvalue rowofnans];
%     Kpeakvalue = [Kpeakvalue rowofnans];
    
    % Anesthetized Chronic
    Npeakamp = []; Npeaklat = []; Ntotalvalue = []; Npeakvalue = [];
    Npeakamp_norm = []; Npeakvalue_norm = [];
    Nnames = fieldnames(ANChronic);
    for i1 = 1:length(Nnames)
        animal = Nnames{i1};
        BF = find((ANChronic.(Nnames{i1}).GS_BF) == (ANChronic.(Nnames{i1}).Frqz'));
        avreclist_BF =  ANChronic.(animal).SingleTrial_AVREC_raw{1,BF};
        try
            avreclist =  ANChronic.(animal).SingleTrial_AVREC_raw{:,(BF+FromBF(ifreq))};
        catch
            avreclist = NaN(1,length(avreclist));
        end
        peakampme = []; peaklat = []; totalvalue = []; peakvalue = [];
        for i2 = 1:50
            if sum(sum(isnan(avreclist))) > 1 || size(avreclist,3)<i2
                peaklat = [peaklat NaN];
                peakampme = [peakampme NaN];
                totalvalue = [totalvalue NaN];
                peakvalue = [peakvalue NaN];
            else
                [latency, ampl, meanOfdata] = iGetPeakData_avrec(avreclist(:,:,i2),[200 300]);
                [~, ~, meanOffulldata] = iGetPeakData_avrec(avreclist(:,:,i2));
                peaklat = [peaklat latency];
                peakampme = [peakampme ampl];
                peakvalue = [peakvalue meanOfdata];
                totalvalue = [totalvalue meanOffulldata];
            end
        end
        
        peakampBF = []; peakvalueBF = [];
        for i2 = 1:50
            if sum(sum(isnan(avreclist))) > 1  || size(avreclist_BF,3)<i2
                peakampBF = [peakampBF NaN];
                peakvalueBF = [peakvalueBF NaN];
            else
                [latency, ampl, meanOfdata] = iGetPeakData_avrec(avreclist_BF(:,:,i2),[200 300]);
                peakampBF = [peakampBF ampl];
                peakvalueBF = [peakvalueBF meanOfdata];
            end
        end
        
        peakamp_norm = peakampme./peakampBF;
        peakvalue_norm = peakvalue./peakvalueBF;
        
        %apply a cutoff threshold
        if sum(isnan(peakampme)) > length(peakampme)*.75 % 25% 
            peakampme = NaN(1,50);
            peaklat = NaN(1,50);
            peakvalue = NaN(1,50);
        end
        
        Npeakamp = [Npeakamp peakampme];
        Npeaklat = [Npeaklat peaklat];
        Ntotalvalue = [Ntotalvalue totalvalue];
        Npeakvalue = [Npeakvalue peakvalue];
        
        Npeakamp_norm = [Npeakamp_norm peakamp_norm];
        Npeakvalue_norm = [Npeakvalue_norm peakvalue_norm];
    end
    
    Npeakamp = [Npeakamp rowofnans rowofnans];
    Npeaklat = [Npeaklat rowofnans rowofnans];
    Ntotalvalue = [Ntotalvalue rowofnans rowofnans];
    Npeakvalue = [Npeakvalue rowofnans rowofnans];
    
    Npeakamp_norm = [Npeakamp_norm rowofnans rowofnans];
    Npeakvalue_norm = [Npeakvalue_norm rowofnans rowofnans];
    
    % Awake Chronic
    Apeakamp = []; Apeaklat = []; Atotalvalue = []; Apeakvalue = [];
    Apeakamp_norm = []; Apeakvalue_norm = [];
    Anames = fieldnames(Awake);
    for i1 = 1:length(Anames)
        animal = Anames{i1};
        BF = find((Awake.(Anames{i1}).GS_BF) == (Awake.(Anames{i1}).Frqz'));
        avreclist_BF =  Awake.(animal).SingleTrial_AVREC_raw{1,BF};
        try
            avreclist =  Awake.(animal).SingleTrial_AVREC_raw{:,(BF+FromBF(ifreq))};
        catch
            avreclist = NaN(1,length(avreclist));
        end
        
        peakampme = []; peaklat = []; totalvalue = []; peakvalue = [];
        for i2 = 1:50
            if sum(sum(isnan(avreclist))) > 1 || size(avreclist,3)<i2
                peaklat = [peaklat NaN];
                peakampme = [peakampme NaN];
                totalvalue = [totalvalue NaN];
                peakvalue = [peakvalue NaN];
            else
                [latency, ampl, meanOfdata] = iGetPeakData_avrec(avreclist(:,:,i2),[200 300]);
                [~, ~, meanOffulldata] = iGetPeakData_avrec(avreclist(:,:,i2));
                peaklat = [peaklat latency];
                peakampme = [peakampme ampl];
                peakvalue = [peakvalue meanOfdata];
                totalvalue = [totalvalue meanOffulldata];
            end
        end
        
        peakampBF = []; peakvalueBF = [];
        for i2 = 1:50
            if sum(sum(isnan(avreclist))) > 1 || size(avreclist_BF,3)<i2
                peakampBF = [peakampBF NaN];
                peakvalueBF = [peakvalueBF NaN];
            else
                [latency, ampl, meanOfdata] = iGetPeakData_avrec(avreclist_BF(:,:,i2),[200 300]);
                peakampBF = [peakampBF ampl];
                peakvalueBF = [peakvalueBF meanOfdata];
            end
        end
        
        peakamp_norm = peakampme./peakampBF;
        peakvalue_norm = peakvalue./peakvalueBF;
        
        %apply a cutoff threshold
        if sum(isnan(peakampme)) > length(peakampme)*.75 % 25% 
            peakampme = NaN(1,50);
            peaklat = NaN(1,50);
            peakvalue = NaN(1,50);
        end
        
        Apeakamp = [Apeakamp peakampme];
        Apeaklat = [Apeaklat peaklat];
        Atotalvalue = [Atotalvalue totalvalue];
        Apeakvalue = [Apeakvalue peakvalue];
        
        Apeakamp_norm = [Apeakamp_norm peakamp_norm];
        Apeakvalue_norm = [Apeakvalue_norm peakvalue_norm];
    end
    
    Apeakamp = [Apeakamp rowofnans rowofnans];
    Apeaklat = [Apeaklat rowofnans rowofnans];
    Atotalvalue = [Atotalvalue rowofnans rowofnans];
    Apeakvalue = [Apeakvalue rowofnans rowofnans];
    
    Apeakamp_norm = [Apeakamp_norm rowofnans rowofnans];
    Apeakvalue_norm = [Apeakvalue_norm rowofnans rowofnans];
    
    % Muscimol Acute
    Mpeakamp = []; Mpeaklat = []; Mtotalvalue = []; Mpeakvalue = [];
    Mpeakamp_norm = []; Mpeakvalue_norm = [];
    Mnames = fieldnames(Muscimol);
    for i1 = 1:length(Mnames)
        animal = Mnames{i1};
        BF = find((Muscimol.(Mnames{i1}).GS_BF) == (Muscimol.(Mnames{i1}).Frqz'));
        avreclist_BF =  Muscimol.(animal).SingleTrial_AVREC_raw{1,BF};
        try
            avreclist =  Muscimol.(animal).SingleTrial_AVREC_raw{:,(BF+FromBF(ifreq))};
        catch
            avreclist = NaN(1,length(avreclist));
        end
        peakampme = []; peaklat = []; totalvalue = []; peakvalue = [];
        for i2 = 1:50
            if isnan(avreclist)
                peaklat = [peaklat NaN];
                peakampme = [peakampme NaN];
                totalvalue = [totalvalue NaN];
                peakvalue = [peakvalue NaN];
            else
                [latency, ampl, meanOfdata] = iGetPeakData_avrec(avreclist(:,:,i2),[200 300]);
                [~, ~, meanOffulldata] = iGetPeakData_avrec(avreclist(:,:,i2));
                peaklat = [peaklat latency];
                peakampme = [peakampme ampl];
                peakvalue = [peakvalue meanOfdata];
                totalvalue = [totalvalue meanOffulldata];
            end
        end
        
        peakampBF = []; peakvalueBF = [];
        for i2 = 1:50
            if isnan(avreclist_BF)
                peakampBF = [peakampBF NaN];
                peakvalueBF = [peakvalueBF NaN];
            else
                [latency, ampl, meanOfdata] = iGetPeakData_avrec(avreclist_BF(:,:,i2),[200 300]);
                peakampBF = [peakampBF ampl];
                peakvalueBF = [peakvalueBF meanOfdata];
            end
        end
        
        peakamp_norm = peakampme./peakampBF;
        peakvalue_norm = peakvalue./peakvalueBF;
        
        %apply a cutoff threshold
        if sum(isnan(peakampme)) > length(peakampme)*.75 % 25% 
            peakampme = NaN(1,50);
            peaklat = NaN(1,50);
            peakvalue = NaN(1,50);
        end
            
        Mpeakamp = [Mpeakamp peakampme];
        Mpeaklat = [Mpeaklat peaklat];
        Mtotalvalue = [Mtotalvalue totalvalue];
        Mpeakvalue = [Mpeakvalue peakvalue];
        
        Mpeakamp_norm = [Mpeakamp_norm peakamp_norm];
        Mpeakvalue_norm = [Mpeakvalue_norm peakvalue_norm];
    end
    
%     Mpeakamp = [Mpeakamp rowofnans];
%     Mpeaklat = [Mpeaklat rowofnans];
%     Mtotalvalue = [Mtotalvalue rowofnans];
%     Mpeakvalue = [Mpeakvalue rowofnans];
     
    Kanimals = horzcat(repmat({Knames{1}},1,50),repmat({Knames{2}},1,50),...
        repmat({Knames{3}},1,50),repmat({Knames{4}},1,50),repmat({Knames{5}},1,50),...
        repmat({Knames{6}},1,50),repmat({Knames{7}},1,50),repmat({Knames{8}},1,50),...
        repmat({Knames{9}},1,50),repmat({Knames{10}},1,50),repmat({Knames{11}},1,50));
    
    Nanimals = horzcat(repmat({Nnames{1}},1,50),repmat({Nnames{2}},1,50),...
        repmat({Nnames{3}},1,50),repmat({Nnames{4}},1,50),repmat({Nnames{5}},1,50),...
        repmat({Nnames{6}},1,50),repmat({Nnames{7}},1,50),repmat({Nnames{8}},1,50),...
        repmat({Nnames{9}},1,50),repmat({'NaN'},1,50),repmat({'NaN'},1,50));
    
    Aanimals = horzcat(repmat({Anames{1}},1,50),repmat({Anames{2}},1,50),...
        repmat({Anames{3}},1,50),repmat({Anames{4}},1,50),repmat({Anames{5}},1,50),...
        repmat({Anames{6}},1,50),repmat({Anames{7}},1,50),repmat({Anames{8}},1,50),...
        repmat({Anames{9}},1,50),repmat({'NaN'},1,50),repmat({'NaN'},1,50));
    
    Manimals = horzcat(repmat({Mnames{1}},1,50),repmat({Mnames{2}},1,50),...
        repmat({Mnames{3}},1,50),repmat({Mnames{4}},1,50),repmat({Mnames{5}},1,50),...
        repmat({Mnames{6}},1,50),repmat({Mnames{7}},1,50),repmat({Mnames{8}},1,50),...
        repmat({Mnames{9}},1,50),repmat({Mnames{10}},1,50),repmat({Mnames{11}},1,50));
    
    %for later tuning curves
%     T = table([Kpeakamp';Apeakamp';Mpeakamp'; repmat({newTicks{ifreq}}, 1, grpsz)';repmat({Group{1}},1,grpsz)']);%,'RowNames',{'PeakAmp','PeakLat','Freq','Group'}
    plottune.peakamp = horzcat(plottune.peakamp, Kpeakamp, Npeakamp, Apeakamp, Mpeakamp);
    plottune.peaklate = horzcat(plottune.peaklate, Kpeaklat, Npeaklat, Apeaklat, Mpeaklat);
    plottune.smallmean = horzcat(plottune.smallmean, Kpeakvalue, Npeakvalue, Apeakvalue, Mpeakvalue);
    plottune.totalmean = horzcat(plottune.totalmean, Ktotalvalue, Ntotalvalue, Atotalvalue, Mtotalvalue);
    
    plottune.peakamp_norm = horzcat(plottune.peakamp_norm, Kpeakamp_norm, Npeakamp_norm, Apeakamp_norm, Mpeakamp_norm);
    plottune.smallmean_norm = horzcat(plottune.smallmean_norm, Kpeakvalue_norm, Npeakvalue_norm, Apeakvalue_norm, Mpeakvalue_norm);
    
    plottune.tgroup = horzcat(plottune.tgroup, repmat({Group{1}},1,grpsz), repmat({Group{4}},1,grpsz), ... 
                                               repmat({Group{2}},1,grpsz), repmat({Group{3}},1,grpsz));
    plottune.animal = horzcat(plottune.animal,Kanimals, Nanimals, Aanimals, Manimals);
    plottune.frqz = horzcat(plottune.frqz, repmat({newTicks{ifreq}}, 1, grpsz), repmat({newTicks{ifreq}}, 1, grpsz), ...
                                           repmat({newTicks{ifreq}}, 1, grpsz), repmat({newTicks{ifreq}}, 1, grpsz));
                                       
end


feature = {'peakamp','peaklate', 'smallmean','totalmean','peakamp_norm','smallmean_norm'};
Name = {'Peak Amplitude', 'Peak Latency', 'RMS between 200 and 300 ms', ...
    'RMS of total Avrec','Normalized Peak Amplitude','Normalized Peak RMS'};


cd(home); cd figs; mkdir('SnglAvRec Tuning'); cd('SnglAvRec Tuning')
PlotAvrecTuning(feature, Name, plottune, home)


plottune.frqzMir = plottune.frqz;
plottune.frqzMir= replace(plottune.frqzMir,"b -2","c +-2");
plottune.frqzMir= replace(plottune.frqzMir,"f +2","c +-2");
plottune.frqzMir= replace(plottune.frqzMir,"c -1","b +-1");
plottune.frqzMir= replace(plottune.frqzMir,"e +1","b +-1");
plottune.frqzMir= replace(plottune.frqzMir,"d BF","a BF");

feature = {'peakamp_norm','smallmean_norm'};
Name = {'Normalized Mirror Peak Amplitude','Normalized Mirror Peak RMS'};

PlotAvrecTuning_Norm(feature, Name, plottune)







