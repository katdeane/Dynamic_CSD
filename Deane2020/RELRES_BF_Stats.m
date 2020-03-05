%% README

%  This script is SPECIFIC for running the RELRES stats for Katrina's
%  Master's thesis. It may be used as a template but please make a copy
%  elsewhere for modicifications

clear all
cd('C:\Users\Katrina\Documents\CortXplorers\MyCode\Dynamic CSD');
warning('OFF');
dbstop if error

home = pwd;
addpath(genpath(home));

%for teg_repeated measures_ANOVA
levels = [2,1];
rowofnans = NaN(1,50);

Ticks = {'a -3','b -2', 'c -1','d BF', 'e +1', 'f +2', 'g +3'};
Ticknames = {'Best Frequency', 'BF - 1','BF - 2','BF - 3','BF + 1', 'BF + 2', 'BF + 3'};
Group = {'Anesthetized','Awake','Muscimol'};
plotdata = struct;
plotdata.group = [];
plotdata.frqz = [];
plotdata.peakamp = [];plotdata.peaklate = [];plotdata.smallmean = [];plotdata.totalmean = [];


%% Load in the appropriate files
cd DATA;
load('AnesthetizedPre_6B_Data.mat')
Anesthetized = Data; clear Data;
load('Awake10dB_6B_Data.mat')
Awake = Data; clear Data;
load('Muscimol_6B_Data.mat')
Muscimol = Data; clear Data;

mkdir('Stats Relres single'); cd('Stats Relres single')

%% BF
disp('***********************STATS FOR BF***********************')
ANpeakamp = []; ANpeaklat = []; ANtotalvalue = []; ANpeakvalue = [];
ANnames = fieldnames(Anesthetized);
for i1 = 1:length(ANnames);
    animal = ANnames{i1};
    BF = find((Anesthetized.(ANnames{i1}).GS_BF) == (Anesthetized.(ANnames{i1}).Frqz'));
    
    relreslist =  Anesthetized.(animal).SingleTrial_RELRES_raw{1,BF};
    peakampme = []; peaklat = []; totalvalue = []; peakvalue = [];
    for i2 = 1:50
        [latency, ampl, meanOfdata] = iGetPeakData_relres(relreslist(:,:,i2),[200 300]);
        [~, ~, meanOffulldata] = iGetPeakData_relres(relreslist(:,:,i2));
        peaklat = [peaklat latency];
        peakampme = [peakampme ampl];
        peakvalue = [peakvalue meanOfdata];
        totalvalue = [totalvalue meanOffulldata];
    end
    
    ANpeakamp = [ANpeakamp peakampme]; 
    ANpeaklat = [ANpeaklat peaklat];
    ANtotalvalue = [ANtotalvalue totalvalue];
    ANpeakvalue = [ANpeakvalue peakvalue];
end

ANpeakamp = [ANpeakamp rowofnans]; 
ANpeaklat = [ANpeaklat rowofnans];
ANtotalvalue = [ANtotalvalue rowofnans];
ANpeakvalue = [ANpeakvalue rowofnans];


AWpeakamp = []; AWpeaklat = []; AWtotalvalue = []; AWpeakvalue = [];
AWnames = fieldnames(Awake);
for i1 = 1:length(AWnames);
    animal = AWnames{i1};
    BF = find((Awake.(AWnames{i1}).GS_BF) == (Awake.(AWnames{i1}).Frqz'));
    
    relreslist =  Awake.(animal).SingleTrial_RELRES_raw{1,BF};
    peakampme = []; peaklat = []; totalvalue = []; peakvalue = [];
    for i2 = 1:50
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
    
    AWpeakamp = [AWpeakamp peakampme];
    AWpeaklat = [AWpeaklat peaklat];
    AWtotalvalue = [AWtotalvalue totalvalue];
    AWpeakvalue = [AWpeakvalue peakvalue];
end

Mpeakamp = []; Mpeaklat = []; Mtotalvalue = []; Mpeakvalue = [];
Mnames = fieldnames(Muscimol);
for i1 = 1:length(Mnames);
    animal = Mnames{i1};
    BF = find((Muscimol.(Mnames{i1}).GS_BF) == (Muscimol.(Mnames{i1}).Frqz'));
    
    relreslist =  Muscimol.(animal).SingleTrial_RELRES_raw{1,BF};
    peakampme = []; peaklat = []; totalvalue = []; peakvalue = [];
    for i2 = 1:50
        [latency, ampl, meanOfdata] = iGetPeakData_relres(relreslist(:,:,i2),[200 300]);
        [~, ~, meanOffulldata] = iGetPeakData_relres(relreslist(:,:,i2));
        peaklat = [peaklat latency];
        peakampme = [peakampme ampl];
        peakvalue = [peakvalue meanOfdata];
        totalvalue = [totalvalue meanOffulldata];
    end
    
    Mpeakamp = [Mpeakamp peakampme];
    Mpeaklat = [Mpeaklat peaklat];
    Mtotalvalue = [Mtotalvalue totalvalue];
    Mpeakvalue = [Mpeakvalue peakvalue];
end

Mpeakamp = [Mpeakamp rowofnans];
Mpeaklat = [Mpeaklat rowofnans];
Mtotalvalue = [Mtotalvalue rowofnans];
Mpeakvalue = [Mpeakvalue rowofnans];

%for later tuning curves
plotdata.peakamp = horzcat(plotdata.peakamp, ANpeakamp, AWpeakamp, Mpeakamp); 
plotdata.peaklate = horzcat(plotdata.peaklate, ANpeaklat, AWpeaklat, Mpeaklat);
plotdata.smallmean = horzcat(plotdata.smallmean, ANpeakvalue, AWpeakvalue, Mpeakvalue);
plotdata.totalmean = horzcat(plotdata.totalmean, ANtotalvalue, AWtotalvalue, Mtotalvalue);
plotdata.group = horzcat(plotdata.group, repmat({Group{1}},1,450), repmat({Group{2}},1,450), repmat({Group{3}},1,450));
plotdata.frqz = horzcat(plotdata.frqz, repmat({Ticks{4}}, 1, 450), repmat({Ticks{4}}, 1, 450), repmat({Ticks{4}}, 1, 450));

% Peak Amp
disp('********Peak Amplitude********')
peakamp_ANvsAW = horzcat(ANpeakamp', AWpeakamp');
peakamp_ANvsM = horzcat(ANpeakamp', Mpeakamp');
var_peakamp = {'Groups','Peak Amplitude'};

disp('**Anesthetized vs Awake**')
AnesthetizedvsAwake = teg_repeated_measures_ANOVA(peakamp_ANvsAW, levels, var_peakamp);
disp('**Anesthetized vs Awake COHEN D**')
ANvsAWcohenD = iMakeCohensD(ANpeakamp, AWpeakamp)
disp('**Anesthetized vs Muscimol**')
AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(peakamp_ANvsM, levels, var_peakamp);
disp('**Anesthetized vs Muscimol COHEN D**')
ANvsMcohenD = iMakeCohensD(ANpeakamp, Mpeakamp)

savefile = 'BF PeakAmp RelresStats.mat';
save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol','ANvsAWcohenD','ANvsMcohenD')

% Peak Lat
disp('********Peak Latency********')
peaklat_ANvsAW = horzcat(ANpeaklat', AWpeaklat');
peaklat_ANvsM = horzcat(ANpeaklat', Mpeaklat');
var_peaklat = {'Groups','Peak Latency'};

disp('**Anesthetized vs Awake**')
AnesthetizedvsAwake = teg_repeated_measures_ANOVA(peaklat_ANvsAW, levels, var_peaklat);
disp('**Anesthetized vs Awake COHEN D**')
ANvsAWcohenD = iMakeCohensD(ANpeaklat, AWpeaklat)
disp('**Anesthetized vs Muscimol**')
AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(peaklat_ANvsM, levels, var_peaklat);
disp('**Anesthetized vs Muscimol COHEN D**')
ANvsMcohenD = iMakeCohensD(ANpeaklat, Mpeaklat)

savefile = 'BF Peaklat RelresStats.mat';
save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol','ANvsAWcohenD','ANvsMcohenD')

% Total Value
disp('********Total Value********')
totalvalue_ANvsAW = horzcat(ANtotalvalue', AWtotalvalue');
totalvalue_ANvsM = horzcat(ANtotalvalue', Mtotalvalue');
var_totalvalue = {'Groups','Total'};

disp('**Anesthetized vs Awake**')
AnesthetizedvsAwake = teg_repeated_measures_ANOVA(totalvalue_ANvsAW, levels, var_totalvalue);
disp('**Anesthetized vs Awake COHEN D**')
ANvsAWcohenD = iMakeCohensD(ANtotalvalue, AWtotalvalue)
disp('**Anesthetized vs Muscimol**')
AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(totalvalue_ANvsM, levels, var_totalvalue);
disp('**Anesthetized vs Muscimol COHEN D**')
ANvsMcohenD = iMakeCohensD(ANtotalvalue, Mtotalvalue)

savefile = 'BF Total RelresStats.mat';
save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol','ANvsAWcohenD','ANvsMcohenD')

% Peak Value
disp('********Peak Value********')
peakvalue_ANvsAW = horzcat(ANpeakvalue', AWpeakvalue');
peakvalue_ANvsM = horzcat(ANpeakvalue', Mpeakvalue');
var_peakvalue = {'Groups','Total'};

disp('**Anesthetized vs Awake**')
AnesthetizedvsAwake = teg_repeated_measures_ANOVA(peakvalue_ANvsAW, levels, var_peakvalue);
disp('**Anesthetized vs Awake COHEN D**')
ANvsAWcohenD = iMakeCohensD(ANpeakvalue, AWpeakvalue)
disp('**Anesthetized vs Muscimol**')
AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(peakvalue_ANvsM, levels, var_peakvalue);
disp('**Anesthetized vs Muscimol COHEN D**')
ANvsMcohenD = iMakeCohensD(ANpeakvalue, Mpeakvalue)

savefile = 'BF PeakValue RelresStats.mat';
save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol','ANvsAWcohenD','ANvsMcohenD')

%% BF - 1
disp('***********************STATS FOR BF - 1***********************')
ANpeakamp = []; ANpeaklat = []; ANtotalvalue = []; ANpeakvalue = [];
ANnames = fieldnames(Anesthetized);
for i1 = 1:length(ANnames);
    animal = ANnames{i1};
    BF = find((Anesthetized.(ANnames{i1}).GS_BF) == (Anesthetized.(ANnames{i1}).Frqz'));
    
    try
        relreslist = Anesthetized.(animal).SingleTrial_RELRES_raw{1,BF-1};
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
    
    ANpeakamp = [ANpeakamp peakampme];
    ANpeaklat = [ANpeaklat peaklat];
    ANtotalvalue = [ANtotalvalue totalvalue];
    ANpeakvalue = [ANpeakvalue peakvalue];
end

ANpeakamp = [ANpeakamp rowofnans];
ANpeaklat = [ANpeaklat rowofnans];
ANtotalvalue = [ANtotalvalue rowofnans];
ANpeakvalue = [ANpeakvalue rowofnans];


AWpeakamp = []; AWpeaklat = []; AWtotalvalue = []; AWpeakvalue = [];
AWnames = fieldnames(Awake);
for i1 = 1:length(AWnames);
    animal = AWnames{i1};
    BF = find((Awake.(AWnames{i1}).GS_BF) == (Awake.(AWnames{i1}).Frqz'));
    
    try
        relreslist = Awake.(animal).SingleTrial_RELRES_raw{1,(BF-1)};
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
    
    AWpeakamp = [AWpeakamp peakampme];
    AWpeaklat = [AWpeaklat peaklat];
    AWtotalvalue = [AWtotalvalue totalvalue];
    AWpeakvalue = [AWpeakvalue peakvalue];
end

Mpeakamp = []; Mpeaklat = []; Mtotalvalue = []; Mpeakvalue = [];
Mnames = fieldnames(Muscimol);
for i1 = 1:length(Mnames);
    animal = Mnames{i1};
    BF = find((Muscimol.(Mnames{i1}).GS_BF) == (Muscimol.(Mnames{i1}).Frqz'));
    
    try
        relreslist = Muscimol.(animal).SingleTrial_RELRES_raw{1,(BF-1)};
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
    
    Mpeakamp = [Mpeakamp peakampme];
    Mpeaklat = [Mpeaklat peaklat];
    Mtotalvalue = [Mtotalvalue totalvalue];
    Mpeakvalue = [Mpeakvalue peakvalue];
end

Mpeakamp = [Mpeakamp rowofnans];
Mpeaklat = [Mpeaklat rowofnans];
Mtotalvalue = [Mtotalvalue rowofnans];
Mpeakvalue = [Mpeakvalue rowofnans];

%for later tuning curves
plotdata.peakamp = horzcat(plotdata.peakamp, ANpeakamp, AWpeakamp, Mpeakamp); 
plotdata.peaklate = horzcat(plotdata.peaklate, ANpeaklat, AWpeaklat, Mpeaklat);
plotdata.smallmean = horzcat(plotdata.smallmean, ANpeakvalue, AWpeakvalue, Mpeakvalue);
plotdata.totalmean = horzcat(plotdata.totalmean, ANtotalvalue, AWtotalvalue, Mtotalvalue);
plotdata.group = horzcat(plotdata.group, repmat({Group{1}},1,450), repmat({Group{2}},1,450), repmat({Group{3}},1,450));
plotdata.frqz = horzcat(plotdata.frqz, repmat({Ticks{3}}, 1, 450), repmat({Ticks{3}}, 1, 450), repmat({Ticks{3}}, 1, 450));

% Peak Amp
disp('********Peak Amplitude********')
peakamp_ANvsAW = horzcat(ANpeakamp', AWpeakamp');
peakamp_ANvsM = horzcat(ANpeakamp', Mpeakamp');
var_peakamp = {'Groups','Peak Amplitude'};

disp('**Anesthetized vs Awake**')
AnesthetizedvsAwake = teg_repeated_measures_ANOVA(peakamp_ANvsAW, levels, var_peakamp);
disp('**Anesthetized vs Awake COHEN D**')
ANvsAWcohenD = iMakeCohensD(ANpeakamp, AWpeakamp)
disp('**Anesthetized vs Muscimol**')
AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(peakamp_ANvsM, levels, var_peakamp);
disp('**Anesthetized vs Muscimol COHEN D**')
ANvsMcohenD = iMakeCohensD(ANpeakamp, Mpeakamp)

savefile = 'BF-1 PeakAmp RelresStats.mat';
save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol','ANvsAWcohenD','ANvsMcohenD')

% Peak Lat
disp('********Peak Latency********')
peaklat_ANvsAW = horzcat(ANpeaklat', AWpeaklat');
peaklat_ANvsM = horzcat(ANpeaklat', Mpeaklat');
var_peaklat = {'Groups','Peak Latency'};

disp('**Anesthetized vs Awake**')
AnesthetizedvsAwake = teg_repeated_measures_ANOVA(peaklat_ANvsAW, levels, var_peaklat);
disp('**Anesthetized vs Awake COHEN D**')
ANvsAWcohenD = iMakeCohensD(ANpeaklat, AWpeaklat)
disp('**Anesthetized vs Muscimol**')
AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(peaklat_ANvsM, levels, var_peaklat);
disp('**Anesthetized vs Muscimol COHEN D**')
ANvsMcohenD = iMakeCohensD(ANpeaklat, Mpeaklat)

savefile = 'BF-1 Peaklat RelresStats.mat';
save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol','ANvsAWcohenD','ANvsMcohenD')

% Total Value
disp('********Total Value********')
totalvalue_ANvsAW = horzcat(ANtotalvalue', AWtotalvalue');
totalvalue_ANvsM = horzcat(ANtotalvalue', Mtotalvalue');
var_totalvalue = {'Groups','Total'};

disp('**Anesthetized vs Awake**')
AnesthetizedvsAwake = teg_repeated_measures_ANOVA(totalvalue_ANvsAW, levels, var_totalvalue);
disp('**Anesthetized vs Awake COHEN D**')
ANvsAWcohenD = iMakeCohensD(ANtotalvalue, AWtotalvalue)
disp('**Anesthetized vs Muscimol**')
AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(totalvalue_ANvsM, levels, var_totalvalue);
disp('**Anesthetized vs Muscimol COHEN D**')
ANvsMcohenD = iMakeCohensD(ANtotalvalue, Mtotalvalue)

savefile = 'BF-1 Total RelresStats.mat';
save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol','ANvsAWcohenD','ANvsMcohenD')

% Peak Value
disp('********Peak Value********')
peakvalue_ANvsAW = horzcat(ANpeakvalue', AWpeakvalue');
peakvalue_ANvsM = horzcat(ANpeakvalue', Mpeakvalue');
var_peakvalue = {'Groups','Total'};

disp('**Anesthetized vs Awake**')
AnesthetizedvsAwake = teg_repeated_measures_ANOVA(peakvalue_ANvsAW, levels, var_peakvalue);
disp('**Anesthetized vs Awake COHEN D**')
ANvsAWcohenD = iMakeCohensD(ANpeakvalue, AWpeakvalue)
disp('**Anesthetized vs Muscimol**')
AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(peakvalue_ANvsM, levels, var_peakvalue);
disp('**Anesthetized vs Muscimol COHEN D**')
ANvsMcohenD = iMakeCohensD(ANpeakvalue, Mpeakvalue)

savefile = 'BF-1 PeakValue RelresStats.mat';
save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol','ANvsAWcohenD','ANvsMcohenD')

%% BF - 2
disp('***********************STATS FOR BF - 2***********************')
ANpeakamp = []; ANpeaklat = []; ANtotalvalue = []; ANpeakvalue = [];
ANnames = fieldnames(Anesthetized);
for i1 = 1:length(ANnames);
    animal = ANnames{i1};
    BF = find((Anesthetized.(ANnames{i1}).GS_BF) == (Anesthetized.(ANnames{i1}).Frqz'));
    
    try
        relreslist = Anesthetized.(animal).SingleTrial_RELRES_raw{1,(BF-2)};
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
    
    ANpeakamp = [ANpeakamp peakampme];
    ANpeaklat = [ANpeaklat peaklat];
    ANtotalvalue = [ANtotalvalue totalvalue];
    ANpeakvalue = [ANpeakvalue peakvalue];
end

ANpeakamp = [ANpeakamp rowofnans];
ANpeaklat = [ANpeaklat rowofnans];
ANtotalvalue = [ANtotalvalue rowofnans];
ANpeakvalue = [ANpeakvalue rowofnans];


AWpeakamp = []; AWpeaklat = []; AWtotalvalue = []; AWpeakvalue = [];
AWnames = fieldnames(Awake);
for i1 = 1:length(AWnames);
    animal = AWnames{i1};
    BF = find((Awake.(AWnames{i1}).GS_BF) == (Awake.(AWnames{i1}).Frqz'));
    
    try
        relreslist = Awake.(animal).SingleTrial_RELRES_raw{1,(BF-2)};
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
    
    AWpeakamp = [AWpeakamp peakampme];
    AWpeaklat = [AWpeaklat peaklat];
    AWtotalvalue = [AWtotalvalue totalvalue];
    AWpeakvalue = [AWpeakvalue peakvalue];
end

Mpeakamp = []; Mpeaklat = []; Mtotalvalue = []; Mpeakvalue = [];
Mnames = fieldnames(Muscimol);
for i1 = 1:length(Mnames);
    animal = Mnames{i1};
    BF = find((Muscimol.(Mnames{i1}).GS_BF) == (Muscimol.(Mnames{i1}).Frqz'));
    
    try
        relreslist = Muscimol.(animal).SingleTrial_RELRES_raw{1,(BF-2)};
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
    
    Mpeakamp = [Mpeakamp peakampme];
    Mpeaklat = [Mpeaklat peaklat];
    Mtotalvalue = [Mtotalvalue totalvalue];
    Mpeakvalue = [Mpeakvalue peakvalue];
end

Mpeakamp = [Mpeakamp rowofnans];
Mpeaklat = [Mpeaklat rowofnans];
Mtotalvalue = [Mtotalvalue rowofnans];
Mpeakvalue = [Mpeakvalue rowofnans];

%for later tuning curves
plotdata.peakamp = horzcat(plotdata.peakamp, ANpeakamp, AWpeakamp, Mpeakamp); 
plotdata.peaklate = horzcat(plotdata.peaklate, ANpeaklat, AWpeaklat, Mpeaklat);
plotdata.smallmean = horzcat(plotdata.smallmean, ANpeakvalue, AWpeakvalue, Mpeakvalue);
plotdata.totalmean = horzcat(plotdata.totalmean, ANtotalvalue, AWtotalvalue, Mtotalvalue);
plotdata.group = horzcat(plotdata.group, repmat({Group{1}},1,450), repmat({Group{2}},1,450), repmat({Group{3}},1,450));
plotdata.frqz = horzcat(plotdata.frqz, repmat({Ticks{2}}, 1, 450), repmat({Ticks{2}}, 1, 450), repmat({Ticks{2}}, 1, 450));

% Peak Amp
disp('********Peak Amplitude********')
peakamp_ANvsAW = horzcat(ANpeakamp', AWpeakamp');
peakamp_ANvsM = horzcat(ANpeakamp', Mpeakamp');
var_peakamp = {'Groups','Peak Amplitude'};

disp('**Anesthetized vs Awake**')
AnesthetizedvsAwake = teg_repeated_measures_ANOVA(peakamp_ANvsAW, levels, var_peakamp);
disp('**Anesthetized vs Awake COHEN D**')
ANvsAWcohenD = iMakeCohensD(ANpeakamp, AWpeakamp)
disp('**Anesthetized vs Muscimol**')
AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(peakamp_ANvsM, levels, var_peakamp);
disp('**Anesthetized vs Muscimol COHEN D**')
ANvsMcohenD = iMakeCohensD(ANpeakamp, Mpeakamp)

savefile = 'BF-2 PeakAmp RelresStats.mat';
save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol','ANvsAWcohenD','ANvsMcohenD')

% Peak Lat
disp('********Peak Latency********')
peaklat_ANvsAW = horzcat(ANpeaklat', AWpeaklat');
peaklat_ANvsM = horzcat(ANpeaklat', Mpeaklat');
var_peaklat = {'Groups','Peak Latency'};

disp('**Anesthetized vs Awake**')
AnesthetizedvsAwake = teg_repeated_measures_ANOVA(peaklat_ANvsAW, levels, var_peaklat);
disp('**Anesthetized vs Awake COHEN D**')
ANvsAWcohenD = iMakeCohensD(ANpeaklat, AWpeaklat)
disp('**Anesthetized vs Muscimol**')
AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(peaklat_ANvsM, levels, var_peaklat);
disp('**Anesthetized vs Muscimol COHEN D**')
ANvsMcohenD = iMakeCohensD(ANpeaklat, Mpeaklat)

savefile = 'BF-2 Peaklat RelresStats.mat';
save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol','ANvsAWcohenD','ANvsMcohenD')

% Total Value
disp('********Total Value********')
totalvalue_ANvsAW = horzcat(ANtotalvalue', AWtotalvalue');
totalvalue_ANvsM = horzcat(ANtotalvalue', Mtotalvalue');
var_totalvalue = {'Groups','Total'};

disp('**Anesthetized vs Awake**')
AnesthetizedvsAwake = teg_repeated_measures_ANOVA(totalvalue_ANvsAW, levels, var_totalvalue);
disp('**Anesthetized vs Awake COHEN D**')
ANvsAWcohenD = iMakeCohensD(ANtotalvalue, AWtotalvalue)
disp('**Anesthetized vs Muscimol**')
AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(totalvalue_ANvsM, levels, var_totalvalue);
disp('**Anesthetized vs Muscimol COHEN D**')
ANvsMcohenD = iMakeCohensD(ANtotalvalue, Mtotalvalue)

savefile = 'BF-2 Total RelresStats.mat';
save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol','ANvsAWcohenD','ANvsMcohenD')

% Peak Value
disp('********Peak Value********')
peakvalue_ANvsAW = horzcat(ANpeakvalue', AWpeakvalue');
peakvalue_ANvsM = horzcat(ANpeakvalue', Mpeakvalue');
var_peakvalue = {'Groups','Total'};

disp('**Anesthetized vs Awake**')
AnesthetizedvsAwake = teg_repeated_measures_ANOVA(peakvalue_ANvsAW, levels, var_peakvalue);
disp('**Anesthetized vs Awake COHEN D**')
ANvsAWcohenD = iMakeCohensD(ANpeakvalue, AWpeakvalue)
disp('**Anesthetized vs Muscimol**')
AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(peakvalue_ANvsM, levels, var_peakvalue);
disp('**Anesthetized vs Muscimol COHEN D**')
ANvsMcohenD = iMakeCohensD(ANpeakvalue, Mpeakvalue)

savefile = 'BF-2 PeakValue RelresStats.mat';
save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol','ANvsAWcohenD','ANvsMcohenD')

%% BF - 3
disp('***********************STATS FOR BF - 3***********************')
ANpeakamp = []; ANpeaklat = []; ANtotalvalue = []; ANpeakvalue = [];
ANnames = fieldnames(Anesthetized);
for i1 = 1:length(ANnames);
    animal = ANnames{i1};
    BF = find((Anesthetized.(ANnames{i1}).GS_BF) == (Anesthetized.(ANnames{i1}).Frqz'));
    
    try
        relreslist = Anesthetized.(animal).SingleTrial_RELRES_raw{1,(BF-3)};
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
    
    ANpeakamp = [ANpeakamp peakampme];
    ANpeaklat = [ANpeaklat peaklat];
    ANtotalvalue = [ANtotalvalue totalvalue];
    ANpeakvalue = [ANpeakvalue peakvalue];
end

ANpeakamp = [ANpeakamp rowofnans];
ANpeaklat = [ANpeaklat rowofnans];
ANtotalvalue = [ANtotalvalue rowofnans];
ANpeakvalue = [ANpeakvalue rowofnans];


AWpeakamp = []; AWpeaklat = []; AWtotalvalue = []; AWpeakvalue = [];
AWnames = fieldnames(Awake);
for i1 = 1:length(AWnames);
    animal = AWnames{i1};
    BF = find((Awake.(AWnames{i1}).GS_BF) == (Awake.(AWnames{i1}).Frqz'));
    
    try
        relreslist = Awake.(animal).SingleTrial_RELRES_raw{1,(BF-3)};
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
    
    AWpeakamp = [AWpeakamp peakampme];
    AWpeaklat = [AWpeaklat peaklat];
    AWtotalvalue = [AWtotalvalue totalvalue];
    AWpeakvalue = [AWpeakvalue peakvalue];
end

Mpeakamp = []; Mpeaklat = []; Mtotalvalue = []; Mpeakvalue = [];
Mnames = fieldnames(Muscimol);
for i1 = 1:length(Mnames);
    animal = Mnames{i1};
    BF = find((Muscimol.(Mnames{i1}).GS_BF) == (Muscimol.(Mnames{i1}).Frqz'));
    
    try
        relreslist = Muscimol.(animal).SingleTrial_RELRES_raw{1,(BF-3)};
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
    
    Mpeakamp = [Mpeakamp peakampme];
    Mpeaklat = [Mpeaklat peaklat];
    Mtotalvalue = [Mtotalvalue totalvalue];
    Mpeakvalue = [Mpeakvalue peakvalue];
end

Mpeakamp = [Mpeakamp rowofnans];
Mpeaklat = [Mpeaklat rowofnans];
Mtotalvalue = [Mtotalvalue rowofnans];
Mpeakvalue = [Mpeakvalue rowofnans];

%for later tuning curves
plotdata.peakamp = horzcat(plotdata.peakamp, ANpeakamp, AWpeakamp, Mpeakamp); 
plotdata.peaklate = horzcat(plotdata.peaklate, ANpeaklat, AWpeaklat, Mpeaklat);
plotdata.smallmean = horzcat(plotdata.smallmean, ANpeakvalue, AWpeakvalue, Mpeakvalue);
plotdata.totalmean = horzcat(plotdata.totalmean, ANtotalvalue, AWtotalvalue, Mtotalvalue);
plotdata.group = horzcat(plotdata.group, repmat({Group{1}},1,450), repmat({Group{2}},1,450), repmat({Group{3}},1,450));
plotdata.frqz = horzcat(plotdata.frqz, repmat({Ticks{1}}, 1, 450), repmat({Ticks{1}}, 1, 450), repmat({Ticks{1}}, 1, 450));

% Peak Amp
disp('********Peak Amplitude********')
peakamp_ANvsAW = horzcat(ANpeakamp', AWpeakamp');
peakamp_ANvsM = horzcat(ANpeakamp', Mpeakamp');
var_peakamp = {'Groups','Peak Amplitude'};

disp('**Anesthetized vs Awake**')
AnesthetizedvsAwake = teg_repeated_measures_ANOVA(peakamp_ANvsAW, levels, var_peakamp);
disp('**Anesthetized vs Awake COHEN D**')
ANvsAWcohenD = iMakeCohensD(ANpeakamp, AWpeakamp)
disp('**Anesthetized vs Muscimol**')
AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(peakamp_ANvsM, levels, var_peakamp);
disp('**Anesthetized vs Muscimol COHEN D**')
ANvsMcohenD = iMakeCohensD(ANpeakamp, Mpeakamp)

savefile = 'BF-3 PeakAmp RelresStats.mat';
save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol','ANvsAWcohenD','ANvsMcohenD')

% Peak Lat
disp('********Peak Latency********')
peaklat_ANvsAW = horzcat(ANpeaklat', AWpeaklat');
peaklat_ANvsM = horzcat(ANpeaklat', Mpeaklat');
var_peaklat = {'Groups','Peak Latency'};

disp('**Anesthetized vs Awake**')
AnesthetizedvsAwake = teg_repeated_measures_ANOVA(peaklat_ANvsAW, levels, var_peaklat);
disp('**Anesthetized vs Awake COHEN D**')
ANvsAWcohenD = iMakeCohensD(ANpeaklat, AWpeaklat)
disp('**Anesthetized vs Muscimol**')
AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(peaklat_ANvsM, levels, var_peaklat);
disp('**Anesthetized vs Muscimol COHEN D**')
ANvsMcohenD = iMakeCohensD(ANpeaklat, Mpeaklat)

savefile = 'BF-3 Peaklat RelresStats.mat';
save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol','ANvsAWcohenD','ANvsMcohenD')

% Total Value
disp('********Total Value********')
totalvalue_ANvsAW = horzcat(ANtotalvalue', AWtotalvalue');
totalvalue_ANvsM = horzcat(ANtotalvalue', Mtotalvalue');
var_totalvalue = {'Groups','Total'};

disp('**Anesthetized vs Awake**')
AnesthetizedvsAwake = teg_repeated_measures_ANOVA(totalvalue_ANvsAW, levels, var_totalvalue);
disp('**Anesthetized vs Awake COHEN D**')
ANvsAWcohenD = iMakeCohensD(ANtotalvalue, AWtotalvalue)
disp('**Anesthetized vs Muscimol**')
AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(totalvalue_ANvsM, levels, var_totalvalue);
disp('**Anesthetized vs Muscimol COHEN D**')
ANvsMcohenD = iMakeCohensD(ANtotalvalue, Mtotalvalue)

savefile = 'BF-3 Total RelresStats.mat';
save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol','ANvsAWcohenD','ANvsMcohenD')

% Peak Value
disp('********Peak Value********')
peakvalue_ANvsAW = horzcat(ANpeakvalue', AWpeakvalue');
peakvalue_ANvsM = horzcat(ANpeakvalue', Mpeakvalue');
var_peakvalue = {'Groups','Total'};

disp('**Anesthetized vs Awake**')
AnesthetizedvsAwake = teg_repeated_measures_ANOVA(peakvalue_ANvsAW, levels, var_peakvalue);
disp('**Anesthetized vs Awake COHEN D**')
ANvsAWcohenD = iMakeCohensD(ANpeakvalue, AWpeakvalue)
disp('**Anesthetized vs Muscimol**')
AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(peakvalue_ANvsM, levels, var_peakvalue);
disp('**Anesthetized vs Muscimol COHEN D**')
ANvsMcohenD = iMakeCohensD(ANpeakvalue, Mpeakvalue)

savefile = 'BF-3 PeakValue RelresStats.mat';
save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol','ANvsAWcohenD','ANvsMcohenD')

%% BF + 1
disp('***********************STATS FOR BF + 1***********************')
ANpeakamp = []; ANpeaklat = []; ANtotalvalue = []; ANpeakvalue = [];
ANnames = fieldnames(Anesthetized);
for i1 = 1:length(ANnames);
    animal = ANnames{i1};
    BF = find((Anesthetized.(ANnames{i1}).GS_BF) == (Anesthetized.(ANnames{i1}).Frqz'));
    
    try
        relreslist = Anesthetized.(animal).SingleTrial_RELRES_raw{1,(BF+1)};
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
    
    ANpeakamp = [ANpeakamp peakampme];
    ANpeaklat = [ANpeaklat peaklat];
    ANtotalvalue = [ANtotalvalue totalvalue];
    ANpeakvalue = [ANpeakvalue peakvalue];
end

ANpeakamp = [ANpeakamp rowofnans];
ANpeaklat = [ANpeaklat rowofnans];
ANtotalvalue = [ANtotalvalue rowofnans];
ANpeakvalue = [ANpeakvalue rowofnans];


AWpeakamp = []; AWpeaklat = []; AWtotalvalue = []; AWpeakvalue = [];
AWnames = fieldnames(Awake);
for i1 = 1:length(AWnames);
    animal = AWnames{i1};
    BF = find((Awake.(AWnames{i1}).GS_BF) == (Awake.(AWnames{i1}).Frqz'));
    
    try
        relreslist = Awake.(animal).SingleTrial_RELRES_raw{1,(BF+1)};
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
    
    AWpeakamp = [AWpeakamp peakampme];
    AWpeaklat = [AWpeaklat peaklat];
    AWtotalvalue = [AWtotalvalue totalvalue];
    AWpeakvalue = [AWpeakvalue peakvalue];
end

Mpeakamp = []; Mpeaklat = []; Mtotalvalue = []; Mpeakvalue = [];
Mnames = fieldnames(Muscimol);
for i1 = 1:length(Mnames);
    animal = Mnames{i1};
    BF = find((Muscimol.(Mnames{i1}).GS_BF) == (Muscimol.(Mnames{i1}).Frqz'));
    
    try
        relreslist = Muscimol.(animal).SingleTrial_RELRES_raw{1,(BF+1)};
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
    
    Mpeakamp = [Mpeakamp peakampme];
    Mpeaklat = [Mpeaklat peaklat];
    Mtotalvalue = [Mtotalvalue totalvalue];
    Mpeakvalue = [Mpeakvalue peakvalue];
end

Mpeakamp = [Mpeakamp rowofnans];
Mpeaklat = [Mpeaklat rowofnans];
Mtotalvalue = [Mtotalvalue rowofnans];
Mpeakvalue = [Mpeakvalue rowofnans];

%for later tuning curves
plotdata.peakamp = horzcat(plotdata.peakamp, ANpeakamp, AWpeakamp, Mpeakamp); 
plotdata.peaklate = horzcat(plotdata.peaklate, ANpeaklat, AWpeaklat, Mpeaklat);
plotdata.smallmean = horzcat(plotdata.smallmean, ANpeakvalue, AWpeakvalue, Mpeakvalue);
plotdata.totalmean = horzcat(plotdata.totalmean, ANtotalvalue, AWtotalvalue, Mtotalvalue);
plotdata.group = horzcat(plotdata.group, repmat({Group{1}},1,450), repmat({Group{2}},1,450), repmat({Group{3}},1,450));
plotdata.frqz = horzcat(plotdata.frqz, repmat({Ticks{5}}, 1, 450), repmat({Ticks{5}}, 1, 450), repmat({Ticks{5}}, 1, 450));

% Peak Amp
disp('********Peak Amplitude********')
peakamp_ANvsAW = horzcat(ANpeakamp', AWpeakamp');
peakamp_ANvsM = horzcat(ANpeakamp', Mpeakamp');
var_peakamp = {'Groups','Peak Amplitude'};

disp('**Anesthetized vs Awake**')
AnesthetizedvsAwake = teg_repeated_measures_ANOVA(peakamp_ANvsAW, levels, var_peakamp);
disp('**Anesthetized vs Awake COHEN D**')
ANvsAWcohenD = iMakeCohensD(ANpeakamp, AWpeakamp)
disp('**Anesthetized vs Muscimol**')
AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(peakamp_ANvsM, levels, var_peakamp);
disp('**Anesthetized vs Muscimol COHEN D**')
ANvsMcohenD = iMakeCohensD(ANpeakamp, Mpeakamp)

savefile = 'BF+1 PeakAmp RelresStats.mat';
save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol','ANvsAWcohenD','ANvsMcohenD')

% Peak Lat
disp('********Peak Latency********')
peaklat_ANvsAW = horzcat(ANpeaklat', AWpeaklat');
peaklat_ANvsM = horzcat(ANpeaklat', Mpeaklat');
var_peaklat = {'Groups','Peak Latency'};

disp('**Anesthetized vs Awake**')
AnesthetizedvsAwake = teg_repeated_measures_ANOVA(peaklat_ANvsAW, levels, var_peaklat);
disp('**Anesthetized vs Awake COHEN D**')
ANvsAWcohenD = iMakeCohensD(ANpeaklat, AWpeaklat)
disp('**Anesthetized vs Muscimol**')
AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(peaklat_ANvsM, levels, var_peaklat);
disp('**Anesthetized vs Muscimol COHEN D**')
ANvsMcohenD = iMakeCohensD(ANpeaklat, Mpeaklat)

savefile = 'BF+1 Peaklat RelresStats.mat';
save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol','ANvsAWcohenD','ANvsMcohenD')

% Total Value
disp('********Total Value********')
totalvalue_ANvsAW = horzcat(ANtotalvalue', AWtotalvalue');
totalvalue_ANvsM = horzcat(ANtotalvalue', Mtotalvalue');
var_totalvalue = {'Groups','Total'};

disp('**Anesthetized vs Awake**')
AnesthetizedvsAwake = teg_repeated_measures_ANOVA(totalvalue_ANvsAW, levels, var_totalvalue);
disp('**Anesthetized vs Awake COHEN D**')
ANvsAWcohenD = iMakeCohensD(ANtotalvalue, AWtotalvalue)
disp('**Anesthetized vs Muscimol**')
AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(totalvalue_ANvsM, levels, var_totalvalue);
disp('**Anesthetized vs Muscimol COHEN D**')
ANvsMcohenD = iMakeCohensD(ANtotalvalue, Mtotalvalue)

savefile = 'BF+1 Total RelresStats.mat';
save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol','ANvsAWcohenD','ANvsMcohenD')

% Peak Value
disp('********Peak Value********')
peakvalue_ANvsAW = horzcat(ANpeakvalue', AWpeakvalue');
peakvalue_ANvsM = horzcat(ANpeakvalue', Mpeakvalue');
var_peakvalue = {'Groups','Total'};

disp('**Anesthetized vs Awake**')
AnesthetizedvsAwake = teg_repeated_measures_ANOVA(peakvalue_ANvsAW, levels, var_peakvalue);
disp('**Anesthetized vs Awake COHEN D**')
ANvsAWcohenD = iMakeCohensD(ANpeakvalue, AWpeakvalue)
disp('**Anesthetized vs Muscimol**')
AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(peakvalue_ANvsM, levels, var_peakvalue);
disp('**Anesthetized vs Muscimol COHEN D**')
ANvsMcohenD = iMakeCohensD(ANpeakvalue, Mpeakvalue)

savefile = 'BF+1 PeakValue RelresStats.mat';
save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol','ANvsAWcohenD','ANvsMcohenD')

%% BF + 2
disp('***********************STATS FOR BF + 2***********************')
ANpeakamp = []; ANpeaklat = []; ANtotalvalue = []; ANpeakvalue = [];
ANnames = fieldnames(Anesthetized);
for i1 = 1:length(ANnames);
    animal = ANnames{i1};
    BF = find((Anesthetized.(ANnames{i1}).GS_BF) == (Anesthetized.(ANnames{i1}).Frqz'));
    
    try
        relreslist = Anesthetized.(animal).SingleTrial_RELRES_raw{1,(BF+2)};
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
    
    ANpeakamp = [ANpeakamp peakampme];
    ANpeaklat = [ANpeaklat peaklat];
    ANtotalvalue = [ANtotalvalue totalvalue];
    ANpeakvalue = [ANpeakvalue peakvalue];
end

ANpeakamp = [ANpeakamp rowofnans];
ANpeaklat = [ANpeaklat rowofnans];
ANtotalvalue = [ANtotalvalue rowofnans];
ANpeakvalue = [ANpeakvalue rowofnans];


AWpeakamp = []; AWpeaklat = []; AWtotalvalue = []; AWpeakvalue = [];
AWnames = fieldnames(Awake);
for i1 = 1:length(AWnames);
    animal = AWnames{i1};
    BF = find((Awake.(AWnames{i1}).GS_BF) == (Awake.(AWnames{i1}).Frqz'));
    
    try
        relreslist = Awake.(animal).SingleTrial_RELRES_raw{1,(BF+2)};
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
    
    AWpeakamp = [AWpeakamp peakampme];
    AWpeaklat = [AWpeaklat peaklat];
    AWtotalvalue = [AWtotalvalue totalvalue];
    AWpeakvalue = [AWpeakvalue peakvalue];
end

Mpeakamp = []; Mpeaklat = []; Mtotalvalue = []; Mpeakvalue = [];
Mnames = fieldnames(Muscimol);
for i1 = 1:length(Mnames);
    animal = Mnames{i1};
    BF = find((Muscimol.(Mnames{i1}).GS_BF) == (Muscimol.(Mnames{i1}).Frqz'));
    
    try
        relreslist = Muscimol.(animal).SingleTrial_RELRES_raw{1,(BF+2)};
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
    
    Mpeakamp = [Mpeakamp peakampme];
    Mpeaklat = [Mpeaklat peaklat];
    Mtotalvalue = [Mtotalvalue totalvalue];
    Mpeakvalue = [Mpeakvalue peakvalue];
end

Mpeakamp = [Mpeakamp rowofnans];
Mpeaklat = [Mpeaklat rowofnans];
Mtotalvalue = [Mtotalvalue rowofnans];
Mpeakvalue = [Mpeakvalue rowofnans];

%for later tuning curves
plotdata.peakamp = horzcat(plotdata.peakamp, ANpeakamp, AWpeakamp, Mpeakamp); 
plotdata.peaklate = horzcat(plotdata.peaklate, ANpeaklat, AWpeaklat, Mpeaklat);
plotdata.smallmean = horzcat(plotdata.smallmean, ANpeakvalue, AWpeakvalue, Mpeakvalue);
plotdata.totalmean = horzcat(plotdata.totalmean, ANtotalvalue, AWtotalvalue, Mtotalvalue);
plotdata.group = horzcat(plotdata.group, repmat({Group{1}},1,450), repmat({Group{2}},1,450), repmat({Group{3}},1,450));
plotdata.frqz = horzcat(plotdata.frqz, repmat({Ticks{6}}, 1, 450), repmat({Ticks{6}}, 1, 450), repmat({Ticks{6}}, 1, 450));

% Peak Amp
disp('********Peak Amplitude********')
peakamp_ANvsAW = horzcat(ANpeakamp', AWpeakamp');
peakamp_ANvsM = horzcat(ANpeakamp', Mpeakamp');
var_peakamp = {'Groups','Peak Amplitude'};

disp('**Anesthetized vs Awake**')
AnesthetizedvsAwake = teg_repeated_measures_ANOVA(peakamp_ANvsAW, levels, var_peakamp);
disp('**Anesthetized vs Awake COHEN D**')
ANvsAWcohenD = iMakeCohensD(ANpeakamp, AWpeakamp)
disp('**Anesthetized vs Muscimol**')
AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(peakamp_ANvsM, levels, var_peakamp);
disp('**Anesthetized vs Muscimol COHEN D**')
ANvsMcohenD = iMakeCohensD(ANpeakamp, Mpeakamp)

savefile = 'BF+2 PeakAmp RelresStats.mat';
save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol','ANvsAWcohenD','ANvsMcohenD')

% Peak Lat
disp('********Peak Latency********')
peaklat_ANvsAW = horzcat(ANpeaklat', AWpeaklat');
peaklat_ANvsM = horzcat(ANpeaklat', Mpeaklat');
var_peaklat = {'Groups','Peak Latency'};

disp('**Anesthetized vs Awake**')
AnesthetizedvsAwake = teg_repeated_measures_ANOVA(peaklat_ANvsAW, levels, var_peaklat);
disp('**Anesthetized vs Awake COHEN D**')
ANvsAWcohenD = iMakeCohensD(ANpeaklat, AWpeaklat)
disp('**Anesthetized vs Muscimol**')
AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(peaklat_ANvsM, levels, var_peaklat);
disp('**Anesthetized vs Muscimol COHEN D**')
ANvsMcohenD = iMakeCohensD(ANpeaklat, Mpeaklat)

savefile = 'BF+2 Peaklat RelresStats.mat';
save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol','ANvsAWcohenD','ANvsMcohenD')

% Total Value
disp('********Total Value********')
totalvalue_ANvsAW = horzcat(ANtotalvalue', AWtotalvalue');
totalvalue_ANvsM = horzcat(ANtotalvalue', Mtotalvalue');
var_totalvalue = {'Groups','Total'};

disp('**Anesthetized vs Awake**')
AnesthetizedvsAwake = teg_repeated_measures_ANOVA(totalvalue_ANvsAW, levels, var_totalvalue);
disp('**Anesthetized vs Awake COHEN D**')
ANvsAWcohenD = iMakeCohensD(ANtotalvalue, AWtotalvalue)
disp('**Anesthetized vs Muscimol**')
AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(totalvalue_ANvsM, levels, var_totalvalue);
disp('**Anesthetized vs Muscimol COHEN D**')
ANvsMcohenD = iMakeCohensD(ANtotalvalue, Mtotalvalue)

savefile = 'BF+2 Total RelresStats.mat';
save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol','ANvsAWcohenD','ANvsMcohenD')

% Peak Value
disp('********Peak Value********')
peakvalue_ANvsAW = horzcat(ANpeakvalue', AWpeakvalue');
peakvalue_ANvsM = horzcat(ANpeakvalue', Mpeakvalue');
var_peakvalue = {'Groups','Total'};

disp('**Anesthetized vs Awake**')
AnesthetizedvsAwake = teg_repeated_measures_ANOVA(peakvalue_ANvsAW, levels, var_peakvalue);
disp('**Anesthetized vs Awake COHEN D**')
ANvsAWcohenD = iMakeCohensD(ANpeakvalue, AWpeakvalue)
disp('**Anesthetized vs Muscimol**')
AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(peakvalue_ANvsM, levels, var_peakvalue);
disp('**Anesthetized vs Muscimol COHEN D**')
ANvsMcohenD = iMakeCohensD(ANpeakvalue, Mpeakvalue)

savefile = 'BF+2 PeakValue RelresStats.mat';
save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol','ANvsAWcohenD','ANvsMcohenD')

%% BF + 3
disp('***********************STATS FOR BF + 3***********************')
ANpeakamp = []; ANpeaklat = []; ANtotalvalue = []; ANpeakvalue = [];
ANnames = fieldnames(Anesthetized);
for i1 = 1:length(ANnames);
    animal = ANnames{i1};
    BF = find((Anesthetized.(ANnames{i1}).GS_BF) == (Anesthetized.(ANnames{i1}).Frqz'));
    
    try
        relreslist = Anesthetized.(animal).SingleTrial_RELRES_raw{1,(BF+3)};
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
    
    ANpeakamp = [ANpeakamp peakampme];
    ANpeaklat = [ANpeaklat peaklat];
    ANtotalvalue = [ANtotalvalue totalvalue];
    ANpeakvalue = [ANpeakvalue peakvalue];
end

ANpeakamp = [ANpeakamp rowofnans];
ANpeaklat = [ANpeaklat rowofnans];
ANtotalvalue = [ANtotalvalue rowofnans];
ANpeakvalue = [ANpeakvalue rowofnans];


AWpeakamp = []; AWpeaklat = []; AWtotalvalue = []; AWpeakvalue = [];
AWnames = fieldnames(Awake);
for i1 = 1:length(AWnames);
    animal = AWnames{i1};
    BF = find((Awake.(AWnames{i1}).GS_BF) == (Awake.(AWnames{i1}).Frqz'));
    
    try
        relreslist = Awake.(animal).SingleTrial_RELRES_raw{1,(BF+3)};
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
    
    AWpeakamp = [AWpeakamp peakampme];
    AWpeaklat = [AWpeaklat peaklat];
    AWtotalvalue = [AWtotalvalue totalvalue];
    AWpeakvalue = [AWpeakvalue peakvalue];
end

Mpeakamp = []; Mpeaklat = []; Mtotalvalue = []; Mpeakvalue = [];
Mnames = fieldnames(Muscimol);
for i1 = 1:length(Mnames);
    animal = Mnames{i1};
    BF = find((Muscimol.(Mnames{i1}).GS_BF) == (Muscimol.(Mnames{i1}).Frqz'));
    
    try
        relreslist = Muscimol.(animal).SingleTrial_RELRES_raw{1,(BF+3)};
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
    
    Mpeakamp = [Mpeakamp peakampme];
    Mpeaklat = [Mpeaklat peaklat];
    Mtotalvalue = [Mtotalvalue totalvalue];
    Mpeakvalue = [Mpeakvalue peakvalue];
end

Mpeakamp = [Mpeakamp rowofnans];
Mpeaklat = [Mpeaklat rowofnans];
Mtotalvalue = [Mtotalvalue rowofnans];
Mpeakvalue = [Mpeakvalue rowofnans];

%for later tuning curves
plotdata.peakamp = horzcat(plotdata.peakamp, ANpeakamp, AWpeakamp, Mpeakamp); 
plotdata.peaklate = horzcat(plotdata.peaklate, ANpeaklat, AWpeaklat, Mpeaklat);
plotdata.smallmean = horzcat(plotdata.smallmean, ANpeakvalue, AWpeakvalue, Mpeakvalue);
plotdata.totalmean = horzcat(plotdata.totalmean, ANtotalvalue, AWtotalvalue, Mtotalvalue);
plotdata.group = horzcat(plotdata.group, repmat({Group{1}},1,450), repmat({Group{2}},1,450), repmat({Group{3}},1,450));
plotdata.frqz = horzcat(plotdata.frqz, repmat({Ticks{7}}, 1, 450), repmat({Ticks{7}}, 1, 450), repmat({Ticks{7}}, 1, 450));

% Peak Amp
disp('********Peak Amplitude********')
peakamp_ANvsAW = horzcat(ANpeakamp', AWpeakamp');
peakamp_ANvsM = horzcat(ANpeakamp', Mpeakamp');
var_peakamp = {'Groups','Peak Amplitude'};

disp('**Anesthetized vs Awake**')
AnesthetizedvsAwake = teg_repeated_measures_ANOVA(peakamp_ANvsAW, levels, var_peakamp);
disp('**Anesthetized vs Awake COHEN D**')
ANvsAWcohenD = iMakeCohensD(ANpeakamp, AWpeakamp)
disp('**Anesthetized vs Muscimol**')
AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(peakamp_ANvsM, levels, var_peakamp);
disp('**Anesthetized vs Muscimol COHEN D**')
ANvsMcohenD = iMakeCohensD(ANpeakamp, Mpeakamp)

savefile = 'BF+3 PeakAmp RelresStats.mat';
save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol','ANvsAWcohenD','ANvsMcohenD')

% Peak Lat
disp('********Peak Latency********')
peaklat_ANvsAW = horzcat(ANpeaklat', AWpeaklat');
peaklat_ANvsM = horzcat(ANpeaklat', Mpeaklat');
var_peaklat = {'Groups','Peak Latency'};

disp('**Anesthetized vs Awake**')
AnesthetizedvsAwake = teg_repeated_measures_ANOVA(peaklat_ANvsAW, levels, var_peaklat);
disp('**Anesthetized vs Awake COHEN D**')
ANvsAWcohenD = iMakeCohensD(ANpeaklat, AWpeaklat)
disp('**Anesthetized vs Muscimol**')
AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(peaklat_ANvsM, levels, var_peaklat);
disp('**Anesthetized vs Muscimol COHEN D**')
ANvsMcohenD = iMakeCohensD(ANpeaklat, Mpeaklat)

savefile = 'BF+3 Peaklat RelresStats.mat';
save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol','ANvsAWcohenD','ANvsMcohenD')

% Total Value
disp('********Total Value********')
totalvalue_ANvsAW = horzcat(ANtotalvalue', AWtotalvalue');
totalvalue_ANvsM = horzcat(ANtotalvalue', Mtotalvalue');
var_totalvalue = {'Groups','Total'};

disp('**Anesthetized vs Awake**')
AnesthetizedvsAwake = teg_repeated_measures_ANOVA(totalvalue_ANvsAW, levels, var_totalvalue);
disp('**Anesthetized vs Awake COHEN D**')
ANvsAWcohenD = iMakeCohensD(ANtotalvalue, AWtotalvalue)
disp('**Anesthetized vs Muscimol**')
AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(totalvalue_ANvsM, levels, var_totalvalue);
disp('**Anesthetized vs Muscimol COHEN D**')
ANvsMcohenD = iMakeCohensD(ANtotalvalue, Mtotalvalue)

savefile = 'BF+3 Total RelresStats.mat';
save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol','ANvsAWcohenD','ANvsMcohenD')

% Peak Value
disp('********Peak Value********')
peakvalue_ANvsAW = horzcat(ANpeakvalue', AWpeakvalue');
peakvalue_ANvsM = horzcat(ANpeakvalue', Mpeakvalue');
var_peakvalue = {'Groups','Total'};

disp('**Anesthetized vs Awake**')
AnesthetizedvsAwake = teg_repeated_measures_ANOVA(peakvalue_ANvsAW, levels, var_peakvalue);
disp('**Anesthetized vs Awake COHEN D**')
ANvsAWcohenD = iMakeCohensD(ANpeakvalue, AWpeakvalue)
disp('**Anesthetized vs Muscimol**')
AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(peakvalue_ANvsM, levels, var_peakvalue);
disp('**Anesthetized vs Muscimol COHEN D**')
ANvsMcohenD = iMakeCohensD(ANpeakvalue, Mpeakvalue)

savefile = 'BF+3 PeakValue RelresStats.mat';
save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol','ANvsAWcohenD','ANvsMcohenD')





%% PEAK AMP Tuning curve
feature = {'peakamp','peaklate', 'smallmean','totalmean'};
Name = {'Maximum Depth', 'Depth latency', 'RMS between 200 and 300 ms', 'RMS of total Relres'};
cd(home); cd figs; mkdir('Group Relres Tuning'); cd('Group Relres Tuning')

for i1 = 1:length(feature);
clear g
figure();

g(1,1)=gramm('x',plotdata.frqz,'y',plotdata.(feature{i1}), 'color', plotdata.group); %,'y',cars.Acceleration,'color',cars.Cylinders,'subset',cars.Cylinders~=3 & cars.Cylinders~=5
g(1,1).stat_summary('type','sem','geom','errorbar'); %mean and sem shown
g(1,1).stat_summary('type','sem','geom','point'); %mean and sem shown
g(1,1).stat_summary('type','sem','geom','area'); %mean and sem shown
g(1,1).set_text_options('base_size',12) 
if strcmp(feature{i1},'peaklate')
    g(1,1).set_names('x','Tone','y','ms','color','Group');
    g(1,1).axe_property('YLim',[230 250]);
elseif strcmp(feature{i1},'peakamp')
    g(1,1).set_names('x','Tone','y','%','color','Group');
    g(1,1).axe_property('YLim',[-0.5 -0.2]);
else
    g(1,1).set_names('x','Tone','y','%','color','Group');
    g(1,1).axe_property('YLim',[0.1 0.4]);
end
g(1,1).set_color_options('map','matlab');
g(1,1).axe_property('XTickLabel',{'-3', '-2', '-1', 'BF', '+1', '+2', '+3'})

g.set_title(Name{i1});
g.draw();
g.export('file_name',Name{i1}, 'file_type','pdf');
g.export('file_name',Name{i1}, 'file_type','png');
close all;
end