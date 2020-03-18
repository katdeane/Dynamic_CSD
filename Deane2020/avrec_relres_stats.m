function avrec_relres_stats(homedir)

%  This script is SPECIFIC for running the AVREC stats for Katrina's
%  Master's thesis. It may be used as a template but please make a copy
%  elsewhere for modicifications

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

%for teg_repeated measures_ANOVA
levels = [2,1];
rowofnans = NaN(1,50); %IMPORTANT: if animal group sizes are equal this isn't needed. Currently these are commented in the code

Ticks = {'a -3','b -2', 'c -1','d BF', 'e +1', 'f +2', 'g +3'};
Ticknames = {'BF-3', 'BF-2','BF-1','BF','BF+1', 'BF+2', 'BF+3'};
Group = {'Anesthetized','Awake','Muscimol','ANChronic'};
FromBF = [-3,-2,-1,0,+1,+2,+3];

%% Load in the appropriate files
cd(homedir);cd DATA;
load('AnesthetizedPre_Data.mat','Data')
Anesthetized = Data; clear Data;
load('Awake10dB_Data.mat','Data')
Awake = Data; clear Data;
load('Muscimol_Data.mat','Data')
Muscimol = Data; clear Data;
load('ANChronic_Data.mat','Data')
ANChronic = Data; clear Data;

mkdir('Stats AvRec Relres single'); cd('Stats AvRec Relres single')
%% Run Stats

plottune = struct;
plottuner = struct;
%for tuning curves
grpsz = 550;
plottune.frqz = [];
plottune.tgroup = [];
plottune.peakamp = [];plottune.peaklate = [];plottune.peakvalue = [];plottune.totalmean = [];
plottuner.frqz = [];
plottuner.tgroup = [];
plottuner.peakamp = [];plottuner.peaklate = [];plottuner.peakvalue = [];plottuner.totalmean = [];

for ifreq = 1:length(Ticks)
    
    disp(['***********************AVREC STATS FOR ' (Ticknames{ifreq}) '***********************'])
    % KETAMINE
    ANpeakamp = nan(1,grpsz); ANpeaklat = nan(1,grpsz); 
    ANtotalvalue = nan(1,grpsz); ANpeakvalue = nan(1,grpsz);
    ThisVec = 1:50;
    ANnames = fieldnames(Anesthetized);
    for i1 = 1:length(ANnames)
        animal = ANnames{i1};
        %Find the best frequency
        BF = find((Anesthetized.(ANnames{i1}).GS_BF) == (Anesthetized.(ANnames{i1}).Frqz'));
        
        try %cut out the section of the matrix necessary (BF, BF-1, etc)
            avreclist =  Anesthetized.(animal).SingleTrial_AVREC_raw{1,(BF+FromBF(ifreq))};
        catch %produce NAN if there isn't an entry here
            avreclist = NaN(1,length(avreclist));
        end
        
        peakampme = nan(1,50); peaklat = nan(1,50); 
        totalvalue = nan(1,50); peakvalue = nan(1,50);
        
        for i2 = 1:50
            if isnan(avreclist)
                continue
            else
                [latency, ampl, meanOfdata] = iGetPeakData_avrec(avreclist(:,:,i2),[200 300]);
                [~, ~, meanOffulldata] = iGetPeakData_avrec(avreclist(:,:,i2));
                peaklat(i2) = latency;
                peakampme(i2) = ampl;
                peakvalue(i2) = meanOfdata;
                totalvalue(i2) = meanOffulldata;
            end
        end
        
        %apply a cutoff threshold
        if sum(isnan(peakampme)) > length(peakampme)*.75 % 25% at least
            peakampme = nan(1,50);
            peaklat = nan(1,50);
            peakvalue = nan(1,50);
        end
        
        ANpeakamp(ThisVec) = peakampme;
        ANpeaklat(ThisVec) = peaklat;
        ANtotalvalue(ThisVec) = totalvalue;
        ANpeakvalue(ThisVec) = peakvalue;
        ThisVec = ThisVec + 50;
    end
        
    % AWAKE
    AWpeakamp = nan(1,grpsz); AWpeaklat = nan(1,grpsz); 
    AWtotalvalue = nan(1,grpsz); AWpeakvalue = nan(1,grpsz);
    ThisVec = 1:50;
    AWnames = fieldnames(Awake);
    for i1 = 1:length(AWnames)
        animal = AWnames{i1};
        BF = find((Awake.(AWnames{i1}).GS_BF) == (Awake.(AWnames{i1}).Frqz'));
        try
            avreclist =  Awake.(animal).SingleTrial_AVREC_raw{:,(BF+FromBF(ifreq))};
        catch
            avreclist = NaN(1,length(avreclist));
        end
        peakampme = nan(1,50); peaklat = nan(1,50); 
        totalvalue = nan(1,50); peakvalue = nan(1,50);
        for i2 = 1:50
            if isnan(avreclist)
                continue
            else
                try
                    [latency, ampl, meanOfdata] = iGetPeakData_avrec(avreclist(:,:,i2),[200 300]);
                    [~, ~, meanOffulldata] = iGetPeakData_avrec(avreclist(:,:,i2));
                    peaklat(i2) = latency;
                    peakampme(i2) = ampl;
                    peakvalue(i2) = meanOfdata;
                    totalvalue(i2) = meanOffulldata;
                catch
                    continue
                end
            end
        end
        
        %apply a cutoff threshold
        if sum(isnan(peakampme)) > length(peakampme)*.75 %&& ifreq ~=5 % 25% at least
            peakampme = nan(1,50);
            peaklat = nan(1,50);
            peakvalue = nan(1,50);
        end
        
        AWpeakamp(ThisVec) = peakampme;
        AWpeaklat(ThisVec) = peaklat;
        AWtotalvalue(ThisVec) = totalvalue;
        AWpeakvalue(ThisVec) = peakvalue;
        ThisVec = ThisVec + 50;
    end
    
    % MUSCIMOL
    Mpeakamp = nan(1,grpsz); Mpeaklat = nan(1,grpsz); 
    Mtotalvalue = nan(1,grpsz); Mpeakvalue = nan(1,grpsz);
    ThisVec = 1:50;
    Mnames = fieldnames(Muscimol);
    for i1 = 1:length(Mnames)
        animal = Mnames{i1};
        BF = find((Muscimol.(Mnames{i1}).GS_BF) == (Muscimol.(Mnames{i1}).Frqz'));
        try
            avreclist =  Muscimol.(animal).SingleTrial_AVREC_raw{:,(BF+FromBF(ifreq))};
        catch
            avreclist = NaN(1,length(avreclist));
        end
        peakampme = nan(1,50); peaklat = nan(1,50); 
        totalvalue = nan(1,50); peakvalue = nan(1,50);
        for i2 = 1:50
            if isnan(avreclist)
                continue
            else
                [latency, ampl, meanOfdata] = iGetPeakData_avrec(avreclist(:,:,i2),[200 300]);
                [~, ~, meanOffulldata] = iGetPeakData_avrec(avreclist(:,:,i2));
                peaklat(i2) = latency;
                peakampme(i2) = ampl;
                peakvalue(i2) = meanOfdata;
                totalvalue(i2) = meanOffulldata;
            end
        end
        
        %apply a cutoff threshold
        if sum(isnan(peakampme)) > length(peakampme)*.75 % 25% at least
            peakampme = nan(1,50);
            peaklat = nan(1,50);
            peakvalue = nan(1,50);
        end
            
        Mpeakamp(ThisVec) = peakampme;
        Mpeaklat(ThisVec) = peaklat;
        Mtotalvalue(ThisVec) = totalvalue;
        Mpeakvalue(ThisVec) = peakvalue;
        ThisVec = ThisVec + 50;
    end
        
    % ANESTHETIZED CHRONIC
    ANCpeakamp = nan(1,grpsz); ANCpeaklat = nan(1,grpsz); 
    ANCtotalvalue = nan(1,grpsz); ANCpeakvalue = nan(1,grpsz);
    ThisVec = 1:50;
    ANCnames = fieldnames(ANChronic);
    for i1 = 1:length(ANCnames)
        animal = ANCnames{i1};
        BF = find((ANChronic.(ANCnames{i1}).GS_BF) == (ANChronic.(ANCnames{i1}).Frqz'));
        try
            avreclist =  ANChronic.(animal).SingleTrial_AVREC_raw{:,(BF+FromBF(ifreq))};
        catch
            avreclist = NaN(1,length(avreclist));
        end
        peakampme = nan(1,50); peaklat = nan(1,50); 
        totalvalue = nan(1,50); peakvalue = nan(1,50);
        for i2 = 1:50
            if isnan(avreclist)
                continue
            else
                try
                    [latency, ampl, meanOfdata] = iGetPeakData_avrec(avreclist(:,:,i2),[200 300]);
                    [~, ~, meanOffulldata] = iGetPeakData_avrec(avreclist(:,:,i2));
                    peaklat(i2) = latency;
                    peakampme(i2) = ampl;
                    peakvalue(i2) = meanOfdata;
                    totalvalue(i2) = meanOffulldata;
                catch
                    continue
                end
            end
        end
        
        %apply a cutoff threshold
        if sum(isnan(peakampme)) > length(peakampme)*.75 % && ifreq ~=5 % 25% at least
            peakampme = NaN(1,50);
            peaklat = NaN(1,50);
            peakvalue = NaN(1,50);
        end
        
        ANCpeakamp(ThisVec) = peakampme;
        ANCpeaklat(ThisVec) = peaklat;
        ANCtotalvalue(ThisVec) = totalvalue;
        ANCpeakvalue(ThisVec) = peakvalue;
        ThisVec = ThisVec + 50;
    end
      
    %for later tuning curves
    plottune.peakamp = horzcat(plottune.peakamp, ANpeakamp, AWpeakamp, Mpeakamp, ANCpeakamp);
    plottune.peaklate = horzcat(plottune.peaklate, ANpeaklat, AWpeaklat, Mpeaklat, ANCpeaklat);
    plottune.peakvalue = horzcat(plottune.peakvalue, ANpeakvalue, AWpeakvalue, Mpeakvalue, ANCpeakvalue);
    plottune.totalmean = horzcat(plottune.totalmean, ANtotalvalue, AWtotalvalue, Mtotalvalue, ANCtotalvalue);
    plottune.tgroup = horzcat(plottune.tgroup, repmat({Group{1}},1,grpsz), ...
        repmat({Group{2}},1,grpsz), repmat({Group{3}},1,grpsz), repmat({Group{4}},1,grpsz)); %#ok<*CCAT1>
    plottune.frqz = horzcat(plottune.frqz, repmat({Ticks{ifreq}}, 1, grpsz), ...
        repmat({Ticks{ifreq}}, 1, grpsz), repmat({Ticks{ifreq}}, 1, grpsz), repmat({Ticks{ifreq}}, 1, grpsz));
    
    % Peak Amp
    % for special case
    ANCpeakamp_noNan = ANCpeakamp(~isnan(ANCpeakamp));
    AWpeakamp_noNan = AWpeakamp(~isnan(AWpeakamp));
    szdif = length(ANCpeakamp_noNan) - length(AWpeakamp_noNan);
    
    if szdif < 0
        for i = 1:abs(szdif)
            ANCpeakamp_noNan = [ANCpeakamp_noNan NaN];
        end
    elseif szdif > 0
        for i = 1:szdif
            AWpeakamp_noNan = [AWpeakamp_noNan NaN];
        end
    elseif isempty(ANCpeakamp_noNan) && isempty(AWpeakamp_noNan) 
        ANCpeakamp_noNan = ANCpeakamp;
        AWpeakamp_noNan = AWpeakamp;
    end

    disp('********Peak Amplitude********')
    peakamp_ANvsAW = horzcat(ANpeakamp', AWpeakamp');
    peakamp_ANCvsAW = horzcat(ANCpeakamp_noNan', AWpeakamp_noNan');
    peakamp_ANvsM = horzcat(ANpeakamp', Mpeakamp');
    var_peakamp = {'Groups','Peak Amplitude'};
    
    disp('**Anesthetized vs Awake**')
    AnesthetizedvsAwake = teg_repeated_measures_ANOVA(peakamp_ANvsAW, levels, var_peakamp);
    disp('**Anesthetized vs Awake COHEN D**')
    ANvsAWcohenD = iMakeCohensD(ANpeakamp, AWpeakamp) %#ok<*NOPRT>
    disp('**ANChronic vs Awake**')
    ANChronicvsAwake = teg_repeated_measures_ANOVA(peakamp_ANCvsAW, levels, var_peakamp);
    disp('**ANChronic vs Awake COHEN D**')
    ANCvsAWcohenD = iMakeCohensD(ANCpeakamp, AWpeakamp)
    disp('**Anesthetized vs Muscimol**')
    AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(peakamp_ANvsM, levels, var_peakamp);
    disp('**Anesthetized vs Muscimol COHEN D**')
    ANvsMcohenD = iMakeCohensD(ANpeakamp, Mpeakamp)
    
    savefile = [Ticknames{ifreq} 'PeakAmp STAvrecStats.mat'];
    save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol','ANChronicvsAwake', ...
        'ANvsAWcohenD', 'ANCvsAWcohenD', 'ANvsMcohenD')
    
    % Peak Lat
    % for special case
    ANCpeaklat_noNan = ANCpeaklat(~isnan(ANCpeaklat));
    AWpeaklat_noNan = AWpeaklat(~isnan(AWpeaklat));
    szdif = length(ANCpeaklat_noNan) - length(AWpeaklat_noNan);
    
    if szdif < 0
        for i = 1:abs(szdif)
            ANCpeaklat_noNan = [ANCpeaklat_noNan NaN];
        end
    elseif szdif > 0
        for i = 1:szdif
            AWpeaklat_noNan = [AWpeaklat_noNan NaN];
        end
    elseif isempty(ANCpeaklat_noNan) && isempty(AWpeaklat_noNan) 
        ANCpeaklat_noNan = ANCpeaklat;
        AWpeaklat_noNan = AWpeaklat;
    end

    disp('********Peak Latency********')
    peaklat_ANvsAW = horzcat(ANpeaklat', AWpeaklat');
    peaklat_ANCvsAW = horzcat(ANCpeaklat_noNan', AWpeaklat_noNan');
    peaklat_ANvsM = horzcat(ANpeaklat', Mpeaklat');
    var_peaklat = {'Groups','Peak Latency'};
    
    disp('**Anesthetized vs Awake**')
    AnesthetizedvsAwake = teg_repeated_measures_ANOVA(peaklat_ANvsAW, levels, var_peaklat);
    disp('**Anesthetized vs Awake COHEN D**')
    ANvsAWcohenD = iMakeCohensD(ANpeaklat, AWpeaklat)
    disp('**ANChronic vs Awake**')
    ANChronicvsAwake = teg_repeated_measures_ANOVA(peaklat_ANCvsAW, levels, var_peaklat);
    disp('**ANChronic vs Awake COHEN D**')
    ANCvsAWcohenD = iMakeCohensD(ANCpeaklat, AWpeaklat)
    disp('**Anesthetized vs Muscimol**')
    AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(peaklat_ANvsM, levels, var_peaklat);
    disp('**Anesthetized vs Muscimol COHEN D**')
    ANvsMcohenD = iMakeCohensD(ANpeaklat, Mpeaklat)
    
    savefile = [Ticknames{ifreq} ' Peaklat STAvrecStats.mat'];
    save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol','ANChronicvsAwake', ...
        'ANvsAWcohenD', 'ANCvsAWcohenD', 'ANvsMcohenD')
    
    % Total Value
    % for special case
    ANCtotalvalue_noNan = ANCtotalvalue(~isnan(ANCtotalvalue));
    AWtotalvalue_noNan = AWtotalvalue(~isnan(AWtotalvalue));
    szdif = length(ANCtotalvalue_noNan) - length(AWtotalvalue_noNan);
    
    if szdif < 0
        for i = 1:abs(szdif)
            ANCtotalvalue_noNan = [ANCtotalvalue_noNan NaN];
        end
    elseif szdif > 0
        for i = 1:szdif
            AWtotalvalue_noNan = [AWtotalvalue_noNan NaN];
        end
    elseif isempty(ANCtotalvalue_noNan) && isempty(AWtotalvalue_noNan) 
        ANCtotalvalue_noNan = ANCtotalvalue;
        AWtotalvalue_noNan = AWtotalvalue;
    end

    disp('********Total Value********')
    totalvalue_ANvsAW = horzcat(ANtotalvalue', AWtotalvalue');
    totalvalue_ANCvsAW = horzcat(ANCtotalvalue_noNan', AWtotalvalue_noNan');
    totalvalue_ANvsM = horzcat(ANtotalvalue', Mtotalvalue');
    var_totalvalue = {'Groups','Total'};
    
    disp('**Anesthetized vs Awake**')
    AnesthetizedvsAwake = teg_repeated_measures_ANOVA(totalvalue_ANvsAW, levels, var_totalvalue);
    disp('**Anesthetized vs Awake COHEN D**')
    ANvsAWcohenD = iMakeCohensD(ANtotalvalue, AWtotalvalue)
    disp('**ANChronic vs Awake**')
    ANChronicvsAwake = teg_repeated_measures_ANOVA(totalvalue_ANCvsAW, levels, var_totalvalue);
    disp('**ANChronic vs Awake COHEN D**')
    ANCvsAWcohenD = iMakeCohensD(ANCtotalvalue, AWtotalvalue)
    disp('**Anesthetized vs Muscimol**')
    AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(totalvalue_ANvsM, levels, var_totalvalue);
    disp('**Anesthetized vs Muscimol COHEN D**')
    ANvsMcohenD = iMakeCohensD(ANtotalvalue, Mtotalvalue)
    
    savefile = [Ticknames{ifreq} ' Total STAvrecStats.mat'];
    save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol','ANChronicvsAwake', ...
        'ANvsAWcohenD', 'ANCvsAWcohenD', 'ANvsMcohenD')
    
    % Peak Value
    % for special case
    ANCpeakvalue_noNan = ANCpeakvalue(~isnan(ANCpeakvalue));
    AWpeakvalue_noNan = AWpeakvalue(~isnan(AWpeakvalue));
    szdif = length(ANCpeakvalue_noNan) - length(AWpeakvalue_noNan);
    
    if szdif < 0
        for i = 1:abs(szdif)
            ANCpeakvalue_noNan = [ANCpeakvalue_noNan NaN];
        end
    elseif szdif > 0
        for i = 1:szdif
            AWpeakvalue_noNan = [AWpeakvalue_noNan NaN];
        end
    elseif isempty(ANCpeakvalue_noNan) && isempty(AWpeakvalue_noNan) 
        ANCpeakvalue_noNan = ANCpeakvalue;
        AWpeakvalue_noNan = AWpeakvalue;
    end

    disp('********Peak Value********')
    peakvalue_ANvsAW = horzcat(ANpeakvalue', AWpeakvalue');
    peakvalue_ANCvsAW = horzcat(ANCpeakvalue_noNan', AWpeakvalue_noNan');
    peakvalue_ANvsM = horzcat(ANpeakvalue', Mpeakvalue');
    var_peakvalue = {'Groups','Total'};
    
    disp('**Anesthetized vs Awake**')
    AnesthetizedvsAwake = teg_repeated_measures_ANOVA(peakvalue_ANvsAW, levels, var_peakvalue);
    disp('**Anesthetized vs Awake COHEN D**')
    ANvsAWcohenD = iMakeCohensD(ANpeakvalue, AWpeakvalue)
    disp('**ANChronic vs Awake**')
    ANChronicvsAwake = teg_repeated_measures_ANOVA(peakvalue_ANCvsAW, levels, var_peakvalue);
    disp('**ANChronic vs Awake COHEN D**')
    ANCvsAWcohenD = iMakeCohensD(ANCpeakvalue, AWpeakvalue)
    disp('**Anesthetized vs Muscimol**')
    AnesthetizedvsMuscimol = teg_repeated_measures_ANOVA(peakvalue_ANvsM, levels, var_peakvalue);
    disp('**Anesthetized vs Muscimol COHEN D**')
    ANvsMcohenD = iMakeCohensD(ANpeakvalue, Mpeakvalue)
    
    savefile = [Ticknames{ifreq} ' PeakValue STAvrecStats.mat'];
    save(savefile, 'AnesthetizedvsAwake', 'AnesthetizedvsMuscimol','ANChronicvsAwake', ...
        'ANvsAWcohenD', 'ANCvsAWcohenD', 'ANvsMcohenD')
end
%%
levels = [2,7];
levels5 = [2,5];
levelsmir = [2,3];
varinames = {'Groups','Frequencies'};

% Set up tables
TAv = table(plottune.frqz,plottune.tgroup,plottune.peakamp,plottune.peaklate,plottune.peakvalue);
TAv.Properties.VariableNames = {'Freq' 'Group' 'PeakAmp' 'PeakLat' 'RMS'};

TRe = table(plottuner.frqz,plottuner.tgroup,plottuner.peakamp,plottuner.peaklate,plottuner.peakvalue);
TRe.Properties.VariableNames = {'Freq' 'Group' 'PeakAmp' 'PeakLat' 'RMS'};

%% AVREC table re-organization

% % AWAKE PEAKAMP
Am3 = TAv.PeakAmp(strcmp(TAv.Group,'Awake')&strcmp(TAv.Freq,'a -3'))';
Am2 = TAv.PeakAmp(strcmp(TAv.Group,'Awake')&strcmp(TAv.Freq,'b -2'))';
Am1 = TAv.PeakAmp(strcmp(TAv.Group,'Awake')&strcmp(TAv.Freq,'c -1'))';
ABF = TAv.PeakAmp(strcmp(TAv.Group,'Awake')&strcmp(TAv.Freq,'d BF'))';
Ap1 = TAv.PeakAmp(strcmp(TAv.Group,'Awake')&strcmp(TAv.Freq,'e +1'))';
Ap2 = TAv.PeakAmp(strcmp(TAv.Group,'Awake')&strcmp(TAv.Freq,'f +2'))';
Ap3 = TAv.PeakAmp(strcmp(TAv.Group,'Awake')&strcmp(TAv.Freq,'g +3'))';

A_PA = horzcat(Am3,Am2,Am1,ABF,Ap1,Ap2,Ap3);
A_PA_5 = horzcat(Am2,Am1,ABF,Ap1,Ap2);
A_PA_5norm = A_PA_5./A_PA_5(:,3);

twoaway = vertcat(Am2,Ap2);
oneaway = vertcat(Am1,Ap1);
central = vertcat(ABF,ABF);
A_PA_5mir = horzcat(central,oneaway,twoaway);
A_PA_5mirnorm = A_PA_5mir./A_PA_5mir(:,1);

% % AWAKE PEAKLAT
Am3 = TAv.PeakLat(strcmp(TAv.Group,'Awake')&strcmp(TAv.Freq,'a -3'))';
Am2 = TAv.PeakLat(strcmp(TAv.Group,'Awake')&strcmp(TAv.Freq,'b -2'))';
Am1 = TAv.PeakLat(strcmp(TAv.Group,'Awake')&strcmp(TAv.Freq,'c -1'))';
ABF = TAv.PeakLat(strcmp(TAv.Group,'Awake')&strcmp(TAv.Freq,'d BF'))';
Ap1 = TAv.PeakLat(strcmp(TAv.Group,'Awake')&strcmp(TAv.Freq,'e +1'))';
Ap2 = TAv.PeakLat(strcmp(TAv.Group,'Awake')&strcmp(TAv.Freq,'f +2'))';
Ap3 = TAv.PeakLat(strcmp(TAv.Group,'Awake')&strcmp(TAv.Freq,'g +3'))';

A_PL = horzcat(Am3,Am2,Am1,ABF,Ap1,Ap2,Ap3);
A_PL_5 = horzcat(Am2,Am1,ABF,Ap1,Ap2);

% % AWAKE RMS
Am3 = TAv.RMS(strcmp(TAv.Group,'Awake')&strcmp(TAv.Freq,'a -3'))';
Am2 = TAv.RMS(strcmp(TAv.Group,'Awake')&strcmp(TAv.Freq,'b -2'))';
Am1 = TAv.RMS(strcmp(TAv.Group,'Awake')&strcmp(TAv.Freq,'c -1'))';
ABF = TAv.RMS(strcmp(TAv.Group,'Awake')&strcmp(TAv.Freq,'d BF'))';
Ap1 = TAv.RMS(strcmp(TAv.Group,'Awake')&strcmp(TAv.Freq,'e +1'))';
Ap2 = TAv.RMS(strcmp(TAv.Group,'Awake')&strcmp(TAv.Freq,'f +2'))';
Ap3 = TAv.RMS(strcmp(TAv.Group,'Awake')&strcmp(TAv.Freq,'g +3'))';

A_RMS = horzcat(Am3,Am2,Am1,ABF,Ap1,Ap2,Ap3);
A_RMS_5 = horzcat(Am2,Am1,ABF,Ap1,Ap2);

% % ANCHRONIC PEAKAMP
Nm3 = TAv.PeakAmp(strcmp(TAv.Group,'ANChronic')&strcmp(TAv.Freq,'a -3'))';
Nm2 = TAv.PeakAmp(strcmp(TAv.Group,'ANChronic')&strcmp(TAv.Freq,'b -2'))';
Nm1 = TAv.PeakAmp(strcmp(TAv.Group,'ANChronic')&strcmp(TAv.Freq,'c -1'))';
NBF = TAv.PeakAmp(strcmp(TAv.Group,'ANChronic')&strcmp(TAv.Freq,'d BF'))';
Np1 = TAv.PeakAmp(strcmp(TAv.Group,'ANChronic')&strcmp(TAv.Freq,'e +1'))';
Np2 = TAv.PeakAmp(strcmp(TAv.Group,'ANChronic')&strcmp(TAv.Freq,'f +2'))';
Np3 = TAv.PeakAmp(strcmp(TAv.Group,'ANChronic')&strcmp(TAv.Freq,'g +3'))';

N_PA = horzcat(Nm3,Nm2,Nm1,NBF,Np1,Np2,Np3);
N_PA_5 = horzcat(Nm2,Nm1,NBF,Np1,Np2);
N_PA_5norm = N_PA_5./N_PA_5(:,3);

twoaway = vertcat(Nm2,Np2);
oneaway = vertcat(Nm1,Np1);
central = vertcat(NBF,NBF);
N_PA_5mir = horzcat(central,oneaway,twoaway);
N_PA_5mirnorm = N_PA_5mir./N_PA_5mir(:,1);

% % ANCHRONIC PEAKLAT
Nm3 = TAv.PeakLat(strcmp(TAv.Group,'ANChronic')&strcmp(TAv.Freq,'a -3'))';
Nm2 = TAv.PeakLat(strcmp(TAv.Group,'ANChronic')&strcmp(TAv.Freq,'b -2'))';
Nm1 = TAv.PeakLat(strcmp(TAv.Group,'ANChronic')&strcmp(TAv.Freq,'c -1'))';
NBF = TAv.PeakLat(strcmp(TAv.Group,'ANChronic')&strcmp(TAv.Freq,'d BF'))';
Np1 = TAv.PeakLat(strcmp(TAv.Group,'ANChronic')&strcmp(TAv.Freq,'e +1'))';
Np2 = TAv.PeakLat(strcmp(TAv.Group,'ANChronic')&strcmp(TAv.Freq,'f +2'))';
Np3 = TAv.PeakLat(strcmp(TAv.Group,'ANChronic')&strcmp(TAv.Freq,'g +3'))';

N_PL = horzcat(Nm3,Nm2,Nm1,NBF,Np1,Np2,Np3);
N_PL_5 = horzcat(Nm2,Nm1,NBF,Np1,Np2);

% % ANCHRONIC RMS
Nm3 = TAv.RMS(strcmp(TAv.Group,'ANChronic')&strcmp(TAv.Freq,'a -3'))';
Nm2 = TAv.RMS(strcmp(TAv.Group,'ANChronic')&strcmp(TAv.Freq,'b -2'))';
Nm1 = TAv.RMS(strcmp(TAv.Group,'ANChronic')&strcmp(TAv.Freq,'c -1'))';
NBF = TAv.RMS(strcmp(TAv.Group,'ANChronic')&strcmp(TAv.Freq,'d BF'))';
Np1 = TAv.RMS(strcmp(TAv.Group,'ANChronic')&strcmp(TAv.Freq,'e +1'))';
Np2 = TAv.RMS(strcmp(TAv.Group,'ANChronic')&strcmp(TAv.Freq,'f +2'))';
Np3 = TAv.RMS(strcmp(TAv.Group,'ANChronic')&strcmp(TAv.Freq,'g +3'))';

N_RMS = horzcat(Nm3,Nm2,Nm1,NBF,Np1,Np2,Np3);
N_RMS_5 = horzcat(Nm2,Nm1,NBF,Np1,Np2);

% % KET PEAKAMP
Km3 = TAv.PeakAmp(strcmp(TAv.Group,'Anesthetized')&strcmp(TAv.Freq,'a -3'))';
Km2 = TAv.PeakAmp(strcmp(TAv.Group,'Anesthetized')&strcmp(TAv.Freq,'b -2'))';
Km1 = TAv.PeakAmp(strcmp(TAv.Group,'Anesthetized')&strcmp(TAv.Freq,'c -1'))';
KBF = TAv.PeakAmp(strcmp(TAv.Group,'Anesthetized')&strcmp(TAv.Freq,'d BF'))';
Kp1 = TAv.PeakAmp(strcmp(TAv.Group,'Anesthetized')&strcmp(TAv.Freq,'e +1'))';
Kp2 = TAv.PeakAmp(strcmp(TAv.Group,'Anesthetized')&strcmp(TAv.Freq,'f +2'))';
Kp3 = TAv.PeakAmp(strcmp(TAv.Group,'Anesthetized')&strcmp(TAv.Freq,'g +3'))';

K_PA = horzcat(Km3,Km2,Km1,KBF,Kp1,Kp2,Kp3);
K_PA_5 = horzcat(Km2,Km1,KBF,Kp1,Kp2);
K_PA_5norm = K_PA_5./K_PA_5(:,3);

twoaway = vertcat(Km2,Kp2);
oneaway = vertcat(Km1,Kp1);
central = vertcat(KBF,KBF);
K_PA_5mir = horzcat(central,oneaway,twoaway);
K_PA_5mirnorm = K_PA_5mir./K_PA_5mir(:,1);


% % KET PEAKLAT
Km3 = TAv.PeakLat(strcmp(TAv.Group,'Anesthetized')&strcmp(TAv.Freq,'a -3'))';
Km2 = TAv.PeakLat(strcmp(TAv.Group,'Anesthetized')&strcmp(TAv.Freq,'b -2'))';
Km1 = TAv.PeakLat(strcmp(TAv.Group,'Anesthetized')&strcmp(TAv.Freq,'c -1'))';
KBF = TAv.PeakLat(strcmp(TAv.Group,'Anesthetized')&strcmp(TAv.Freq,'d BF'))';
Kp1 = TAv.PeakLat(strcmp(TAv.Group,'Anesthetized')&strcmp(TAv.Freq,'e +1'))';
Kp2 = TAv.PeakLat(strcmp(TAv.Group,'Anesthetized')&strcmp(TAv.Freq,'f +2'))';
Kp3 = TAv.PeakLat(strcmp(TAv.Group,'Anesthetized')&strcmp(TAv.Freq,'g +3'))';

K_PL = horzcat(Km3,Km2,Km1,KBF,Kp1,Kp2,Kp3);
K_PL_5 = horzcat(Km2,Km1,KBF,Kp1,Kp2);

% % KET RMS
Km3 = TAv.RMS(strcmp(TAv.Group,'Anesthetized')&strcmp(TAv.Freq,'a -3'))';
Km2 = TAv.RMS(strcmp(TAv.Group,'Anesthetized')&strcmp(TAv.Freq,'b -2'))';
Km1 = TAv.RMS(strcmp(TAv.Group,'Anesthetized')&strcmp(TAv.Freq,'c -1'))';
KBF = TAv.RMS(strcmp(TAv.Group,'Anesthetized')&strcmp(TAv.Freq,'d BF'))';
Kp1 = TAv.RMS(strcmp(TAv.Group,'Anesthetized')&strcmp(TAv.Freq,'e +1'))';
Kp2 = TAv.RMS(strcmp(TAv.Group,'Anesthetized')&strcmp(TAv.Freq,'f +2'))';
Kp3 = TAv.RMS(strcmp(TAv.Group,'Anesthetized')&strcmp(TAv.Freq,'g +3'))';

K_RMS = horzcat(Km3,Km2,Km1,KBF,Kp1,Kp2,Kp3);
K_RMS_5 = horzcat(Km2,Km1,KBF,Kp1,Kp2);

% % MUSC PEAKAMP
Mm3 = TAv.PeakAmp(strcmp(TAv.Group,'Muscimol')&strcmp(TAv.Freq,'a -3'))';
Mm2 = TAv.PeakAmp(strcmp(TAv.Group,'Muscimol')&strcmp(TAv.Freq,'b -2'))';
Mm1 = TAv.PeakAmp(strcmp(TAv.Group,'Muscimol')&strcmp(TAv.Freq,'c -1'))';
MBF = TAv.PeakAmp(strcmp(TAv.Group,'Muscimol')&strcmp(TAv.Freq,'d BF'))';
Mp1 = TAv.PeakAmp(strcmp(TAv.Group,'Muscimol')&strcmp(TAv.Freq,'e +1'))';
Mp2 = TAv.PeakAmp(strcmp(TAv.Group,'Muscimol')&strcmp(TAv.Freq,'f +2'))';
Mp3 = TAv.PeakAmp(strcmp(TAv.Group,'Muscimol')&strcmp(TAv.Freq,'g +3'))';

M_PA = horzcat(Mm3,Mm2,Mm1,MBF,Mp1,Mp2,Mp3);
M_PA_5 = horzcat(Mm2,Mm1,MBF,Mp1,Mp2);
M_PA_5norm = M_PA_5./M_PA_5(:,3);

twoaway = vertcat(Mm2,Mp2);
oneaway = vertcat(Mm1,Mp1);
central = vertcat(MBF,MBF);
M_PA_5mir = horzcat(central,oneaway,twoaway);
M_PA_5mirnorm = M_PA_5mir./M_PA_5mir(:,1);

% % MUSC PEAKLAT
Mm3 = TAv.PeakLat(strcmp(TAv.Group,'Muscimol')&strcmp(TAv.Freq,'a -3'))';
Mm2 = TAv.PeakLat(strcmp(TAv.Group,'Muscimol')&strcmp(TAv.Freq,'b -2'))';
Mm1 = TAv.PeakLat(strcmp(TAv.Group,'Muscimol')&strcmp(TAv.Freq,'c -1'))';
MBF = TAv.PeakLat(strcmp(TAv.Group,'Muscimol')&strcmp(TAv.Freq,'d BF'))';
Mp1 = TAv.PeakLat(strcmp(TAv.Group,'Muscimol')&strcmp(TAv.Freq,'e +1'))';
Mp2 = TAv.PeakLat(strcmp(TAv.Group,'Muscimol')&strcmp(TAv.Freq,'f +2'))';
Mp3 = TAv.PeakLat(strcmp(TAv.Group,'Muscimol')&strcmp(TAv.Freq,'g +3'))';

M_PL = horzcat(Mm3,Mm2,Mm1,MBF,Mp1,Mp2,Mp3);
M_PL_5 = horzcat(Mm2,Mm1,MBF,Mp1,Mp2);

% % MUSC RMS
Mm3 = TAv.RMS(strcmp(TAv.Group,'Muscimol')&strcmp(TAv.Freq,'a -3'))';
Mm2 = TAv.RMS(strcmp(TAv.Group,'Muscimol')&strcmp(TAv.Freq,'b -2'))';
Mm1 = TAv.RMS(strcmp(TAv.Group,'Muscimol')&strcmp(TAv.Freq,'c -1'))';
MBF = TAv.RMS(strcmp(TAv.Group,'Muscimol')&strcmp(TAv.Freq,'d BF'))';
Mp1 = TAv.RMS(strcmp(TAv.Group,'Muscimol')&strcmp(TAv.Freq,'e +1'))';
Mp2 = TAv.RMS(strcmp(TAv.Group,'Muscimol')&strcmp(TAv.Freq,'f +2'))';
Mp3 = TAv.RMS(strcmp(TAv.Group,'Muscimol')&strcmp(TAv.Freq,'g +3'))';

M_RMS = horzcat(Mm3,Mm2,Mm1,MBF,Mp1,Mp2,Mp3);
M_RMS_5 = horzcat(Mm2,Mm1,MBF,Mp1,Mp2);


h=figure; 
K_meanline = nanmean(K_PA_5norm,1);
K_stdbar = nanstd(K_PA_5norm,1);
N_meanline = nanmean(N_PA_5norm,1);
N_stdbar = nanstd(N_PA_5norm,1);
A_meanline = nanmean(A_PA_5norm,1);
A_stdbar = nanstd(A_PA_5norm,1);
M_meanline = nanmean(M_PA_5norm,1);
M_stdbar = nanstd(M_PA_5norm,1);
errorbar(K_meanline, K_stdbar,'LineWidth',1.5)
hold on
errorbar(N_meanline, N_stdbar,'LineWidth',1.5)
errorbar(A_meanline, A_stdbar,'LineWidth',1.5)
errorbar(M_meanline, M_stdbar,'LineWidth',1.5)
legend('Ketamine','AnChronic','Awake','Muscimol')
xlim([0,6])

savefig(h,'Avrec tuning Peak Amp Normalized'); close all;


%% AVREC STATS
disp('*********************** AVREC STATS FULL CURVE ***********************')
disp('************** PEAKAMP ***************')
com_KvA = horzcat(A_PA, K_PA);
com_NvA = horzcat(A_PA, N_PA);
com_KvM = horzcat(K_PA, M_PA);

disp('**Anesthetized vs Awake**')
AvK_AvrecPeakAmp = teg_repeated_measures_ANOVA(com_KvA, levels, varinames);
AvK_AvrecPeakAmpCD = iMakeCohensD(K_PA, A_PA)

disp('**ANChronic vs Awake**')
AvN_AvrecPeakAmp = teg_repeated_measures_ANOVA(com_NvA, levels, varinames);
AvN_AvrecPeakAmpCD = iMakeCohensD(N_PA, A_PA)

disp('**Anesthetized vs Muscimol**')
MvK_AvrecPeakAmp = teg_repeated_measures_ANOVA(com_KvM, levels, varinames);
MvK_AvrecPeakAmpCD = iMakeCohensD(K_PA, M_PA)

disp('************** PEAKLAT ***************')
com_KvsAw = horzcat(A_PL, K_PL);
com_NvsAw = horzcat(A_PL, N_PL);
com_KvM = horzcat(K_PL, M_PL);

disp('**Anesthetized vs Awake**')
AvK_AvrecPeakLat = teg_repeated_measures_ANOVA(com_KvsAw, levels, varinames);
AvK_AvrecPeakLatCD = iMakeCohensD(K_PL, A_PL)

disp('**ANChronic vs Awake**')
AvN_AvrecPeakLat = teg_repeated_measures_ANOVA(com_NvsAw, levels, varinames);
AvN_AvrecPeakLatCD = iMakeCohensD(N_PL, A_PL)

disp('**Anesthetized vs Muscimol**')
MvK_AvrecPeakLat = teg_repeated_measures_ANOVA(com_KvM, levels, varinames);
MvK_AvrecPeakLatCD = iMakeCohensD(K_PL, M_PL)

disp('**************   RMS   ***************')
com_KvsAw = horzcat(A_RMS, K_RMS);
com_NvsAw = horzcat(A_RMS, N_RMS);
com_KvM = horzcat(K_RMS, M_RMS);

disp('**Anesthetized vs Awake**')
AvK_AvrecRMS = teg_repeated_measures_ANOVA(com_KvsAw, levels, varinames);
AvK_AvrecRMSCD = iMakeCohensD(K_RMS, A_RMS)

disp('**ANChronic vs Awake**')
AvN_AvrecRMS = teg_repeated_measures_ANOVA(com_NvsAw, levels, varinames);
AvN_AvrecRMSCD = iMakeCohensD(N_RMS, A_RMS)

disp('**Anesthetized vs Muscimol**')
MvK_AvrecRMS = teg_repeated_measures_ANOVA(com_KvM, levels, varinames);
MvK_AvrecRMSCD = iMakeCohensD(K_RMS, M_RMS)

disp('*********************** AVREC STATS 5 octaves ***********************')
disp('************** PEAKAMP ***************')
com_KvsAw = horzcat(A_PA_5, K_PA_5);
com_NvsAw = horzcat(A_PA_5, N_PA_5);
com_KvM = horzcat(K_PA_5, M_PA_5);

disp('**Anesthetized vs Awake**')
AvK_AvrecPeakAmp_mini = teg_repeated_measures_ANOVA(com_KvsAw, levels5, varinames);
AvK_AvrecPeakAmp_miniCD = iMakeCohensD(K_PA_5, A_PA_5)

disp('**ANChronic vs Awake**')
AvN_AvrecPeakAmp_mini = teg_repeated_measures_ANOVA(com_NvsAw, levels5, varinames);
AvN_AvrecPeakAmp_miniCD = iMakeCohensD(N_PA_5, A_PA_5)

disp('**Anesthetized vs Muscimol**')
MvK_AvrecPeakAmp_mini = teg_repeated_measures_ANOVA(com_KvM, levels5, varinames);
MvK_AvrecPeakAmp_miniCD = iMakeCohensD(K_PA_5, M_PA_5)

disp('************** Normalizied PEAKAMP ***************')
com_KvsAw = horzcat(A_PA_5norm, K_PA_5norm);
com_NvsAw = horzcat(A_PA_5norm, N_PA_5norm);
com_KvM = horzcat(K_PA_5norm, M_PA_5norm);

disp('**Anesthetized vs Awake**')
AvK_AvrecPeakAmpNorm_mini = teg_repeated_measures_ANOVA(com_KvsAw, levels5, varinames);
AvK_AvrecPeakAmpNorm_miniCD = iMakeCohensD(K_PA_5norm, A_PA_5norm)

disp('**ANChronic vs Awake**')
AvN_AvrecPeakAmpNorm_mini = teg_repeated_measures_ANOVA(com_NvsAw, levels5, varinames);
AvN_AvrecPeakAmpNorm_miniCD = iMakeCohensD(N_PA_5norm, A_PA_5norm)

disp('**Anesthetized vs Muscimol**')
MvK_AvrecPeakAmpNorm_mini = teg_repeated_measures_ANOVA(com_KvM, levels5, varinames);
MvK_AvrecPeakAmpNorm_miniCD = iMakeCohensD(K_PA_5norm, M_PA_5norm)

disp('************** Normalizied Mirror PEAKAMP ***************')
com_KvsAw = horzcat(A_PA_5mirnorm, K_PA_5mirnorm);
com_NvsAw = horzcat(A_PA_5mirnorm, N_PA_5mirnorm);
com_KvM = horzcat(K_PA_5mirnorm, M_PA_5mirnorm);

disp('**Anesthetized vs Awake**')
AvK_AvrecPeakAmpMirNorm_mini = teg_repeated_measures_ANOVA(com_KvsAw, levelsmir, varinames);
AvK_AvrecPeakAmpMirNorm_miniCD = iMakeCohensD(K_PA_5mirnorm, A_PA_5mirnorm)

disp('**ANChronic vs Awake**')
AvN_AvrecPeakAmpMirNorm_mini = teg_repeated_measures_ANOVA(com_NvsAw, levelsmir, varinames);
AvN_AvrecPeakAmpMirNorm_miniCD = iMakeCohensD(N_PA_5mirnorm, A_PA_5mirnorm)

disp('**Anesthetized vs Muscimol**')
MvK_AvrecPeakAmpMirNorm_mini = teg_repeated_measures_ANOVA(com_KvM, levelsmir, varinames);
MvK_AvrecPeakAmpMirNorm_miniCD = iMakeCohensD(K_PA_5mirnorm, M_PA_5mirnorm)


disp('************** PEAKLAT ***************')
com_KvsAw = horzcat(A_PL_5, K_PL_5);
com_NvsAw = horzcat(A_PL_5, N_PL_5);
com_KvM = horzcat(K_PL_5, M_PL_5);

disp('**Anesthetized vs Awake**')
AvK_AvrecPeakLat_mini = teg_repeated_measures_ANOVA(com_KvsAw, levels5, varinames);
AvK_AvrecPeakLat_miniCD = iMakeCohensD(K_PL_5, A_PL_5)

disp('**ANChronic vs Awake**')
AvN_AvrecPeakLat_mini = teg_repeated_measures_ANOVA(com_NvsAw, levels5, varinames);
AvN_AvrecPeakLat_miniCD = iMakeCohensD(N_PL_5, A_PL_5)

disp('**Anesthetized vs Muscimol**')
MvK_AvrecPeakLat_mini = teg_repeated_measures_ANOVA(com_KvM, levels5, varinames);
MvK_AvrecPeakLat_miniCD = iMakeCohensD(K_PL_5, M_PL_5)

disp('**************   RMS   ***************')
com_KvsAw = horzcat(A_RMS_5, K_RMS_5);
com_NvsAw = horzcat(A_RMS_5, N_RMS_5);
com_KvM = horzcat(K_RMS_5, M_RMS_5);

disp('**Anesthetized vs Awake**')
AvK_AvrecRMS_mini = teg_repeated_measures_ANOVA(com_KvsAw, levels5, varinames);
AvK_AvrecRMS_miniCD = iMakeCohensD(K_RMS_5, A_RMS_5)

disp('**ANChronic vs Awake**')
AvN_AvrecRMS_mini = teg_repeated_measures_ANOVA(com_NvsAw, levels5, varinames);
AvN_AvrecRMS_miniCD = iMakeCohensD(N_RMS_5, A_RMS_5)

disp('**Anesthetized vs Muscimol**')
MvK_AvrecRMS_mini = teg_repeated_measures_ANOVA(com_KvM, levels5, varinames);
MvK_AvrecRMS_miniCD = iMakeCohensD(K_RMS_5, M_RMS_5)


%% Save
savefile = 'Full Tuning curve Avrec.mat';
save(savefile, 'AvK_AvrecPeakAmp', 'AvK_AvrecPeakAmpCD', 'AvN_AvrecPeakAmp', ...
    'AvN_AvrecPeakAmpCD', 'MvK_AvrecPeakAmp', 'MvK_AvrecPeakAmpCD', 'AvK_AvrecPeakLat', ...
    'AvK_AvrecPeakLatCD', 'AvN_AvrecPeakLat', 'AvN_AvrecPeakLatCD', 'MvK_AvrecPeakLat', ...
    'MvK_AvrecPeakLatCD', 'AvK_AvrecRMS', 'AvK_AvrecRMSCD', 'AvN_AvrecRMS', ...
    'AvN_AvrecRMSCD', 'MvK_AvrecRMS', 'MvK_AvrecRMSCD')

savefile = '5 Octave Tuning curve Avrec.mat';
save(savefile, 'AvK_AvrecPeakAmp_mini', 'AvK_AvrecPeakAmp_miniCD', ...
    'AvK_AvrecPeakAmpNorm_mini', 'AvK_AvrecPeakAmpNorm_miniCD', 'AvK_AvrecPeakAmpMirNorm_mini', ... 
    'AvK_AvrecPeakAmpMirNorm_miniCD','AvN_AvrecPeakAmp_mini', ...
    'AvN_AvrecPeakAmp_miniCD', 'AvN_AvrecPeakAmpNorm_mini', 'AvN_AvrecPeakAmpNorm_miniCD', ...
    'AvN_AvrecPeakAmpMirNorm_mini', 'AvN_AvrecPeakAmpMirNorm_miniCD', ...
    'MvK_AvrecPeakAmp_mini', 'MvK_AvrecPeakAmp_miniCD', 'MvK_AvrecPeakAmpNorm_mini', ...
    'MvK_AvrecPeakAmpNorm_miniCD', 'MvK_AvrecPeakAmpMirNorm_mini', 'MvK_AvrecPeakAmpMirNorm_miniCD', ...
    'AvK_AvrecPeakLat_mini', 'AvK_AvrecPeakLat_miniCD', ...
    'AvN_AvrecPeakLat_mini', 'AvN_AvrecPeakLat_miniCD', 'MvK_AvrecPeakLat_mini', ...
    'MvK_AvrecPeakLat_miniCD', 'AvK_AvrecRMS_mini', 'AvK_AvrecRMS_miniCD', ...
    'AvN_AvrecRMS_mini', 'AvN_AvrecRMS_miniCD', 'MvK_AvrecRMS_mini', 'MvK_AvrecRMS_miniCD')










