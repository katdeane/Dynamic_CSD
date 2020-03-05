%% 10ms Bins for Avrec

% this code will eventually take Avrec data, bin it into 10ms chunks, take 
% the rms of that and plot the progression.


%% Load and initialize
clear
warning('OFF'); 

% Change directory to your working folder
if exist('D:\MyCode\Dynamic_CSD_Analysis','dir') == 7
    cd('D:\MyCode\Dynamic_CSD_Analysis');
elseif exist('D:\Dynamic_CSD_Analysis','dir') == 7
    cd('D:\Dynamic_CSD_Analysis');
elseif exist('C:\Users\kedea\Documents\Dynamic_CSD_Analysis','dir') == 7
    cd('C:\Users\kedea\Documents\Dynamic_CSD_Analysis')
end

home = pwd;
addpath(genpath(home));

cd DATA;
load('AnesthetizedPre_Data.mat')
Ketamine = Data; clear Data;
load('Awake10dB_Data.mat')
Awake = Data; clear Data;
load('Muscimol_Data.mat')
Muscimol = Data; clear Data;

chunksize = 10; %10 ms chunks through CSD
% if you don't want the BF, what do you want?
BFdis = 0;      % = 0 if BF, -1 for 1 below BF, etc.

%% Ketamine RMS and avrec/avrec data

% preallocate lists
K_Animal = fieldnames(Ketamine);
RMS_PeakAmp = cell(size(K_Animal));
RMS_PeakLat = cell(size(K_Animal));
RMS_Mean = cell(size(K_Animal));
AV_PeakAmp = cell(size(K_Animal));
AV_PeakLat = cell(size(K_Animal));
AV_Mean = cell(size(K_Animal));
RMS = cell(size(K_Animal));
AV = cell(size(K_Animal));
RMS_PeakAmp_Avg = cell(size(K_Animal));
RMS_PeakLat_Avg = cell(size(K_Animal));
RMS_Mean_Avg = cell(size(K_Animal));
AV_PeakAmp_Avg = cell(size(K_Animal));
AV_PeakLat_Avg = cell(size(K_Animal));
AV_Mean_Avg = cell(size(K_Animal));
RMS_Avg = cell(size(K_Animal));
AV_Avg = cell(size(K_Animal));

% loop through each animal in the group
for iAn = 1:length(K_Animal)
    BF = find(Ketamine.(K_Animal{iAn}).Frqz == Ketamine.(K_Animal{iAn}).GS_BF);
    
    
    % % AVERAGE % %
    % get out AvRec
    if BF+BFdis <= 0
        curAR = NaN([1 600]);
    else
        curAR = Ketamine.(K_Animal{iAn}).AVREC_raw{1,BF+BFdis}';
    end
    
    % preallocate lists to store from loop
    chunking = 1:chunksize:size(curAR,2); %numbers to iterate (every 10ms)
    avrecRMS_Avg = zeros(size(chunking)); %list for RMS output (1:60)
    avrecCSD_Avg = curAR(1:600); %list for avrec output (1:600)
    
    % loop through each trial in steps of 10ms to calc avrec/rms of avrec
    for ichunk = chunking
        chunkAR = curAR(:,ichunk:ichunk+9);
        avrecRMS_Avg(chunking == ichunk) = rms(chunkAR);
    end
    
    [pks, locs, ~, p] = findpeaks(avrecRMS_Avg(20:40));
    if isempty(p)
        RMSpeakamp_Avg(itrial) = NaN;
        RMSpeaklat_Avg(itrial) = NaN;
        RMSmean_Avg(itrial) = NaN;
        AVpeakamp_Avg(itrial) = NaN;
        AVpeaklat_Avg(itrial) = NaN;
        AVmean_Avg(itrial) = NaN;
    else
        [~, maxInd]=max(p); %which peak is most prominent
        RMSpeakamp_Avg = pks(maxInd);
        RMSpeaklat_Avg = locs(maxInd);
        RMSmean_Avg = nanmean(avrecRMS_Avg(20:40));
        
        [pks, locs, ~, p] = findpeaks(avrecCSD_Avg(200:400)*-1);
        [~, maxInd]=max(p); %which peak is most prominent
        AVpeakamp_Avg = pks(maxInd);
        AVpeaklat_Avg = locs(maxInd);
        AVmean_Avg = nanmean(avrecRMS_Avg(20:40));
    end
    
    % % SINGLE TRIAL % %
    % get out list of single trial CSD's
    if BF+BFdis <= 0
        cursnglAR = NaN([600 1 50]);
    else
        cursnglAR = Ketamine.(K_Animal{iAn}).SingleTrial_AVREC_raw{1,BF+BFdis};
    end
    
    % preallocate lists
    trials = 1:size(cursnglAR,3);
    RMSpeakamp = zeros(size(trials));
    RMSpeaklat = zeros(size(trials));
    RMSmean = zeros(size(trials));
    AVpeakamp = zeros(size(trials));
    AVpeaklat = zeros(size(trials));
    AVmean = zeros(size(trials));
    animalRMS = cell(size(trials));
    animalAV = cell(size(trials));
    
    % loop through each trial of the current CSD
    for itrial = trials
        
        trialAR = cursnglAR(:,:,itrial)';
        
        % preallocate lists to store from loop
        chunking = 1:chunksize:size(cursnglAR,1); %numbers to iterate (every 10ms)
        avrecRMS = zeros(size(chunking)); %list for RMS output (1:60)
        avrecCSD = trialAR(1:600); 
        
        % loop through each trial in steps of 10ms to calc avrec/rms of avrec
        for ichunk = chunking
            chunkAR = trialAR(:,ichunk:ichunk+9);
            avrecRMS(chunking == ichunk) = rms(chunkAR);
        end
        
        [pks, locs, ~, p] = findpeaks(avrecRMS(20:40));
        if isempty(p)
            RMSpeakamp(itrial) = NaN;
            RMSpeaklat(itrial) = NaN;
            RMSmean(itrial) = NaN;
            AVpeakamp(itrial) = NaN;
            AVpeaklat(itrial) = NaN;
            AVmean(itrial) = NaN;
        else
            [~, maxInd]=max(p); %which peak is most prominent
            RMSpeakamp(itrial) = pks(maxInd);
            RMSpeaklat(itrial) = locs(maxInd);
            RMSmean(itrial) = nanmean(avrecRMS(20:40));
            
            [pks, locs, ~, p] = findpeaks(avrecCSD(200:400)*-1);
            [~, maxInd]=max(p); %which peak is most prominent
            AVpeakamp(itrial) = pks(maxInd);
            AVpeaklat(itrial) = locs(maxInd);
            AVmean(itrial) = nanmean(avrecRMS(20:40));
        end
        
        animalRMS{itrial} = avrecRMS;
        animalAV{itrial} = avrecCSD;
    end
    

    RMS_PeakAmp{iAn} = RMSpeakamp;
    RMS_PeakLat{iAn} = RMSpeaklat;
    RMS_Mean{iAn} = RMSmean;
    AV_PeakAmp{iAn} = AVpeakamp;
    AV_PeakLat{iAn} = AVpeaklat;
    AV_Mean{iAn} = AVmean;
    RMS{iAn} = animalRMS;
    AV{iAn} = animalAV;
    RMS_PeakAmp_Avg{iAn} = RMSpeakamp_Avg;
    RMS_PeakLat_Avg{iAn} = RMSpeaklat_Avg;
    RMS_Mean_Avg{iAn} = RMSmean_Avg;
    AV_PeakAmp_Avg{iAn} = AVpeakamp_Avg;
    AV_PeakLat_Avg{iAn} = AVpeaklat_Avg;
    AV_Mean_Avg{iAn} = AVmean_Avg;
    RMS_Avg{iAn} = avrecRMS_Avg;
    AV_Avg{iAn} = avrecCSD_Avg;
end

Animal = K_Animal;
Ket_T = table(Animal, RMS, AV, RMS_PeakAmp, RMS_PeakLat, RMS_Mean, AV_PeakAmp,...
    AV_PeakLat, AV_Mean, RMS_Avg, AV_Avg, RMS_PeakAmp_Avg, RMS_PeakLat_Avg, ...
    RMS_Mean_Avg, AV_PeakAmp_Avg, AV_PeakLat_Avg, AV_Mean_Avg);
  
%% Awake RMS and avrec/avrec data

% preallocate lists
A_Animal = fieldnames(Awake);
RMS_PeakAmp = cell(size(A_Animal));
RMS_PeakLat = cell(size(A_Animal));
RMS_Mean = cell(size(A_Animal));
AV_PeakAmp = cell(size(A_Animal));
AV_PeakLat = cell(size(A_Animal));
AV_Mean = cell(size(A_Animal));
RMS = cell(size(A_Animal));
AV = cell(size(A_Animal));
RMS_PeakAmp_Avg = cell(size(A_Animal));
RMS_PeakLat_Avg = cell(size(A_Animal));
RMS_Mean_Avg = cell(size(A_Animal));
AV_PeakAmp_Avg = cell(size(A_Animal));
AV_PeakLat_Avg = cell(size(A_Animal));
AV_Mean_Avg = cell(size(A_Animal));
RMS_Avg = cell(size(A_Animal));
AV_Avg = cell(size(A_Animal));

% loop through each animal in the group
for iAn = 1:length(A_Animal)
    
    % get out list of single trial CSD's
    BF = find(Awake.(A_Animal{iAn}).Frqz == Awake.(A_Animal{iAn}).GS_BF);   
    
    % % AVERAGE % %
    % get out AvRec
    if BF+BFdis <= 0
        curAR = NaN([1 600]);
    else
        curAR = Awake.(A_Animal{iAn}).AVREC_raw{1,BF+BFdis}';
    end
    
    % preallocate lists to store from loop
    chunking = 1:chunksize:size(curAR,2); %numbers to iterate (every 10ms)
    avrecRMS_Avg = zeros(size(chunking)); %list for RMS output (1:60)
    avrecCSD_Avg = curAR(1:600); %list for avrec output (1:600)
    
    % loop through each trial in steps of 10ms to calc avrec/rms of avrec
    for ichunk = chunking
        chunkAR = curAR(:,ichunk:ichunk+9);
        avrecRMS_Avg(chunking == ichunk) = rms(chunkAR);
    end
    
    [pks, locs, ~, p] = findpeaks(avrecRMS_Avg(20:40));
    if isempty(p)
        RMSpeakamp_Avg(itrial) = NaN;
        RMSpeaklat_Avg(itrial) = NaN;
        RMSmean_Avg(itrial) = NaN;
        AVpeakamp_Avg(itrial) = NaN;
        AVpeaklat_Avg(itrial) = NaN;
        AVmean_Avg(itrial) = NaN;
    else
        [~, maxInd]=max(p); %which peak is most prominent
        RMSpeakamp_Avg = pks(maxInd);
        RMSpeaklat_Avg = locs(maxInd);
        RMSmean_Avg = nanmean(avrecRMS_Avg(20:40));
        
        [pks, locs, ~, p] = findpeaks(avrecCSD_Avg(200:400)*-1);
        [~, maxInd]=max(p); %which peak is most prominent
        AVpeakamp_Avg = pks(maxInd);
        AVpeaklat_Avg = locs(maxInd);
        AVmean_Avg = nanmean(avrecRMS_Avg(20:40));
    end
    
    % % SINGLE TRIAL % %
    % get out list of single trial CSD's
    if BF+BFdis <= 0
        cursnglAR = NaN([600 1 50]);
    else
        cursnglAR = Awake.(A_Animal{iAn}).SingleTrial_AVREC_raw{1,BF+BFdis};
    end
    
    % preallocate lists
    trials = 1:size(cursnglAR,3);
    RMSpeakamp = zeros(size(trials));
    RMSpeaklat = zeros(size(trials));
    RMSmean = zeros(size(trials));
    AVpeakamp = zeros(size(trials));
    AVpeaklat = zeros(size(trials));
    AVmean = zeros(size(trials));
    animalRMS = cell(size(trials));
    animalAV = cell(size(trials));
    
    % loop through each trial of the current CSD
    for itrial = trials
        
        trialAR = cursnglAR(:,:,itrial)';
        
        % preallocate lists to store from loop
        chunking = 1:chunksize:size(cursnglAR,1); %numbers to iterate (every 10ms)
        avrecRMS = zeros(size(chunking)); %list for RMS output (1:60)
        avrecCSD = trialAR(1:600); 
        
        % loop through each trial in steps of 10ms to calc avrec/rms of avrec
        for ichunk = chunking
            chunkAR = trialAR(:,ichunk:ichunk+9);
            avrecRMS(chunking == ichunk) = rms(chunkAR);
        end
        
        [pks, locs, ~, p] = findpeaks(avrecRMS(20:40));
        if isempty(p)
            RMSpeakamp(itrial) = NaN;
            RMSpeaklat(itrial) = NaN;
            RMSmean(itrial) = NaN;
            AVpeakamp(itrial) = NaN;
            AVpeaklat(itrial) = NaN;
            AVmean(itrial) = NaN;
        else
            [~, maxInd]=max(p); %which peak is most prominent
            RMSpeakamp(itrial) = pks(maxInd);
            RMSpeaklat(itrial) = locs(maxInd);
            RMSmean(itrial) = nanmean(avrecRMS(20:40));
            
            [pks, locs, ~, p] = findpeaks(avrecCSD(200:400)*-1);
            [~, maxInd]=max(p); %which peak is most prominent
            AVpeakamp(itrial) = pks(maxInd);
            AVpeaklat(itrial) = locs(maxInd);
            AVmean(itrial) = nanmean(avrecRMS(20:40));
        end
        animalRMS{itrial} = avrecRMS;
        animalAV{itrial} = avrecCSD;
        
    end
    

    RMS_PeakAmp{iAn} = RMSpeakamp;
    RMS_PeakLat{iAn} = RMSpeaklat;
    RMS_Mean{iAn} = RMSmean;
    AV_PeakAmp{iAn} = AVpeakamp;
    AV_PeakLat{iAn} = AVpeaklat;
    AV_Mean{iAn} = AVmean;
    RMS{iAn} = animalRMS;
    AV{iAn} = animalAV;
    RMS_PeakAmp_Avg{iAn} = RMSpeakamp_Avg;
    RMS_PeakLat_Avg{iAn} = RMSpeaklat_Avg;
    RMS_Mean_Avg{iAn} = RMSmean_Avg;
    AV_PeakAmp_Avg{iAn} = AVpeakamp_Avg;
    AV_PeakLat_Avg{iAn} = AVpeaklat_Avg;
    AV_Mean_Avg{iAn} = AVmean_Avg;
    RMS_Avg{iAn} = avrecRMS_Avg;
    AV_Avg{iAn} = avrecCSD_Avg;
end
Animal = A_Animal;
Aw_T = table(Animal, RMS, AV, RMS_PeakAmp, RMS_PeakLat, RMS_Mean, AV_PeakAmp,...
    AV_PeakLat, AV_Mean, RMS_Avg, AV_Avg, RMS_PeakAmp_Avg, RMS_PeakLat_Avg, ...
    RMS_Mean_Avg, AV_PeakAmp_Avg, AV_PeakLat_Avg, AV_Mean_Avg);

%% Muscimol RMS and avrec/avrec data

% preallocate lists
M_Animal = fieldnames(Muscimol);
RMS_PeakAmp = cell(size(M_Animal));
RMS_PeakLat = cell(size(M_Animal));
RMS_Mean = cell(size(M_Animal));
AV_PeakAmp = cell(size(M_Animal));
AV_PeakLat = cell(size(M_Animal));
AV_Mean = cell(size(M_Animal));
RMS = cell(size(M_Animal));
AV = cell(size(M_Animal));
RMS_PeakAmp_Avg = cell(size(M_Animal));
RMS_PeakLat_Avg = cell(size(M_Animal));
RMS_Mean_Avg = cell(size(M_Animal));
AV_PeakAmp_Avg = cell(size(M_Animal));
AV_PeakLat_Avg = cell(size(M_Animal));
AV_Mean_Avg = cell(size(M_Animal));
RMS_Avg = cell(size(M_Animal));
AV_Avg = cell(size(M_Animal));

% loop through each animal in the group
for iAn = 1:length(M_Animal)
    
    % get out list of single trial CSD's
    BF = find(Muscimol.(M_Animal{iAn}).Frqz == Muscimol.(M_Animal{iAn}).GS_BF);   
    
    % % AVERAGE % %
    % get out AvRec
    if BF+BFdis <= 0
        curAR = NaN([1 600]);
    else
        curAR = Muscimol.(M_Animal{iAn}).AVREC_raw{1,BF+BFdis}';
    end
    
    % preallocate lists to store from loop
    chunking = 1:chunksize:size(curAR,2); %numbers to iterate (every 10ms)
    avrecRMS_Avg = zeros(size(chunking)); %list for RMS output (1:60)
    avrecCSD_Avg = curAR(1:600); %list for avrec output (1:600)
    
    % loop through each trial in steps of 10ms to calc avrec/rms of avrec
    for ichunk = chunking
        chunkAR = curAR(:,ichunk:ichunk+9);
        avrecRMS_Avg(chunking == ichunk) = rms(chunkAR);
    end
    
    [pks, locs, ~, p] = findpeaks(avrecRMS_Avg(20:40));
    if isempty(p)
        RMSpeakamp_Avg(itrial) = NaN;
        RMSpeaklat_Avg(itrial) = NaN;
        RMSmean_Avg(itrial) = NaN;
        AVpeakamp_Avg(itrial) = NaN;
        AVpeaklat_Avg(itrial) = NaN;
        AVmean_Avg(itrial) = NaN;
    else
        [~, maxInd]=max(p); %which peak is most prominent
        RMSpeakamp_Avg = pks(maxInd);
        RMSpeaklat_Avg = locs(maxInd);
        RMSmean_Avg = nanmean(avrecRMS_Avg(20:40));
        
        [pks, locs, ~, p] = findpeaks(avrecCSD_Avg(200:400)*-1);
        [~, maxInd]=max(p); %which peak is most prominent
        AVpeakamp_Avg = pks(maxInd);
        AVpeaklat_Avg = locs(maxInd);
        AVmean_Avg = nanmean(avrecRMS_Avg(20:40));
    end
    
    % % SINGLE TRIAL % %
    % get out list of single trial CSD's
    if BF+BFdis <= 0
        cursnglAR = NaN([600 1 50]);
    else
        cursnglAR = Muscimol.(M_Animal{iAn}).SingleTrial_AVREC_raw{1,BF+BFdis};
    end
    
    % preallocate lists
    trials = 1:size(cursnglAR,3);
    RMSpeakamp = zeros(size(trials));
    RMSpeaklat = zeros(size(trials));
    RMSmean = zeros(size(trials));
    AVpeakamp = zeros(size(trials));
    AVpeaklat = zeros(size(trials));
    AVmean = zeros(size(trials));
    animalRMS = cell(size(trials));
    animalAV = cell(size(trials));
    
    % loop through each trial of the current CSD
    for itrial = trials
        
        trialAR = cursnglAR(:,:,itrial)';
        
        % preallocate lists to store from loop
        chunking = 1:chunksize:size(cursnglAR,1); %numbers to iterate (every 10ms)
        avrecRMS = zeros(size(chunking)); %list for RMS output (1:60)
        avrecCSD = trialAR(1:600); 
        
        % loop through each trial in steps of 10ms to calc avrec/rms of avrec
        for ichunk = chunking
            chunkAR = trialAR(:,ichunk:ichunk+9);
            avrecRMS(chunking == ichunk) = rms(chunkAR);
        end
        
        [pks, locs, ~, p] = findpeaks(avrecRMS(20:40));
        if isempty(p)
            RMSpeakamp(itrial) = NaN;
            RMSpeaklat(itrial) = NaN;
            RMSmean(itrial) = NaN;
            AVpeakamp(itrial) = NaN;
            AVpeaklat(itrial) = NaN;
            AVmean(itrial) = NaN;
        else
            [~, maxInd]=max(p); %which peak is most prominent
            RMSpeakamp(itrial) = pks(maxInd);
            RMSpeaklat(itrial) = locs(maxInd);
            RMSmean(itrial) = nanmean(avrecRMS(20:40));
            
            [pks, locs, ~, p] = findpeaks(avrecCSD(200:400)*-1);
            [~, maxInd]=max(p); %which peak is most prominent
            AVpeakamp(itrial) = pks(maxInd);
            AVpeaklat(itrial) = locs(maxInd);
            AVmean(itrial) = nanmean(avrecRMS(20:40));
        end
        
        animalRMS{itrial} = avrecRMS;
        animalAV{itrial} = avrecCSD;
    end
    

    RMS_PeakAmp{iAn} = RMSpeakamp;
    RMS_PeakLat{iAn} = RMSpeaklat;
    RMS_Mean{iAn} = RMSmean;
    AV_PeakAmp{iAn} = AVpeakamp;
    AV_PeakLat{iAn} = AVpeaklat;
    AV_Mean{iAn} = AVmean;
    RMS{iAn} = animalRMS;
    AV{iAn} = animalAV;
    RMS_PeakAmp_Avg{iAn} = RMSpeakamp_Avg;
    RMS_PeakLat_Avg{iAn} = RMSpeaklat_Avg;
    RMS_Mean_Avg{iAn} = RMSmean_Avg;
    AV_PeakAmp_Avg{iAn} = AVpeakamp_Avg;
    AV_PeakLat_Avg{iAn} = AVpeaklat_Avg;
    AV_Mean_Avg{iAn} = AVmean_Avg;
    RMS_Avg{iAn} = avrecRMS_Avg;
    AV_Avg{iAn} = avrecCSD_Avg;
end

Animal = M_Animal;
Mus_T = table(Animal, RMS, AV, RMS_PeakAmp, RMS_PeakLat, RMS_Mean, AV_PeakAmp,...
    AV_PeakLat, AV_Mean, RMS_Avg, AV_Avg, RMS_PeakAmp_Avg, RMS_PeakLat_Avg, ...
    RMS_Mean_Avg, AV_PeakAmp_Avg, AV_PeakLat_Avg, AV_Mean_Avg);

%% concatonate full table
Klist = {'Ket' 'Ket' 'Ket' 'Ket' 'Ket' 'Ket' 'Ket' 'Ket' 'Ket' 'Ket' 'Ket'}'; 
Kextra = {'Ket' 'Ket' 'Ket' 'Ket' 'Ket'};
Alist = {'Aw' 'Aw' 'Aw' 'Aw' 'Aw' 'Aw' 'Aw' 'Aw' 'Aw'}'; 
Aextra = {'Aw' 'Aw' 'Aw' 'Aw' 'Aw' 'Aw'};
Mlist = {'Mus' 'Mus' 'Mus' 'Mus' 'Mus' 'Mus' 'Mus' 'Mus' 'Mus' 'Mus' 'Mus'}'; 
Mextra = {'Mus' 'Mus' 'Mus' 'Mus' 'Mus'};
Groups = vertcat(Klist, Alist, Mlist);

% correct group sizes. There has to be a better way to do this...
Kfor60 = horzcat(Klist', Klist', Klist', Klist', Klist', Kextra);
Afor60 = horzcat(Alist', Alist', Alist', Alist', Alist', Alist', Aextra);
Mfor60 = horzcat(Mlist', Mlist', Mlist', Mlist', Mlist', Mextra);
Kfor600 = horzcat(Kfor60,Kfor60,Kfor60,Kfor60,Kfor60,Kfor60,Kfor60,Kfor60,Kfor60,Kfor60);
Afor600 = horzcat(Afor60,Afor60,Afor60,Afor60,Afor60,Afor60,Afor60,Afor60,Afor60,Afor60);
Mfor600 = horzcat(Mfor60,Mfor60,Mfor60,Mfor60,Mfor60,Mfor60,Mfor60,Mfor60,Mfor60,Mfor60);
Kfor3k = horzcat(Kfor600, Kfor600, Kfor600, Kfor600, Kfor600);
Afor3k = horzcat(Afor600, Afor600, Afor600, Afor600, Afor600);
Mfor3k = horzcat(Mfor600, Mfor600, Mfor600, Mfor600, Mfor600);
Kfor30k = horzcat(Kfor3k, Kfor3k, Kfor3k, Kfor3k, Kfor3k, Kfor3k, Kfor3k, Kfor3k, Kfor3k, Kfor3k);
Afor30k = horzcat(Afor3k, Afor3k, Afor3k, Afor3k, Afor3k, Afor3k, Afor3k, Afor3k, Afor3k, Afor3k);
Mfor30k = horzcat(Mfor3k, Mfor3k, Mfor3k, Mfor3k, Mfor3k, Mfor3k, Mfor3k, Mfor3k, Mfor3k, Mfor3k);

K_Groups60 = {Kfor60 Kfor60 Kfor60 Kfor60 Kfor60 Kfor60 Kfor60 Kfor60 Kfor60 Kfor60 Kfor60}';
A_Groups60 = {Afor60 Afor60 Afor60 Afor60 Afor60 Afor60 Afor60 Afor60 Afor60}';
M_Groups60 = {Mfor60 Mfor60 Mfor60 Mfor60 Mfor60 Mfor60 Mfor60 Mfor60 Mfor60 Mfor60 Mfor60}';
K_Groups600 = {Kfor600 Kfor600 Kfor600 Kfor600 Kfor600 Kfor600 Kfor600 Kfor600 Kfor600 Kfor600 Kfor600}';
A_Groups600 = {Afor600 Afor600 Afor600 Afor600 Afor600 Afor600 Afor600 Afor600 Afor600}';
M_Groups600 = {Mfor600 Mfor600 Mfor600 Mfor600 Mfor600 Mfor600 Mfor600 Mfor600 Mfor600 Mfor600 Mfor600}';
K_Groups3k = {Kfor3k Kfor3k Kfor3k Kfor3k Kfor3k Kfor3k Kfor3k Kfor3k Kfor3k Kfor3k Kfor3k}';
A_Groups3k = {Afor3k Afor3k Afor3k Afor3k Afor3k Afor3k Afor3k Afor3k Afor3k}';
M_Groups3k = {Mfor3k Mfor3k Mfor3k Mfor3k Mfor3k Mfor3k Mfor3k Mfor3k Mfor3k Mfor3k Mfor3k}';
K_Groups30k = {Kfor30k Kfor30k Kfor30k Kfor30k Kfor30k Kfor30k Kfor30k Kfor30k Kfor30k Kfor30k Kfor30k}';
A_Groups30k = {Afor30k Afor30k Afor30k Afor30k Afor30k Afor30k Afor30k Afor30k Afor30k}';
M_Groups30k = {Mfor30k Mfor30k Mfor30k Mfor30k Mfor30k Mfor30k Mfor30k Mfor30k Mfor30k Mfor30k Mfor30k}';

Groups60 = vertcat(K_Groups60, A_Groups60, M_Groups60); 
Groups600 = vertcat(K_Groups600, A_Groups600, M_Groups600);
Groups3k = vertcat(K_Groups3k, A_Groups3k, M_Groups3k);
Groups30k = vertcat(K_Groups30k, A_Groups30k, M_Groups30k);
T60 = {1:60}; Time60 = repmat(T60,[31 1]);
T600 = {1:600}; Time600 = repmat(T600,[31 1]);
dT60 = 1:60; T3k = {repmat(dT60, [1 50])}; Time3k = repmat(T3k, [31 1]);
dT600 = 1:600; T30k = {repmat(dT600, [1 50])}; Time30k = repmat(T30k, [31 1]);

T = vertcat(Ket_T,Aw_T,Mus_T);
T.Group = Groups;
T.Groups60 = Groups60;
T.Groups600 = Groups600;
T.Groups3k = Groups3k;
T.Groups30k = Groups30k;
T.Time60 = Time60;
T.Time600 = Time600;
T.Time3k = Time3k;
T.Time30k = Time30k;

% sizing to fit the cookie cutter
for i = 1:size(T,1)
    
    for ii = 1:50
        
        if size(T.RMS{i,1},2)<50 %these animals have 49 trials and need padding to fit labels
            T.RMS{i,1}{1,ii} = NaN([1,60]);
            T.AV{i,1}{1,ii} = NaN([1,600]);
        end
        
    T.RMS{i,1}{1,ii} = T.RMS{i,1}{1,ii}(1:60);
    T.AV{i,1}{1,ii} = T.AV{i,1}{1,ii}(1:600);
    end
    
    RMS = horzcat(T.RMS{i,1}{1,:});  
    T.RMS{i,1} = RMS;
    AV = horzcat(T.AV{i,1}{1,:});  
    T.AV{i,1} = AV;
    
    % cutting the averages
    T.RMS_Avg{i,1} = T.RMS_Avg{i,1}(1:60);
    T.AV_Avg{i,1} = T.AV_Avg{i,1}(1:600);
end

%% PLOT RMS progression of AVG CSDs
cd(home); cd figs; mkdir('AvRec_BINS'); cd AvRec_BINS 

plotRMS_Avg = horzcat(T.RMS_Avg{:,1});
plotT60 = horzcat(T.Time60{:,1});
plotG60 = horzcat(T.Groups60{:,1});

clear g 

g = gramm('x', plotT60, 'y', plotRMS_Avg, 'color', plotG60);
% g.geom_point()
g.stat_summary() %'type', 'sem', 'geom', 'bar'
g.set_title(['AVG RMS of Avrec from 10ms bins ' num2str(BFdis)]);
g.set_color_options('map','matlab');
g.set_names('x','Binned Time','y','RMS AvRec','color','Group');
g.axe_property('YLim',[0 0.002]);
g.draw();
g.export('file_name',['AVG RMS of Avrec from 10ms bins ' num2str(BFdis)], 'file_type','png');
close all

%% PLOT Avrec progression of AVG CSDs

plotAV_Avg = horzcat(T.AV_Avg{:,1});
plotT600 = horzcat(T.Time600{:,1});
plotG600 = horzcat(T.Groups600{:,1});

clear g 

g = gramm('x', plotT600, 'y', plotAV_Avg, 'color', plotG600);
g.stat_summary() %'type', 'sem', 'geom', 'line'
g.set_title(['AVG Avrec from 10ms bins ' num2str(BFdis)]);
g.set_color_options('map','matlab');
g.set_names('x','Time','y','AvRec','color','Group');
g.axe_property('YLim',[0 0.002]);
g.draw();
g.export('file_name',['AVG Avrec from 10ms bins ' num2str(BFdis)], 'file_type','png');
close all

%% PLOT RMS progression of SINGLE CSDs

plotRMS = horzcat(T.RMS{:,1});
plotT3k = horzcat(T.Time3k{:,1});
plotG3k = horzcat(T.Groups3k{:,1});

clear g 

g = gramm('x', plotT3k, 'y', plotRMS, 'color', plotG3k);
g.stat_summary() %'type', 'sem', 'geom', 'bar'
g.set_title(['Single RMS of Avrec from 10ms bins ' num2str(BFdis)]);
g.set_color_options('map','matlab');
g.set_names('x','Binned Time','y','RMS AvRec','color','Group');
g.axe_property('YLim',[0 0.002]);
g.draw();
g.export('file_name',['Single RMS of Avrec from 10ms bins ' num2str(BFdis)], 'file_type','png');
close all

%% PLOT Avrec progression of SINGLE CSDs

plotAV = horzcat(T.AV{:,1});
plotT30k = horzcat(T.Time30k{:,1});
plotG30k = horzcat(T.Groups30k{:,1});

clear g 

g = gramm('x', plotT30k, 'y', plotAV, 'color', plotG30k);
g.stat_summary() %'type', 'sem', 'geom', 'bar'
g.set_title(['Single Avrec from 10ms bins ' num2str(BFdis)]);
g.set_color_options('map','matlab');
g.set_names('x','Time','y','AvRec','color','Group');
g.axe_property('YLim',[0 0.002]);
g.draw();
g.export('file_name',['Single Avrec from 10ms bins ' num2str(BFdis)], 'file_type','png');
close all

%%
% SINGLE TRIAL STUFF
clear g
g = gramm('x',T.RMS_PeakAmp,'y',T.RMS_Mean,'color',T.Animal, 'column', T.Group);
g.geom_point();
g.set_title(['sngle RMS Peak Amplitude vs Mean_' num2str(BFdis)]);
g.set_names('x','Highest Peak Amplitude','y','Mean during tone','color','Animal');
g.draw();
g.export('file_name',['sngle RMS Peak Amplitude vs Mean_' num2str(BFdis)], 'file_type','png');
close all

% clear g %want to just visualize mean but I'm really loosing it
% g = gramm('x',T.RMS_Mean,'color',T.Animal, 'column', T.Group);
% g.geom_jitter();
% g.set_title(['sngle RMS Peak Amplitude vs Mean_' num2str(BFdis)]);
% g.set_names('x','Highest Peak Amplitude','y','Mean during tone','color','Animal');
% g.draw();
% g.export('file_name',['sngle RMS Peak Amplitude vs Mean_' num2str(BFdis)], 'file_type','png');
% close all

clear g
g = gramm('x',T.RMS_PeakLat,'y',T.RMS_PeakAmp,'color',T.Animal, 'column', T.Group);
g.geom_point();
g.set_title(['sngle RMS Peak Latency vs Amplitude_' num2str(BFdis)]);
g.draw();
g.export('file_name',['sngle RMS Peak Latency vs Amplitude_' num2str(BFdis)], 'file_type','png');
close all

clear g
g = gramm('x',T.AV_PeakLat,'y',T.AV_PeakAmp,'color',T.Animal, 'column', T.Group);
g.geom_point();
g.set_title(['sngle AV Peak Latency vs Amplitude_' num2str(BFdis)]);
g.draw();
g.export('file_name',['sngle AV Peak Latency vs Amplitude_' num2str(BFdis)], 'file_type','png');
close all


clear g
g = gramm('x',T.AV_PeakLat,'y',T.RMS_PeakLat,'color',T.Animal, 'column', T.Group); 
g.geom_point();
g.set_title(['sngle AV vs RMS Peak Latency_' num2str(BFdis)]);
g.draw();
g.export('file_name',['sngle AV vs RMS Peak Latency_' num2str(BFdis)], 'file_type','png');
close all


clear g
g = gramm('x',T.AV_PeakAmp,'y',T.RMS_PeakAmp,'color',T.Animal, 'column', T.Group); 
g.geom_point();
g.set_title(['sngle AV vs RMS Peak Amplitude_' num2str(BFdis)]);
g.draw();
g.export('file_name',['sngle AV vs RMS Peak Amplitude_' num2str(BFdis)], 'file_type','png');
close all

% % AVERAGE STUFF
clear g
g = gramm('x',T.RMS_PeakAmp_Avg,'y',T.RMS_Mean_Avg,'color',T.Animal, 'column', T.Group);
g.geom_point();
g.set_title(['AVG RMS Peak Amplitude vs Mean_' num2str(BFdis)]);
g.set_names('x','Highest Peak Amplitude','y','Mean during tone','color','Animal');
g.draw();
g.export('file_name',['AVG RMS Peak Amplitude vs Mean_' num2str(BFdis)], 'file_type','png');
close all


clear g
g = gramm('x',T.RMS_PeakLat_Avg,'y',T.RMS_PeakAmp_Avg,'color',T.Animal, 'column', T.Group); 
g.geom_point();
g.set_title(['RMS Peak Latency vs Amplitude_' num2str(BFdis)]);
g.draw();
g.export('file_name',['RMS Peak Latency vs Amplitude_' num2str(BFdis)], 'file_type','png');
close all


clear g
g = gramm('x',T.AV_PeakLat_Avg,'y',T.AV_PeakAmp_Avg,'color',T.Animal, 'column', T.Group); 
g.geom_point();
g.set_title(['AV Peak Latency vs Amplitude_' num2str(BFdis)]);
g.draw();
g.export('file_name',['AV Peak Latency vs Amplitude_' num2str(BFdis)], 'file_type','png');
close all


clear g
g = gramm('x',T.AV_PeakLat_Avg,'y',T.RMS_PeakLat_Avg,'color',T.Animal, 'column', T.Group);
g.geom_point();
g.set_title(['AV vs RMS Peak Latency_' num2str(BFdis)]);
g.draw();
g.export('file_name',['AV vs RMS Peak Latency_' num2str(BFdis)], 'file_type','png');
close all


clear g
g = gramm('x',T.AV_PeakAmp_Avg,'y',T.RMS_PeakAmp_Avg,'color',T.Animal, 'column', T.Group); 
g.geom_point();
g.set_title(['AV vs RMS Peak Amplitude_' num2str(BFdis)]);
g.draw();
g.export('file_name',['AV vs RMS Peak Amplitude_' num2str(BFdis)], 'file_type','png');
close all

%% Stats of all the things
figlist = {'Musc_Box','Musc_Stat','Aw_Box','Aw_Stat'};
Fluff = {NaN([1,50])};

ARmsMean = T.RMS_Mean(strcmp('Aw', T.Group));
ARmsMean = vertcat(ARmsMean, Fluff, Fluff);
ARmsMean = horzcat(ARmsMean{:,1}, NaN, NaN)'; %bah. throws an error. can't be bothered atm

KRmsMean = T.RMS_Mean(strcmp('Ket', T.Group));
KRmsMean = horzcat(KRmsMean{:,1})';

MRmsMean = T.RMS_Mean(strcmp('Mus', T.Group));
MRmsMean = horzcat(MRmsMean{:,1})';

AvK_RmsMean = horzcat(ARmsMean, KRmsMean);
MvK_RmsMean = horzcat(MRmsMean, KRmsMean);

[AP,ASTATS] = vartestn(AvK_RmsMean, 'testtype','BrownForsythe');
[MP,MSTATS] = vartestn(MvK_RmsMean, 'testtype','BrownForsythe');


for ifig = 1:4
h = gcf;
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(h, ['RMS Mean Single ' figlist{ifig} '_BF'])
saveas(gcf, ['RMS Mean Single ' figlist{ifig} '_BF.pdf'])
close (h)
end