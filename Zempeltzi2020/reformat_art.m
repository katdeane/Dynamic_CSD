function [BL, tone, Fs, SWEEP, frqz]=reformat_art(DATA)
%This takes a .mat from the chronic electrode recordings with just a DATA struct and formats it to
%be able to provide necessary input for SingleTrialCSD function

%% set variables %%

BL = 200; %ms before tone played
tone = 200; %duration of tone
Fs = (DATA.SR); %sampling rate

%% Create SWEEP %%

% squeeze out artifact ridden trials
% list of artifact trials:
xout_list = unique(DATA.ART.intchanepoch(:,4));
% get the list that you want still in (without artifact)
triallist = 1:size(DATA.LFP,3);
% for seperating out 10dB, use just up to the current trial number
xout_list = xout_list(xout_list <= length(triallist));
triallist(xout_list) = NaN;
triallist_in = triallist(~isnan(triallist));
% Take just those and overwrite the LFP
data_lfp = DATA.LFP(:,:,triallist_in);
% epochs need to match the LFP list
epochs_code =  [DATA.EPOCHS.code]';
epochlist = nan(1,size(data_lfp,3));
count = 0;

% fill epochlist with correct trials
for iE = 1:length(epochs_code)
    if find(triallist_in == iE)
        count = count +1;
        epochlist(count) = epochs_code(iE);
    end
end

% get epochlist into structure
SWEEP = struct('Header',{zeros},'AD',{zeros}); %otherwise use the EPOCH.code list

if Fs == 1000 %this is already correct
    for i = 1:size(data_lfp,3)
        SWEEP(i).Header = epochlist(i);
        SWEEP(i).AD = data_lfp(:,:,i)';
    end
else
    for i = 1:size(data_lfp,3)
        downscale1 = data_lfp(:,:,i);
        downscale2 = (downscale1(:,1:2:end)+downscale1(:,2:2:end))./2;        
%         try strcmp(DATA.FILE(1:3),'H:\'); %These need to be re-scaled -
%         currently solved this by taking new conversion method from awake
%         group!
%             SWEEP(i).AD = (downscale2')/1000;
%         catch
            SWEEP(i).AD = downscale2'; %replaces what's in field for that iteration with the 2D LFP that corresponds to that code
            SWEEP(i).Header = epochlist(i);
%         end
        
    end
end

% fix frqz length
frqz = [250 500 1000 2000 4000 8000 16000 250 500 1000 2000 4000 8000 16000 250 500 1000 2000 4000 8000 16000]; 
% MAY NEED to create longer list for full list of codes. if so use code
% below to dynamically allow for list size change

numstim = unique([DATA.EPOCHS.code]);
if length(frqz)>length(numstim)
   frqz = frqz(1:length(numstim)); 
end
