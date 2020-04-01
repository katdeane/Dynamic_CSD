% First run script to extract the relevant data from the huge files
% initially sent to me

%% Init
clear
sourceFolder = 'E:\Andrew Curran\Collaborations\Katrina Spectral Analysis\Dynamic_CSD_Analysis Folders\DATA';

d = dir([sourceFolder '\*Data.mat']);

%% Option 1
% for ii = 1:length(d)
%     sourceFile = [d(ii).folder '\' d(ii).name];
%     load(sourceFile)
%     %Structures within always called 'Data'
%     fn = fieldnames(Data);
%     for jj = 1:length(fn)
%         lfpDat.(d(ii).name(1:end-9)).lfp = Data.(fn{jj}).singletrialLFP;
%         lfpDat.(d(ii).name(1:end-9)).bf = Data.(fn{jj}).GS_BF;
%         %S.(fn{jj})
%         lfpDat.(d(ii).name(1:end-9)).tones = Data.(fn{jj}).Frqz;
%         saveFile = [sourceFolder '\' fn{jj} '_' d(ii).name(1:end-9)];
%         save(saveFile,'lfpDat')
%         clear lfpDat
%     end
% 
% end


%% Other option
for ii = 1:length(d)
    sourceFile = [d(ii).folder '\' d(ii).name];
    load(sourceFile)
    %Structures within always called 'Data'
    fn = fieldnames(Data);
    for jj = 1:length(fn)
        if contains(fn{jj},'10dB')
            newFold = [d(ii).folder '\' fn{jj}(1:end-5) '\']; 
        else
            newFold = [d(ii).folder '\' fn{jj} '\'];
        end
        mkdir(newFold)
        tones = Data.(fn{jj}).Frqz;
        for kk = 1:length(tones)
            nextNewFold = [newFold d(ii).name(1:end-9) '\LFP\'];
            mkdir(nextNewFold)
            lfpDat.lfpTrials = Data.(fn{jj}).singletrialLFP{kk};
            lfpDat.bf = Data.(fn{jj}).GS_BF;
            lfpDat.sinkOn = Data.(fn{jj}).Sinkonset;
            lfpDat.sinkOff = Data.(fn{jj}).Sinkoffset;
            csdDat.csdTrials = Data.(fn{jj}).singletrialCSD{kk};
            csdDat.bf = Data.(fn{jj}).GS_BF;
            csdDat.sinkOn = Data.(fn{jj}).Sinkonset;
            csdDat.sinkOff = Data.(fn{jj}).Sinkoffset;
            if lfpDat.bf==tones(kk)
                toneStr = 'BFtoneNumber';
            else
                toneStr = 'toneNumber';
            end
            saveFile = [nextNewFold toneStr num2str(kk) '_' num2str(tones(kk)) 'Hz'];
            save(saveFile,'lfpDat')
            nextNewFold = [newFold d(ii).name(1:end-9) '\CSD\'];
            mkdir(nextNewFold)
            saveFile = [nextNewFold toneStr num2str(kk) '_' num2str(tones(kk)) 'Hz'];
            save(saveFile,'csdDat')
            clear lfpDat
            clear csdDat
        end
        
    end

end
