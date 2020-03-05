function [DATA] = Groupsorting_single(Data,names,singletrialpara, DimDat,SINKs,NumFreq,mode)
% This function structures the data for group analysis, storing data points
% per parameter per sink. Binning (1 or 0) and Mirroring (1 or 0) is handled here as well as
% dismissing stray data points off BF.

% From this code, tuning curves may be more easily plotted 

for ip = 1:length (singletrialpara)
    for id = 1:DimDat
        Entries = length(SINKs);  
        %% Generating Data per Parameter per Sink per Single Trial (*50)
        for is = 1:Entries
            clear dataholder
            dataholder = cell(length(names),NumFreq); %holding cells to be filled in by data entry
            BF_Pos = (NumFreq+1)/2; %so that the BF is center
            
            for iN = 1:length(names)
                
                if length([Data(id).(names{iN}).(singletrialpara{ip}).(SINKs{is})]) > ((NumFreq+1)/2)*50
                    curLayer = [Data(id).(names{iN}).(singletrialpara{ip}).(SINKs{is})]';
                    curLayer = reshape(curLayer, 2,length(curLayer)/2);
                elseif length([Data(id).(names{iN}).(singletrialpara{ip}).(SINKs{is})]) <= ((NumFreq+1)/2)*50
                    curLayer = [Data(id).(names{iN}).(singletrialpara{ip}).(SINKs{is})]';
                    GSLayer = [Data(id).(names{iN}).(singletrialpara{ip}).(SINKs{is})]';
                    %GSBF =find(Data(id).(names{iN}).Frqz == Data(id).(names{iN}).GS_BF); %already stored based on average data
                    
                    %fing freqavg for getting BF of GS sorted
                    GSAvg = [];
                    GSLaycell = {};
                    for iavg = 1:50:length(GSLayer)
                        gsavg = nanmean(GSLayer(iavg:iavg+49));
                        GSAvg = [GSAvg gsavg];
                    end
                    
                    %find freqavg for getting BF of ST sorted 
                    FreqAvg = [];
                    curLaycell = {};
                    for iavg = 1:50:length(curLayer)
                        freqavg = nanmean(curLayer(iavg:iavg+49));
                        FreqAvg = [FreqAvg freqavg];
                        
                        curLaycell = [curLaycell {(curLayer(iavg:iavg+49))'}];
                    end
                    
                    switch mode
                        case 'GS'
                            CurBF = find(GSAvg == nanmax(GSAvg));
                        case 'ST'
                            CurBF = find(FreqAvg == nanmax(FreqAvg));
                    end
                   
                    if isempty(CurBF)
                        CurBF =1;
                    end
                    
                    dataholder(iN,(BF_Pos+1-CurBF):(BF_Pos+1-CurBF)-1+(length(Data(id).(names{iN}).Frqz))) = curLaycell;
                end
            end
            
            dataholder = cellfun(@fillwnans, dataholder, 'UniformOutput', false);
            
            %% Save and Exit
            
            if isstruct(Data(1).(names{1}).(singletrialpara{ip})) % Check whether input is sink based
                DATA(id).(singletrialpara{ip}).(SINKs{is}) = dataholder;

            else
                DATA(id).(singletrialpara{ip}) = dataholder;

            end
        end
    end
end
end