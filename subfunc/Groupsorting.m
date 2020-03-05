function [DATA,SHIFT] = Groupsorting(Data,names,para, DimDat,SINKs,NumFreq,mode,mirror,nThresh,bin)
% This function structures the data for group analysis, storing data points
% per parameter per sink. Binning (1 or 0) and Mirroring (1 or 0) is handled here as well as
% dismissing stray data points off BF.

% From this code, tuning curves may be more easily plotted 

for i2 = 1:length (para)
    for i3 = 1:DimDat
        if isstruct(Data(1).(names{1}).(para{i2})) % Check wether input is sink based
            Entries = length(SINKs);
        else
            Entries = 1;
        end
        
        %% Generating Data per Parameter/Sink
        for i4 = 1:Entries
            clear dataholder
            dataholder = nan(length(names),NumFreq); %nan holding places to be filled in by data entry
            BF_Pos = (NumFreq+1)/2;
            ALL = [];
            
            for i5 = 1:length(names)
                BFshift_pos =[];
                if ~isempty(Data(i3).(names{i5}))
                    
                    try %if sink-based (change to if-statement pls)
                        if length([Data(i3).(names{i5}).(para{i2}).(SINKs{i4})]) > (NumFreq+1)/2
                            curLayer = [Data(i3).(names{i5}).(para{i2}).(SINKs{i4})];
                            curLayer = reshape(curLayer, 2,length(curLayer)/2);
                        elseif length([Data(i3).(names{i5}).(para{i2}).(SINKs{i4})]) <= (NumFreq+1)/2
                            curLayer = [Data(i3).(names{i5}).(para{i2}).(SINKs{i4})];
                            %                                 if strcmp(para{i2},'tempSinkRMS')
                            %                                 GSBF =find(([Data(i3).(names{i5}).tempSinkRMS.IVE]) == nanmax([Data(i3).(names{i5}).tempSinkRMS.IVE]));
                            %                                 else
                            %                                 GSBF =find(([Data(i3).(names{i5}).SinkRMS.IVE]) == nanmax([Data(i3).(names{i5}).SinkRMS.IVE]));
                            %                                 end
                            GSBF = find([Data(i3).(names{i5}).SinkRMS.IVE] == nanmax([Data(i3).(names{i5}).SinkRMS.IVE])); %find(Data(i3).(names{i5}).Frqz == Data(i3).(names{i5}).GS_BF);
                            switch mode
                                case 'GS'
                                    CurBF = GSBF;
                                case 'ST' %in this case, timing should still be held to a strength parameter
                                    if strcmp('SinkPeakLate',para(i2)) %strcmp('Sinkonset',para(i2)) ||
                                        CurBF =find(([Data(i3).(names{i5}).SinkRMS.(SINKs{i4})]) == nanmax([Data(i3).(names{i5}).SinkRMS.(SINKs{i4})]));
                                    elseif strcmp('Sinkonset',para(i2))
                                        CurBF =find(([Data(i3).(names{i5}).(para{i2}).(SINKs{i4})]) == nanmin([Data(i3).(names{i5}).(para{i2}).(SINKs{i4})]));
                                    else
                                        CurBF =find(([Data(i3).(names{i5}).(para{i2}).(SINKs{i4})]) == nanmax([Data(i3).(names{i5}).(para{i2}).(SINKs{i4})]));
                                    end
                            end
                                    
                            if isempty(GSBF)
                                BFshift_pos = find([Data(i3).(names{i5}).(para{i2}).(SINKs{i4})] == nanmax([Data(i3).(names{i5}).(para{i2}).(SINKs{i4})]))-0;
                            else
                                BFshift_pos = find([Data(i3).(names{i5}).(para{i2}).(SINKs{i4})] == nanmax([Data(i3).(names{i5}).(para{i2}).(SINKs{i4})]))- GSBF;
                            end
                            if isempty(BFshift_pos)
                                BFshift_pos =NaN;
                            end
                            
                            if isempty(CurBF)
                                CurBF =1;
                            end
                            
                            dataholder(i5,(BF_Pos+1-CurBF):(BF_Pos+1-CurBF)-1+length(Data(i3).(names{i5}).Frqz)) =curLayer;
                        end
                    catch
                        if length([Data(i3).(names{i5}).(para{i2})]) == 1
                            curLayer = [Data(i3).(names{i5}).(para{i2})];
                            dataholder(i5,BF_Pos) =curLayer;
                            
                        elseif length([Data(i3).(names{i5}).(para{i2})]) > (NumFreq+1)/2
                            curLayer = Data(i3).(names{i5}).(para{i2})(1:length(Data(i3).(names{i5}).Frqz));
                            if size(curLayer,1) ~= 1
                                curLayer = reshape(curLayer, 2,length(curLayer)/2);
                            end
                            
                        elseif length([Data(i3).(names{i5}).(para{i2})]) <= (NumFreq+1)/2
                            curLayer = [Data(i3).(names{i5}).(para{i2})];
                            
                            if strcmp(para{i2},'tempSinkRMS')
                                GSBF =find(([Data(i3).(names{i5}).tempSinkRMS.IVE]) == nanmax([Data(i3).(names{i5}).tempSinkRMS.IVE]));
                            else
                                GSBF =find(([Data(i3).(names{i5}).SinkRMS.IVE]) == nanmax([Data(i3).(names{i5}).SinkRMS.IVE]));
                            end
                            
                            switch mode
                                case 'GS'
                                    CurBF = GSBF;
                                case 'ST'
                                    CurBF =find(([Data(i3).(names{i5}).(para{i2})]) == nanmax([Data(i3).(names{i5}).(para{i2})]));
                            end
                            
                            if isempty(CurBF)
                                CurBF =1;
                            end
                        end
                        dataholder(i5,(BF_Pos+1-CurBF):(BF_Pos+1-CurBF)-1+length(Data(i3).(names{i5}).Frqz)) = curLayer(1:length(Data(i3).(names{i5}).Frqz));
                    end
                else
                    BFshift_pos = NaN;
                end
                
                ALL = [ALL BFshift_pos];
            end
            
            %% Binning and Mirroring
            clear data_bin_mirror
            
            if mirror == 1 && bin == 0
                data_bin_mirror = nan(size(dataholder));
                
                for i6 = 1:BF_Pos
                    if i6 == 1
                        data_bin_mirror(:,BF_Pos)= dataholder(:,BF_Pos);
                    else
                        data_bin_mirror(:,BF_Pos+i6-1)= nanmean([dataholder(:,BF_Pos-i6+1), dataholder(:,BF_Pos+i6-1)],2);
                    end
                end
                dataholder = data_bin_mirror;
                
            elseif mirror == 1 && bin == 1
                data_bin_mirror = nan(size(dataholder));
                
                for i6 = 1:BF_Pos
                    if i6 == 1
                        data_bin_mirror(:,BF_Pos)= dataholder(:,BF_Pos);
                    else
                        clear neg pos BIN
                        try
                            pos(:,1) = dataholder(:,BF_Pos+i6*2-3);
                            pos(:,2) = dataholder(:,BF_Pos+i6*2-2);
                            neg(:,2)= dataholder(:,BF_Pos-i6*2+3);
                            neg(:,1)= dataholder(:,BF_Pos-i6*2+2);
                            BIN = horzcat(neg,pos);
                            BIN = nanmean(BIN,2);
                            data_bin_mirror(:,BF_Pos+i6-1)= BIN;
                        end
                    end
                end
                dataholder = data_bin_mirror;
                
            elseif mirror == 0 && bin == 1
                data_bin_mirror = nan(size(dataholder));
                for i6 = 1:BF_Pos
                    if i6 == 1
                        data_bin_mirror(:,BF_Pos)= dataholder(:,BF_Pos);
                    else
                        clear neg pos BIN
                        try
                            pos(:,1) = dataholder(:,BF_Pos+i6*2-3);
                            pos(:,2) = dataholder(:,BF_Pos+i6*2-2);
                            neg(:,2)= dataholder(:,BF_Pos-i6*2+3);
                            neg(:,1)= dataholder(:,BF_Pos-i6*2+2);
                            neg = nanmean(neg,2);
                            pos = nanmean(pos,2);
                            
                            data_bin_mirror(:,BF_Pos+i6-1)= pos;
                            data_bin_mirror(:,BF_Pos-i6+1)= neg;
                        end
                        
                    end
                end
                dataholder = data_bin_mirror;  
            end
            
            %% Save and Exit
            
            BinSum = sum(~isnan(dataholder));
            dis = find(BinSum< (nThresh*(length(names))));
            dataholder(:,dis) = nan;
            if isstruct(Data(1).(names{1}).(para{i2})) % Check whether input is sink based
                DATA(i3).(para{i2}).(SINKs{i4}) = dataholder;
                SHIFT(i3).(para{i2}).(SINKs{i4}) = ALL;
            else
                DATA(i3).(para{i2}) = dataholder;
                SHIFT(i3).(para{i2}) = ALL;
            end
        end
    end
end
end

