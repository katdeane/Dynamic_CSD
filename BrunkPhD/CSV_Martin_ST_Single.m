dbstop if error
clear

h = load('K:\CSD_dynamic_analysis\DATA\Input_HighP_7post_Data.mat');
c = load('K:\CSD_dynamic_analysis\DATA\Input_Control_7post_n7_Data.mat');
y = load('K:\CSD_dynamic_analysis\DATA\Input_YFP_7post_Data.mat');

order={'h' 'c' 'y'};
% order={'h' 'c'};
para={'Full_Single_RMS_AVREC' 'Early_Single_RMS_AVREC' 'Late_Single_RMS_AVREC'...
    'Full_Single_RMS_RELRES' 'Early_Single_RMS_RELRES' 'Late_Single_RMS_RELRES'...
    'Full_Single_RMS_ABSRES' 'Early_Single_RMS_ABSRES' 'Late_Single_RMS_ABSRES'...
    'SingleSinkPeakAmp' 'SingleSinkRMS'};
sinks={'I_IIL' 'IVE' 'VaE' 'VbE' 'VIE' 'InfE' 'VIL'};
Tone = {'BF' 'nearBF' 'nonBF'};


Matrix{1,1} = 'Full RMS AVREC'; 
Matrix{1,2} = 'Early RMS AVREC';
Matrix{1,3} = 'Late RMS AVREC';
Matrix{1,4} = 'Full RMS RelRes';
Matrix{1,5} = 'Early RMS RelRes';
Matrix{1,6} = 'Late RMS RelRes';
Matrix{1,7} = 'Full RMS AbsRes';
Matrix{1,8} = 'Early RMS AbsRes';
Matrix{1,9} = 'Late RMS AbsRes';

Matrix{1,10} = 'SinkPeak L I/II';
Matrix{1,11} = 'SinkPeak E III/IV';
Matrix{1,12} = 'SinkPeak E Va';
Matrix{1,13} = 'SinkPeak E Vb';
Matrix{1,14} = 'SinkPeak E VIa';
Matrix{1,15} = 'SinkPeak E Inf';
Matrix{1,16} = 'SinkPeak L VI';

Matrix{1,17} = 'SinkRMS L I/II';
Matrix{1,18} = 'SinkRMS E III/IV';
Matrix{1,19} = 'SinkRMS E Va';
Matrix{1,20} = 'SinkRMS E Vb';
Matrix{1,21} = 'SinkRMS E VIa';
Matrix{1,22} = 'SinkRMS E Inf';
Matrix{1,23} = 'SinkRMS L VI';

Matrix{1,24} ='Subject';	
Matrix{1,25} ='Group';	
Matrix{1,26} ='Measurement';	
Matrix{1,27} ='Frequency';


Index = 2;
for i1 = 1:length(order)
    IDs = fieldnames(eval([order{i1} '.Data']));
    DATA = eval([order{i1} '.Data']);
    
    for i2 = 1:length(IDs)
        for i3 = 1:length([DATA.(IDs{i2})])
            
            for i4 = 1:length(Tone)            
                 
                for i5 = 1:length(para)
%                     keyboard
%                 BF = DATA(i3).(IDs{i2}).GS_BF;
%                 Frqz = DATA(i3).(IDs{i2}).Frqz;
%                 BF_Pos = find(Frqz == BF); 
                                
%                 BF = nanmax(DATA(i3).(IDs{i2}).(para{i5}));
%                 BF_Pos = find(DATA(i3).(IDs{i2}).(para{i5}) == BF);
                    if i5 == 10
                        S=0;
                    elseif i5 == 11
                        S=6;
                    end
                    if i5 < 10
%                         keyboard
                BF = max(nanmean(DATA(i3).(IDs{i2}).(para{i5}),2)');
                BF_Pos = find(nanmean(DATA(i3).(IDs{i2}).(para{i5}),2)' == BF);
                    if i4 == 1
                        clear ALL
                        ALL = DATA(i3).(IDs{i2}).(para{i5})(BF_Pos,:)';
                        
                        if isempty(ALL), keyboard, end;
%                         if length(ALL) < 50, keyboard; end
                        if length(ALL) < 50 
                            diff2 = 50-length(ALL);
                            ALL = vertcat(ALL, nan(diff2,1));
                        end
                        Matrix(Index:Index+(length(ALL))-1,i5) = cellstr(num2str(ALL));
                        Matrix(Index:Index+(length(ALL))-1,24) =cellstr(IDs{i2});
                        Matrix(Index:Index+(length(ALL))-1,25) =cellstr(order{i1});
                        Matrix(Index:Index+(length(ALL))-1,26) =cellstr(DATA(i3).(IDs{i2}).Condition);
                        Matrix(Index:Index+(length(ALL))-1,27) =cellstr(Tone{i4});
                        
                    elseif i4 == 2

                        clear ALL neg pos
                        try
                            neg(:,1) =DATA(i3).(IDs{i2}).(para{i5})(BF_Pos-1,:)';
                        end
                        try
                            neg(:,2) =DATA(i3).(IDs{i2}).(para{i5})(BF_Pos-2,:)';
                        end
                        try
                            pos(:,1) =DATA(i3).(IDs{i2}).(para{i5})(BF_Pos+1,:)';
                        end
                        try
                            pos(:,2) =DATA(i3).(IDs{i2}).(para{i5})(BF_Pos+2,:)';
                        end
                        if exist('neg','var') == 0
                            neg = [];
                        end
                        if exist('pos','var') == 0
                            pos = [];
                        end
%                         neg = nanmean(neg,2); pos = nanmean(pos,2);
                        ALL = nanmean([neg, pos],2);
                        if isempty(ALL), keyboard, end;
                        if length(ALL) < 50 
                            diff2 = 50-length(ALL);
                            ALL = vertcat(ALL, nan(diff2,1));
                        end
                        Matrix(Index:Index+(length(ALL))-1,i5) = cellstr(num2str(ALL));
                        Matrix(Index:Index+(length(ALL))-1,24) =cellstr(IDs{i2});
                        Matrix(Index:Index+(length(ALL))-1,25) =cellstr(order{i1});
                        Matrix(Index:Index+(length(ALL))-1,26) =cellstr(DATA(i3).(IDs{i2}).Condition);
                        Matrix(Index:Index+(length(ALL))-1,27) =cellstr(Tone{i4});
                        
                    else
                        clear ALL neg pos
                        try
                            neg(:,1) =DATA(i3).(IDs{i2}).(para{i5})(BF_Pos-3,:)';
                        end
                        try
                            neg(:,2) =DATA(i3).(IDs{i2}).(para{i5})(BF_Pos-4,:)';
                        end
                        try
                            pos(:,1) =DATA(i3).(IDs{i2}).(para{i5})(BF_Pos+3,:)';
                        end
                        try
                            pos(:,2) =DATA(i3).(IDs{i2}).(para{i5})(BF_Pos+4,:)';
                        end
                        if exist('neg','var') == 0
                            neg = nan(50,2);
                        end
                        if exist('pos','var') == 0
                            pos = nan(50,2);
                        end
%                         neg = nanmean(neg,2); pos = nanmean(pos,2);
                        if length(neg) < 50
                          diff2 = 50-length(neg);
                          neg = vertcat(neg, nan(diff2, size(neg,2)));
                        end
                        if length(pos) < 50
                          diff2 = 50-length(pos);
                          pos = vertcat(pos, nan(diff2, size(pos,2)));
                        end
                        ALL = nanmean([neg, pos],2);
                        if isempty(ALL), keyboard, end;
                        if length(ALL) < 50, keyboard; end
                        Matrix(Index:Index+(length(ALL))-1,i5) = cellstr(num2str(ALL));
                        Matrix(Index:Index+(length(ALL))-1,24) =cellstr(IDs{i2});
                        Matrix(Index:Index+(length(ALL))-1,25) =cellstr(order{i1});
                        Matrix(Index:Index+(length(ALL))-1,26) =cellstr(DATA(i3).(IDs{i2}).Condition);
                        Matrix(Index:Index+(length(ALL))-1,27) =cellstr(Tone{i4});
                    end 
                    else

                        
                        for i6 = 1:length(sinks)
%                              keyboard

                            entries =[DATA(i3).(IDs{i2}).(para{i5}).(sinks{i6})];
                            
                            tones = length(DATA(i3).(IDs{i2}).Frqz);
                            if length(entries) == tones*50
                            elseif length(entries) < tones*50
                                diff2 = tones*50 - length(entries);
                                entries = [entries nan(1, diff2)];
                            else
                                keyboard
                            end
                            BF = max(nanmean(reshape(entries,tones,[]),2)');
                            BF_Pos = find(nanmean(reshape(entries,tones,[]),2)' == BF);
%                 BF = nanmax([DATA(i3).(IDs{i2}).(para{i5}).(sinks{i6})]);
%                 BF_Pos = find([DATA(i3).(IDs{i2}).(para{i5}).(sinks{i6})] == BF);
                           if i4 == 1
                                clear ALL
                                ALL = DATA(i3).(IDs{i2}).(para{i5})(BF_Pos).(sinks{i6})';
                                if isempty(ALL), keyboard, end;
%                                 if length(ALL) < 50, keyboard; end
                                if length(ALL) < 50 
                                    diff2 = 50-length(ALL);
                                    ALL = vertcat(ALL, nan(diff2,1));
                                end
                                
%                                 keyboard
                                Matrix(Index:Index+(length(ALL))-1,(i5-1+i6+S)) = cellstr(num2str(ALL));
                                Matrix(Index:Index+(length(ALL))-1,24) =cellstr(IDs{i2});
                                Matrix(Index:Index+(length(ALL))-1,25) =cellstr(order{i1});
                                Matrix(Index:Index+(length(ALL))-1,26) =cellstr(DATA(i3).(IDs{i2}).Condition);
                                Matrix(Index:Index+(length(ALL))-1,27) =cellstr(Tone{i4});

                            elseif i4 == 2

                                clear ALL neg pos
                                try
                                    neg(:,1) =DATA(i3).(IDs{i2}).(para{i5})(BF_Pos-1).(sinks{i6})';
                                    if length(neg) <= 50
                                        diff2 = 50 - length(neg) ;
                                        neg = vertcat(neg, nan(1,diff2));
                                    end
                                end
                                try
                                    neg(:,2) =DATA(i3).(IDs{i2}).(para{i5})(BF_Pos-2).(sinks{i6})';
                                end
                                try
                                    pos(:,1) =DATA(i3).(IDs{i2}).(para{i5})(BF_Pos+1).(sinks{i6})';
                                end
                                try
                                    pos(:,2) =DATA(i3).(IDs{i2}).(para{i5})(BF_Pos+2).(sinks{i6})';
                                end
                                if exist('neg','var') == 0
                                    neg = nan(50,2);
                                end
                                if exist('pos','var') == 0
                                    pos = nan(50,2);
                                end
%                                 neg = nanmean(neg,2); pos = nanmean(pos,2);
                                if length(neg)~= length(pos)
                                    if length(neg)> length(pos)
                                      diff =  length(neg)-length(pos);
                                      pos = vertcat(pos ,nan(diff,size(pos,2)));
                                    else
                                      diff = length(pos)-length(neg); 
                                      neg = vertcat(neg ,nan(diff,size(neg,2)));
                                    end
                                    if size(neg,1) == 1
                                        neg = neg';
                                    end
                                    if size(pos,1) == 1
                                        pos = pos';
                                    end    
                                end
                                if size(neg,1) < 50
                                    neg = neg';
                                end
                                if size(pos,1) < 50
                                    pos = pos';
                                end   
                                try
                                ALL = nanmean([neg, pos],2);
                                catch
                                ALL = nanmean([neg', pos'],2);   
                                end
                                if isempty(ALL), keyboard, end;
                                if length(ALL) < 50 
                                 diff2 = 50-length(ALL);
                                 ALL = vertcat(ALL, nan(diff2,1));
                                end
                                Matrix(Index:Index+(length(ALL))-1,(i5-1+i6+S)) = cellstr(num2str(ALL));
                                Matrix(Index:Index+(length(ALL))-1,24) =cellstr(IDs{i2});
                                Matrix(Index:Index+(length(ALL))-1,25) =cellstr(order{i1});
                                Matrix(Index:Index+(length(ALL))-1,26) =cellstr(DATA(i3).(IDs{i2}).Condition);
                                Matrix(Index:Index+(length(ALL))-1,27) =cellstr(Tone{i4});

                           else
%                                 keyboard
                                clear ALL neg pos
                                try
                                    neg(:,1) =DATA(i3).(IDs{i2}).(para{i5})(BF_Pos-3).(sinks{i6})';
                                end
                                try
                                    neg(:,2) =DATA(i3).(IDs{i2}).(para{i5})(BF_Pos-4).(sinks{i6})';
                                end
                                try
                                    pos(:,1) =DATA(i3).(IDs{i2}).(para{i5})(BF_Pos+3).(sinks{i6})';
                                end
                                try
                                    pos(:,2) =DATA(i3).(IDs{i2}).(para{i5})(BF_Pos+4).(sinks{i6})';
                                end
                                if exist('neg','var') == 0
                                    neg = nan(50,2);
                                end
                                if exist('pos','var') == 0
                                    pos = nan(50,2);
                                end
%                                 neg = nanmean(neg,2); pos = nanmean(pos,2);
                                if length(neg)~= length(pos)
                                    if length(neg)> length(pos)
                                      diff =  length(neg)-length(pos);
                                      pos = vertcat(pos ,nan(diff,size(pos,2)))';
                                    else
                                      diff = length(pos)-length(neg); 
                                      neg = vertcat(neg ,nan(diff,size(neg,2)))';
                                    end
                                    if size(neg,1) == 1
                                        neg = neg';
                                    end
                                    if size(pos,1) == 1
                                        pos = pos';
                                    end   
                                end
                                if size(neg,2) ~= size(pos,2)
                                neg = neg';
                                end
                                try
                                ALL = nanmean([neg, pos],2);
                                catch
                                ALL = nanmean([neg', pos],2);
                                end
                                if isempty(ALL), keyboard, end;
                                if length(ALL) < 50 
                                 diff2 = 50-length(ALL);
                                 ALL = vertcat(ALL, nan(diff2,1));
                                end
                                Matrix(Index:Index+(length(ALL))-1,(i5-1+i6+S)) = cellstr(num2str(ALL));
                                Matrix(Index:Index+(length(ALL))-1,24) =cellstr(IDs{i2});
                                Matrix(Index:Index+(length(ALL))-1,25) =cellstr(order{i1});
                                Matrix(Index:Index+(length(ALL))-1,26) =cellstr(DATA(i3).(IDs{i2}).Condition);
                                Matrix(Index:Index+(length(ALL))-1,27) =cellstr(Tone{i4});
                            end  
                        end                         
                    end
                end
                Index = Index+length(ALL);
            end           
        end        
    end    
end