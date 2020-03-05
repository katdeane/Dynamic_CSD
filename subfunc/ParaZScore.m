function [Zscrd_DATA, DATA_BU] = ParaZScore(DATA,Para,Sinks, PRE, nThresh,method)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

DATA_BU = DATA;

for i1 = 1:length(Para)
    
    if isstruct(DATA_BU(1).(Para{i1}))
        entries = length(Sinks);
    else
        entries = 1;
    end
    
    for i2 = 1: length(DATA_BU)

        for i3=1:entries
            try
            dummy = [DATA_BU.(Para{i1})];
            dummy = [dummy.(Sinks{i3})];
            dummy2 = sum(~isnan(dummy));
            dummy(:,dummy2 <= round(size(dummy,1)*nThresh)) = NaN;
            
            DIM = size(DATA_BU(1).(Para{i1}).(Sinks{i3}));
            catch
            dummy = [DATA_BU.(Para{i1})]; 
            dummy2 = sum(~isnan(dummy));
            dummy(:,dummy2 <= round(size(dummy,1)*nThresh)) = NaN;
            DIM = size(DATA_BU(1).(Para{i1}));            
            end
            
            if strcmp(method,'bin')
                
            M_all = []; STD_all = [];
                for i4 = 1:DIM(1) 
                    SA = [];
                    RAW = dummy(i4,:);
                    
                   for i5 = 1:length(DATA)
                       SA =vertcat(SA, RAW(1:DIM(2)));                      
                       RAW(1:DIM(2)) = [];                       
                   end
                   
                 M = nanmean(SA);
                 STD = nanstd(SA);

                    try
                     DATA(i2).(Para{i1}).(Sinks{i3})(i4,:)=(DATA_BU(i2).(Para{i1}).(Sinks{i3})(i4,:)-M)./STD;
                    catch
                     DATA(i2).(Para{i1})(i4,:)=(DATA_BU(i2).(Para{i1})(i4,:)-M)./STD;  
                    end
                end
                
                
            else
                M_all = []; STD_all = [];
                for i4 = 1:DIM(1)                
                 M = nanmean(dummy(i4,:));
                 STD = nanstd(dummy(i4,:));  

                 M_all = vertcat(M_all,M); STD_all = vertcat(STD_all,STD);             
                end
                
                try
                 DATA(i2).(Para{i1}).(Sinks{i3})=(DATA_BU(i2).(Para{i1}).(Sinks{i3})-M_all)./STD_all;
                catch
                 DATA(i2).(Para{i1})=(DATA_BU(i2).(Para{i1})-M_all)./STD_all;  
                end
            end           
        end 
    end  

end
Zscrd_DATA = DATA;

end