function [latency, ampl, meanOfdata] = iGetPeakData_avrec(data,range)

if sum(data) == 0
    latency = NaN;
    ampl = NaN;
    meanOfdata = NaN;
    
else
    if nargin == 1
        
        [pks, locs, ~, p] = findpeaks(data);
        
        [str, maxInd]=max(p); %which peak is most prominent 
        
        if str < 0.00008 %arbitrary threshold for me to try it out
            ampl = NaN;
            latency = NaN;
            meanOfdata = NaN;
        else
            ampl = pks(maxInd);
            latency = locs(maxInd);
            meanOfdata = nanmean(abs(data));
            
        end
        
        
                
    else
        
        limDat = data(range(1):range(2));
        limDatplus = data(range(1):(range(2)+100));
        
        [pks, locs, ~, p] = findpeaks(limDat);
        
        [str, maxInd]=max(p);
        
        if str < 0.0005
            ampl = NaN;
            latency = NaN;
            meanOfdata = NaN;
        else
            ampl = pks(maxInd);
            latency = locs(maxInd)+range(1)-1;
            meanOfdata = nanmean(abs(limDatplus));
        end
        
        
        
    end
end