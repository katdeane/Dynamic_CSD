function [latency, ampl, meanOfdata] = iGetPeakData_relres(data,range)

if nargin == 1
    
    [pks, locs, ~, p] = findpeaks((data*-1));
    
    [~, maxInd]=max(p); %which peak is most prominent
    
    ampl = pks(maxInd); %what's the amplitude of the peak there
    ampl = ampl*-1;
    latency = locs(maxInd); %what's the latency of that peak there
    
    meanOfdata = nanmean(abs(data));
    
    
else
    
    limDat = data(range(1):range(2))*-1;
    limDatplus = data(range(1):(range(2)+100));
    
    [pks, locs, ~, p] = findpeaks(limDat);
    
    [~, maxInd]=max(p);
    
    ampl = pks(maxInd);
    ampl = ampl*-1;
    latency = locs(maxInd)+range(1)-1;
    
    meanOfdata = nanmean(abs(limDatplus));
    
    
end