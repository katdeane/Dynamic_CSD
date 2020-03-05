function [DUR,RMS,PAMP,PLAT,INT] = sink_dura_group(Layer,AvgCSD,BL,tone,SWEEP,ChanOrder, Baseline)
%This function produces the duration of sinks within pre-specified
%layers.

dbstop if error

%Strictness of detection can be adjusted by changing 'std_lev'. The early input
%onset is currently cut to 300 sampling points, which can be adjusted with
%'onset_lev'. Offset of the early sink in each layer specifies the detection point from which
%the onset of the late sink can start. In layers 4, 5b, and 6, strictness
%of detection of pre-sink is set dynamically within the code from the thresold
%level pre_lev which is the standard deviation away from the mean and the
%point at which we take the ratio of fwhm (the pre-sink will only
%be detected if there is already a normal sink in the area and if so, it
%will combine the sink and pre-sink into one signal box).

onset_lev = BL+65;
std_detect = 1.1;
std_lev = 1.5;
pre_std_lev = .8;

Order = {'VbE','IVE','VIaE','VIbE','VaE','I_IIE','VbL','IVL','VIaL','VIbL','VaL','I_IIL','all_chan'};
threshold_std = 2;
threshold_dur = 0.005;
Latency_HighCutoff = 0.2;
Latency_LowCutoff = 0.015;

for i1= 1:length(AvgCSD) %length of stimuli
    
    %create std and mean line from previously generated mat or from current
    %CSDs' BL horizontally stacking
    if ~isempty(Baseline)
        std_BL = Baseline.threshstd;
        mean_BL = Baseline.threshmean;
    else
        global_BL = cell2mat(cellfun(@(x) x(:,1:BL),AvgCSD, 'UniformOutput',0));
        std_BL = nanstd(global_BL(:)); %the standard deviation of the summed baseline
        mean_BL = nanmean(global_BL(:)); %the mean of the summed baseline
    end
    
    for i2 = 1:length(Order) %Order defined just above
        
        %to take correct channels for layers pre-set in group .m's
        if strcmp(Order{i2}, 'all_chan')
            Chan = (1:size(AvgCSD{i1},1)); % all layers/channels
        else
            Chan = Layer.(Order{i2});
        end        
        
        %current stimulus avgcsd only in the current layer
        rawCSD = (nanmean(AvgCSD{i1}(Chan,:)))*-1; %to calculate latency and amplitude from unaltered signal
        
        %zero all source info and shape the data for sink detection
        holdAvgCSD = AvgCSD{i1}(Chan,:);
        zerosource =find(holdAvgCSD(:,:)>=0);
        holdAvgCSD(zerosource)= 0; %equates all positive values (denoting sources) to zero
        zeroCSD_layer = (nanmean(holdAvgCSD))*-1; %flips negative values to positive 
        
        g = gausswin(10); %generates a gausswin distribution of 15 points from 1
        g = g/sum(g); %turns it into a percentage distribution
        AvgCSD_layer = conv(zeroCSD_layer, g, 'same'); %normalizes CSD to distribution
        AvgCSD_layer(:,end) = mean_BL - std_BL; %so that it definitely has an end
        
        %define thresholds for signal detection
        T = ones(1,length(AvgCSD_layer));
        thresh_mean = (mean_BL + (std_BL*std_detect))*T; 
        thresh_std = (mean_BL+(std_BL*std_lev))*T; 
        thresh_fwhm = mean_BL+(std_BL*pre_std_lev); 
        
        %% Early Input %%
        if i2 <= (length(Order)/2)

            AvgCSD_layer(:,1:BL) = mean_BL - std_BL; %so that the first point is counted as an actual intercept
                
            %find intercept points
            P = InterX([1:length(AvgCSD_layer);thresh_mean],[1:length(AvgCSD_layer);AvgCSD_layer]);
            P = P(1,:);

            rmslist = NaN(1,length(P)-1);
            for i3 = 1:length(P)-1
                %if the first point is 300 or less and if there's a peak following it, calculate rms
                if P(i3) <= (onset_lev) && nanmax(AvgCSD_layer(:,P(i3):P(i3+1))) > mean(thresh_std)
                    rmslist(i3) = rms(AvgCSD_layer(:,P(i3):P(i3+1))); 
                end
            end
            Sink_time=[P(find(rmslist == max(rmslist))) P(find(rmslist == max(rmslist))+1)];
            
            %% Pre-sink %%
            if i2 <= 4 && sum(~isnan(rmslist)) ~= 0 %if it's layer 4, 5b, or 6 early and there is a sink
                sinkmax = nanmax(AvgCSD_layer(:,Sink_time(1):Sink_time(2)));
                ratio = round(sinkmax/thresh_fwhm);
                fwhm = (((sinkmax-thresh_fwhm)/ratio)+thresh_fwhm); 
                fwhm = fwhm*T;
                
                pre_rmslist = NaN(1,length(P)-1); %new list to compare only the smaller early sink
                for i3 = 1:length(P)-1
                    %if onset is before 250 and there is a the peak needed
                    if P(i3) <= (250) && nanmax(AvgCSD_layer(:,P(i3):P(i3+1))) > mean(fwhm)
                        pre_rmslist(i3) = rms(AvgCSD_layer(:,P(i3):P(i3+1)));
                    end
                end
                
                if sum(~isnan(pre_rmslist)) ~= 0 %did we find a sink here?
                    pre_rms = find(~isnan(pre_rmslist)); 
                    %take the onset of the first sink below 250 and the offset of the biggest sink below 300
                    Sink_time=[P(pre_rms(1)) P(find(rmslist == max(rmslist))+1)]; 
                end
            end
            
        %% Late Input 
        else
            if ~isnan(DUR(i1).(Order{i2-(floor(length(Order)/2))})(2)) %set starting point as end of early sink if exists or based on onset level chosen
                early = DUR(i1).(Order{i2-(floor(length(Order)/2))})(2);
            else 
                early = onset_lev;
            end
             
            AvgCSD_layer(:,1:early) = mean_BL - std_BL;
            
            P = InterX([1:length(AvgCSD_layer);thresh_mean],[1:length(AvgCSD_layer);AvgCSD_layer]);
            P =P(1,:);
            
            rmslist = NaN(1,length(P)-1);
            for i3 = 1:length(P)-1
                if nanmax(AvgCSD_layer(:,P(i3):P(i3+1))) > mean(thresh_std)
                    rmslist(i3) = rms(AvgCSD_layer(:,P(i3):P(i3+1)));
                end
            end
            Sink_time=[P(find(rmslist == max(rmslist))) P(find(rmslist == max(rmslist))+1)];
            
        end %end of sink detection loop
        
       %% Fill data into structures
       
        %to fill the DUR with the correct sink time point (or NaN's if no time points)
        if isempty(Sink_time)
            Sink_time = [NaN NaN];
        end
        
        %to determine peak amplitude and latency 
        if isnan(Sink_time(1))
            peakamp = NaN; peaklat = NaN; 
        else
            peakamp = nanmax(rawCSD(:,Sink_time(1):Sink_time(2)));
            peaklat = (find(rawCSD(:,Sink_time(1):Sink_time(2)) == nanmax(rawCSD(:,Sink_time(1):Sink_time(2)))))+Sink_time(1);
        end
        
        PAMP(i1).(Order{i2}) = peakamp;
        PLAT(i1).(Order{i2}) = round(peaklat);
        DUR(i1).(Order{i2}) = round(Sink_time);
        
        state = round((Sink_time(1):Sink_time(end))-BL);
        if isnan(state)
            state = 1:2;
        end
        
%         LayerCSD = cellfun(@(x)mean(x(Chan,:)),AvgCSD, 'UniformOutput', false);
        RMSINT = sqrt((rawCSD*-1).^2)';
        try % all of above to generate state
            [~, ~, ~, RMS_AvgRecCSD,~, ~,~, ~] =...
                ExtractCSDBasedPar(1, zeroCSD_layer,BL,...
                1000, 1, 0, threshold_std, threshold_dur,... 
                Latency_HighCutoff, Latency_LowCutoff,state);
        catch
            fprintf('You are getting full NaNs! /n')
            RMS_AvgRecCSD = NaN;
        end
        
        if length(state) == 2
            RMS_AvgRecCSD = NaN;
        end
        
        RMS(i1).(Order{i2})= RMS_AvgRecCSD;
        INT(i1).(Order{i2})= trapz(RMSINT(state+BL));
    end
end