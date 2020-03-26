function [DUR,RMS,SINGLE_RMS,SINT,PAMP,SINGLE_SinkPeak,PLAT,SINGLE_PeakLat,INT] = ...
    sink_dura_single(Layer,AvgCSD,SingleTrialCSD,BL, Baseline)
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

onset_lev = BL+65; %NOTE: Awake animals have a presink and then an early sink which starts a bit later than 50ms
std_detect = 1.1;
std_lev = 1.5;
pre_std_lev = .8;

Order = {'VbE','IVE','VIaE','VIbE','VaE','I_IIE','InfE',...
    'VbL','IVL','VIaL','VIbL','VaL','I_IIL','InfL','all_chan'};
PAMP = struct;
PLAT = struct;
DUR = struct;
RMS = struct;
SINT = struct;
INT = struct;
SINGLE_SinkPeak = struct;
SINGLE_PeakLat = struct;
SINGLE_RMS = struct;

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
        rawCSD_single = (nanmean(SingleTrialCSD{i1}(Chan,:,:))) *-1;
        
        %zero all source info and shape the data for sink detection
        holdAvgCSDzero = AvgCSD{i1}(Chan,:);
        zerosource = find(holdAvgCSDzero(:,:)>=0);
        holdAvgCSDzero(zerosource) = 0; %#ok<*FNDSB> %equates all positive values (denoting sources) to zero
        zeroCSD_layer = (nanmean(holdAvgCSDzero))*-1; %flips negative values to positive 
        
        %nan all source for INT calculation
        holdAvgCSDnan = AvgCSD{i1}(Chan,:);
        holdAvgCSDnan(zerosource) = NaN;
        nanCSD = (nanmean(holdAvgCSDnan))*-1;
        
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
                %if the first point less than the onset level and if there's a peak following it, calculate rms
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
        
        %% Avg Data Structures
        
        %to fill the DUR with the correct sink time point (or NaN's if no time points)
        if isempty(Sink_time)
            Sink_time = [NaN NaN];
        end
        
        %to determine peak amplitude and latency
        if isnan(Sink_time(1))
            peakamp = NaN; peaklat = NaN; sinkrms = NaN; sinkint = NaN;
        else
            sinkint = nanmean(nanCSD(:,Sink_time(1):Sink_time(2)));
            sinkrms = rms(nanCSD(:,Sink_time(1):Sink_time(2)));
            peakamp = nanmax(nanCSD(:,Sink_time(1):Sink_time(2)));
            peaklat = (find(nanCSD(:,Sink_time(1):Sink_time(2)) == nanmax(rawCSD(:,Sink_time(1):Sink_time(2)))))+Sink_time(1);
        end
        
        PAMP(i1).(Order{i2}) = peakamp; % peak amplitude
        PLAT(i1).(Order{i2}) = round(peaklat); % peak latency
        DUR(i1).(Order{i2}) = round(Sink_time); % sink duration
        RMS(i1).(Order{i2})= sinkrms; % sink rms
        SINT(i1).(Order{i2})= sinkint; % sink integral (mean)
        
        state = round((Sink_time(1):Sink_time(end))-BL);
        if isnan(state)
            state = 1:2;
        end
        
        RMSINT = sqrt((rawCSD*-1).^2)';
        
        INT(i1).(Order{i2})= trapz(RMSINT(state+BL)); % integral of sink curve to get V*s/mm?
        
        %% Single Trial Data Structures
        
        if isnan(sinkrms)
            SINGLE_SinkPeak(i1).(Order{i2})= nan(50,1);
            SINGLE_PeakLat(i1).(Order{i2})= nan(50,1);
            SINGLE_RMS(i1).(Order{i2})= nan(50,1);
            continue % skip and go to next loop
        end
         
        for i4 = 1:size(rawCSD_single,3)
            curRun = rawCSD_single(:,Sink_time(1):Sink_time(2),i4);

            SINGLE_RMS(i1).(Order{i2})(i4) = rms(curRun);
            SINGLE_SinkPeak(i1).(Order{i2})(i4)= max(curRun);
            SINGLE_PeakLat(i1).(Order{i2})(i4)= find(curRun == max(curRun),1) + 200;
        end
        
        % if only 49 trils (case with some chronic animals)
        if length(SINGLE_SinkPeak(i1).(Order{i2})) == 49 
            SINGLE_SinkPeak(i1).(Order{i2})(50)= NaN;
            SINGLE_PeakLat(i1).(Order{i2})(50)= NaN;
            SINGLE_RMS(i1).(Order{i2})(50)= NaN;
        end
        % if more than 50 trials (%case with other chronic animals)
        if length(SINGLE_SinkPeak(i1).(Order{i2})) > 50 
            SINGLE_SinkPeak(i1).(Order{i2}) = SINGLE_SinkPeak(i1).(Order{i2})(1:50);
            SINGLE_PeakLat(i1).(Order{i2}) = SINGLE_PeakLat(i1).(Order{i2})(1:50);
            SINGLE_RMS(i1).(Order{i2}) = SINGLE_RMS(i1).(Order{i2})(1:50);
        end
    end
end