function structout = consec_sinks(structin, sinkonset_struct, num_sinks, dur_stim, start_time)

if ~exist('start_time','var')
    start_time = 1; % 3 ms to avoid that the sink starts directly at 0 ms
end
if ~exist('num_sinks','var')
    num_sinks  = 2;
end
if ~exist('dur_stim','var')
    dur_stim   = 1000;
end

% structin   = Data(imeas).(para{ipar})(istim).(layer{ilay});
% sinkonset_struct = Data(imeas).Sinkonset(istim).(layer{ilay});

%preallocation of onset detection window containers
structout  = nan(1,num_sinks);
det_on     = nan(1,num_sinks);
det_off    = nan(1,num_sinks);
det_jump   = dur_stim/num_sinks;

if det_jump > 90
    det_dur = 90; % so that onset detection window is 3:65 ms
else
    det_dur = det_jump - 1; % 1 ms space between detection windows 
end

% fill detection window containers
for idet = 1:num_sinks
    if idet == 1
        det_on(idet) = start_time;
    else
        det_on(idet) = det_on(idet-1) + det_jump;
    end
    det_off(idet) = det_on(idet) + det_dur;
end


%% Take features

% this runs through all sinks for each detection window opportunity
for itake = 1:num_sinks
    for isink = 1:length(structin)
        
        if isnan(structin(isink))
            continue
        end
        % if sink onset is within detection window
        if sinkonset_struct(isink) > det_on(itake) && sinkonset_struct(isink) < det_off(itake)
            
            if isnan(structout(itake))
                % pull out sink feature and place it into appropriate bin
                structout(itake) =  structin(isink);
            else
                % if there is already a sink detected here, take the bigger
                % one
                if structin(isink) > structout(itake)
                    structout(itake) =  structin(isink);
                end
            end
             
        end
    end
end