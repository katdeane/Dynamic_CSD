function layer_spikes = StimLayerCat(Stim_times, Layer, CurLay, AnNum, CortDepth)

%% Sweep Spike Channel sorting 
% author: Katrina Deane 8/30/12

% input:    the stimulus' spike time points arranged by layer buffered with nans
%           the current animal's layer definitions
%           the current layer being run
%           the current animal being run
%  ouput:   all spike time points concatonated to single vector for layer

% Used in code "dynamic_spike" to pre-sort spikes. 

% NOTE: using rounded values from SweepSpikeChan function

%%
if strcmp(CurLay,'All')
    layer_chan = 1:CortDepth;
else
    layer_chan = str2num(Layer.(CurLay){1,AnNum});
end

for iChan = 1:length(layer_chan)
    if iChan == 1
        layer_spikes = Stim_times(layer_chan(iChan),:);
    else
        layer_spikes = horzcat(layer_spikes,Stim_times(layer_chan(iChan),:));
    end
end
