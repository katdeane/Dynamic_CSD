function [SinkRelRes] = Sink_RelRes(AvgLayerRelRes,DUR,L)
%UNTITLED Summary of this function goes here
%   fuck bla bla shit...
keyboard
SINKS = fieldnames(L);
Container = DUR; % simply overright existing stuff
for i1 = 1:length(AvgLayerRelRes)
    
    dummy_RelRes = AvgLayerRelRes{i1};
    
    for i2 = 1:length(SINKS)
        try
        Sink_window = dummy_RelRes([L.(SINKS{i2})],[DUR(i1).(SINKS{i2})(1):DUR(i1).(SINKS{i2})(2)]);
        keyboard
%         [~, ~, RMS_AvgRecCSD,~, ~,~, ~] =...
%                 ExtractCSDBasedPar_01(AvgLayerRelRes,BL,...
%                 tone, 1000, 1, 0, threshold_std, threshold_dur,...
%                 Latency_HighCutoff, Latency_LowCutoff,state);
        catch
        end
    end
end



x=1;
SinkRelRes = x;

end

