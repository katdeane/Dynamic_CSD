%% DOCUMENT YOUR CODE PLS THANKS 
%%
Session = 11;
animal = 1;
Sink = 'IVE';


IDs = fieldnames(Data);

Frqz = Data(Session).(IDs{animal}).Frqz;
AVG_Peak_BF = Data(Session).(IDs{animal}).GS_BF;

PB_Pos = find(Frqz ==AVG_Peak_BF);

AVG_Single = [Data(Session).(IDs{animal}).SingleSinkPeakAmp.(Sink)];
B = reshape(AVG_Single,[50,length([Data(Session).(IDs{animal}).SinkPeakAmp.(Sink)])]);

B = nanmean(B);


Lippi_Single = [Data(Session).(IDs{animal}).SingleSinkPeakAmp_AVG.(Sink)];
C = reshape(Lippi_Single,[50,length([Data(Session).(IDs{animal}).SinkPeakAmp.(Sink)])]);

C = nanmean(C);

figure,
plot(1:length([Data(Session).(IDs{animal}).SinkPeakAmp.(Sink)]),...
    [Data(Session).(IDs{animal}).SinkPeakAmp.(Sink)],'LineWidth',2)
hold on
plot(1:length([Data(Session).(IDs{animal}).SinkPeakAmp.(Sink)]),...
    B,'LineWidth',2)
plot(1:length([Data(Session).(IDs{animal}).SinkPeakAmp.(Sink)]),...
    C,'LineWidth',2)
plot(PB_Pos,nanmax([Data(Session).(IDs{animal}).SinkPeakAmp.(Sink)]),'r*')


xticks(1:1:length([Data(Session).(IDs{animal}).SinkPeakAmp.(Sink)]))
xticklabels(strsplit(num2str(Data(Session).(IDs{animal}).Frqz)))

title([IDs{animal} ' ' Data(Session).(IDs{animal}).Condition ' ' Sink])

legend("evoked","induced","Lippilyse")