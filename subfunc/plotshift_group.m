function plotshift_group(Data,names,SINKs,ST_Shift,home,CurAn)
% Called by GroupAnalysis_fnc_temporal.m to generate figures. Figures saved
% in generated group folder in figs folder

for i2 = 1:length(fields(ST_Shift))
    FNames=fields(ST_Shift);
    
    for i3 = 1: length(fields(Data))
        if i3 == 1
            X = [Data.(names{i3})];
            X = {X.Condition};
        elseif length([Data.(names{i3})]) > length(X)
            X = [Data.(names{i3})];
            X = {X.Condition};
        end
    end
    
    for i3 = 1: length(ST_Shift)
        if isstruct(ST_Shift(i3).(FNames{i2})) && i3 ==1
            Container = struct;
            Entries = length(SINKs);
        elseif i3 ==1
            Container = struct;
            Entries = 1;
        end
        
        for i4 = 1:Entries
            if Entries == 1
                Container.Mean(i3) = nanmean(ST_Shift(i3).(FNames{i2}));
                Container.RMS_Mean(i3) = nanmean(sqrt(ST_Shift(i3).(FNames{i2}).^2));
                Container.RelShift(i3) = nansum(ST_Shift(i3).(FNames{i2}))./nansum(abs(ST_Shift(i3).(FNames{i2})));
            else
                Container.Mean.(SINKs{i4})(i3) = nanmean(ST_Shift(i3).(FNames{i2}).(SINKs{i4}));
                Container.RMS_Mean.(SINKs{i4})(i3) = nanmean(sqrt(ST_Shift(i3).(FNames{i2}).(SINKs{i4}).^2));
                Container.RelShift.(SINKs{i4})(i3) = nansum(ST_Shift(i3).(FNames{i2}).(SINKs{i4}))./nansum(abs(ST_Shift(i3).(FNames{i2}).(SINKs{i4})));
            end
        end
    end
    
    h=figure('units','normalized','outerposition',[0 0 1 1], 'Name',['rel. Sink Tuning shifts towards granular layer ', FNames{i2}]);
    for i3 = 1:Entries
        subplot(Entries,1,i3)
        if Entries == 1
            bar(Container.RMS_Mean,0.5)
            hold on
            bar(Container.Mean,0.25)
            bar(Container.RelShift,0.25)
            title([FNames{i2}])
            xticks (1:1:length(X)); xticklabels(X); xtickangle(-45); set(gca, 'Fontsize',8);
            legend('Abs','mean','RelShift')
        else
            bar(Container.RMS_Mean.(SINKs{i3}),0.5)
            hold on
            bar(Container.Mean.(SINKs{i3}),0.25)
            bar(Container.RelShift.(SINKs{i3}),0.25)
            title([FNames{i2} ' ' SINKs{i3}])
            xticks (1:1:length(X)); xticklabels(X); xtickangle(-45); set(gca, 'Fontsize',8);
            legend('Abs','mean','RelShift')
            
        end
    end
    cd (home), cd figs, cd (['Group_' CurAn(1:end-2)]);
    savefig(h,[CurAn(1:end-4), ' rel Sink Tuning shifts towards granular layer ', FNames{i2}]); close all
end