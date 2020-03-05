function plotshift_single(Data,names,SINKs,ST_Shift, NumAnimals,home,CurAn,SParam)
% Called by GroupAnalysis_fnc_temporal.m to generate figures. Figures saved
% in generated group folder in figs folder
for i2 = 1:length(SParam)
  h=figure('units','normalized','outerposition',[0 0 1 1], 'Name',[CurAn(1:end-4), ' Single animals GS shift', SParam{i2}]);      
        for i3 = 1: length(fields(Data))
            if i3 == 1
            X = [Data.(names{i3})];
            X = {X.Condition};
            elseif length([Data.(names{i3})]) > length(X)
            X = [Data.(names{i3})];
            X = {X.Condition};
            end
        end
        Labels = X;
        X = [ST_Shift.(SParam{i2})];

        for i3 = 1:length(SINKs)
             Y = [X.(SINKs{i3})];
             Y = reshape(Y,NumAnimals,[]);
             subplot(length(SINKs),1,i3)
             if size(Y,1) == length(names)
             plot(Y')
             else
             plot(Y)
             end
             ylabel(SINKs{i3})
             xticks (1:1:length(Labels)); xticklabels(Labels); xtickangle(-45); set(gca, 'Fontsize',8);
             if i3 ==1
                 legend(names)
             end
        end
        cd (home), cd figs, cd (['Group_' CurAn(1:end-2)]);
        savefig(h,[CurAn(1:end-4), ' Single animals GS shift ', SParam{i2}]); close all
end    