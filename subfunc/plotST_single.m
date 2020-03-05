function plotST_single(Data,Data_STBased,names,SINKs,para,zscore,nThresh,bin,mirror,home,CurAn,BF_Pos)
% Called by GroupAnalysis_fnc_temporal.m to generate figures. Figures saved
% in generated group folder in figs folder
SA =[];
    for i2 = 1:length(para)
        for i3 =1: length(names)
            dummy2 = [];
            if isstruct(Data(1).(names{1}).(para{i2})) % Check whether input is sink based
                Entries = length(SINKs);
            else
                Entries = 1;
            end
            for i4 = 1:Entries
                dummy =[];
                for i5 =1:length(Data_STBased)
                    if ~isempty(Data(i5).(names{i3}))
                        if Entries == 1
                            dummy = vertcat(dummy, Data_STBased(i5).(para{i2})(i3,:));
                        else
                            dummy = vertcat(dummy, Data_STBased(i5).(para{i2}).(SINKs{i4})(i3,:));
                        end
                    end
                end
                dummy2(:,:,i4) = dummy;
            end
            SA.(names{i3}).(para{i2}).measurements =dummy2;
        end
    end
    
    for i2 = 1:length(names)
        for i3 = 1:length(para)
            h=figure('units','normalized','outerposition',[0 0 1 1],'Name',[CurAn(1:end-4) '_' names{i2}...
                '_ST based tempAnal_'  para{i3} '_Thresh_' num2str(nThresh) '_zscored_' num2str(zscore) '_bin_' num2str(bin) '_mirror_' num2str(mirror)]);
            if isstruct(Data(1).(names{1}).(para{i3})) % Check whether input is sink based
                Entries = length(SINKs);
            else
                Entries = 1;
            end
            if zscore == 0
                if strcmp(Data(3).(names{1}).Condition ,'Pre_3')
                    setPre = 3;
                    Pre = SA.(names{i2}).(para{i3}).measurements(1:setPre,:,:);
                    PreBL = nanmean(Pre,1);
                else
                    setPre = 1;
                    Pre = SA.(names{i2}).(para{i3}).measurements(1,:,:);
                    PreBL = Pre;
                end
                Pre = Pre./PreBL;
            else
                if strcmp(Data(3).(names{1}).Condition ,'Pre_3')
                    setPre = 3;
                else
                    setPre = 1;
                end
                Pre = ones(size(SA.(names{i2}).(para{i3}).measurements(1,:,:))); PreBL = Pre;
                
            end
            PreSTD = 0.25;
            
            dummyBF = SA.(names{i2}).(para{i3}).measurements(:,BF_Pos,:)./PreBL(:,BF_Pos,:);
            dummyLNearBF = SA.(names{i2}).(para{i3}).measurements(:,BF_Pos-1,:)./PreBL(:,BF_Pos-1,:);
            dummyLNonBF = SA.(names{i2}).(para{i3}).measurements(:,BF_Pos-2,:)./PreBL(:,BF_Pos-2,:);
            dummyLoffBF = SA.(names{i2}).(para{i3}).measurements(:,BF_Pos-3,:)./PreBL(:,BF_Pos-3,:);
            dummyHNearBF = SA.(names{i2}).(para{i3}).measurements(:,BF_Pos+1,:)./PreBL(:,BF_Pos+1,:);
            dummyHNonBF = SA.(names{i2}).(para{i3}).measurements(:,BF_Pos+2,:)./PreBL(:,BF_Pos+2,:);
            dummyHoffBF = SA.(names{i2}).(para{i3}).measurements(:,BF_Pos+3,:)./PreBL(:,BF_Pos+3,:);
            
            for i4 =1:size (SA.(names{i2}).(para{i3}).measurements,3)
                X = [Data.(names{i2})];
                X = {X.Condition};
                
                subplot(size (SA.(names{i2}).(para{i3}).measurements,3),7,1+((i4-1)*7))
                plot(dummyLoffBF(:,:,i4))
                
                xticks (1:1:length(X)); xticklabels(X); xtickangle(-45); set(gca, 'Fontsize',8);
                hold on
                try
                    if Entries == 1
                        STD2 = 0.25; % STD2 = 2*nanstd(Pre(:,BF_Pos-3,i4));
                        M = nanmean(Pre(:,BF_Pos-3,i4));
                        plot([1 length(dummyBF)],[M+STD2 M+STD2],'LineStyle',':','Color','black')
                    else
                        plot([1 length(dummyBF)],[1+PreSTD 1+PreSTD],'LineStyle',':','Color','black')
                    end
                catch
                end
                if zscore == 0; ylim([0.5 2]); end
                title(['LowoffBF ' SINKs{i4}])
                
                subplot(size (SA.(names{i2}).(para{i3}).measurements,3),7,2+((i4-1)*7))
                plot(dummyLNonBF(:,:,i4))
                xticks (1:1:length(X)); xticklabels(X); xtickangle(-45); set(gca, 'Fontsize',8);
                hold on
                try
                    if Entries == 1
                        STD2 = 0.25; % STD2 = 2*nanstd(Pre(:,BF_Pos-3,i4));
                        M = nanmean(Pre(:,BF_Pos-2,i4));
                        plot([1 length(dummyBF)],[M+STD2 M+STD2],'LineStyle',':','Color','black')
                    else
                        plot([1 length(dummyBF)],[1+PreSTD 1+PreSTD],'LineStyle',':','Color','black')
                    end
                catch
                end
                if zscore == 0; ylim([0.5 2]); end
                title(['LownonBF ' SINKs{i4}])
                
                subplot(size (SA.(names{i2}).(para{i3}).measurements,3),7,3+((i4-1)*7))
                plot(dummyLNearBF(:,:,i4))
                xticks (1:1:length(X)); xticklabels(X); xtickangle(-45); set(gca, 'Fontsize',8);
                hold on
                try
                    if Entries == 1
                        STD2 = 0.25; % STD2 = 2*nanstd(Pre(:,BF_Pos-3,i4));
                        M = nanmean(Pre(:,BF_Pos-1,i4));
                        plot([1 length(dummyBF)],[M+STD2 M+STD2],'LineStyle',':','Color','black')
                    else
                        plot([1 length(dummyBF)],[1+PreSTD 1+PreSTD],'LineStyle',':','Color','black')
                    end
                catch
                end
                if zscore == 0; ylim([0.5 2]); end
                title(['LownearBF ' SINKs{i4}])
                
                subplot(size (SA.(names{i2}).(para{i3}).measurements,3),7,4+((i4-1)*7))
                plot(dummyBF(:,:,i4))
                xticks (1:1:length(X)); xticklabels(X); xtickangle(-45); set(gca, 'Fontsize',8);
                hold on
                try if Entries == 1
                        STD2 = 0.25; % STD2 = 2*nanstd(Pre(:,BF_Pos-3,i4));
                        M = nanmean(Pre(:,BF_Pos,i4));
                        plot([1 length(dummyBF)],[M+STD2 M+STD2],'LineStyle',':','Color','black')
                    else
                        plot([1 length(dummyBF)],[1+PreSTD 1+PreSTD],'LineStyle',':','Color','black')
                    end
                catch
                end
                if zscore == 0; ylim([0.5 2]); end
                title(['BF ' SINKs{i4}])
                
                subplot(size (SA.(names{i2}).(para{i3}).measurements,3),7,5+((i4-1)*7))
                plot(dummyHNearBF(:,:,i4))
                xticks (1:1:length(X)); xticklabels(X); xtickangle(-45); set(gca, 'Fontsize',8);
                hold on
                try if Entries == 1
                        STD2 = 0.25; % STD2 = 2*nanstd(Pre(:,BF_Pos-3,i4));
                        M = nanmean(Pre(:,BF_Pos+1,i4));
                        plot([1 length(dummyBF)],[M+STD2 M+STD2],'LineStyle',':','Color','black')
                    else
                        plot([1 length(dummyBF)],[1+PreSTD 1+PreSTD],'LineStyle',':','Color','black')
                    end
                catch
                end
                if zscore == 0; ylim([0.5 2]); end
                title(['HighnearBF ' SINKs{i4}])
                
                subplot(size (SA.(names{i2}).(para{i3}).measurements,3),7,6+((i4-1)*7))
                plot(dummyHNonBF(:,:,i4))
                xticks (1:1:length(X)); xticklabels(X); xtickangle(-45); set(gca, 'Fontsize',8);
                hold on
                try if Entries == 1
                        STD2 = 0.25; % STD2 = 2*nanstd(Pre(:,BF_Pos-3,i4));
                        M = nanmean(Pre(:,BF_Pos+2,i4));
                        plot([1 length(dummyBF)],[M+STD2 M+STD2],'LineStyle',':','Color','black')
                    else
                        plot([1 length(dummyBF)],[1+PreSTD 1+PreSTD],'LineStyle',':','Color','black')
                    end
                catch
                end
                if zscore == 0; ylim([0.5 2]); end
                title(['HighnonBF ' SINKs{i4}])
                
                subplot(size (SA.(names{i2}).(para{i3}).measurements,3),7,7+((i4-1)*7))
                plot(dummyHoffBF(:,:,i4))
                xticks (1:1:length(X)); xticklabels(X); xtickangle(-45); set(gca, 'Fontsize',8);
                hold on
                try if Entries == 1
                        STD2 = 0.25; % STD2 = 2*nanstd(Pre(:,BF_Pos-3,i4));
                        M = nanmean(Pre(:,BF_Pos+3,i4));
                        plot([1 length(dummyBF)],[M+STD2 M+STD2],'LineStyle',':','Color','black')
                    else
                        plot([1 length(dummyBF)],[1+PreSTD 1+PreSTD],'LineStyle',':','Color','black')
                    end
                catch
                end
                if zscore == 0; ylim([0.5 2]); end
                title(['HighoffBF ' SINKs{i4}])
                
            end
            cd (home), cd figs, cd (['Group_' CurAn(1:end-2)]);
            savefig(h,[CurAn(1:end-4) '_' names{i2} '_ST_tempAnal_' para{i3}...
                '_Thresh_' num2str(nThresh) '_zscored_' num2str(zscore) '_bin_' num2str(bin) '_mirror_' num2str(mirror) ...
                '.fig']); close all
        end
    end