function plotGS_group(Data,Data_GSBased,names,SINKs,para,zscore,nThresh,bin,mirror,home,CurAn)
% Called by GroupAnalysis_fnc_temporal.m to generate figures. Figures saved
% in generated group folder in figs folder

for i3 = 1:length(para)
    WORKgs.(para{i3})=[];
    if ~isstruct(Data_GSBased(1).(para{i2}))
        WORKgs.(para{i3})=[];
    else
        for i4 =1:length(SINKs)
            WORKgs.(para{i3}).(SINKs{i4})=[];
        end
    end
end

clear  PreDummy
for i2 = 1:length(para)
    
    if isstruct(Data(1).(names{1}).(para{i2})) % Check wether input is sink based
        Entries = length(SINKs);
    else
        Entries = 1;
    end
    
    for i3 = 1:Entries
        % calculate Pre means for normalization
        clear PreDummy
        for i4 = 1:setPre
            if Entries == 1
                PreDummy(:,:,i4) = Data_GSBased(i4).(para{i2});
            else
                PreDummy(:,:,i4) = Data_GSBased(i4).(para{i2}).(SINKs{i3});
            end
        end
        for i4 = 1:size(PreDummy,1)
            X = PreDummy(i4,:,:);
            Y = reshape(X,size(X,1)*size(X,2),size(X,3));
            Y = Y';
            Y = nanmean(Y);
            PreDummy2(i4,:) =Y;
        end
        PreDummy = PreDummy2; clear X PreDummy2 dummyBF dummyLNearBF  dummyLNonBF dummyHNearBF dummyHNonBF dummyLoffBF dummyHoffBF
        if zscore == 1
            PreDummy = ones(size(PreDummy));
        end
        
        if Entries == 1
            X =Data_GSBased;
            PreDummylength = length(X);
        else
            X =[Data_GSBased.(para{i2})];
            PreDummylength = length(X);
        end
        
        %Calculate normalized temp curves
        DummyMean =[]; DummySEM =[];
        for i4 = 1: PreDummylength
            if Entries == 1
                Dummy = X(i4).(para{i2});
                DummyMean = vertcat(DummyMean, nanmean(Dummy./PreDummy));
                DummySEM = vertcat(DummySEM, nanstd(Dummy./PreDummy)./sqrt(sum(~isnan(Dummy))));
            else
                Dummy = X(i4).(SINKs{i3});
                try
                    DummyMean = vertcat(DummyMean, nanmean(Dummy./PreDummy));
                    DummySEM = vertcat(DummySEM, nanstd(Dummy./PreDummy)./sqrt(sum(~isnan(Dummy))));
                catch
                    %                     DummyMean = vertcat(DummyMean, nanmean(Dummy./PreDummy));
                    %                     DummySEM = vertcat(DummySEM, nanstd(Dummy./PreDummy)/sqrt(size(Dummy,1)));
                end
            end
        end
        if Entries == 1
            WORKgs.(para{i2}).Mean=DummyMean;
            WORKgs.(para{i2}).SEM=DummySEM;
            WORKgs.(para{i2}).posPre2STD =1+2*nanstd(DummyMean(1:setPre,:));
            WORKgs.(para{i2}).negPre2STD =1-2*nanstd(DummyMean(1:setPre,:));
        else
            WORKgs.(para{i2}).(SINKs{i3}).Mean=DummyMean;
            WORKgs.(para{i2}).(SINKs{i3}).SEM=DummySEM;
            WORKgs.(para{i2}).(SINKs{i3}).posPre2STD =1+2*nanstd(DummyMean(1:setPre,:));
            WORKgs.(para{i2}).(SINKs{i3}).negPre2STD =1-2*nanstd(DummyMean(1:setPre,:));
        end
    end
end

for i3 = 1:length(para)
    
    h=figure('units','normalized','outerposition',[0 0 1 1],'Name',[ 'Group_analysis '...
        CurAn(1:end-4) '_GS _tempAnal_' para{i3} '_Thresh_' num2str(nThresh)...
        '_zscored_' num2str(zscore) '_bin_' num2str(bin) '_mirror_' num2str(mirror) '.fig']);
    
    for i4 = 1:length(names)
        X=[Data.(names{i4})];
        X={X.Condition};
        if i4 ==1
            XLab = X;
        elseif length(X) > length(XLab)
            XLab = X;
        end
    end
    
    if isstruct(Data(1).(names{1}).(para{i3})) % Check wether input is sink based
        Entries = length(SINKs);
    else
        Entries = 1;
        clear dummyBF dummyLNearBF dummyLNonBF dummyHNearBF dummyHNonBF dummyLoffBF dummyHoffBF
    end
    
    for i4 = 1:Entries
        try
            if Entries == 1
                try
                    dummyBF(:,1) = WORKgs.(para{i3}).(SINKs{i4}).Mean(:,BF_Pos,:);dummyBF(:,2) = WORKgs.(para{i3}).(SINKs{i4}).SEM(:,BF_Pos,:);
                    dummyLNearBF(:,1) =  WORKgs.(para{i3}).(SINKs{i4}).Mean(:,BF_Pos-1,:);dummyLNearBF(:,2) = WORKgs.(para{i3}).(SINKs{i4}).SEM(:,BF_Pos-1,:);
                    dummyLNonBF(:,1) =  WORKgs.(para{i3}).(SINKs{i4}).Mean(:,BF_Pos-2,:);dummyLNonBF(:,2) = WORKgs.(para{i3}).(SINKs{i4}).SEM(:,BF_Pos-2,:);
                    dummyLoffBF(:,1) =  WORKgs.(para{i3}).(SINKs{i4}).Mean(:,BF_Pos-3,:);dummyLoffBF(:,2) = WORKgs.(para{i3}).(SINKs{i4}).SEM(:,BF_Pos-3,:);
                    dummyHNearBF(:,1) =  WORKgs.(para{i3}).(SINKs{i4}).Mean(:,BF_Pos+1,:);dummyHNearBF(:,2) = WORKgs.(para{i3}).(SINKs{i4}).SEM(:,BF_Pos+1,:);
                    dummyHNonBF(:,1) =  WORKgs.(para{i3}).(SINKs{i4}).Mean(:,BF_Pos+2,:);dummyHNonBF(:,2) = WORKgs.(para{i3}).(SINKs{i4}).SEM(:,BF_Pos+2,:);
                    dummyHoffBF(:,1) =  WORKgs.(para{i3}).(SINKs{i4}).Mean(:,BF_Pos+3,:);dummyHoffBF(:,2) = WORKgs.(para{i3}).(SINKs{i4}).SEM(:,BF_Pos+3,:);
                    Y =WORKgs.(para{i3}).(SINKs{i4}).posPre2STD; Y2=WORKgs.(para{i3}).(SINKs{i4}).negPre2STD;
                catch
                    dummyBF(:,1) = WORKgs.(para{i3}).Mean(:,BF_Pos,:);dummyBF(:,2) = WORKgs.(para{i3}).SEM(:,BF_Pos,:);
                    dummyLNearBF(:,1) =  WORKgs.(para{i3}).Mean(:,BF_Pos-1,:);dummyLNearBF(:,2) = WORKgs.(para{i3}).SEM(:,BF_Pos-1,:);
                    dummyLNonBF(:,1) =  WORKgs.(para{i3}).Mean(:,BF_Pos-2,:);dummyLNonBF(:,2) = WORKgs.(para{i3}).SEM(:,BF_Pos-2,:);
                    dummyLoffBF(:,1) =  WORKgs.(para{i3}).Mean(:,BF_Pos-3,:);dummyLoffBF(:,2) = WORKgs.(para{i3}).SEM(:,BF_Pos-3,:);
                    dummyHNearBF(:,1) =  WORKgs.(para{i3}).Mean(:,BF_Pos+1,:);dummyHNearBF(:,2) = WORKgs.(para{i3}).SEM(:,BF_Pos+1,:);
                    dummyHNonBF(:,1) =  WORKgs.(para{i3}).Mean(:,BF_Pos+2,:);dummyHNonBF(:,2) = WORKgs.(para{i3}).SEM(:,BF_Pos+2,:);
                    dummyHoffBF(:,1) =  WORKgs.(para{i3}).Mean(:,BF_Pos+3,:);dummyHoffBF(:,2) = WORKgs.(para{i3}).SEM(:,BF_Pos+3,:);
                    Y =WORKgs.(para{i3}).posPre2STD; Y2=WORKgs.(para{i3}).negPre2STD;
                end
                
                if Entries == 1
                    subplot(Entries,7,1+((i4-1)*7))
                else
                    subplot(Entries,7,1+((i4-1)*7))
                end
                shadedErrorBar([],dummyLoffBF(:,1),dummyLoffBF(:,2),'lineprops',{'r','markerfacecolor','r'})
                hold on
                try
                    if zscore == 0
                        if Entries == 1
                            M = nanmean(dummyLoffBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyLoffBF(1:3,1));
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyLoffBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                        else
                            M = nanmean(dummyLoffBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyLoffBF(1:3,1));
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyLoffBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                            %                 plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[1+PreSTD 1+PreSTD],'k.:');
                            %                 plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[1-PreSTD 1-PreSTD],'k.:');
                            %                 plot(1:length(dummyLoffBF(:,1)) ,(dummyLoffBF(:,1) >= (1+ PreSTD))*1.8 ,'r*');
                        end
                    else
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos-3)+Y2(BF_Pos-3) Y(BF_Pos-3)+Y2(BF_Pos-3)],'k.:');
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos-3)-Y2(BF_Pos-3) Y(BF_Pos-3)-Y2(BF_Pos-3)],'k.:');
                    end
                catch
                    if zscore == 0
                        if Entries == 1
                            M = nanmean(dummyLNonBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyLNonBF(1:3,1));
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyLNonBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                        else
                            M = nanmean(dummyLNonBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyLNonBF(1:3,1));
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyLNonBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                            %                 plot([1 size(WORKgs.(para{i3}).Mean,1)],[1+PreSTD 1+PreSTD],'k.:');
                            %                 plot([1 size(WORKgs.(para{i3}).Mean,1)],[1-PreSTD 1-PreSTD],'k.:');
                            %                  plot(1:length(dummyLoffBF(:,1)) ,(dummyLNonBF(:,1) >= (1+ PreSTD))*1.8 ,'r*');
                        end
                    else
                        plot([1 size(WORKgs.(para{i3}).Mean,1)],[Y(BF_Pos-3)+Y2(BF_Pos-3) Y(BF_Pos-3)+Y2(BF_Pos-3)],'k.:');
                        plot([1 size(WORKgs.(para{i3}).Mean,1)],[Y(BF_Pos-3)-Y2(BF_Pos-3) Y(BF_Pos-3)-Y2(BF_Pos-3)],'k.:');
                    end
                end
                if zscore == 0; ylim([0.5 2]); end
                hold off
                title(['LowOffBF ' SINKs{i4}])
                xticks (1:1:length(XLab)); xticklabels(XLab); xtickangle(-45); set(gca, 'Fontsize',8);
                
                if Entries == 1
                    subplot(Entries,7,2+((i4-1)*7))
                else
                    subplot(Entries,7,2+((i4-1)*7))
                end
                shadedErrorBar([],dummyLNonBF(:,1),dummyLNonBF(:,2),'lineprops',{'r','markerfacecolor','r'})
                hold on
                try
                    if zscore == 0
                        if Entries == 1
                            M = nanmean(dummyLNonBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyLNonBF(1:3,1));
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyLNonBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                        else
                            M = nanmean(dummyLNonBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyLNonBF(1:3,1));
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyLNonBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                            %                 plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[1+PreSTD 1+PreSTD],'k.:');
                            %                 plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[1-PreSTD 1-PreSTD],'k.:');
                            %                  plot(1:length(dummyLoffBF(:,1)) ,(dummyLNonBF(:,1) >= (1+ PreSTD))*1.8 ,'r*');
                        end
                    else
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos-2)+Y2(BF_Pos-2) Y(BF_Pos-2)+Y2(BF_Pos-2)],'k.:');
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos-2)-Y2(BF_Pos-2) Y(BF_Pos-2)-Y2(BF_Pos-2)],'k.:');
                        plot(1:length(dummyLoffBF(:,1)) ,(dummyLNonBF(:,1) >= (1+ PreSTD))*1.8 ,'r*');
                    end
                catch
                    if zscore == 0
                        if Entries == 1
                            M = nanmean(dummyLNonBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyLNonBF(1:3,1));
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyLNonBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                        else
                            M = nanmean(dummyLNonBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyLNonBF(1:3,1));
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyLNonBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                            %                 plot([1 size(WORKgs.(para{i3}).Mean,1)],[1+PreSTD 1+PreSTD],'k.:');
                            %                 plot([1 size(WORKgs.(para{i3}).Mean,1)],[1-PreSTD 1-PreSTD],'k.:');
                            %                  plot(1:length(dummyLoffBF(:,1)) ,(dummyLNonBF(:,1) >= (1+ PreSTD))*1.8 ,'r*');
                        end
                    else
                        plot([1 size(WORKgs.(para{i3}).Mean,1)],[Y(BF_Pos-2)+Y2(BF_Pos-2) Y(BF_Pos-2)+Y2(BF_Pos-2)],'k.:');
                        plot([1 size(WORKgs.(para{i3}).Mean,1)],[Y(BF_Pos-2)-Y2(BF_Pos-2) Y(BF_Pos-2)-Y2(BF_Pos-2)],'k.:');
                    end
                end
                
                if zscore == 0; ylim([0.5 2]); end
                hold off
                title(['LowNonBF ' SINKs{i4}])
                xticks (1:1:length(XLab)); xticklabels(XLab); xtickangle(-45); set(gca, 'Fontsize',8);
                if Entries == 1
                    subplot(Entries,7,3+((i4-1)*7))
                else
                    subplot(Entries,7,3+((i4-1)*7))
                end
                shadedErrorBar([],dummyLNearBF(:,1),dummyLNearBF(:,2),'lineprops',{'r','markerfacecolor','r'})
                hold on
                try
                    if zscore == 0
                        if Entries == 1
                            M = nanmean(dummyLNearBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyLNearBF(1:3,1));
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyLNearBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                        else
                            M = nanmean(dummyLNearBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyLNearBF(1:3,1));
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyLNearBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                            %                 plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[1+PreSTD 1+PreSTD],'k.:');
                            %                 plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[1-PreSTD 1-PreSTD],'k.:');
                            %                  plot(1:length(dummyLoffBF(:,1)) ,(dummyLNearBF(:,1) >= (1+ PreSTD))*1.8 ,'r*');
                        end
                    else
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos-1)+Y2(BF_Pos-1) Y(BF_Pos-1)+Y2(BF_Pos-1)],'k.:');
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos-1)-Y2(BF_Pos-1) Y(BF_Pos-1)-Y2(BF_Pos-1)],'k.:');
                    end
                catch
                    if zscore == 0
                        if Entries == 1
                            M = nanmean(dummyLNearBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyLNearBF(1:3,1));
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyLNearBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                        else
                            M = nanmean(dummyLNearBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyLNearBF(1:3,1));
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyLNearBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                            %                 plot([1 size(WORKgs.(para{i3}).Mean,1)],[1+PreSTD 1+PreSTD],'k.:');
                            %                 plot([1 size(WORKgs.(para{i3}).Mean,1)],[1-PreSTD 1-PreSTD],'k.:');
                            %                 plot(1:length(dummyLoffBF(:,1)) ,(dummyLNearBF(:,1) >= (1+ PreSTD))*1.8 ,'r*');
                        end
                    else
                        plot([1 size(WORKgs.(para{i3}).Mean,1)],[Y(BF_Pos-1)+Y2(BF_Pos-1) Y(BF_Pos-1)+Y2(BF_Pos-1)],'k.:');
                        plot([1 size(WORKgs.(para{i3}).Mean,1)],[Y(BF_Pos-1)-Y2(BF_Pos-1) Y(BF_Pos-1)-Y2(BF_Pos-1)],'k.:');
                    end
                end
                
                if zscore == 0; ylim([0.5 2]); end
                hold off
                title(['LowNearBF ' SINKs{i4}])
                xticks (1:1:length(XLab)); xticklabels(XLab); xtickangle(-45); set(gca, 'Fontsize',8);
                if Entries == 1
                    subplot(Entries,7,4+((i4-1)*7))
                else
                    subplot(Entries,7,4+((i4-1)*7))
                end
                shadedErrorBar([],dummyBF(:,1),dummyBF(:,2),'lineprops',{'r','markerfacecolor','r'})
                hold on
                try
                    if zscore == 0
                        if Entries == 1
                            M = nanmean(dummyBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyBF(1:3,1));
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                        else
                            M = nanmean(dummyBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyBF(1:3,1));
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                            %                 plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[1+PreSTD 1+PreSTD],'k.:');
                            %                 plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[1-PreSTD 1-PreSTD],'k.:');
                            %                 plot(1:length(dummyLoffBF(:,1)) ,(dummyBF(:,1) >= (1+ PreSTD))*1.8 ,'r*');
                        end
                    else
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos)+Y2(BF_Pos) Y(BF_Pos)+Y2(BF_Pos)],'k.:');
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos)-Y2(BF_Pos) Y(BF_Pos)-Y2(BF_Pos)],'k.:');
                    end
                catch
                    if zscore == 0
                        if Entries == 1
                            M = nanmean(dummyBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyBF(1:3,1));
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                        else
                            M = nanmean(dummyBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyBF(1:3,1));
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                            %                 plot([1 size(WORKgs.(para{i3}).Mean,1)],[1+PreSTD 1+PreSTD],'k.:');
                            %                 plot([1 size(WORKgs.(para{i3}).Mean,1)],[1-PreSTD 1-PreSTD],'k.:');
                            %                 plot(1:length(dummyLoffBF(:,1)) ,(dummyBF(:,1) >= (1+ PreSTD))*1.8 ,'r*');
                        end
                    else
                        plot([1 size(WORKgs.(para{i3}).Mean,1)],[Y(BF_Pos)+Y2(BF_Pos) Y(BF_Pos)+Y2(BF_Pos)],'k.:');
                        plot([1 size(WORKgs.(para{i3}).Mean,1)],[Y(BF_Pos)-Y2(BF_Pos) Y(BF_Pos)-Y2(BF_Pos)],'k.:');
                    end
                end
                
                if zscore == 0; ylim([0.5 2]); end
                hold off
                title(['BF ' SINKs{i4}])
                xticks (1:1:length(XLab)); xticklabels(XLab); xtickangle(-45); set(gca, 'Fontsize',8);
                if Entries == 1
                    subplot(Entries,7,5+((i4-1)*7))
                else
                    subplot(Entries,7,5+((i4-1)*7))
                end
                shadedErrorBar([],dummyHNearBF(:,1),dummyHNearBF(:,2),'lineprops',{'r','markerfacecolor','r'})
                hold on
                try
                    if zscore == 0
                        if Entries == 1
                            M = nanmean(dummyHNearBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyHNearBF(1:3,1));
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyHNearBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                        else
                            M = nanmean(dummyHNearBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyHNearBF(1:3,1));
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyHNearBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                            %                 plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[1+PreSTD 1+PreSTD],'k.:');
                            %                 plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[1-PreSTD 1-PreSTD],'k.:');
                            %                 plot(1:length(dummyLoffBF(:,1)) ,(dummyHNearBF(:,1) >= (1+ PreSTD))*1.8 ,'r*');
                        end
                    else
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos+1)+Y2(BF_Pos+1) Y(BF_Pos+1)+Y2(BF_Pos+1)],'k.:');
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos+1)-Y2(BF_Pos+1) Y(BF_Pos+1)-Y2(BF_Pos+1)],'k.:');
                    end
                catch
                    if zscore == 0
                        if Entries == 1
                            M = nanmean(dummyHNearBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyHNearBF(1:3,1));
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyHNearBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                        else
                            M = nanmean(dummyHNearBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyHNearBF(1:3,1));
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyHNearBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                            %                 plot([1 size(WORKgs.(para{i3}).Mean,1)],[1+PreSTD 1+PreSTD],'k.:');
                            %                 plot([1 size(WORKgs.(para{i3}).Mean,1)],[1-PreSTD 1-PreSTD],'k.:');
                            %                 plot(1:length(dummyLoffBF(:,1)) ,(dummyHNearBF(:,1) >= (1+ PreSTD))*1.8 ,'r*');
                        end
                    else
                        plot([1 size(WORKgs.(para{i3}).Mean,1)],[Y(BF_Pos+1)+Y2(BF_Pos+1) Y(BF_Pos+1)+Y2(BF_Pos+1)],'k.:');
                        plot([1 size(WORKgs.(para{i3}).Mean,1)],[Y(BF_Pos+1)-Y2(BF_Pos+1) Y(BF_Pos+1)-Y2(BF_Pos+1)],'k.:');
                    end
                end
                
                if zscore == 0; ylim([0.5 2]); end
                hold off
                title(['HighNearBF ' SINKs{i4}])
                xticks (1:1:length(XLab)); xticklabels(XLab); xtickangle(-45); set(gca, 'Fontsize',8);
                if Entries == 1
                    subplot(Entries,7,6+((i4-1)*7))
                else
                    subplot(Entries,7,6+((i4-1)*7))
                end
                shadedErrorBar([],dummyHNonBF(:,1),dummyHNonBF(:,2),'lineprops',{'r','markerfacecolor','r'})
                hold on
                try
                    if zscore == 0
                        if Entries == 1
                            M = nanmean(dummyHNonBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyHNonBF(1:3,1));
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyHNonBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                        else
                            M = nanmean(dummyHNonBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyHNonBF(1:3,1));
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyHNonBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                            %                 plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[1+PreSTD 1+PreSTD],'k.:');
                            %                 plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[1-PreSTD 1-PreSTD],'k.:');
                            %                 plot(1:length(dummyLoffBF(:,1)) ,(dummyHNonBF(:,1) >= (1+ PreSTD))*1.8 ,'r*');
                        end
                    else
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos+2)+Y2(BF_Pos+2) Y(BF_Pos+2)+Y2(BF_Pos+2)],'k.:');
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos+2)-Y2(BF_Pos+2) Y(BF_Pos+2)-Y2(BF_Pos+2)],'k.:');
                    end
                catch
                    if zscore == 0
                        if Entries == 1
                            M = nanmean(dummyHNonBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyHNonBF(1:3,1));
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyHNonBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                        else
                            M = nanmean(dummyHNonBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyHNonBF(1:3,1));
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyHNonBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                            %                 plot([1 size(WORKgs.(para{i3}).Mean,1)],[1+PreSTD 1+PreSTD],'k.:');
                            %                 plot([1 size(WORKgs.(para{i3}).Mean,1)],[1-PreSTD 1-PreSTD],'k.:');
                            %                 plot(1:length(dummyLoffBF(:,1)) ,(dummyHNonBF(:,1) >= (1+ PreSTD))*1.8 ,'r*');
                        end
                    else
                        plot([1 size(WORKgs.(para{i3}).Mean,1)],[Y(BF_Pos+2)+Y2(BF_Pos+2) Y(BF_Pos+2)+Y2(BF_Pos+2)],'k.:');
                        plot([1 size(WORKgs.(para{i3}).Mean,1)],[Y(BF_Pos+2)-Y2(BF_Pos+2) Y(BF_Pos+2)-Y2(BF_Pos+2)],'k.:');
                    end
                end
                
                if zscore == 0; ylim([0.5 2]); end
                hold off
                title(['HighNonBF ' SINKs{i4}])
                xticks (1:1:length(XLab)); xticklabels(XLab); xtickangle(-45); set(gca, 'Fontsize',8);
                if Entries == 1
                    subplot(Entries,7,7+((i4-1)*7))
                else
                    subplot(Entries,7,7+((i4-1)*7))
                end
                shadedErrorBar([],dummyHoffBF(:,1),dummyHoffBF(:,2),'lineprops',{'r','markerfacecolor','r'})
                hold on
                try
                    if zscore == 0
                        if Entries == 1
                            M = nanmean(dummyHoffBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyHoffBF(1:3,1));
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyHoffBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                        else
                            M = nanmean(dummyHoffBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyHoffBF(1:3,1));
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyHoffBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                            %                 plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[1+PreSTD 1+PreSTD],'k.:');
                            %                 plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[1-PreSTD 1-PreSTD],'k.:');
                            %                 plot(1:length(dummyLoffBF(:,1)) ,(dummyHoffBF(:,1) >= (1+ PreSTD))*1.8 ,'r*');
                        end
                    else
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos+3)+Y2(BF_Pos+3) Y(BF_Pos+3)+Y2(BF_Pos+3)],'k.:');
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos+3)-Y2(BF_Pos+3) Y(BF_Pos+3)-Y2(BF_Pos+3)],'k.:');
                    end
                catch
                    if zscore == 0
                        if Entries == 1
                            M = nanmean(dummyHoffBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyHoffBF(1:3,1));
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyHoffBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                        else
                            M = nanmean(dummyHoffBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyHoffBF(1:3,1));
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKgs.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyHoffBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                            %                 plot([1 size(WORKgs.(para{i3}).Mean,1)],[1+PreSTD 1+PreSTD],'k.:');
                            %                 plot([1 size(WORKgs.(para{i3}).Mean,1)],[1-PreSTD 1-PreSTD],'k.:');
                            %                 plot(1:length(dummyLoffBF(:,1)) ,(dummyHoffBF(:,1) >= (1+ PreSTD))*1.8 ,'r*');
                        end
                    else
                        plot([1 size(WORKgs.(para{i3}).Mean,1)],[Y(BF_Pos+3)+Y2(BF_Pos+3) Y(BF_Pos+3)+Y2(BF_Pos+3)],'k.:');
                        plot([1 size(WORKgs.(para{i3}).Mean,1)],[Y(BF_Pos+3)-Y2(BF_Pos+3) Y(BF_Pos+3)-Y2(BF_Pos+3)],'k.:');
                    end
                end
                
                if zscore == 0; ylim([0.5 2]); end
                hold off
                title(['HighOffBF ' SINKs{i4}])
                xticks (1:1:length(XLab)); xticklabels(XLab); xtickangle(-45); set(gca, 'Fontsize',8)
                
            else
                clear dummyBF dummyLNearBF dummyLNonBF dummyHNearBF dummyHNonBF dummyLoffBF dummyHoffBF
                dummyBF(:,1) = WORKgs.(para{i3}).(SINKs{i4}).Mean(:,BF_Pos,:);dummyBF(:,2) = WORKgs.(para{i3}).(SINKs{i4}).SEM(:,BF_Pos,:);
                dummyLNearBF(:,1) =  WORKgs.(para{i3}).(SINKs{i4}).Mean(:,BF_Pos-1,:);dummyLNearBF(:,2) = WORKgs.(para{i3}).(SINKs{i4}).SEM(:,BF_Pos-1,:);
                dummyLNonBF(:,1) =  WORKgs.(para{i3}).(SINKs{i4}).Mean(:,BF_Pos-2,:);dummyLNonBF(:,2) = WORKgs.(para{i3}).(SINKs{i4}).SEM(:,BF_Pos-2,:);
                dummyLoffBF(:,1) =  WORKgs.(para{i3}).(SINKs{i4}).Mean(:,BF_Pos-3,:);dummyLoffBF(:,2) = WORKgs.(para{i3}).(SINKs{i4}).SEM(:,BF_Pos-3,:);
                dummyHNearBF(:,1) =  WORKgs.(para{i3}).(SINKs{i4}).Mean(:,BF_Pos+1,:);dummyHNearBF(:,2) = WORKgs.(para{i3}).(SINKs{i4}).SEM(:,BF_Pos+1,:);
                dummyHNonBF(:,1) =  WORKgs.(para{i3}).(SINKs{i4}).Mean(:,BF_Pos+2,:);dummyHNonBF(:,2) = WORKgs.(para{i3}).(SINKs{i4}).SEM(:,BF_Pos+2,:);
                dummyHoffBF(:,1) =  WORKgs.(para{i3}).(SINKs{i4}).Mean(:,BF_Pos+3,:);dummyHoffBF(:,2) = WORKgs.(para{i3}).(SINKs{i4}).SEM(:,BF_Pos+3,:);
                
                subplot(Entries,7,1+((i4-1)*7))
                shadedErrorBar([],dummyLoffBF(:,1),dummyLoffBF(:,2),'lineprops',{'r','markerfacecolor','r'})
                hold on
                try
                    if zscore == 0
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[1+PreSTD 1+PreSTD],'k.:');
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[1-PreSTD 1-PreSTD],'k.:');
                        plot(1:length(dummyLoffBF(:,1)) ,(dummyLoffBF(:,1) >= (1+ PreSTD))*1.8 ,'r*');
                    else
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos-3)+Y2(BF_Pos-3) Y(BF_Pos-3)+Y2(BF_Pos-3)],'k.:');
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos-3)-Y2(BF_Pos-3) Y(BF_Pos-3)-Y2(BF_Pos-3)],'k.:');
                    end
                catch
                end
                if zscore == 0; ylim([0.5 2]); end
                hold off
                title(['LowOffBF ' SINKs{i4}])
                xticks (1:1:length(XLab)); xticklabels(XLab); xtickangle(-45); set(gca, 'Fontsize',8);
                
                subplot(Entries,7,2+((i4-1)*7))
                shadedErrorBar([],dummyLNonBF(:,1),dummyLNonBF(:,2),'lineprops',{'r','markerfacecolor','r'})
                hold on
                try
                    if zscore == 0
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[1+PreSTD 1+PreSTD],'k.:');
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[1-PreSTD 1-PreSTD],'k.:');
                        plot(1:length(dummyLoffBF(:,1)) ,(dummyLNonBF(:,1) >= (1+ PreSTD))*1.8 ,'r*');
                    else
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos-2)+Y2(BF_Pos-2) Y(BF_Pos-2)+Y2(BF_Pos-2)],'k.:');
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos-2)-Y2(BF_Pos-2) Y(BF_Pos-2)-Y2(BF_Pos-2)],'k.:');
                    end
                catch
                end
                if zscore == 0; ylim([0.5 2]); end
                hold off
                title(['LowNonBF ' SINKs{i4}])
                xticks (1:1:length(XLab)); xticklabels(XLab); xtickangle(-45); set(gca, 'Fontsize',8);
                
                subplot(Entries,7,3+((i4-1)*7))
                shadedErrorBar([],dummyLNearBF(:,1),dummyLNearBF(:,2),'lineprops',{'r','markerfacecolor','r'})
                hold on
                try
                    if zscore == 0
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[1+PreSTD 1+PreSTD],'k.:');
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[1-PreSTD 1-PreSTD],'k.:');
                        plot(1:length(dummyLoffBF(:,1)) ,(dummyLNearBF(:,1) >= (1+ PreSTD))*1.8 ,'r*');
                    else
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos-1)+Y2(BF_Pos-1) Y(BF_Pos-1)+Y2(BF_Pos-1)],'k.:');
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos-1)-Y2(BF_Pos-1) Y(BF_Pos-1)-Y2(BF_Pos-1)],'k.:');
                    end
                catch
                end
                if zscore == 0; ylim([0.5 2]); end
                hold off
                title(['LowNearBF ' SINKs{i4}])
                xticks (1:1:length(XLab)); xticklabels(XLab); xtickangle(-45); set(gca, 'Fontsize',8);
                
                subplot(Entries,7,4+((i4-1)*7))
                shadedErrorBar([],dummyBF(:,1),dummyBF(:,2),'lineprops',{'r','markerfacecolor','r'})
                hold on
                try
                    if zscore == 0
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[1+PreSTD 1+PreSTD],'k.:');
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[1-PreSTD 1-PreSTD],'k.:');
                        plot(1:length(dummyLoffBF(:,1)) ,(dummyBF(:,1) >= (1+ PreSTD))*1.8 ,'r*');
                    else
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos)+Y2(BF_Pos) Y(BF_Pos)+Y2(BF_Pos)],'k.:');
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos)-Y2(BF_Pos) Y(BF_Pos)-Y2(BF_Pos)],'k.:');
                    end
                catch
                end
                if zscore == 0; ylim([0.5 2]); end
                hold off
                title(['BF ' SINKs{i4}])
                xticks (1:1:length(XLab)); xticklabels(XLab); xtickangle(-45); set(gca, 'Fontsize',8);
                
                subplot(Entries,7,5+((i4-1)*7))
                shadedErrorBar([],dummyHNearBF(:,1),dummyHNearBF(:,2),'lineprops',{'r','markerfacecolor','r'})
                hold on
                try
                    if zscore == 0
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[1+PreSTD 1+PreSTD],'k.:');
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[1-PreSTD 1-PreSTD],'k.:');
                        plot(1:length(dummyLoffBF(:,1)) ,(dummyHNearBF(:,1) >= (1+ PreSTD))*1.8 ,'r*');
                    else
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos+1)+Y2(BF_Pos+1) Y(BF_Pos+1)+Y2(BF_Pos+1)],'k.:');
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos+1)-Y2(BF_Pos+1) Y(BF_Pos+1)-Y2(BF_Pos+1)],'k.:');
                    end
                catch
                end
                if zscore == 0; ylim([0.5 2]); end
                hold off
                title(['HighNearBF ' SINKs{i4}])
                xticks (1:1:length(XLab)); xticklabels(XLab); xtickangle(-45); set(gca, 'Fontsize',8);
                
                subplot(Entries,7,6+((i4-1)*7))
                shadedErrorBar([],dummyHNonBF(:,1),dummyHNonBF(:,2),'lineprops',{'r','markerfacecolor','r'})
                hold on
                try
                    if zscore == 0
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[1+PreSTD 1+PreSTD],'k.:');
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[1-PreSTD 1-PreSTD],'k.:');
                        plot(1:length(dummyLoffBF(:,1)) ,(dummyHNonBF(:,1) >= (1+ PreSTD))*1.8 ,'r*');
                    else
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos+2)+Y2(BF_Pos+2) Y(BF_Pos+2)+Y2(BF_Pos+2)],'k.:');
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos+2)-Y2(BF_Pos+2) Y(BF_Pos+2)-Y2(BF_Pos+2)],'k.:');
                    end
                catch
                end
                if zscore == 0; ylim([0.5 2]); end
                hold off
                title(['HighNonBF ' SINKs{i4}])
                xticks (1:1:length(XLab)); xticklabels(XLab); xtickangle(-45); set(gca, 'Fontsize',8);
                
                subplot(Entries,7,7+((i4-1)*7))
                shadedErrorBar([],dummyHoffBF(:,1),dummyHoffBF(:,2),'lineprops',{'r','markerfacecolor','r'})
                hold on
                try
                    if zscore == 0
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[1+PreSTD 1+PreSTD],'k.:');
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[1-PreSTD 1-PreSTD],'k.:');
                        plot(1:length(dummyLoffBF(:,1)) ,(dummyHoffBF(:,1) >= (1+ PreSTD))*1.8 ,'r*');
                    else
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos+3)+Y2(BF_Pos+3) Y(BF_Pos+3)+Y2(BF_Pos+3)],'k.:');
                        plot([1 size(WORKgs.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos+3)-Y2(BF_Pos+3) Y(BF_Pos+3)-Y2(BF_Pos+3)],'k.:');
                    end
                catch
                end
                if zscore == 0; ylim([0.5 2]); end
                hold off
                title(['HighOffBF ' SINKs{i4}])
                xticks (1:1:length(XLab)); xticklabels(XLab); xtickangle(-45); set(gca, 'Fontsize',8);
            end
        end
    end
    cd (home), cd figs, cd (['Group_' CurAn(1:end-2)]);
    savefig(h,['Group_analysis ' CurAn(1:end-4)...
        '_GS based tempAnal_' para{i3} '_Thresh_' num2str(nThresh)...
        '_zscored_' num2str(zscore) '_bin_' num2str(bin) '_mirror_' num2str(mirror)...
        '.fig']); close all
end