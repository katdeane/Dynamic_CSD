function plotST_group(Data,Data_STBased,names,SINKs,para,zscore,nThresh,bin,mirror,home,CurAn)
% Called by GroupAnalysis_fnc_temporal.m to generate figures. Figures saved
% in generated group folder in figs folder

% generate empty structure skeleton to fill with data at each sink per parameter
for i3 = 1:length(para)
    WORKst.(para{i3})=[];
    if ~isstruct(Data_STBased(1).(para{i3}))
        WORKst.(para{i3})=[];
    else
        for i4 =1:length(SINKs)
            WORKst.(para{i3}).(SINKs{i4})=[];
        end
    end
end

for ip = 1:length(para)
    if isstruct(Data(1).(names{1}).(para{ip})) % Check whether input is sink based
        Entries = length(SINKs);
    else
        Entries = 1;
    end
    for i3 = 1:Entries
        % calculate Pre means for normalization
        clear PreDummy
        for i4 = 1:setPre
            if Entries == 1
                PreDummy(:,:,i4) = Data_STBased(i4).(para{ip});
            else
                PreDummy(:,:,i4) = Data_STBased(i4).(para{ip}).(SINKs{i3});
            end
        end
        for i4 = 1:size(PreDummy,1)
            X = PreDummy(i4,:,:);
            Y = reshape(X,size(X,1)*size(X,2),size(X,3));
            Y = Y';
            Y2 =nanstd(Y);
            Y = nanmean(Y);
            PreDummy2(i4,:) =Y;
        end
        
        PreDummy = PreDummy2; clear X PreDummy2 dummyBF dummyLNearBF  dummyLNonBF dummyHNearBF dummyHNonBF dummyLoffBF dummyHoffBF
        if zscore == 1
            PreDummy = ones(size(PreDummy));
        end
        
        if Entries == 1
            X =Data_STBased;
            PreDummylength = length(X);
        else
            X =[Data_STBased.(para{ip})];
            PreDummylength = length(X);
        end
        
        %Calculate normalized temp curves
        DummyMean =[]; DummySEM =[];
        for i4 = 1: PreDummylength
            if Entries == 1
                Dummy = X(i4).(para{ip});
                DummyMean = vertcat(DummyMean, nanmean(Dummy./PreDummy));
                DummySEM = vertcat(DummySEM, nanstd(Dummy./PreDummy)./sqrt(sum(~isnan(Dummy))));           %  sqrt(size(Dummy,1)));
            else
                Dummy = X(i4).(SINKs{i3});
                DummyMean = vertcat(DummyMean, nanmean(Dummy./PreDummy));
                DummySEM = vertcat(DummySEM, nanstd(Dummy./PreDummy)./sqrt(sum(~isnan(Dummy))));%/sqrt(size(Dummy,1)));
            end
        end
        
        if Entries == 1
            WORKst.(para{ip}).Mean=DummyMean;
            WORKst.(para{ip}).SEM=DummySEM;
            WORKst.(para{ip}).posPre2STD =1+2*nanstd(DummyMean(1:setPre,:));
            WORKst.(para{ip}).negPre2STD =1-2*nanstd(DummyMean(1:setPre,:));
        else
            WORKst.(para{ip}).(SINKs{i3}).Mean=DummyMean;
            WORKst.(para{ip}).(SINKs{i3}).SEM=DummySEM;
            WORKst.(para{ip}).(SINKs{i3}).posPre2STD =1+2*nanstd(DummyMean(1:setPre,:));
            WORKst.(para{ip}).(SINKs{i3}).negPre2STD =1-2*nanstd(DummyMean(1:setPre,:));
        end
    end
end

for i3 = 1:length(para)
    
    h=figure('units','normalized','outerposition',[0 0 1 1],'Name',[ 'Group_analysis '...
        CurAn(1:end-4) '_ST_tempAnal_' para{i3} '_zscored_'...
        num2str(zscore) '_bin_' num2str(bin) '_mirror_' num2str(mirror)]);
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
        if Entries == 1
            try
                dummyBF(:,1) = WORKst.(para{i3}).Mean(:,BF_Pos,:);dummyBF(:,2) = WORKst.(para{i3}).SEM(:,BF_Pos,:);
                dummyLNearBF(:,1) =  WORKst.(para{i3}).Mean(:,BF_Pos-1,:);dummyLNearBF(:,2) = WORKst.(para{i3}).SEM(:,BF_Pos-1,:);
                dummyLNonBF(:,1) =  WORKst.(para{i3}).Mean(:,BF_Pos-2,:);dummyLNonBF(:,2) = WORKst.(para{i3}).SEM(:,BF_Pos-2,:);
                dummyLoffBF(:,1) =  WORKst.(para{i3}).Mean(:,BF_Pos-3,:);dummyLoffBF(:,2) = WORKst.(para{i3}).SEM(:,BF_Pos-3,:);
                dummyHNearBF(:,1) =  WORKst.(para{i3}).Mean(:,BF_Pos+1,:);dummyHNearBF(:,2) = WORKst.(para{i3}).SEM(:,BF_Pos+1,:);
                dummyHNonBF(:,1) =  WORKst.(para{i3}).Mean(:,BF_Pos+2,:);dummyHNonBF(:,2) = WORKst.(para{i3}).SEM(:,BF_Pos+2,:);
                dummyHoffBF(:,1) =  WORKst.(para{i3}).Mean(:,BF_Pos+2,:);dummyHoffBF(:,2) = WORKst.(para{i3}).SEM(:,BF_Pos+2,:);
                
                subplot(Entries,7,1+((i4-1)*7))
                shadedErrorBar([],dummyLoffBF(:,1),dummyLoffBF(:,2),'lineprops',{'r','markerfacecolor','r'})
                hold on
                %             keyboard
                try
                    if zscore == 0
                        if Entries == 1
                            M = nanmean(dummyLoffBF(1:3,1));
                            %                        STD2 = 3*nanmean(dummyLoffBF(1:3,1));
                            STD2 = 0.25; %STD2 = 3*nanmean(dummyLoffBF(1:3,2));
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyLoffBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                            
                        else
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[1+PreSTD 1+PreSTD],'k.:');
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[1-PreSTD 1-PreSTD],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyLoffBF(:,1) >= (1+ PreSTD))*1.8 ,'r*');
                        end
                    else
                        plot([1 size(WORKst.(para{i3}).Mean,1)],[Y(BF_Pos-3)+Y2(BF_Pos-3) Y(BF_Pos-3)+Y2(BF_Pos-3)],'k.:');
                        plot([1 size(WORKst.(para{i3}).Mean,1)],[Y(BF_Pos-3)-Y2(BF_Pos-3) Y(BF_Pos-3)-Y2(BF_Pos-3)],'k.:');
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
                        if Entries == 1
                            M = nanmean(dummyLNonBF(1:3,1));
                            STD2 = 0.25; %STD2 = 3*nanmean(dummyLNonBF(1:3,2));
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyLNonBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                        else
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[1+PreSTD 1+PreSTD],'k.:');
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[1-PreSTD 1-PreSTD],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyLNonBF(:,1) >= (1+ PreSTD))*1.8 ,'r*');
                        end
                    else
                        plot([1 size(WORKst.(para{i3}).Mean,1)],[Y(BF_Pos-2)+Y2(BF_Pos-2) Y(BF_Pos-2)+Y2(BF_Pos-2)],'k.:');
                        plot([1 size(WORKst.(para{i3}).Mean,1)],[Y(BF_Pos-2)-Y2(BF_Pos-2) Y(BF_Pos-2)-Y2(BF_Pos-2)],'k.:');
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
                        if Entries == 1
                            M = nanmean(dummyLNearBF(1:3,1));
                            STD2 = 0.25; %STD2 = 3*nanmean(dummyLNearBF(1:3,2));
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyLNearBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                        else
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[1+PreSTD 1+PreSTD],'k.:');
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[1-PreSTD 1-PreSTD],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyLNearBF(:,1) >= (1+ PreSTD))*1.8 ,'r*');
                        end
                    else
                        plot([1 size(WORKst.(para{i3}).Mean,1)],[Y(BF_Pos-1)+Y2(BF_Pos-1) Y(BF_Pos-1)+Y2(BF_Pos-1)],'k.:');
                        plot([1 size(WORKst.(para{i3}).Mean,1)],[Y(BF_Pos-1)-Y2(BF_Pos-1) Y(BF_Pos-1)-Y2(BF_Pos-1)],'k.:');
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
                        if Entries == 1
                            M = nanmean(dummyBF(1:3,1));
                            STD2 = 0.25; %STD2 = 3*nanmean(dummyBF(1:3,2));
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                        else
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[1+PreSTD 1+PreSTD],'k.:');
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[1-PreSTD 1-PreSTD],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyBF(:,1) >= (1+ PreSTD))*1.8 ,'r*');
                        end
                    else
                        plot([1 size(WORKst.(para{i3}).Mean,1)],[Y(BF_Pos)+Y2(BF_Pos) Y(BF_Pos)+Y2(BF_Pos)],'k.:');
                        plot([1 size(WORKst.(para{i3}).Mean,1)],[Y(BF_Pos)-Y2(BF_Pos) Y(BF_Pos)-Y2(BF_Pos)],'k.:');
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
                        if Entries == 1
                            M = nanmean(dummyHNearBF(1:3,1));
                            STD2 = 4*nanmean(dummyHNearBF(1:3,2));
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyHNearBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                        else
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[1+PreSTD 1+PreSTD],'k.:');
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[1-PreSTD 1-PreSTD],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyHNearBF(:,1) >= (1+ PreSTD))*1.8 ,'r*');
                        end
                    else
                        plot([1 size(WORKst.(para{i3}).Mean,1)],[Y(BF_Pos+1)+Y2(BF_Pos+1) Y(BF_Pos+1)+Y2(BF_Pos+1)],'k.:');
                        plot([1 size(WORKst.(para{i3}).Mean,1)],[Y(BF_Pos+1)-Y2(BF_Pos+1) Y(BF_Pos+1)-Y2(BF_Pos+1)],'k.:');
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
                        if Entries == 1
                            M = nanmean(dummyHNonBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyHNonBF(1:3,2));
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyHNonBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                        else
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[1+PreSTD 1+PreSTD],'k.:');
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[1-PreSTD 1-PreSTD],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyHNonBF(:,1) >= (1+ PreSTD))*1.8 ,'r*');
                        end
                    else
                        plot([1 size(WORKst.(para{i3}).Mean,1)],[Y(BF_Pos+2)+Y2(BF_Pos+2) Y(BF_Pos+2)+Y2(BF_Pos+2)],'k.:');
                        plot([1 size(WORKst.(para{i3}).Mean,1)],[Y(BF_Pos+2)-Y2(BF_Pos+2) Y(BF_Pos+2)-Y2(BF_Pos+2)],'k.:');
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
                        if Entries == 1
                            M = nanmean(dummyHoffBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyHoffBF(1:3,2));
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyHoffBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                        else
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[1+PreSTD 1+PreSTD],'k.:');
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[1-PreSTD 1-PreSTD],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyHoffBF(:,1) >= (1+ PreSTD))*1.8 ,'r*');
                        end
                    else
                        plot([1 size(WORKst.(para{i3}).Mean,1)],[Y(BF_Pos+3)+Y2(BF_Pos+3) Y(BF_Pos+3)+Y2(BF_Pos+3)],'k.:');
                        plot([1 size(WORKst.(para{i3}).Mean,1)],[Y(BF_Pos+3)-Y2(BF_Pos+3) Y(BF_Pos+3)-Y2(BF_Pos+3)],'k.:');
                    end
                catch
                end
                if zscore == 0; ylim([0.5 2]); end
                hold off
                title(['HighOffBF ' SINKs{i4}])
                xticks (1:1:length(XLab)); xticklabels(XLab); xtickangle(-45); set(gca, 'Fontsize',8);
            catch
            end
            
            %keyboard
        else
            clear dummyBF dummyLNearBF dummyLNonBF dummyHNearBF dummyHNonBF dummyLoffBF dummyHoffBF
            dummyBF(:,1) = WORKst.(para{i3}).(SINKs{i4}).Mean(:,BF_Pos,:);dummyBF(:,2) = WORKst.(para{i3}).(SINKs{i4}).SEM(:,BF_Pos,:);
            dummyLNearBF(:,1) =  WORKst.(para{i3}).(SINKs{i4}).Mean(:,BF_Pos-1,:);dummyLNearBF(:,2) = WORKst.(para{i3}).(SINKs{i4}).SEM(:,BF_Pos-1,:);
            dummyLNonBF(:,1) =  WORKst.(para{i3}).(SINKs{i4}).Mean(:,BF_Pos-2,:);dummyLNonBF(:,2) = WORKst.(para{i3}).(SINKs{i4}).SEM(:,BF_Pos-2,:);
            dummyLoffBF(:,1) =  WORKst.(para{i3}).(SINKs{i4}).Mean(:,BF_Pos-3,:);dummyLoffBF(:,2) = WORKst.(para{i3}).(SINKs{i4}).SEM(:,BF_Pos-3,:);
            dummyHNearBF(:,1) =  WORKst.(para{i3}).(SINKs{i4}).Mean(:,BF_Pos+1,:);dummyHNearBF(:,2) = WORKst.(para{i3}).(SINKs{i4}).SEM(:,BF_Pos+1,:);
            dummyHNonBF(:,1) =  WORKst.(para{i3}).(SINKs{i4}).Mean(:,BF_Pos+2,:);dummyHNonBF(:,2) = WORKst.(para{i3}).(SINKs{i4}).SEM(:,BF_Pos+2,:);
            dummyHoffBF(:,1) =  WORKst.(para{i3}).(SINKs{i4}).Mean(:,BF_Pos+3,:);dummyHoffBF(:,2) = WORKst.(para{i3}).(SINKs{i4}).SEM(:,BF_Pos+3,:);
            
            try    %MB what is going on here with this giant try thing?
                subplot(Entries,7,1+((i4-1)*7))
                shadedErrorBar([],dummyLoffBF(:,1),dummyLoffBF(:,2),'lineprops',{'r','markerfacecolor','r'})
                hold on
                try
                    if zscore == 0
                        if Entries == 1
                            M = nanmean(dummyLoffBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyLoffBF(1:3,2));
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyLoffBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                        else
                            M = nanmean(dummyLoffBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyLoffBF(1:3,2));
                            plot([1 size(WORKst.(para{i3}).(SINKs{i4}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKst.(para{i3}).(SINKs{i4}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyLoffBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                            %                 plot([1 size(WORK.(para{i3}).(SINKs{i4}).Mean,1)],[1+PreSTD 1+PreSTD],'k.:');
                            %                 plot([1 size(WORK.(para{i3}).(SINKs{i4}).Mean,1)],[1-PreSTD 1-PreSTD],'k.:');
                            %                 plot(1:length(dummyLoffBF(:,1)) ,(dummyLoffBF(:,1) >= (1+ PreSTD))*1.8 ,'r*');
                        end
                    else
                        plot([1 size(WORKst.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos-3)+Y2(BF_Pos-3) Y(BF_Pos-3)+Y2(BF_Pos-3)],'k.:');
                        plot([1 size(WORKst.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos-3)-Y2(BF_Pos-3) Y(BF_Pos-3)-Y2(BF_Pos-3)],'k.:');
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
                        if Entries == 1
                            M = nanmean(dummyLNonBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyLNonBF(1:3,2));
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyLNonBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                        else
                            M = nanmean(dummyLNonBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyLNonBF(1:3,2));
                            plot([1 size(WORKst.(para{i3}).(SINKs{i4}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKst.(para{i3}).(SINKs{i4}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyLNonBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                            %                 plot([1 size(WORK.(para{i3}).(SINKs{i4}).Mean,1)],[1+PreSTD 1+PreSTD],'k.:');
                            %                 plot([1 size(WORK.(para{i3}).(SINKs{i4}).Mean,1)],[1-PreSTD 1-PreSTD],'k.:');
                            %                 plot(1:length(dummyLoffBF(:,1)) ,(dummyLNonBF(:,1) >= (1+ PreSTD))*1.8 ,'r*');
                        end
                    else
                        plot([1 size(WORKst.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos-2)+Y2(BF_Pos-2) Y(BF_Pos-2)+Y2(BF_Pos-2)],'k.:');
                        plot([1 size(WORKst.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos-2)-Y2(BF_Pos-2) Y(BF_Pos-2)-Y2(BF_Pos-2)],'k.:');
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
                        if Entries == 1
                            M = nanmean(dummyLNearBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyLNearBF(1:3,2));
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyLNearBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                        else
                            M = nanmean(dummyLNearBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyLNearBF(1:3,2));
                            plot([1 size(WORKst.(para{i3}).(SINKs{i4}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKst.(para{i3}).(SINKs{i4}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyLNearBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                            %                 plot([1 size(WORK.(para{i3}).(SINKs{i4}).Mean,1)],[1+PreSTD 1+PreSTD],'k.:');
                            %                 plot([1 size(WORK.(para{i3}).(SINKs{i4}).Mean,1)],[1-PreSTD 1-PreSTD],'k.:');
                            %                 plot(1:length(dummyLoffBF(:,1)) ,(dummyLNearBF(:,1) >= (1+ PreSTD))*1.8 ,'r*');
                        end
                    else
                        plot([1 size(WORKst.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos-1)+Y2(BF_Pos-1) Y(BF_Pos-1)+Y2(BF_Pos-1)],'k.:');
                        plot([1 size(WORKst.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos-1)-Y2(BF_Pos-1) Y(BF_Pos-1)-Y2(BF_Pos-1)],'k.:');
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
                        if Entries == 1
                            M = nanmean(dummyBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyBF(1:3,2));
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                        else
                            M = nanmean(dummyBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyBF(1:3,2));
                            plot([1 size(WORKst.(para{i3}).(SINKs{i4}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKst.(para{i3}).(SINKs{i4}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                            %                 plot([1 size(WORK.(para{i3}).(SINKs{i4}).Mean,1)],[1+PreSTD 1+PreSTD],'k.:');
                            %                 plot([1 size(WORK.(para{i3}).(SINKs{i4}).Mean,1)],[1-PreSTD 1-PreSTD],'k.:');
                            %                 plot(1:length(dummyLoffBF(:,1)) ,(dummyBF(:,1) >= (1+ PreSTD))*1.8 ,'r*');
                        end
                    else
                        plot([1 size(WORKst.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos)+Y2(BF_Pos) Y(BF_Pos)+Y2(BF_Pos)],'k.:');
                        plot([1 size(WORKst.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos)-Y2(BF_Pos) Y(BF_Pos)-Y2(BF_Pos)],'k.:');
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
                        if Entries == 1
                            M = nanmean(dummyHNearBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyHNearBF(1:3,2));
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyHNearBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                        else
                            M = nanmean(dummyHNearBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyHNearBF(1:3,2));
                            plot([1 size(WORKst.(para{i3}).(SINKs{i4}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKst.(para{i3}).(SINKs{i4}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyHNearBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                            %                 plot([1 size(WORK.(para{i3}).(SINKs{i4}).Mean,1)],[1+PreSTD 1+PreSTD],'k.:');
                            %                 plot([1 size(WORK.(para{i3}).(SINKs{i4}).Mean,1)],[1-PreSTD 1-PreSTD],'k.:');
                            %                 plot(1:length(dummyLoffBF(:,1)) ,(dummyHNearBF(:,1) >= (1+ PreSTD))*1.8 ,'r*');
                        end
                    else
                        plot([1 size(WORKst.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos+1)+Y2(BF_Pos+1) Y(BF_Pos+1)+Y2(BF_Pos+1)],'k.:');
                        plot([1 size(WORKst.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos+1)-Y2(BF_Pos+1) Y(BF_Pos+1)-Y2(BF_Pos+1)],'k.:');
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
                        if Entries == 1
                            M = nanmean(dummyHNonBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyHNonBF(1:3,2));
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyHNonBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                        else
                            M = nanmean(dummyHNonBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyHNonBF(1:3,2));
                            plot([1 size(WORKst.(para{i3}).(SINKs{i4}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKst.(para{i3}).(SINKs{i4}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyHNonBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                            %                 plot([1 size(WORK.(para{i3}).(SINKs{i4}).Mean,1)],[1+PreSTD 1+PreSTD],'k.:');
                            %                 plot([1 size(WORK.(para{i3}).(SINKs{i4}).Mean,1)],[1-PreSTD 1-PreSTD],'k.:');
                            %                 plot(1:length(dummyLoffBF(:,1)) ,(dummyHNonBF(:,1) >= (1+ PreSTD))*1.8 ,'r*');
                        end
                    else
                        plot([1 size(WORKst.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos+2)+Y2(BF_Pos+2) Y(BF_Pos+2)+Y2(BF_Pos+2)],'k.:');
                        plot([1 size(WORKst.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos+2)-Y2(BF_Pos+2) Y(BF_Pos+2)-Y2(BF_Pos+2)],'k.:');
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
                        if Entries == 1
                            M = nanmean(dummyHoffBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyHoffBF(1:3,2));
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKst.(para{i3}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyHoffBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                        else
                            M = nanmean(dummyHoffBF(1:3,1));
                            STD2 = 0.25; %STD2 = 4*nanmean(dummyHoffBF(1:3,2));
                            plot([1 size(WORKst.(para{i3}).(SINKs{i4}).Mean,1)],[M+STD2 M+STD2],'k.:');
                            plot([1 size(WORKst.(para{i3}).(SINKs{i4}).Mean,1)],[M-STD2 M-STD2],'k.:');
                            plot(1:length(dummyLoffBF(:,1)) ,(dummyHoffBF(:,1) >= (M+ STD2))*1.8 ,'r*');
                            %                 plot([1 size(WORK.(para{i3}).(SINKs{i4}).Mean,1)],[1+PreSTD 1+PreSTD],'k.:');
                            %                 plot([1 size(WORK.(para{i3}).(SINKs{i4}).Mean,1)],[1-PreSTD 1-PreSTD],'k.:');
                            %                 plot(1:length(dummyLoffBF(:,1)) ,(dummyHoffBF(:,1) >= (1+ PreSTD))*1.8 ,'r*');
                        end
                    else
                        plot([1 size(WORKst.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos+3)+Y2(BF_Pos+3) Y(BF_Pos+3)+Y2(BF_Pos+3)],'k.:');
                        plot([1 size(WORKst.(para{i3}).(SINKs{i4}).Mean,1)],[Y(BF_Pos+3)-Y2(BF_Pos+3) Y(BF_Pos+3)-Y2(BF_Pos+3)],'k.:');
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
        '_ST_tempAnal_' para{i3} '_zscored_' num2str(zscore) '_bin_' ...
        num2str(bin) '_mirror_' num2str(mirror) '.fig']); close all
end