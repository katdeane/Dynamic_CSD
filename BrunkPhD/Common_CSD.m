%% DOCUMENT YOUR CODE PLS THANKS
%%

%  if ~exist('Data','var')
    load('K:\CSD_dynamic_analysis\DATA\Input_Control_7post_n7_Data.mat')
    run('K:\CSD_dynamic_analysis\groups\Input_Control_7post_n7.m')
%   end
%  load('K:\CSD_dynamic_analysis\DATA\Input_YFP_7post_Data.mat')
%  run('K:\CSD_dynamic_analysis\groups\Input_YFP_7post.m')
%  
%  load('K:\CSD_dynamic_analysis\DATA\Input_HighP_7post_Data.mat')
%  run('K:\CSD_dynamic_analysis\groups\Input_HighP_7post.m')
 Group = 'Control ';
 warning ('off','all');
 clear common_layer
 dbstop if error
 addpath('K:\CSD_dynamic_analysis\subfunc')
 
 Threshold=0.0;
 
 FRQ = [125, 250, 500, 1000, 2000, 4000, 8000, 16000, 32000];
% SinkIDs =   {'I_IIE' 'I_IIL' 'IVE' 'IVL' 'VaE' 'VaL' 'VbE' 'VbL' 'VIE' 'VIL'}; % perhaps InfE/L instead of Vb/VIa
  SinkIDs =   {'I_IIL' 'IVE' 'VaE'  'VbE'  'VIE' 'VIL'};
 
 Sinks = fieldnames(Layer);
 
 %%% Find all means & STDs
 for i1 = 1:length(Sinks)
     Chan = char(Layer.(Sinks{i1}));
     
     for i2 = 1:length(Chan)
         if i2 == 1
         common_layer.(Sinks{i1}).mean = nanmean(str2num(Chan(i2,:)));
         common_layer.(Sinks{i1}).STD = nanstd(str2num(Chan(i2,:)));
         else
         common_layer.(Sinks{i1}).mean = [common_layer.(Sinks{i1}).mean, nanmean(str2num(Chan(i2,:)))];
         common_layer.(Sinks{i1}).STD  = [common_layer.(Sinks{i1}).STD,nanstd(str2num(Chan(i2,:)))];
         end
     end
     
 end
 
 %%% make common mean & STD
 for i1 = 1:length(Sinks)
     common_layer.(Sinks{i1}).mean = nanmean(common_layer.(Sinks{i1}).mean);
     common_layer.(Sinks{i1}).STD = nanmean(common_layer.(Sinks{i1}).STD);     
 end
 
 % GS-based sorting 
    MATRIX = nan(length(Data)*length(animals),length(FRQ)*2-1,length(FRQ),length(SinkIDs),3); 
 % all entries, possible frequencies along granular BF, all possible BFs,
 % all sinks, 1= onset latencies, 2= peak amp latencies, 3= peak amp
 
 % ST-based sorting according to Peak Amps
    ST_MATRIX = nan(length(Data)*length(animals),length(FRQ)*2-1,length(FRQ),length(SinkIDs),3,2); 
 % all entries, possible frequencies along granular BF, all possible BFs,
 % all sinks, 1= onset latencies, 2= peak amp latencies, 3= peak amp, shift
 % towards GS
 
 index = 1;
 
 for i1 = 1:length(Data)
     try
         for i2 = 1:length(animals)
             Cur_Frq = Data(i1).(animals{i2}).Frqz;
%              keyboard
%              Cur_BF = Data(i1).(animals{i2}).GS_BF;
             Cur_BFPos = find([Data(i1).(animals{i2}).SinkPeakAmp.IVE] ==...
                 nanmax([Data(i1).(animals{i2}).SinkPeakAmp.IVE]));
%              Cur_BFPos = find(Cur_Frq == Cur_BF);
             START =find(FRQ == Cur_Frq(1));             
             
             for  i3 = 1:length(SinkIDs)
                 for i4 = 1:length(Cur_Frq)  
                     
                     %GS-based
                     % peak onset latency
                 MATRIX(index,length(FRQ)-Cur_BFPos+i4,START-1+i4,i3,1) =...
                     Data(i1).(animals{i2}).Sinkonset.(SinkIDs{i3})(i4); 
                 
                     % peak latency
                 MATRIX(index,length(FRQ)-Cur_BFPos+i4,START-1+i4,i3,2) =...
                     Data(i1).(animals{i2}). SinkPeakLate(i4).(SinkIDs{i3})-200; % dismiss baseline
                 
                    % peak amp
                 MATRIX(index,length(FRQ)-Cur_BFPos+i4,START-1+i4,i3,3) =...
                     Data(i1).(animals{i2}). SinkPeakAmp(i4).(SinkIDs{i3});
                 
                 % ST-based
                     ST = [Data(i1).(animals{i2}). SinkPeakAmp.(SinkIDs{i3})];
                     ST_Pos = find(ST == nanmax(ST));
                     SHIFT = Cur_BFPos-ST_Pos;
% if SHIFT ~= 0
%     keyboard
% end
                     % peak onset latency
                 ST_MATRIX(index,length(FRQ)-Cur_BFPos+i4+SHIFT,START-1+i4,i3,1,1) =...
                     Data(i1).(animals{i2}).Sinkonset.(SinkIDs{i3})(i4); 
                 
                     % peak latency
                 ST_MATRIX(index,length(FRQ)-Cur_BFPos+i4+SHIFT,START-1+i4,i3,2,1) =...
                     Data(i1).(animals{i2}). SinkPeakLate(i4).(SinkIDs{i3})-200; % dismiss baseline
                 
                    % peak amp
                 ST_MATRIX(index,length(FRQ)-Cur_BFPos+i4+SHIFT,START-1+i4,i3,3,1) =...
                     Data(i1).(animals{i2}). SinkPeakAmp(i4).(SinkIDs{i3});
                    
                 ST_MATRIX(index,:,START-1+i4,i3,:,2) =...
                     SHIFT;
                    
                 end
             end
          index = index+1;            
         end
     catch
     end
     
     clear Cur_Frq Cur_BF Cur_BFPos START 
 end

 
 for i1 = 1:size(MATRIX,5)
     
      for i2 = 1:size(MATRIX,4)
          clear Dummy Dummy_SEM p clear Test
         X = sum(~isnan(MATRIX(:,length(FRQ),:,i2,i1)));
         X = X(:); Y = X/sum(X); Z = find(X/sum(X) >= Threshold);         
        
         if i1 == 1
         Name = [SinkIDs{i2} ' Sink onset GS-based ' ];
         elseif i1 == 2
             Name = [SinkIDs{i2} ' Peak amp latency GS-based ' ];
         else
             Name = [SinkIDs{i2} ' Peak amp GS-based ' ];
         end
       
        figure('Name',Name)

        subplot(2,1,1)
        Tit = title(['BFs ' Name]);
        set(Tit,'Interpreter', 'none')

        Dummy = nanmean(MATRIX(:,length(FRQ),:,i2,i1)); Dummy = Dummy (:);
        Dummy_SEM = nanstd(MATRIX(:,length(FRQ),:,i2,i1)); Dummy_SEM = Dummy_SEM (:);Dummy_SEM = Dummy_SEM./sqrt(X);
        shadedErrorBar(Z,Dummy(Z),Dummy_SEM(Z),'lineprops','r')

        Test = MATRIX(:,length(FRQ),Z,i2,i1);
        Test = reshape(Test,[],size(Test,3),1);        
        [p,tbl,stats] =kruskalwallis(Test,[],'off');
        
        if p >= 0.05
        STD =nanstd(MATRIX(:,length(FRQ),:,i2,i1));STD =STD(:);
        D =[FRQ', Dummy, STD];
        D(D == 0) = NaN; Results=weightedfit(D); y = Results.slope*FRQ+Results.Intercept;
        hold on, plot(y,'b:','LineWidth',1)
         txt = ['weighted linear fit y=' num2str(Results.slope) '*x+' num2str(Results.Intercept)];
         text(0.25,nanmax(y)-nanmax(y)*.1,txt)
        else

        STD =nanstd(MATRIX(:,length(FRQ),:,i2,i1));STD =STD(:);
        D =[FRQ', Dummy];
        D(isnan(D)) = 0;
        [~,~,slope, intercept,MSE, R2, S] = logfit(FRQ',Dummy,'logx');
        y = (intercept)+(slope)*log10(FRQ);
         hold on, plot(y,'b:','LineWidth',1)
         txt = ['logarithmic fit y=' num2str(slope) '*log10(x)+' num2str(intercept)];
         text(0.25,nanmax(y)-nanmax(y)*.1,txt)
        end
%         [results,means] = multcompare(stats,'CType','bonferroni');
        for e = 1:length(Z)
            txt = [num2str(round(Y(Z(e))*100,2)) ' %'];
            text(Z(e)-.25,Dummy(Z(e))+Dummy_SEM(Z(e)),txt)
        end
        ylim([0 max(Dummy(Z))+max(Dummy_SEM(Z))]), xlim([0 10])
        xticks(1:1:length(FRQ));
        XLAB = num2str(FRQ); XLAB = split(XLAB);xticklabels(XLAB)
        xlabel('BF in Hz')
        
        txt = ['p-Values from Kruskal Wallis Test = ' num2str(p)];
        text (0.25, max(Dummy(Z))*.1, txt)
        
        subplot(2,1,2)
        clear Dummy Dummy_SEM p clear Test
        hold on 
        Dummy2 = [];
        for e = 1:length(Z)
            X2 = sum(~isnan(MATRIX(:,:,Z(e),i2,i1)));
            Y2 = X2/max(X2);
            Z2 = find(X2/max(X2) >= Threshold); 
            Centre = find(Z2 == length(FRQ));
            
            Dummy = nanmean(MATRIX(:,:,Z(e),i2,i1));
            Dummy_SEM = nanstd(MATRIX(:,:,Z(e),i2,i1));
            Dummy_SEM = Dummy_SEM./sqrt(X2);

            if ~isnan(Dummy(length(FRQ)))
            errorbar(Z(e)+1-Centre:Z(e)-Centre+length(Z2),Dummy(Z2),Dummy_SEM(Z2),'LineWidth',2)
            plot(Z(e),Dummy(length(FRQ)),'r--o','LineWidth',2)
            Dummy2(:,e) = Dummy;
            end
        end      
        [p,tbl,stats] =kruskalwallis(Dummy2,[],'off');
        
        ylim([0 nanmax(nanmax(nanmax(MATRIX(:,:,Z,i2,i1))))]), xlim([0 10])
        xticks(1:1:length(FRQ));
        XLAB = num2str(FRQ); XLAB = split(XLAB);xticklabels(XLAB)
        xlabel('Frequency [Hz]')
        Tit = title(['BF tuning' Name]);
        set(Tit,'Interpreter', 'none')
        
        txt = ['p-Values from Kruskal Wallis Test = ' num2str(p)];
        text (0.25, nanmean(nanmax(nanmean(MATRIX(:,:,Z,i2,i1))))*.1, txt)
%          keyboard
        try
         cd ('figs')
        end
        set(gcf, 'PaperType', 'A4');
        set(gcf, 'PaperOrientation', 'landscape');
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPositionMode', 'auto');
        set(gcf, 'PaperPosition', [0.2 0.1 29 20 ]); 
        savefig(gcf,[Group Name ' GS-based Threshold = ' num2str(Threshold) '.fig'])
        saveas(gcf,[Group Name ' GS-based Threshold = ' num2str(Threshold) '.pdf'])
        close all
      end

%         savefig(gcf,filename)
%         close all
 end
 
  for i1 = 1:size(ST_MATRIX,5)
     
      for i2 = 1:size(ST_MATRIX,4)

          clear Dummy Dummy_SEM p clear Test
         X = sum(~isnan(ST_MATRIX(:,length(FRQ),:,i2,i1)));
         X = X(:); Y = X/sum(X); Z = find(X/sum(X) >= Threshold);         
        
         if i1 == 1
         Name = [SinkIDs{i2} ' Sink onset ST-based ' ];
         elseif i1 == 2
             Name = [SinkIDs{i2} ' Peak amp latency ST-based ' ];
         else
             Name = [SinkIDs{i2} ' Peak amp ST-based ' ];
         end
       
        figure('Name',Name)
%         set(0,'DefaultFigureWindowStyle','docked')
        subplot(2,1,1)
        Tit = title(['BFs ' Name]);
        set(Tit,'Interpreter', 'none')

        Dummy = nanmean(ST_MATRIX(:,length(FRQ),:,i2,i1)); Dummy = Dummy (:);
        Dummy_SEM = nanstd(ST_MATRIX(:,length(FRQ),:,i2,i1)); Dummy_SEM = Dummy_SEM (:);Dummy_SEM = Dummy_SEM./sqrt(X);
%         try
        shadedErrorBar(Z,Dummy(Z),Dummy_SEM(Z),'lineprops','r')

        Test = ST_MATRIX(:,length(FRQ),Z,i2,i1);
        Test = reshape(Test,[],size(Test,3),1);        
        [p,tbl,stats] =kruskalwallis(Test,[],'off');
        
        if p >= 0.05
        D =[FRQ', Dummy, Dummy_SEM];
        D(D == 0) = NaN; Results=weightedfit(D); y = Results.slope*FRQ+Results.Intercept;
        hold on, plot(y,'b:','LineWidth',1)
         txt = ['weighted linear fit y=' num2str(Results.slope) '*x+' num2str(Results.Intercept)];
         text(0.25,nanmax(y)-nanmax(y)*.1,txt)
        else
        STD =nanstd(MATRIX(:,length(FRQ),:,i2,i1));STD =STD(:);
        D =[FRQ', Dummy];
        D(isnan(D)) = 0;
        [~,~,slope, intercept,MSE, R2, S] = logfit(FRQ',Dummy,'logx');
        y = (intercept)+(slope)*log10(FRQ);
         hold on, plot(y,'b:','LineWidth',1)
         txt = ['logarithmic fit y=' num2str(slope) '*log10(x)+' num2str(intercept)];
         text(0.25,nanmax(y)-nanmax(y)*.1,txt)
        end
        
        
%         [results,means] = multcompare(stats,'CType','bonferroni');
        for e = 1:length(Z)
            txt = [num2str(round(Y(Z(e))*100,2)) ' %'];
            text(Z(e)-.25,Dummy(Z(e))+Dummy_SEM(Z(e)),txt)
        end
        ylim([0 max(Dummy(Z))+max(Dummy_SEM(Z))]), xlim([0 10])
        xticks(1:1:length(FRQ));
        XLAB = num2str(FRQ); XLAB = split(XLAB);xticklabels(XLAB)
        xlabel('BF in Hz')
        
        txt = ['p-Values from Kruskal Wallis Test = ' num2str(p)];
        text (0.25, max(Dummy(Z))*.1, txt)
%         end
        subplot(2,1,2)
        clear Dummy Dummy_SEM p clear Test
        hold on 
        Dummy2 = [];
        for e = 1:length(Z)
            X2 = sum(~isnan(ST_MATRIX(:,:,Z(e),i2,i1)));
            Y2 = X2/max(X2);
            Z2 = find(X2/max(X2) >= Threshold); 
            Centre = find(Z2 == length(FRQ));
            
            Dummy = nanmean(ST_MATRIX(:,:,Z(e),i2,i1));
            Dummy_SEM = nanstd(ST_MATRIX(:,:,Z(e),i2,i1));
            Dummy_SEM = Dummy_SEM./sqrt(X2);

            if ~isnan(Dummy(length(FRQ)))
                try
            errorbar(Z(e)+1-Centre:Z(e)-Centre+length(Z2),Dummy(Z2),Dummy_SEM(Z2),'LineWidth',2)
            plot(Z(e),Dummy(Z2(length(FRQ))),'r--o','LineWidth',2)
            Dummy2(:,e) = Dummy;
                end
            end
        end      
        [p,tbl,stats] =kruskalwallis(Dummy2,[],'off');

        ylim([0 nanmax(nanmax(nanmax(ST_MATRIX(:,:,Z,i2,i1))))]), xlim([0 10])
        xticks(1:1:length(FRQ));
        XLAB = num2str(FRQ); XLAB = split(XLAB);xticklabels(XLAB)
        xlabel('Frequency [Hz]')
        Tit = title(['BF tuning' Name]);
        set(Tit,'Interpreter', 'none')
        
        txt = ['p-Values from Kruskal Wallis Test = ' num2str(p)];
        text (0.25, nanmean(nanmax(nanmean(ST_MATRIX(:,:,Z,i2,i1))))*.1, txt)

%         keyboard
        try
         cd ('figs')
        end
        set(gcf, 'PaperType', 'A4');
        set(gcf, 'PaperOrientation', 'landscape');
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPositionMode', 'auto');
        set(gcf, 'PaperPosition', [0.2 0.1 29 20 ]); 

        savefig(gcf,[Group Name ' ST-based Threshold = ' num2str(Threshold) '.fig'])
        saveas(gcf,[Group Name ' ST-based Threshold = ' num2str(Threshold) '.pdf'])
        close all
      end

%         savefig(gcf,filename)
%         close all
 end
 
%  keyboard
 I_II = (common_layer.I_IIL.mean+common_layer.I_IIL.STD +...
     common_layer.IVE.mean-common_layer.IVE.STD)/2;
 IV = (common_layer.IVE.mean+common_layer.IVE.STD +...
     common_layer.VaE.mean-common_layer.VaE.STD)/2;
 Va = (common_layer.VaE.mean+common_layer.VaE.STD +...
     common_layer.VbE.mean-common_layer.VbE.STD)/2;
 Vb = (common_layer.VbE.mean+common_layer.VbE.STD +...
     common_layer.VIE.mean-common_layer.VIE.STD)/2;
 VIa = (common_layer.VIE.mean+common_layer.VIE.STD +...
     common_layer.VIL.mean-common_layer.VIL.STD)/2;
for i4 =1:3 
 figure('Name','BFs GS sorted')
 set( gca, 'Color', [0.5,0.5,0.5] )
 hold on
 CM = {'black' 'red' 'blue' 'green' 'cyan' 'magenta' 'yellow'};
 index =1;
% if i4 < 3
for i1=2:8   

    for i2=9 %BF

    for i3 = 1:length(SinkIDs)
        Dummy = nanmean((MATRIX(:,:,i1,i3,i4)));
        Dummy_SEM = nanstd((MATRIX(:,:,i1,i3,i4)))./sqrt(sum(~isnan(MATRIX(:,:,i1,i3,i4))));        
        errorbarxy(Dummy(i2),common_layer.(SinkIDs{i3}).mean,Dummy_SEM(i2),common_layer.(SinkIDs{i3}).STD,CM{index},'LineWidth',2);
    end
    set(gca,'ydir','reverse')
    ylim([1,32])
    xlim([0,400])
    index = index +1;
    end
    
    
end
% else
% %     keyboard
% end
plot ([1 400],[I_II I_II],'k--'),t =text(380,mean([1 ,I_II]), 'I/II'); t.FontSize = 20; t.FontWeight = 'bold';
plot ([1 400],[IV IV],'k--'),t =text(380,mean([I_II ,IV]), 'III/IV'); t.FontSize = 20; t.FontWeight = 'bold';
plot ([1 400],[Va Va],'k--'),t =text(380,mean([IV ,Va]), 'Va'); t.FontSize = 20; t.FontWeight = 'bold';
plot ([1 400],[Vb Vb],'k--'),t =text(380,mean([Va ,Vb]), 'Vb'); t.FontSize = 20; t.FontWeight = 'bold';
plot ([1 400],[VIa VIa],'k--'),t =text(380,mean([Vb ,VIa]), 'VIa'); t.FontSize = 20; t.FontWeight = 'bold';
t =text(380,mean([VIa ,32]), 'VIb'); t.FontSize = 20; t.FontWeight = 'bold';

t = text(350,2, '250 Hz'); t.Color = 'black';
t = text(350,3, '500 Hz'); t.Color = 'red';
t = text(350,4, '1000 Hz'); t.Color = 'blue';
t = text(350,5, '2000 Hz'); t.Color = 'green';
t = text(350,6, '4000 Hz'); t.Color = 'cyan';
t = text(350,7, '8000 Hz'); t.Color = 'magenta';
t = text(350,8, '16000 Hz'); t.Color = 'yellow';
if i4 == 3
     xlim([0,0.01])
end
end

for i4 =1:3 
 figure('Name','BFs ST sorted')
 set( gca, 'Color', [0.5,0.5,0.5] )
 hold on
 CM = {'black' 'red' 'blue' 'green' 'cyan' 'magenta' 'yellow'};
 index =1; 
 
for i1=2:8   

    for i2=9

    for i3 = 1:length(SinkIDs)
        Dummy = nanmean((ST_MATRIX(:,:,i1,i3,i4)));
        Dummy_SEM = nanstd((ST_MATRIX(:,:,i1,i3,i4)))./sqrt(sum(~isnan(ST_MATRIX(:,:,i1,i3,i4))));        
        errorbarxy(Dummy(i2),common_layer.(SinkIDs{i3}).mean,Dummy_SEM(i2),common_layer.(SinkIDs{i3}).STD,CM{index},'LineWidth',2);
    end
    set(gca,'ydir','reverse')
    ylim([1,32])
    xlim([0,400])
    index = index +1;
    end
end
plot ([1 400],[I_II I_II],'k--'),t =text(380,mean([1 ,I_II]), 'I/II'); t.FontSize = 20; t.FontWeight = 'bold';
plot ([1 400],[IV IV],'k--'),t =text(380,mean([I_II ,IV]), 'III/IV'); t.FontSize = 20; t.FontWeight = 'bold';
plot ([1 400],[Va Va],'k--'),t =text(380,mean([IV ,Va]), 'Va'); t.FontSize = 20; t.FontWeight = 'bold';
plot ([1 400],[Vb Vb],'k--'),t =text(380,mean([Va ,Vb]), 'Vb'); t.FontSize = 20; t.FontWeight = 'bold';
plot ([1 400],[VIa VIa],'k--'),t =text(380,mean([Vb ,VIa]), 'VIa'); t.FontSize = 20; t.FontWeight = 'bold';
t =text(380,mean([VIa ,32]), 'VIb'); t.FontSize = 20; t.FontWeight = 'bold';

t = text(350,2, '250 Hz'); t.Color = 'black';
t = text(350,3, '500 Hz'); t.Color = 'red';
t = text(350,4, '1000 Hz'); t.Color = 'blue';
t = text(350,5, '2000 Hz'); t.Color = 'green';
t = text(350,6, '4000 Hz'); t.Color = 'cyan';
t = text(350,7, '8000 Hz'); t.Color = 'magenta';
t = text(350,8, '16000 Hz'); t.Color = 'yellow';
if i4 == 3
    camroll(90)
     xlim([0,0.01])
t = text(0.009,2, '250 Hz'); t.Color = 'black';
t = text(0.009,4, '500 Hz'); t.Color = 'red';
t = text(0.009,6, '1000 Hz'); t.Color = 'blue';
t = text(0.009,8, '2000 Hz'); t.Color = 'green';
t = text(0.009,10, '4000 Hz'); t.Color = 'cyan';
t = text(0.009,12, '8000 Hz'); t.Color = 'magenta';
t = text(0.009,14, '16000 Hz'); t.Color = 'yellow';
end

end


% 
% 
% 
% 
% Y =[nanmean(MATRIX2(:,9,:,3))];
% Y2 =[nanstd(MATRIX2(:,9,:,3))./sqrt(sum(~isnan(MATRIX2(:,9,:,3))))];
% 
% figure
% set(0,'DefaultFigureWindowStyle','docked')
% plot(ST),figure
% set(0,'DefaultFigureWindowStyle','docked')
% plot(FRQ)
