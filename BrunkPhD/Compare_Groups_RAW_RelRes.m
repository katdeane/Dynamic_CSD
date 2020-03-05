clear,dbstop if error,warning ('off','all');
mkdir Figs_Group_comparison
Home = pwd; FigFolder = [Home '\Figs_Group_comparison'];
cd ('K:\CSD_dynamic_analysis\DATA\Output')
addpath('K:\CSD_dynamic_analysis\subfunc')

AllOpto = load('Output_Input_ALL_OPTO_7post_Data_Threshold_25_Zscore_0_binned_1.mat');
HighP = load('Output_Input_HighP_7post_Data_Threshold_25_Zscore_0_binned_1.mat');
Control = load('Output_Input_Control_7post_Data_Threshold_25_Zscore_0_binned_1.mat');
YFP = load('Output_Input_YFP_7post_Data_Threshold_25_Zscore_0_binned_1.mat');

Groups = {'Control','YFP','AllOpto','HighP'};

para = {'Full_RMS_RELRES','Early_RMS_RELRES','Late_RMS_RELRES',...
    };

Cond = [3, 4, 10]; %Pre3, Combi, Post6
Cond2 = {'Pre','Combi','Post'};
Sinks = AllOpto.Data.Sinks;

for i1 = 1:length(para)
    
    if isstruct(AllOpto.Data.GS_based(1).(para{i1}))
        Entries = length(AllOpto.Data.Sinks);
    else
        Entries = 1;
    end
    sorting = {'GS_based'};
    for i2 = 1:length(sorting)    
        for i3 = 1:Entries
            if Entries == 1            
            CondName = para{i1};
            else
            CondName = [para{i1} ' Sink ' Sinks{i3}];
            end

            h= figure('Name',['RAW_' para{i1} '_' sorting{i2} '_' CondName ],...
                'PaperType','a4letter',...
                'PaperOrientation','Portrait',...
                'PaperUnits','centimeters', ...
                'PaperPosition',[0.63 0.63 19.72 28.41],...    
                'PaperSize',[20.98 29.68]);
    
            for i4 = 1:5 %Low Non, Low Near, BF, High Near, High Non 
                subplot(5,1,i4)
                clear X XL, X = []; XL = {};
                for i5 = 1:length(Groups)% number of groups
% keyboard
if Entries == 1
    clear PreNorm Pre Combi Post
   PreNorm(:,:,1) = eval([(Groups{i5}) '.Data.' (sorting{i2}) '(1).' (para{i1})]);
   PreNorm(:,:,2) = eval([(Groups{i5}) '.Data.' (sorting{i2}) '(2).' (para{i1})]);
   PreNorm(:,:,3) = eval([(Groups{i5}) '.Data.' (sorting{i2}) '(3).' (para{i1})]);
   Pre(:,:) = eval([(Groups{i5}) '.Data.' (sorting{i2}) '(' num2str(Cond(1)) ').' (para{i1})]);
   Combi(:,:) = eval([(Groups{i5}) '.Data.' (sorting{i2}) '(' num2str(Cond(2)) ').' (para{i1})]);
   Post(:,:) = eval([(Groups{i5}) '.Data.' (sorting{i2}) '(' num2str(Cond(3)) ').' (para{i1})]);
else
    clear PreNorm Pre Combi Post
   PreNorm(:,:,1) = eval([(Groups{i5}) '.Data.' (sorting{i2}) '(1).' (para{i1}) '.' (Sinks{i3})]);
   PreNorm(:,:,2) = eval([(Groups{i5}) '.Data.' (sorting{i2}) '(2).' (para{i1}) '.' (Sinks{i3})]);
   PreNorm(:,:,3) = eval([(Groups{i5}) '.Data.' (sorting{i2}) '(3).' (para{i1}) '.' (Sinks{i3})]);
   Pre(:,:) = eval([(Groups{i5}) '.Data.' (sorting{i2}) '(' num2str(Cond(1)) ').' (para{i1}) '.' (Sinks{i3})]);
   Combi(:,:) = eval([(Groups{i5}) '.Data.' (sorting{i2}) '(' num2str(Cond(2)) ').' (para{i1}) '.' (Sinks{i3})]);
   Post(:,:) = eval([(Groups{i5}) '.Data.' (sorting{i2}) '(' num2str(Cond(3)) ').' (para{i1}) '.' (Sinks{i3})]);
end

PreNorm = PreNorm (:,eval([Groups{i5} '.Data.BF_Pos'])-2:eval([Groups{i5} '.Data.BF_Pos'])+2,:);
PreNorm = nanmean(PreNorm,3);
Pre = Pre(:,eval([Groups{i5} '.Data.BF_Pos'])-2:eval([Groups{i5} '.Data.BF_Pos'])+2);
Combi = Combi(:,eval([Groups{i5} '.Data.BF_Pos'])-2:eval([Groups{i5} '.Data.BF_Pos'])+2);
Post = Post(:,eval([Groups{i5} '.Data.BF_Pos'])-2:eval([Groups{i5} '.Data.BF_Pos'])+2);
PreNorm = 1;

Pre = Pre./PreNorm*100;
Combi = Combi./PreNorm*100;
Post = Post./PreNorm*100;

y = [Pre(:,i4), Combi(:,i4),Post(:,i4)];
sem =(nanstd(y)./sqrt(sum(~isnan(y))));
y = nanmean(y);
x = [i5-.2 i5 i5+.2];
X = [X x];
XL = [XL Cond2];

hold on
bar(x,y)
% errorbar(x,y,nanstd([Pre(:,i4), Combi(:,i4),Post(:,i4)]))
errorb(x,y,sem)
xticks(X)
xticklabels(XL)
title (['BF shift' num2str(i4-3)])
% ylim([0.25 3.5])

hold off

lgd =legend(Groups,'Location','best','AutoUpdate','off');
lgd.FontSize = 6;
                end
                for i5 = 1:length(Groups)% number of groups
% keyboard
if Entries == 1
    clear PreNorm Pre Combi Post
   PreNorm(:,:,1) = eval([(Groups{i5}) '.Data.' (sorting{i2}) '(1).' (para{i1})]);
   PreNorm(:,:,2) = eval([(Groups{i5}) '.Data.' (sorting{i2}) '(2).' (para{i1})]);
   PreNorm(:,:,3) = eval([(Groups{i5}) '.Data.' (sorting{i2}) '(3).' (para{i1})]);
   Pre(:,:) = eval([(Groups{i5}) '.Data.' (sorting{i2}) '(' num2str(Cond(1)) ').' (para{i1})]);
   Combi(:,:) = eval([(Groups{i5}) '.Data.' (sorting{i2}) '(' num2str(Cond(2)) ').' (para{i1})]);
   Post(:,:) = eval([(Groups{i5}) '.Data.' (sorting{i2}) '(' num2str(Cond(3)) ').' (para{i1})]);
else
    clear PreNorm Pre Combi Post
   PreNorm(:,:,1) = eval([(Groups{i5}) '.Data.' (sorting{i2}) '(1).' (para{i1}) '.' (Sinks{i3})]);
   PreNorm(:,:,2) = eval([(Groups{i5}) '.Data.' (sorting{i2}) '(2).' (para{i1}) '.' (Sinks{i3})]);
   PreNorm(:,:,3) = eval([(Groups{i5}) '.Data.' (sorting{i2}) '(3).' (para{i1}) '.' (Sinks{i3})]);
   Pre(:,:) = eval([(Groups{i5}) '.Data.' (sorting{i2}) '(' num2str(Cond(1)) ').' (para{i1}) '.' (Sinks{i3})]);
   Combi(:,:) = eval([(Groups{i5}) '.Data.' (sorting{i2}) '(' num2str(Cond(2)) ').' (para{i1}) '.' (Sinks{i3})]);
   Post(:,:) = eval([(Groups{i5}) '.Data.' (sorting{i2}) '(' num2str(Cond(3)) ').' (para{i1}) '.' (Sinks{i3})]);
end

PreNorm = PreNorm (:,eval([Groups{i5} '.Data.BF_Pos'])-2:eval([Groups{i5} '.Data.BF_Pos'])+2,:);
PreNorm = nanmean(PreNorm,3);
Pre = Pre(:,eval([Groups{i5} '.Data.BF_Pos'])-2:eval([Groups{i5} '.Data.BF_Pos'])+2);
Combi = Combi(:,eval([Groups{i5} '.Data.BF_Pos'])-2:eval([Groups{i5} '.Data.BF_Pos'])+2);
Post = Post(:,eval([Groups{i5} '.Data.BF_Pos'])-2:eval([Groups{i5} '.Data.BF_Pos'])+2);
PreNorm = 1;

Pre = Pre./PreNorm*100;
Combi = Combi./PreNorm*100;
Post = Post./PreNorm*100;

y = [Pre(:,i4), Combi(:,i4),Post(:,i4)];
sem =(nanstd(y)./sqrt(sum(~isnan(y))));
y = nanmean(y);
x = [i5-.2 i5 i5+.2];
X = [X x];
XL = [XL Cond2];
% keyboard
hold on
errorbar(x,y,nanstd([Pre(:,i4), Combi(:,i4),Post(:,i4)]),'LineWidth',2,'Color','red')
[h,p,ci,stats] = ttest2(Pre(:,i4),Combi(:,i4));
if p <= 0.01
     
    plot ([i5-.2 i5],[2 2],'Color','black','LineWidth',2)
    plot([(i5-.2 + i5)/2 (i5-.2 + i5)/2+.05],[1.8 1.8],'r*')
    txt = ['p = ' num2str(p)];
    t=text(i5-.2,1.5,txt);
    t.FontSize = 6;
elseif p <= 0.05
    plot ([i5-.2 i5],[2 2],'Color','black','LineWidth',2)
    plot((i5-.2 + i5)/2,1.8,'r*')
    txt = ['p = ' num2str(p)];
    t=text(i5-.2, 1.5,txt);
    t.FontSize = 6;
end

[h,p] = ttest2(Pre(:,i4),Post(:,i4));
if p <= 0.01
%      keyboard
    plot ([i5-.2 i5+.2],[2.1 2.1],'Color','black','LineWidth',2)
    plot([i5 i5+.05],[2.3 2.3],'r*')
    txt = ['p = ' num2str(p)];
    t=text(i5,2.4,txt);
    t.FontSize = 6;
    
elseif p <= 0.05
    plot ([i5-.2 i5+.2],[2.1 2.1],'Color','black','LineWidth',2)
    plot(i5,2.3,'r*')
    txt = ['p = ' num2str(p)];
    t=text(i5,2.4,txt);
    t.FontSize = 6;
end
hold off
                end               
            end
%             keyboard
            t = text(1,23,['RAW_' para{i1} ' ' sorting{i2} ' ' CondName]);
            s = t.FontSize;
            t.FontSize = 12;
            
            cd (FigFolder)
            saveas(gcf,['RAW_' para{i1} '_' sorting{i2} '_' CondName],'pdf');
            close all
            
%   keyboard          
%                 h= figure('Name',[para{i1} '_' sorting{i2} '_' CondName '_BandWidth'],...
%                 'PaperType','a4letter',...
%                 'PaperOrientation','Portrait',...
%                 'PaperUnits','centimeters', ...
%                 'PaperPosition',[0.63 0.63 19.72 28.41],...    
%                 'PaperSize',[20.98 29.68]);        
%             
            
            
            
        end
    end
end

