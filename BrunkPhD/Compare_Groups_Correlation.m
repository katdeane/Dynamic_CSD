clear,dbstop if error,warning ('off','all');
mkdir Figs_Group_comparison
Home = pwd; FigFolder = [Home '\Figs_Group_comparison'];
cd ('K:\CSD_dynamic_analysis\DATA\Output')
addpath('K:\CSD_dynamic_analysis\subfunc')

% HighP = load('Output_Input_ALL_OPTO_7post_Data_Threshold_25_Zscore_0_binned_1.mat');
HighP = load('Output_Input_HighP_7post_Data_Threshold_25_Zscore_0_binned_1.mat');
Control = load('Output_Input_Control_7post_Data_Threshold_25_Zscore_0_binned_1.mat');
% Control = load('Output_Input_Control_7post_Data_Threshold_25_Zscore_0_binned_1.mat');

% Groups = {'Control','Control','HighP','HighP'};
Groups = {'Control','Control','HighP','HighP'};

% para1 = {'SinkRMS','tempSinkRMS','SinkPeakAmp','Full_RMS_AVREC','Early_RMS_AVREC',...
%     'Late_RMS_AVREC','Full_RMS_RELRES','Early_RMS_RELRES','Late_RMS_RELRES'};
% 
% para2 = {'Full_RMS_RELRES','Early_RMS_RELRES','Late_RMS_RELRES',...
%     'Full_RMS_ABSRES','Early_RMS_ABSRES','Late_RMS_ABSRES','LPs'};

para1 = {'Full_RMS_AVREC','Early_RMS_AVREC',...
    'Late_RMS_AVREC'};

para2 = {'Full_RMS_RELRES','Early_RMS_RELRES','Late_RMS_RELRES',...
    };

% Combis

Cond = [3, 4, 10]; %Pre3, Combi, Post6
Cond2 = {'Pre','Combi','Post'};
sorting = {'GS_based','ST_based'};
Sinks = HighP.Data.Sinks;

%% Do correlations para1 = yaxis, para2 = xaxis
Output = struct([Groups{1} '_1'],[],Groups{2},[],[Groups{3} '_1'],[],Groups{4},[]);


for i0 = 1:length(sorting)
    
for i1 = 1:length(para1)
    
    if isstruct(HighP.Data.GS_based(1).(para1{i1}))
        Entries = length(HighP.Data.Sinks);
    else
        Entries = 1;
    end
    
    clear PreP1 CombiP1 PostP1
    PreP1 = NaN(size(HighP.Data.GS_based(1).Full_RMS_AVREC,1),5,4,Entries);
    CombiP1 = NaN(size(HighP.Data.GS_based(1).Full_RMS_AVREC,1),5,4,Entries);
    PostP1 = NaN(size(HighP.Data.GS_based(1).Full_RMS_AVREC,1),5,4,Entries);
    
    for i2 = 1:Entries % create Pre, Combi, Post Matrix for Para1 (animals,Frqzbins,Groups,Sinks)        
      
         try 
            % Control 
            PreP1(1:size(Control.Data.GS_based(1).Full_RMS_AVREC,1),:,1,i2) =...
                Control.Data.(sorting{i0})(Cond(1)).(para1{i1}).(Sinks{i2})(:,Control.Data.BF_Pos-2:Control.Data.BF_Pos+2);
            CombiP1(1:size(Control.Data.GS_based(1).Full_RMS_AVREC,1),:,1,i2)=...
                Control.Data.(sorting{i0})(Cond(2)).(para1{i1}).(Sinks{i2})(:,Control.Data.BF_Pos-2:Control.Data.BF_Pos+2);
            PostP1(1:size(Control.Data.GS_based(1).Full_RMS_AVREC,1),:,1,i2)=...
                Control.Data.(sorting{i0})(Cond(3)).(para1{i1}).(Sinks{i2})(:,Control.Data.BF_Pos-2:Control.Data.BF_Pos+2);

            % Control
            PreP1(1:size(Control.Data.GS_based(1).Full_RMS_AVREC,1),:,2,i2) =...
                Control.Data.(sorting{i0})(Cond(1)).(para1{i1}).(Sinks{i2})(:,Control.Data.BF_Pos-2:Control.Data.BF_Pos+2);
            CombiP1(1:size(Control.Data.GS_based(1).Full_RMS_AVREC,1),:,2,i2)=...
                Control.Data.(sorting{i0})(Cond(2)).(para1{i1}).(Sinks{i2})(:,Control.Data.BF_Pos-2:Control.Data.BF_Pos+2);
            PostP1(1:size(Control.Data.GS_based(1).Full_RMS_AVREC,1),:,2,i2)=...
                Control.Data.(sorting{i0})(Cond(3)).(para1{i1}).(Sinks{i2})(:,Control.Data.BF_Pos-2:Control.Data.BF_Pos+2);

            % All_Opto
            PreP1(1:size(HighP.Data.GS_based(1).Full_RMS_AVREC,1),:,3,i2) =...
                HighP.Data.(sorting{i0})(Cond(1)).(para1{i1}).(Sinks{i2})(:,HighP.Data.BF_Pos-2:HighP.Data.BF_Pos+2);
            CombiP1(1:size(HighP.Data.GS_based(1).Full_RMS_AVREC,1),:,3,i2)=...
                HighP.Data.(sorting{i0})(Cond(2)).(para1{i1}).(Sinks{i2})(:,HighP.Data.BF_Pos-2:HighP.Data.BF_Pos+2);
            PostP1(1:size(HighP.Data.GS_based(1).Full_RMS_AVREC,1),:,3,i2)=...
                HighP.Data.(sorting{i0})(Cond(3)).(para1{i1}).(Sinks{i2})(:,HighP.Data.BF_Pos-2:HighP.Data.BF_Pos+2);

            % HighP
            PreP1(1:size(HighP.Data.GS_based(1).Full_RMS_AVREC,1),:,4,i2) =...
                HighP .Data.(sorting{i0})(Cond(1)).(para1{i1}).(Sinks{i2})(:,HighP .Data.BF_Pos-2:HighP .Data.BF_Pos+2);
            CombiP1(1:size(HighP.Data.GS_based(1).Full_RMS_AVREC,1),:,4,i2)=...
                HighP .Data.(sorting{i0})(Cond(2)).(para1{i1}).(Sinks{i2})(:,HighP .Data.BF_Pos-2:HighP .Data.BF_Pos+2);
            PostP1(1:size(HighP.Data.GS_based(1).Full_RMS_AVREC,1),:,4,i2)=...
                HighP .Data.(sorting{i0})(Cond(3)).(para1{i1}).(Sinks{i2})(:,HighP .Data.BF_Pos-2:HighP .Data.BF_Pos+2);
         catch
            % Control
            PreP1(1:size(Control.Data.GS_based(1).Full_RMS_AVREC,1),:,1,i2) =...
                Control.Data.(sorting{i0})(Cond(1)).(para1{i1})(:,Control.Data.BF_Pos-2:Control.Data.BF_Pos+2);
            CombiP1(1:size(Control.Data.GS_based(1).Full_RMS_AVREC,1),:,1,i2)= ...
                Control.Data.(sorting{i0})(Cond(2)).(para1{i1})(:,Control.Data.BF_Pos-2:Control.Data.BF_Pos+2);
            PostP1(1:size(Control.Data.GS_based(1).Full_RMS_AVREC,1),:,1,i2)=...
                Control.Data.(sorting{i0})(Cond(3)).(para1{i1})(:,Control.Data.BF_Pos-2:Control.Data.BF_Pos+2);  

            % Control
            PreP1(1:size(Control.Data.GS_based(1).Full_RMS_AVREC,1),:,2,i2) =...
                Control.Data.(sorting{i0})(Cond(1)).(para1{i1})(:,Control.Data.BF_Pos-2:Control.Data.BF_Pos+2);
            CombiP1(1:size(Control.Data.GS_based(1).Full_RMS_AVREC,1),:,2,i2)=...
                Control.Data.(sorting{i0})(Cond(2)).(para1{i1})(:,Control.Data.BF_Pos-2:Control.Data.BF_Pos+2);
            PostP1(1:size(Control.Data.GS_based(1).Full_RMS_AVREC,1),:,2,i2)=...
                Control.Data.(sorting{i0})(Cond(3)).(para1{i1})(:,Control.Data.BF_Pos-2:Control.Data.BF_Pos+2);

            % All_Opto
            PreP1(1:size(HighP.Data.GS_based(1).Full_RMS_AVREC,1),:,3,i2) =...
                HighP.Data.(sorting{i0})(Cond(1)).(para1{i1})(:,HighP.Data.BF_Pos-2:HighP.Data.BF_Pos+2);
            CombiP1(1:size(HighP.Data.GS_based(1).Full_RMS_AVREC,1),:,3,i2)=...
                HighP.Data.(sorting{i0})(Cond(2)).(para1{i1})(:,HighP.Data.BF_Pos-2:HighP.Data.BF_Pos+2);
            PostP1(1:size(HighP.Data.GS_based(1).Full_RMS_AVREC,1),:,3,i2)=...
                HighP.Data.(sorting{i0})(Cond(3)).(para1{i1})(:,HighP.Data.BF_Pos-2:HighP.Data.BF_Pos+2);

            % HighP
            PreP1(1:size(HighP.Data.GS_based(1).Full_RMS_AVREC,1),:,4,i2) =...
                HighP .Data.(sorting{i0})(Cond(1)).(para1{i1})(:,HighP .Data.BF_Pos-2:HighP .Data.BF_Pos+2);
            CombiP1(1:size(HighP.Data.GS_based(1).Full_RMS_AVREC,1),:,4,i2)=...
                HighP .Data.(sorting{i0})(Cond(2)).(para1{i1})(:,HighP .Data.BF_Pos-2:HighP .Data.BF_Pos+2);
            PostP1(1:size(HighP.Data.GS_based(1).Full_RMS_AVREC,1),:,4,i2)=...
                HighP .Data.(sorting{i0})(Cond(3)).(para1{i1})(:,HighP .Data.BF_Pos-2:HighP .Data.BF_Pos+2);
         end

    end

        for i2 = 1:length(para2) % create Pre, Combi, Post Matrix for Para1 (animals,Frqzbins,Groups,Sinks)  
        if ~strcmp(para1{i1},para2{i2}) % only compare different parameters  
%               keyboard
             if ~strcmp(para2{i2},'LPs')
                if isstruct(HighP.Data.GS_based(1).(para2{i2}))
                    Entries = length(HighP.Data.Sinks);
                else
                    Entries = 1;
                end
                Entries = 1;
             end
            
            if strcmp(para2{i2},'LPs')% get Lever pressing max if entries exist

                ID_Con =Control.Data.names;
                ID_Control =Control.Data.names;
                ID_HighP =HighP.Data.names;
                ID_HighP = HighP.Data.names;
                
                ID_Groups = {'ID_Con','ID_Control','ID_HighP','ID_HighP'};
                LP_Groups = {'LP_Con','LP_Control','LP_HighP','LP_HighP'};
                
                LP_Con = []; LP_Control = []; LP_HighP = []; LP_HighP = [];

                for i3 = 1:length(ID_Groups)
                    for i4 = 1:length(eval(ID_Groups{i3}))
                        CA = eval(ID_Groups{i3});
                        LP = GOTs_LPs(CA{i4});
                        if i3 ==1
                            LP_Con = [eval(LP_Groups{i3}) LP];
                        elseif i3 ==2
                            LP_Control = [eval(LP_Groups{i3}) LP];
                        elseif i3 ==3
                            LP_HighP = [eval(LP_Groups{i3}) LP];
                        else
                            LP_HighP = [eval(LP_Groups{i3}) LP];
                        end
                    end
                end
 
            else
% keyboard
            clear PreP2 CombiP2 PostP2
            PreP2 = NaN(size(HighP.Data.GS_based(1).Full_RMS_AVREC,1),5,4,Entries);
            CombiP2 = NaN(size(HighP.Data.GS_based(1).Full_RMS_AVREC,1),5,4,Entries);
            PostP2 = NaN(size(HighP.Data.GS_based(1).Full_RMS_AVREC,1),5,4,Entries);
            
            for i3 = 1:Entries

                 try 
                    % Control 
                    PreP2(1:size(Control.Data.GS_based(1).Full_RMS_AVREC,1),:,1,i3) =...
                        Control.Data.(sorting{i0})(Cond(1)).(para2{i2}).(Sinks{i3})(:,Control.Data.BF_Pos-2:Control.Data.BF_Pos+2);
                    CombiP2(1:size(Control.Data.GS_based(1).Full_RMS_AVREC,1),:,1,i3)=...
                        Control.Data.(sorting{i0})(Cond(2)).(para22{i2}).(Sinks{i3})(:,Control.Data.BF_Pos-2:Control.Data.BF_Pos+2);
                    PostP2(1:size(Control.Data.GS_based(1).Full_RMS_AVREC,1),:,1,i3)=...
                        Control.Data.(sorting{i0})(Cond(3)).(para2{i2}).(Sinks{i3})(:,Control.Data.BF_Pos-2:Control.Data.BF_Pos+2);

                    % Control
                    PreP2(1:size(Control.Data.GS_based(1).Full_RMS_AVREC,1),:,2,i3) =...
                        Control.Data.(sorting{i0})(Cond(1)).(para2{i2}).(Sinks{i3})(:,Control.Data.BF_Pos-2:Control.Data.BF_Pos+2);
                    CombiP2(1:size(Control.Data.GS_based(1).Full_RMS_AVREC,1),:,2,i3)=...
                        Control.Data.(sorting{i0})(Cond(2)).(para2{i2}).(Sinks{i3})(:,Control.Data.BF_Pos-2:Control.Data.BF_Pos+2);
                    PostP2(1:size(Control.Data.GS_based(1).Full_RMS_AVREC,1),:,2,i3)=...
                        Control.Data.(sorting{i0})(Cond(3)).(para2{i2}).(Sinks{i3})(:,Control.Data.BF_Pos-2:Control.Data.BF_Pos+2);

                    % All_Opto
                    PreP2(1:size(HighP.Data.GS_based(1).Full_RMS_AVREC,1),:,3,i3) =...
                        HighP.Data.(sorting{i0})(Cond(1)).(para2{i2}).(Sinks{i3})(:,HighP.Data.BF_Pos-2:HighP.Data.BF_Pos+2);
                    CombiP2(1:size(HighP.Data.GS_based(1).Full_RMS_AVREC,1),:,3,i3)=...
                        HighP.Data.(sorting{i0})(Cond(2)).(para2{i2}).(Sinks{i3})(:,HighP.Data.BF_Pos-2:HighP.Data.BF_Pos+2);
                    PostP2(1:size(HighP.Data.GS_based(1).Full_RMS_AVREC,1),:,3,i3)=...
                        HighP.Data.(sorting{i0})(Cond(3)).(para2{i2}).(Sinks{i3})(:,HighP.Data.BF_Pos-2:HighP.Data.BF_Pos+2);

                    % HighP
                    PreP2(1:size(HighP.Data.GS_based(1).Full_RMS_AVREC,1),:,4,i3) =...
                        HighP .Data.(sorting{i0})(Cond(1)).(para2{i2}).(Sinks{i3})(:,HighP .Data.BF_Pos-2:HighP .Data.BF_Pos+2);
                    CombiP2(1:size(HighP.Data.GS_based(1).Full_RMS_AVREC,1),:,4,i3)=...
                        HighP .Data.(sorting{i0})(Cond(2)).(para2{i2}).(Sinks{i3})(:,HighP .Data.BF_Pos-2:HighP .Data.BF_Pos+2);
                    PostP2(1:size(HighP.Data.GS_based(1).Full_RMS_AVREC,1),:,4,i3)=...
                        HighP .Data.(sorting{i0})(Cond(3)).(para2{i2}).(Sinks{i3})(:,HighP .Data.BF_Pos-2:HighP .Data.BF_Pos+2);
                 catch
                    % Control
                    PreP2(1:size(Control.Data.GS_based(1).Full_RMS_AVREC,1),:,1,i3) =...
                        Control.Data.(sorting{i0})(Cond(1)).(para2{i2})(:,Control.Data.BF_Pos-2:Control.Data.BF_Pos+2);
                    CombiP2(1:size(Control.Data.GS_based(1).Full_RMS_AVREC,1),:,1,i3)= ...
                        Control.Data.(sorting{i0})(Cond(2)).(para2{i2})(:,Control.Data.BF_Pos-2:Control.Data.BF_Pos+2);
                    PostP2(1:size(Control.Data.GS_based(1).Full_RMS_AVREC,1),:,1,i3)=...
                        Control.Data.(sorting{i0})(Cond(3)).(para2{i2})(:,Control.Data.BF_Pos-2:Control.Data.BF_Pos+2);  

                    % Control
                    PreP2(1:size(Control.Data.GS_based(1).Full_RMS_AVREC,1),:,2,i3) =...
                        Control.Data.(sorting{i0})(Cond(1)).(para2{i2})(:,Control.Data.BF_Pos-2:Control.Data.BF_Pos+2);
                    CombiP2(1:size(Control.Data.GS_based(1).Full_RMS_AVREC,1),:,2,i3)=...
                        Control.Data.(sorting{i0})(Cond(2)).(para2{i2})(:,Control.Data.BF_Pos-2:Control.Data.BF_Pos+2);
                    PostP2(1:size(Control.Data.GS_based(1).Full_RMS_AVREC,1),:,2,i3)=...
                        Control.Data.(sorting{i0})(Cond(3)).(para2{i2})(:,Control.Data.BF_Pos-2:Control.Data.BF_Pos+2);

                    % All_Opto
                    PreP2(1:size(HighP.Data.GS_based(1).Full_RMS_AVREC,1),:,3,i3) =...
                        HighP.Data.(sorting{i0})(Cond(1)).(para2{i2})(:,HighP.Data.BF_Pos-2:HighP.Data.BF_Pos+2);
                    CombiP2(1:size(HighP.Data.GS_based(1).Full_RMS_AVREC,1),:,3,i3)=...
                        HighP.Data.(sorting{i0})(Cond(2)).(para2{i2})(:,HighP.Data.BF_Pos-2:HighP.Data.BF_Pos+2);
                    PostP2(1:size(HighP.Data.GS_based(1).Full_RMS_AVREC,1),:,3,i3)=...
                        HighP.Data.(sorting{i0})(Cond(3)).(para2{i2})(:,HighP.Data.BF_Pos-2:HighP.Data.BF_Pos+2);

                    % HighP
                    PreP2(1:size(HighP.Data.GS_based(1).Full_RMS_AVREC,1),:,4,i3) =...
                        HighP .Data.(sorting{i0})(Cond(1)).(para2{i2})(:,HighP .Data.BF_Pos-2:HighP .Data.BF_Pos+2);
                    CombiP2(1:size(HighP.Data.GS_based(1).Full_RMS_AVREC,1),:,4,i3)=...
                        HighP .Data.(sorting{i0})(Cond(2)).(para2{i2})(:,HighP .Data.BF_Pos-2:HighP .Data.BF_Pos+2);
                    PostP2(1:size(HighP.Data.GS_based(1).Full_RMS_AVREC,1),:,4,i3)=...
                        HighP .Data.(sorting{i0})(Cond(3)).(para2{i2})(:,HighP .Data.BF_Pos-2:HighP .Data.BF_Pos+2);
                 end                 
            end

            end         
            clear X
%           keyboard
            for i3 = 1:length(Groups)
                for i4 = 1:length(Cond2)
                    if i4 == 1
                        X = PreP2;
                        Y = PreP1;
                    elseif i4 == 2
                        X = CombiP2;
                        Y = CombiP1;
                    else
                        X = PostP2;
                        Y = PostP1;
                    end
                    
                    for i5 = 1:size(PreP1,4)
                        if size(PreP1,4)== 1
                            CondName = Cond2{i4};
                        else
                            CondName = [Cond2{i4} '_' Sinks{i5}];
                        end
                
                h= figure('Name',['Correlation ' Groups{i3} '_' para1{i1} '_Vs_' para2{i2} ' ' sorting{i0} '_' CondName ],...
                    'PaperType','a4letter',...
                    'PaperOrientation','Portrait',...
                    'PaperUnits','centimeters', ...
                    'PaperPosition',[0.63 0.63 19.72 28.41],...    
                    'PaperSize',[20.98 29.68]); 
%                  keyboard
                try
                    clear Container
%                     BF_Pos = eval([Groups{i3} '.Data.BF_Pos']);
                        BF_Pos =3;
                        
                    subplot(5,1,1)
                    scatter(X(:,BF_Pos-2,i5),Y(:,BF_Pos-2,i3,i5))
                    test = [X(:,BF_Pos-2,i5), Y(:,BF_Pos-2,i3,i5)];
                    test = sum(~isnan(test),2);
                    
                    Xnan = X((test(:,1)== 2),BF_Pos-2,i5);
                    Ynan = Y((test(:,1)== 2),BF_Pos-2,i3,i5);
                    try
                    [rho,pval] = corr(Xnan,Ynan,'Type','Pearson');
%                     [h,p,ks2stat] =kstest2(Xnan,Ynan); 
                    
                    fitvars = polyfit(Xnan , Ynan , 1);
                    m = fitvars(1);
                    c = fitvars(2);
                    
                    hold on
                    Xnan2 =vertcat(0, Xnan);
                    Y_Fit = polyval(fitvars, Xnan2);
                    plot(Xnan2,Y_Fit);  
                    message = sprintf(['y = ' num2str(m) '*x +' num2str(c) ...
                        '\nlinear Correlation (Pearson) R = ' num2str(rho) ' p = ' num2str(pval)]);
                    text(min(Xnan),max(Y_Fit),message);
                    hold off
                    catch
                     m = nan; c = nan;    
                    end
                    title('Low non BF')
                    xlabel(para2{i2}, 'Interpreter', 'none'); ylabel([para1{i2} '_' CondName], 'Interpreter', 'none');
                    container.m(1) = m; container.c(1) = c;
                    

                    subplot(5,1,2)
                    scatter(X(:,BF_Pos-1,i5),Y(:,BF_Pos-1,i3,i5))
                    test = [X(:,BF_Pos-1,i5), Y(:,BF_Pos-1,i3,i5)];
                    test = sum(~isnan(test),2);
                    
                    Xnan = X((test(:,1)== 2),BF_Pos-1,i5);
                    Ynan = Y((test(:,1)== 2),BF_Pos-1,i3,i5);
                    try
                    [rho,pval] = corr(Xnan,Ynan,'Type','Pearson');
                    
                    fitvars = polyfit(Xnan , Ynan , 1);
                    m = fitvars(1);
                    c = fitvars(2);
                    
                    hold on
                    Xnan2 =vertcat(0, Xnan);
                    Y_Fit = polyval(fitvars, Xnan2);
                    plot(Xnan2,Y_Fit);  
                    message = sprintf(['y = ' num2str(m) '*x +' num2str(c) ...
                        '\nlinear Correlation (Pearson) R = ' num2str(rho) ' p = ' num2str(pval)]);
                    text(min(Xnan),max(Y_Fit),message);
                    hold off
                    catch
                     m = nan; c = nan;    
                    end
                    title('Low near BF')
                    xlabel(para2{i2}, 'Interpreter', 'none'); ylabel([para1{i2} '_' CondName], 'Interpreter', 'none');
                    container.m(2) = m; container.c(2) = c;
                    
                    
                    subplot(5,1,3)
                    scatter(X(:,BF_Pos,i5),Y(:,BF_Pos,i3,i5))
                    test = [X(:,BF_Pos,i5), Y(:,BF_Pos,i3,i5)];
                    test = sum(~isnan(test),2);
                    
                    Xnan = X((test(:,1)== 2),BF_Pos,i5);
                    Ynan = Y((test(:,1)== 2),BF_Pos,i3,i5);
                    try
                    [rho,pval] = corr(Xnan,Ynan,'Type','Pearson');
                    
                    fitvars = polyfit(Xnan , Ynan , 1);
                    m = fitvars(1);
                    c = fitvars(2);
                    
                    hold on
                    Xnan2 =vertcat(0, Xnan);
                    Y_Fit = polyval(fitvars, Xnan2);
                    plot(Xnan2,Y_Fit);  
                    message = sprintf(['y = ' num2str(m) '*x +' num2str(c) ...
                        '\nlinear Correlation (Pearson) R = ' num2str(rho) ' p = ' num2str(pval)]);
                    text(min(Xnan),max(Y_Fit),message);
                    hold off                    
                    catch
                     m = nan; c = nan;    
                    end
                    title('BF')
                    xlabel(para2{i2}, 'Interpreter', 'none'); ylabel([para1{i2} '_' CondName], 'Interpreter', 'none');

                    container.m(3) = m; container.c(3) = c;
                    
                    
                    subplot(5,1,4)
                    scatter(X(:,BF_Pos+1,i5),Y(:,BF_Pos+1,i3,i5))
                    test = [X(:,BF_Pos+1,i5), Y(:,BF_Pos+1,i3,i5)];
                    test = sum(~isnan(test),2);
                    
                    Xnan = X((test(:,1)== 2),BF_Pos+1,i5);
                    Ynan = Y((test(:,1)== 2),BF_Pos+1,i3,i5);
                    try
                    [rho,pval] = corr(Xnan,Ynan,'Type','Pearson');
                    
                    fitvars = polyfit(Xnan , Ynan , 1);
                    m = fitvars(1);
                    c = fitvars(2);
                    
                    
                    hold on
                    Xnan2 =vertcat(0, Xnan);
                    Y_Fit = polyval(fitvars, Xnan2);
                    plot(Xnan2,Y_Fit);  
                    message = sprintf(['y = ' num2str(m) '*x +' num2str(c) ...
                        '\nlinear Correlation (Pearson) R = ' num2str(rho) ' p = ' num2str(pval)]);
                    text(min(Xnan),max(Y_Fit),message);
                    hold off
                    catch
                     m = nan; c = nan;    
                    end
                    title('High near BF')
                    xlabel(para2{i2}, 'Interpreter', 'none'); ylabel([para1{i2} '_' CondName], 'Interpreter', 'none');
                    container.m(4) = m; container.c(4) = c;
                    
                    
                    subplot(5,1,5)
                    scatter(X(:,BF_Pos+2,i5),Y(:,BF_Pos+2,i3,i5))
                    test = [X(:,BF_Pos+2,i5), Y(:,BF_Pos+2,i3,i5)];
                    test = sum(~isnan(test),2);
                    
                    Xnan = X((test(:,1)== 2),BF_Pos+2,i5);
                    Ynan = Y((test(:,1)== 2),BF_Pos+2,i3,i5);
                    try
                    [rho,pval] = corr(Xnan,Ynan,'Type','Pearson');
                    
                    fitvars = polyfit(Xnan , Ynan , 1);
                    m = fitvars(1);
                    c = fitvars(2);
                    
                    hold on
                    Xnan2 =vertcat(0, Xnan);
                    Y_Fit = polyval(fitvars, Xnan2);
                    plot(Xnan2,Y_Fit);  
                    message = sprintf(['y = ' num2str(m) '*x +' num2str(c) ...
                        '\nlinear Correlation (Pearson) R = ' num2str(rho) ' p = ' num2str(pval)]);
                    text(min(Xnan),max(Y_Fit),message);
                    hold off
                    catch
                     m = nan; c = nan;    
                    end
                    title('High non BF')
                    xlabel(para2{i2}, 'Interpreter', 'none'); ylabel([para1{i2} '_' CondName], 'Interpreter', 'none');
                    container.m(5) = m; container.c(5) = c;
%                     keyboard
                    %% output
%                     Output.(Groups{i3}).(sorting{i0}) =[];
                    Output.(Groups{i3}).(sorting{i0}).([para1{i1} '_Vs_' para2{i2} ]).(CondName)=[];
                    Output.(Groups{i3}).(sorting{i0}).([para1{i1} '_Vs_' para2{i2} ]).(CondName).m = container.m';
                    Output.(Groups{i3}).(sorting{i0}).([para1{i1} '_Vs_' para2{i2} ]).(CondName).c = container.c';
                    cd (FigFolder)
                    
                    savefig(gcf,['Correlation ' Groups{i3} '_' para1{i1} '_Vs_' para2{i2} ' ' sorting{i0} '_' CondName])
                    saveas(gcf,['Correlation ' Groups{i3} '_' para1{i1} '_Vs_' para2{i2} ' ' sorting{i0} '_' CondName],'pdf');
                    close all

                catch
%                     keyboard
                    
                end
%                 keyboard
                  h= figure('Name',['Correlation all ' Groups{i3} '_' para1{i1} '_Vs_' para2{i2} ' ' sorting{i0} '_' CondName ],...
                    'PaperType','a4letter',...
                    'PaperOrientation','Portrait',...
                    'PaperUnits','centimeters', ...
                    'PaperPosition',[0.63 0.63 19.72 28.41],...    
                    'PaperSize',[20.98 29.68]);
                
                X1 = X(:,:,i3); % groups in 3rd dimension sinks in 4th dimension
                X1 = X1(:);
                
                Y1 = Y(:,:,:,i5);
                Y1 = Y1(:,:,i3);
                Y1 = Y1(:);
                test = [X1, Y1];
                test = sum(~isnan(test),2);
                Xnan = X1(test(:,1)== 2); Ynan = Y1(test(:,1)== 2);
                scatter(Xnan,Ynan); 
                
                    try
                    [rho,pval] = corr(Xnan,Ynan,'Type','Pearson');
                    
                    fitvars = polyfit(Xnan , Ynan , 1);
                    m = fitvars(1);
                    c = fitvars(2);
                    
                    hold on
                    Xnan2 =vertcat(0, Xnan);
                    Y_Fit = polyval(fitvars, Xnan2);
                    plot(Xnan2,Y_Fit);  
                    message = sprintf(['y = ' num2str(m) '*x +' num2str(c) ...
                        '\nlinear Correlation (Pearson) R = ' num2str(rho) ' p = ' num2str(pval)]);
                    text(min(Xnan),max(Y_Fit),message);
                    hold off
                    catch
                     m = nan; c = nan;    
                    end
                    title(['All  Correlation ' Groups{i3} '_' para1{i1} '_Vs_' para2{i2} ' ' sorting{i0} '_' CondName])
                    xlabel(para2{i2}, 'Interpreter', 'none'); ylabel([para1{i2} '_' CondName], 'Interpreter', 'none');
                    container.m(1) = m; container.c(1) = c;
%                     Output.(Groups{i3}).(sorting{i0}) =[];
                    Output.(Groups{i3}).(sorting{i0}).([para1{i1} '_Vs_' para2{i2} ]).(CondName)=[];
                    Output.(Groups{i3}).(sorting{i0}).([para1{i1} '_Vs_' para2{i2} ]).(CondName).m = container.m';
                    Output.(Groups{i3}).(sorting{i0}).([para1{i1} '_Vs_' para2{i2} ]).(CondName).c = container.c';
                    cd (FigFolder)

%                     pause(0.25)
                    savefig(gcf,['All  Correlation ' Groups{i3} '_' para1{i1} '_Vs_' para2{i2} ' ' sorting{i0} '_' CondName])
                    saveas(gcf,['All Correlation ' Groups{i3} '_' para1{i1} '_Vs_' para2{i2} ' ' sorting{i0} '_' CondName],'pdf');
                    close all                
                    Output.(Groups{i3}).(sorting{i0}).([para1{i1} '_Vs_' para2{i2} ]).(CondName).data = [Xnan,Ynan];
                    end
                end
            end

        else
    %                 h= figure('Name',[para1{i1} '_' sorting{i3} '_' CondName '_BandWidth'],...
    %                 'PaperType','a4letter',...
    %                 'PaperOrientation','Portrait',...
    %                 'PaperUnits','centimeters', ...
    %                 'PaperPosition',[0.63 0.63 19.72 28.41],...    
    %                 'PaperSize',[20.98 29.68]); 
        end
        end
    
end
end
keyboard

cd (Home)
cd DATA
mkdir ('Correlation')
cd ('Correlation')
save(['Correlation_analysis_' date],'Output');
