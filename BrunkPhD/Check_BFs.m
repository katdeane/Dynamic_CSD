%% DOCUMENT YOUR CODE PLS THANKS
%%

if ~exist('C','var')
C = load('K:\CSD_dynamic_analysis\DATA\Input_Control_7post_n7_Data.mat');
O = load('K:\CSD_dynamic_analysis\DATA\Input_HighP_7post_Data.mat');
Y = load('K:\CSD_dynamic_analysis\DATA\Input_YFP_7post_Data.mat');
end
C_Names = fieldnames(C.Data);
O_Names = fieldnames(O.Data);
Y_Names = fieldnames(Y.Data);

BF_C = nan(length(C_Names),3,2);
BF_O = nan(length(O_Names),3,2);
BF_Y = nan(length(Y_Names),3,2);

for i1 = 1:length(C_Names)
    
    for i2 = 1:3
       BF_C(i1,i2,1) = C.Data(i2).(C_Names{i1}).GS_BF; 
       vec(1,i2) = C.Data(i2).(C_Names{i1}).GS_BF; 
    end

    val = unique(vec);
    cnt = histc(vec,val);
    BF = find(cnt == max(cnt));
    if length (BF) == 1
    BF_C(i1,1,2) = val(BF);
    else
    BF_C(i1,1,2) = val(2);
    end

end

for i1 = 1:length(O_Names)
    
    for i2 = 1:3
       BF_O(i1,i2) = O.Data(i2).(O_Names{i1}).GS_BF; 
       vec(1,i2) = O.Data(i2).(O_Names{i1}).GS_BF; 
    end
    
    val = unique(vec);
    cnt = histc(vec,val);
    BF = find(cnt == max(cnt));
    if length (BF) == 1
    BF_O(i1,1,2) = val(BF);
    else
    BF_O(i1,1,2) = val(2);
    end
    
end

for i1 = 1:length(Y_Names)
    
    for i2 = 1:3
       BF_Y(i1,i2) = Y.Data(i2).(Y_Names{i1}).GS_BF; 
       vec(1,i2) = Y.Data(i2).(Y_Names{i1}).GS_BF; 
    end
    
    val = unique(vec);
    cnt = histc(vec,val);
    BF = find(cnt == max(cnt));
    if length (BF) == 1
    BF_Y(i1,1,2) = val(BF);
    else
    BF_Y(i1,1,2) = val(2);
    end
    
end

frqz = [125, 250, 500, 1000, 2000, 4000, 8000, 16000, 32000];
Idx = 1:length(frqz);


vec_1 = BF_C(:,1,2); 
val_1 = unique(vec_1);
cnt_1 = histc(vec_1,val_1);

vec_2 = BF_O(:,1,2); 
val_2 = unique(vec_2);
cnt_2 = histc(vec_2,val_2);

vec_3 = BF_Y(:,1,2); 
val_3 = unique(vec_3);
cnt_3 = histc(vec_3,val_3);

Hit_1 = frqz == val_1;
Hit_2 = frqz == val_2;
Hit_3 = frqz == val_3;

Hit_1 = sum(Hit_1);
Hit_2 = sum(Hit_2);
Hit_3 = sum(Hit_3);

figure
hold on
bar(Idx(logical(Hit_1))-0.25,cnt_1',0.25)
bar(Idx(logical(Hit_2)),cnt_2',0.25)
bar(Idx(logical(Hit_3))+0.25,cnt_3',0.25)

yticks(0:1:7)
xlim ([1-0.5 length(frqz)+0.5])
xticks = Idx;
xticklabels(strsplit(num2str(frqz)));
set(gca,'FontSize',16)

xlabel('Initial BFs in Hz','Fontsize',20,'FontWeight','bold')
ylabel('Animal count','Fontsize',20,'FontWeight','bold')

lgd = legend('Control','C1V1', 'YFP');
lgd.FontSize = 16;

mean (vec_1)/1000,std (vec_1)/sqrt(7)/1000
mean (vec_2)/1000,std (vec_2)/sqrt(12)/1000
mean (vec_3)/1000,std (vec_3)/sqrt(7)/1000


x =[vertcat(vec_1, nan(5,1)),vec_2,vertcat(vec_3, nan(5,1))]
[p,tbl,stats] = kruskalwallis(x)
[results,means] = multcompare(stats,'CType','bonferroni')





