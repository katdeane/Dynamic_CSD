% Max Paper ideas:
dbstop if error
clear  
close all
% addpath ('\\linsamba\BFdata\CortXplorer\_MATLAB CODE\Routines')
addpath('K:\CSD_dynamic_analysis\subfunc')
%%% Fig 2 Grand AVREC/ RELRES 
load('K:\CSD_dynamic_analysis\DATA\Input_HighP_7post_Data.mat')

% get names
ID = fieldnames(Data);
M = [1, 4, 9]
%Fig2  Pre 3 Full AVREC/RelRes GS BF / Laser 4
AVREC_1 = nan(800,12); RelRes_1 = nan(800,12);
AVREC_2 = nan(800,12); RelRes_2 = nan(800,12);
AVREC_3 = nan(800,12); RelRes_3 = nan(800,12);
for i1 = 1:length(ID)
    BF = find(Data(3).(ID{i1}).Frqz == Data(3).(ID{i1}).GS_BF);
    AVREC_1(1:length(Data(M(1)).(ID{i1}).AVREC_raw {BF}),i1)=Data(M(1)).(ID{i1}).AVREC_raw {BF};
    RelRes_1 (1:length(Data(M(1)).(ID{i1}).RELRES_raw {BF}),i1)=Data(M(1)).(ID{i1}).RELRES_raw {BF};
    AVREC_2(1:length(Data(M(2)).(ID{i1}).AVREC_raw {BF}),i1)=Data(M(2)).(ID{i1}).AVREC_raw {BF};
    RelRes_2 (1:length(Data(M(2)).(ID{i1}).RELRES_raw {BF}),i1)=Data(M(2)).(ID{i1}).RELRES_raw {BF};
    AVREC_3(1:length(Data(M(3)).(ID{i1}).AVREC_raw {BF}),i1)=Data(M(3)).(ID{i1}).AVREC_raw {BF};
    RelRes_3 (1:length(Data(M(3)).(ID{i1}).RELRES_raw {BF}),i1)=Data(M(3)).(ID{i1}).RELRES_raw {BF};
end

figure
subplot(2,1,1)

hold on
 area([200 250],[0.002 0.002],'FaceColor','y','FaceAlpha',0.05)
area([280 500],[0.002 0.002],'FaceColor','r','FaceAlpha',0.05)
plot(nanmean(AVREC_1'),'LineWidth',2)
plot(nanmean(AVREC_2'),'LineWidth',2)
plot(nanmean(AVREC_3'),'LineWidth',2)
plot([200 200],[0 0.00250],'--w')
plot([400 400],[0 0.00250],'--w')
xlim([100 500])
xticks([100 150 200 250 300 350 400 450 500])
xticklabels({'-100','-50', '0','50', '100','150', '200','250','300'})
set(gca,'fontsize',12,'FontWeight','bold')
xlabel('Time [ms]','FontSize',16,'FontWeight','bold')
ylabel('AVREC [mV/mm²]','FontSize',16,'FontWeight','bold')
title('AVREC High Performer n = 12','FontSize',20,'FontWeight','bold')
legend('Early','Late','Pre','Opto','Post','Location','best')

subplot(2,1,2)
hold on
area([200  200 250 250],[-30 20 20 -30],'FaceColor','y','FaceAlpha',0.05)
area([280 280 500 500],[-30 20 20 -30],'FaceColor','r','FaceAlpha',0.05)
plot(nanmean(RelRes_1'*100),'LineWidth',2)

plot(nanmean(RelRes_2'*100),'LineWidth',2)
plot(nanmean(RelRes_3'*100),'LineWidth',2)
plot([200 200],[-30 15],'--w')
plot([400 400],[-30 15],'--w')
xlim([100 500])
xticks([100 150 200 250 300 350 400 450 500])
xticklabels({'-100','-50', '0','50', '100','150', '200','250','300'})
set(gca,'fontsize',12,'FontWeight','bold')
xlabel('Time [ms]','FontSize',16,'FontWeight','bold')
ylabel('RelResCSD [%]','FontSize',16,'FontWeight','bold')
title('RelResCSD High Performer n = 12','FontSize',20,'FontWeight','bold')
legend('Early','Late','Pre','Opto','Post','Location','best')

% keyboard
AVREC_1 = nan(800,12,3); RelRes_1 = nan(800,12,3);
AVREC_2 = nan(800,12); RelRes_2 = nan(800,12);
AVREC_3 = nan(800,12,3); RelRes_3 = nan(800,12,3);
AVREC_4 = nan(800,12,2); RelRes_4 = nan(800,12,2);
AVREC_5 = nan(800,12,2); RelRes_5 = nan(800,12,2);
for i1 = 1:length(ID)
%     BF = find(Data(3).(ID{i1}).Frqz == Data(3).(ID{i1}).GS_BF);
    BF =find(Data(1).(ID{i1}).Early_RMS_AVREC == max(Data(1).(ID{i1}).Early_RMS_AVREC));
    AVREC_1(1:length(Data(1).(ID{i1}).AVREC_raw {BF}),i1,1)=Data(M(1)).(ID{i1}).AVREC_raw {BF};
    RelRes_1 (1:length(Data(1).(ID{i1}).RELRES_raw {BF}),i1,1)=Data(M(1)).(ID{i1}).RELRES_raw {BF};
    
    BF =find(Data(2).(ID{i1}).Early_RMS_AVREC == max(Data(2).(ID{i1}).Early_RMS_AVREC));
    AVREC_1(1:length(Data(2).(ID{i1}).AVREC_raw {BF}),i1,2)=Data(2).(ID{i1}).AVREC_raw {BF};
    RelRes_1 (1:length(Data(2).(ID{i1}).RELRES_raw {BF}),i1,2)=Data(2).(ID{i1}).RELRES_raw {BF};
    
    BF =find(Data(3).(ID{i1}).Early_RMS_AVREC == max(Data(3).(ID{i1}).Early_RMS_AVREC));
    AVREC_1(1:length(Data(3).(ID{i1}).AVREC_raw {BF}),i1,3)=Data(3).(ID{i1}).AVREC_raw {BF};
    RelRes_1 (1:length(Data(3).(ID{i1}).RELRES_raw {BF}),i1,3)=Data(3).(ID{i1}).RELRES_raw {BF};
    
    BF =find(Data(4).(ID{i1}).Early_RMS_AVREC == max(Data(4).(ID{i1}).Early_RMS_AVREC));
    AVREC_2(1:length(Data(4).(ID{i1}).AVREC_raw {BF}),i1)=Data(4).(ID{i1}).AVREC_raw {BF};
    RelRes_2 (1:length(Data(4).(ID{i1}).RELRES_raw {BF}),i1)=Data(4).(ID{i1}).RELRES_raw {BF};
    
    BF =find(Data(5).(ID{i1}).Early_RMS_AVREC == max(Data(5).(ID{i1}).Early_RMS_AVREC));
    AVREC_3(1:length(Data(5).(ID{i1}).AVREC_raw {BF}),i1,1)=Data(5).(ID{i1}).AVREC_raw {BF};
    RelRes_3 (1:length(Data(5).(ID{i1}).RELRES_raw {BF}),i1,1)=Data(5).(ID{i1}).RELRES_raw {BF};
    
    BF =find(Data(6).(ID{i1}).Early_RMS_AVREC == max(Data(6).(ID{i1}).Early_RMS_AVREC));
    AVREC_3(1:length(Data(6).(ID{i1}).AVREC_raw {BF}),i1,2)=Data(6).(ID{i1}).AVREC_raw {BF};
    RelRes_3 (1:length(Data(6).(ID{i1}).RELRES_raw {BF}),i1,2)=Data(6).(ID{i1}).RELRES_raw {BF};
    
    BF =find(Data(7).(ID{i1}).Early_RMS_AVREC == max(Data(7).(ID{i1}).Early_RMS_AVREC));
    AVREC_3(1:length(Data(7).(ID{i1}).AVREC_raw {BF}),i1,3)=Data(7).(ID{i1}).AVREC_raw {BF};
    RelRes_3 (1:length(Data(7).(ID{i1}).RELRES_raw {BF}),i1,3)=Data(7).(ID{i1}).RELRES_raw {BF};
    
    BF =find(Data(8).(ID{i1}).Early_RMS_AVREC == max(Data(8).(ID{i1}).Early_RMS_AVREC));
    AVREC_4(1:length(Data(8).(ID{i1}).AVREC_raw {BF}),i1,1)=Data(8).(ID{i1}).AVREC_raw {BF};
    RelRes_4 (1:length(Data(8).(ID{i1}).RELRES_raw {BF}),i1,1)=Data(8).(ID{i1}).RELRES_raw {BF};
    
    BF =find(Data(9).(ID{i1}).Early_RMS_AVREC == max(Data(9).(ID{i1}).Early_RMS_AVREC));
    AVREC_5(1:length(Data(9).(ID{i1}).AVREC_raw {BF}),i1,2)=Data(9).(ID{i1}).AVREC_raw {BF};
    RelRes_5 (1:length(Data(9).(ID{i1}).RELRES_raw {BF}),i1,2)=Data(9).(ID{i1}).RELRES_raw {BF};
    try
        BF =find(Data(10).(ID{i1}).Early_RMS_AVREC == max(Data(10).(ID{i1}).Early_RMS_AVREC));
    AVREC_5(1:length(Data(10).(ID{i1}).AVREC_raw {BF}),i1,1)=Data(10).(ID{i1}).AVREC_raw {BF};
    RelRes_5 (1:length(Data(10).(ID{i1}).RELRES_raw {BF}),i1,1)=Data(10).(ID{i1}).RELRES_raw {BF};
    catch,end
    try
        BF =find(Data(11).(ID{i1}).Early_RMS_AVREC == max(Data(11).(ID{i1}).Early_RMS_AVREC));
    AVREC_5(1:length(Data(11).(ID{i1}).AVREC_raw {BF}),i1,2)=Data(11).(ID{i1}).AVREC_raw {BF};
    RelRes_5 (1:length(Data(11).(ID{i1}).RELRES_raw {BF}),i1,2)=Data(11).(ID{i1}).RELRES_raw_raw {BF};
    catch,end
end
AVREC_1 = nanmean(AVREC_1,3);
AVREC_3 = nanmean(AVREC_3,3);
AVREC_4 = nanmean(AVREC_4,3);
 AVREC_5 = nanmean(AVREC_5,3);
RelRes_1 = nanmean(RelRes_1,3);
RelRes_3 = nanmean(RelRes_3,3);
RelRes_4 = nanmean(RelRes_4,3);
 RelRes_5 = nanmean(RelRes_5,3);

figure
subplot(2,1,1)
hold on
area([200 250 ],[0.0025 0.0025 ],'FaceColor','y','FaceAlpha',0.05)
area([280 500 ],[0.0025 0.0025 ],'FaceColor','r','FaceAlpha',0.05)
plot(nanmean(AVREC_1'),'LineWidth',2)
plot(nanmean(AVREC_2'),'LineWidth',2)
plot(nanmean(AVREC_3'),'LineWidth',2)
plot(nanmean(AVREC_4'),'LineWidth',2)
 plot(nanmean(AVREC_5'),'LineWidth',2)
plot([200 200],[0 0.00250],'--w','LineWidth',2)
plot([400 400],[0 0.00250],'--w','LineWidth',2)
xlim([100 500])
xticks([100 150 200 250 300 350 400 450 500])
xticklabels({'-100','-50', '0','50', '100','150', '200','250','300'})
set(gca,'fontsize',12,'FontWeight','bold')
xlabel('Time [ms]','FontSize',16,'FontWeight','bold')
ylabel('AVREC [mV/mm²]','FontSize',16,'FontWeight','bold')
title('AVREC High Performer n = 12','FontSize',20,'FontWeight','bold')
legend('Early','Late','Pre 1-3','Opto','Post 1-3','Post 4-5','Post 6-7','Location','best')

subplot(2,1,2)
hold on
area([200  200 250 250],[-30 15 15 -30],'FaceColor','y','FaceAlpha',0.05)
area([280 280 500 500],[-30 15 15 -30],'FaceColor','r','FaceAlpha',0.05)
plot(nanmean(RelRes_1'*100),'LineWidth',2)
plot(nanmean(RelRes_2'*100),'LineWidth',2)
plot(nanmean(RelRes_3'*100),'LineWidth',2)
plot(nanmean(RelRes_4'*100),'LineWidth',2)
plot(nanmean(RelRes_5'*100),'LineWidth',2)
plot([200 200],[-30 15],'--w','LineWidth',2)
plot([400 400],[-30 15],'--w','LineWidth',2)
xlim([100 500])
xticks([100 150 200 250 300 350 400 450 500])
xticklabels({'-100','-50', '0','50', '100','150', '200','250','300'})
set(gca,'fontsize',12,'FontWeight','bold')
xlabel('Time [ms]','FontSize',16,'FontWeight','bold')
ylabel('RelResCSD [%]','FontSize',16,'FontWeight','bold')
title('RELRESCSD High Performer n = 12','FontSize',20,'FontWeight','bold')
legend('Early','Late','Pre 1-3','Opto','Post 1-3','Post 4-5','Post 6-7','Location','best')


load('K:\CSD_dynamic_analysis\DATA\Input_YFP_7post_Data.mat')

% get names
ID = fieldnames(Data);
M = [1, 4, 9]
%Fig2  Pre 3 Full AVREC/RelRes GS BF / Laser 4
AVREC_1 = nan(800,12); RelRes_1 = nan(800,12);
AVREC_2 = nan(800,12); RelRes_2 = nan(800,12);
AVREC_3 = nan(800,12); RelRes_3 = nan(800,12);
for i1 = 1:length(ID)
    BF = find(Data(3).(ID{i1}).Frqz == Data(3).(ID{i1}).GS_BF);
    AVREC_1(1:length(Data(M(1)).(ID{i1}).AVREC_raw {BF}),i1)=Data(M(1)).(ID{i1}).AVREC_raw {BF};
    RelRes_1 (1:length(Data(M(1)).(ID{i1}).RELRES_raw {BF}),i1)=Data(M(1)).(ID{i1}).RELRES_raw {BF};
    AVREC_2(1:length(Data(M(2)).(ID{i1}).AVREC_raw {BF}),i1)=Data(M(2)).(ID{i1}).AVREC_raw {BF};
    RelRes_2 (1:length(Data(M(2)).(ID{i1}).RELRES_raw {BF}),i1)=Data(M(2)).(ID{i1}).RELRES_raw {BF};
    AVREC_3(1:length(Data(M(3)).(ID{i1}).AVREC_raw {BF}),i1)=Data(M(3)).(ID{i1}).AVREC_raw {BF};
    RelRes_3 (1:length(Data(M(3)).(ID{i1}).RELRES_raw {BF}),i1)=Data(M(3)).(ID{i1}).RELRES_raw {BF};
end
% keyboard
figure
subplot(2,1,1)
hold on
area([200 250 ],[0.0025 0.0025 ],'FaceColor','y','FaceAlpha',0.05)
area([280 500 ],[0.0025 0.0025 ],'FaceColor','r','FaceAlpha',0.05)
plot(nanmean(AVREC_1'),'LineWidth',2)
plot(nanmean(AVREC_2'),'LineWidth',2)
plot(nanmean(AVREC_3'),'LineWidth',2)
plot([200 200],[0 0.00250],'--w','LineWidth',2)
plot([400 400],[0 0.00250],'--w','LineWidth',2)
xlim([100 500])
xticks([100 150 200 250 300 350 400 450 500])
xticklabels({'-100','-50', '0','50', '100','150', '200','250','300'})
set(gca,'fontsize',12,'FontWeight','bold')
xlabel('Time [ms]','FontSize',16,'FontWeight','bold')
ylabel('AVREC [mV/mm²]','FontSize',16,'FontWeight','bold')
title('AVREC YFP n = 7','FontSize',20,'FontWeight','bold')
legend('Early','Late','Pre','Opto','Post','Location','best')

subplot(2,1,2)
hold on
area([200  200 250 250],[-30 15 15 -30],'FaceColor','y','FaceAlpha',0.05)
area([280 280 500 500],[-30 15 15 -30],'FaceColor','r','FaceAlpha',0.05)
plot(nanmean(RelRes_1'*100),'LineWidth',2)
plot(nanmean(RelRes_2'*100),'LineWidth',2)
plot(nanmean(RelRes_3'*100),'LineWidth',2)
plot([200 200],[-30 15],'--w','LineWidth',2)
plot([400 400],[-30 15],'--w','LineWidth',2)
xlim([100 500])
xticks([100 150 200 250 300 350 400 450 500])
xticklabels({'-100','-50', '0','50', '100','150', '200','250','300'})
set(gca,'fontsize',12,'FontWeight','bold')
xlabel('Time [ms]','FontSize',16,'FontWeight','bold')
ylabel('RelResCSD [%]','FontSize',16,'FontWeight','bold')
title('RelResCSD YFP n = 4','FontSize',20,'FontWeight','bold')
legend('Early','Late','Pre','Opto','Post','Location','best')
% keyboard

AVREC_1 = nan(800,12,3); RelRes_1 = nan(800,12,3);
AVREC_2 = nan(800,12); RelRes_2 = nan(800,12);
AVREC_3 = nan(800,12,3); RelRes_3 = nan(800,12,3);
AVREC_4 = nan(800,12,2); RelRes_4 = nan(800,12,2);
AVREC_5 = nan(800,12,2); RelRes_5 = nan(800,12,2);
for i1 = 1:length(ID)
%     BF = find(Data(3).(ID{i1}).Frqz == Data(3).(ID{i1}).GS_BF);
    BF =find(Data(1).(ID{i1}).Early_RMS_AVREC == max(Data(1).(ID{i1}).Early_RMS_AVREC));
    AVREC_1(1:length(Data(1).(ID{i1}).AVREC_raw {BF}),i1,1)=Data(M(1)).(ID{i1}).AVREC_raw {BF};
    RelRes_1 (1:length(Data(1).(ID{i1}).RELRES_raw {BF}),i1,1)=Data(M(1)).(ID{i1}).RELRES_raw {BF};
    
    BF =find(Data(2).(ID{i1}).Early_RMS_AVREC == max(Data(2).(ID{i1}).Early_RMS_AVREC));
    AVREC_1(1:length(Data(2).(ID{i1}).AVREC_raw {BF}),i1,2)=Data(2).(ID{i1}).AVREC_raw {BF};
    RelRes_1 (1:length(Data(2).(ID{i1}).RELRES_raw {BF}),i1,2)=Data(2).(ID{i1}).RELRES_raw {BF};
    
    BF =find(Data(3).(ID{i1}).Early_RMS_AVREC == max(Data(3).(ID{i1}).Early_RMS_AVREC));
    AVREC_1(1:length(Data(3).(ID{i1}).AVREC_raw {BF}),i1,3)=Data(3).(ID{i1}).AVREC_raw {BF};
    RelRes_1 (1:length(Data(3).(ID{i1}).RELRES_raw {BF}),i1,3)=Data(3).(ID{i1}).RELRES_raw {BF};
    
    BF =find(Data(4).(ID{i1}).Early_RMS_AVREC == max(Data(4).(ID{i1}).Early_RMS_AVREC));
    AVREC_2(1:length(Data(4).(ID{i1}).AVREC_raw {BF}),i1)=Data(4).(ID{i1}).AVREC_raw {BF};
    RelRes_2 (1:length(Data(4).(ID{i1}).RELRES_raw {BF}),i1)=Data(4).(ID{i1}).RELRES_raw {BF};
    
    BF =find(Data(5).(ID{i1}).Early_RMS_AVREC == max(Data(5).(ID{i1}).Early_RMS_AVREC));
    AVREC_3(1:length(Data(5).(ID{i1}).AVREC_raw {BF}),i1,1)=Data(5).(ID{i1}).AVREC_raw {BF};
    RelRes_3 (1:length(Data(5).(ID{i1}).RELRES_raw {BF}),i1,1)=Data(5).(ID{i1}).RELRES_raw {BF};
    
    BF =find(Data(6).(ID{i1}).Early_RMS_AVREC == max(Data(6).(ID{i1}).Early_RMS_AVREC));
    AVREC_3(1:length(Data(6).(ID{i1}).AVREC_raw {BF}),i1,2)=Data(6).(ID{i1}).AVREC_raw {BF};
    RelRes_3 (1:length(Data(6).(ID{i1}).RELRES_raw {BF}),i1,2)=Data(6).(ID{i1}).RELRES_raw {BF};
    
    BF =find(Data(7).(ID{i1}).Early_RMS_AVREC == max(Data(7).(ID{i1}).Early_RMS_AVREC));
    AVREC_3(1:length(Data(7).(ID{i1}).AVREC_raw {BF}),i1,3)=Data(7).(ID{i1}).AVREC_raw {BF};
    RelRes_3 (1:length(Data(7).(ID{i1}).RELRES_raw {BF}),i1,3)=Data(7).(ID{i1}).RELRES_raw {BF};
    
    BF =find(Data(8).(ID{i1}).Early_RMS_AVREC == max(Data(8).(ID{i1}).Early_RMS_AVREC));
    AVREC_4(1:length(Data(8).(ID{i1}).AVREC_raw {BF}),i1,1)=Data(8).(ID{i1}).AVREC_raw {BF};
    RelRes_4 (1:length(Data(8).(ID{i1}).RELRES_raw {BF}),i1,1)=Data(8).(ID{i1}).RELRES_raw {BF};
    
    BF =find(Data(9).(ID{i1}).Early_RMS_AVREC == max(Data(9).(ID{i1}).Early_RMS_AVREC));
    AVREC_5(1:length(Data(9).(ID{i1}).AVREC_raw {BF}),i1,2)=Data(9).(ID{i1}).AVREC_raw {BF};
    RelRes_5 (1:length(Data(9).(ID{i1}).RELRES_raw {BF}),i1,2)=Data(9).(ID{i1}).RELRES_raw {BF};
    try
        BF =find(Data(10).(ID{i1}).Early_RMS_AVREC == max(Data(10).(ID{i1}).Early_RMS_AVREC));
    AVREC_5(1:length(Data(10).(ID{i1}).AVREC_raw {BF}),i1,1)=Data(10).(ID{i1}).AVREC_raw {BF};
    RelRes_5 (1:length(Data(10).(ID{i1}).RELRES_raw {BF}),i1,1)=Data(10).(ID{i1}).RELRES_raw {BF};
    catch,end
    try
        BF =find(Data(11).(ID{i1}).Early_RMS_AVREC == max(Data(11).(ID{i1}).Early_RMS_AVREC));
    AVREC_5(1:length(Data(11).(ID{i1}).AVREC_raw {BF}),i1,2)=Data(11).(ID{i1}).AVREC_raw {BF};
    RelRes_5 (1:length(Data(11).(ID{i1}).RELRES_raw {BF}),i1,2)=Data(11).(ID{i1}).RELRES_raw {BF};
    catch,end
end
AVREC_1 = nanmean(AVREC_1,3);
AVREC_3 = nanmean(AVREC_3,3);
AVREC_4 = nanmean(AVREC_4,3);
AVREC_5 = nanmean(AVREC_5,3);
RelRes_1 = nanmean(RelRes_1,3);
RelRes_3 = nanmean(RelRes_3,3);
RelRes_4 = nanmean(RelRes_4,3);
RelRes_5 = nanmean(RelRes_5,3);

figure
subplot(2,1,1)
hold on
area([200 250 ],[0.0025 0.0025 ],'FaceColor','y','FaceAlpha',0.05)
area([280 500 ],[0.0025 0.0025 ],'FaceColor','r','FaceAlpha',0.05)
plot(nanmean(AVREC_1'),'LineWidth',2)
plot(nanmean(AVREC_2'),'LineWidth',2)
plot(nanmean(AVREC_3'),'LineWidth',2)
plot(nanmean(AVREC_4'),'LineWidth',2)
plot(nanmean(AVREC_5'),'LineWidth',2)
plot([200 200],[0 0.00250],'--w','LineWidth',2)
plot([400 400],[0 0.00250],'--w','LineWidth',2)
xlim([100 500])
xticks([100 150 200 250 300 350 400 450 500])
xticklabels({'-100','-50', '0','50', '100','150', '200','250','300'})
set(gca,'fontsize',12,'FontWeight','bold')
xlabel('Time [ms]','FontSize',16,'FontWeight','bold')
ylabel('AVREC [mV/mm²]','FontSize',16,'FontWeight','bold')
title('AVREC YFP n = 7','FontSize',20,'FontWeight','bold')
legend('Early','Late','Pre 1-3','Opto','Post 1-3','Post 4-5','Post 6-7','Location','best')

subplot(2,1,2)
hold on
area([200  200 250 250],[-30 15 15 -30],'FaceColor','y','FaceAlpha',0.05)
area([280 280 500 500],[-30 15 15 -30],'FaceColor','r','FaceAlpha',0.05)
plot(nanmean(RelRes_1'*100),'LineWidth',2)
plot(nanmean(RelRes_2'*100),'LineWidth',2)
plot(nanmean(RelRes_3'*100),'LineWidth',2)
plot(nanmean(RelRes_4'*100),'LineWidth',2)
plot(nanmean(RelRes_5'*100),'LineWidth',2)
plot([200 200],[-30 15],'--w','LineWidth',2)
plot([400 400],[-30 15],'--w','LineWidth',2)
xlim([100 500])
xticks([100 150 200 250 300 350 400 450 500])
xticklabels({'-100','-50', '0','50', '100','150', '200','250','300'})
set(gca,'fontsize',12,'FontWeight','bold')
xlabel('Time [ms]','FontSize',16,'FontWeight','bold')
ylabel('RelResCSD [%]','FontSize',16,'FontWeight','bold')
title('RelResCSD YFP n = 4','FontSize',20,'FontWeight','bold')
legend('Early','Late','Pre 1-3','Opto','Post 1-3','Post 4-5','Post 6-7','Location','best')




load('K:\CSD_dynamic_analysis\DATA\Input_Control_7post_n7_Data.mat')

% get names
ID = fieldnames(Data);
M = [1, 4, 9]
%Fig2  Pre 3 Full AVREC/RelRes GS BF / Laser 4
AVREC_1 = nan(800,12); RelRes_1 = nan(800,12);
AVREC_2 = nan(800,12); RelRes_2 = nan(800,12);
AVREC_3 = nan(800,12); RelRes_3 = nan(800,12);
for i1 = 1:length(ID)
    BF = find(Data(3).(ID{i1}).Frqz == Data(3).(ID{i1}).GS_BF);
    AVREC_1(1:length(Data(M(1)).(ID{i1}).AVREC_raw {BF}),i1)=Data(M(1)).(ID{i1}).AVREC_raw {BF};
    RelRes_1 (1:length(Data(M(1)).(ID{i1}).RELRES_raw {BF}),i1)=Data(M(1)).(ID{i1}).RELRES_raw {BF};
    AVREC_2(1:length(Data(M(2)).(ID{i1}).AVREC_raw {BF}),i1)=Data(M(2)).(ID{i1}).AVREC_raw {BF};
    RelRes_2 (1:length(Data(M(2)).(ID{i1}).RELRES_raw {BF}),i1)=Data(M(2)).(ID{i1}).RELRES_raw {BF};
    AVREC_3(1:length(Data(M(3)).(ID{i1}).AVREC_raw {BF}),i1)=Data(M(3)).(ID{i1}).AVREC_raw {BF};
    RelRes_3 (1:length(Data(M(3)).(ID{i1}).RELRES_raw {BF}),i1)=Data(M(3)).(ID{i1}).RELRES_raw {BF};
end

figure
subplot(2,1,1)
hold on
area([200 250 ],[0.0025 0.0025 ],'FaceColor','y','FaceAlpha',0.05)
area([280 500 ],[0.0025 0.0025 ],'FaceColor','r','FaceAlpha',0.05)
plot(nanmean(AVREC_1'),'LineWidth',2)
plot(nanmean(AVREC_2'),'LineWidth',2)
plot(nanmean(AVREC_3'),'LineWidth',2)
plot([200 200],[0 0.00250],'--w')
plot([400 400],[0 0.0025],'--w')
xlim([100 500])
xticks([100 150 200 250 300 350 400 450 500])
xticklabels({'-100','-50', '0','50', '100','150', '200','250','300'})
set(gca,'fontsize',12,'FontWeight','bold')
xlabel('Time [ms]','FontSize',16,'FontWeight','bold')
ylabel('AVREC [mV/mm²]','FontSize',16,'FontWeight','bold')
title('AVREC Control n = 7','FontSize',20,'FontWeight','bold')
legend('Early','Late','Pre','Opto','Post','Location','best')

subplot(2,1,2)
hold on
area([200  200 250 250],[-30 15 15 -30],'FaceColor','y','FaceAlpha',0.05)
area([280 280 500 500],[-30 15 15 -30],'FaceColor','r','FaceAlpha',0.05)
plot(nanmean(RelRes_1'*100),'LineWidth',2)
plot(nanmean(RelRes_2'*100),'LineWidth',2)
plot(nanmean(RelRes_3'*100),'LineWidth',2)
plot([200 200],[-30 15],'--w')
plot([400 400],[-30 15],'--w')
xlim([100 500])
xticks([100 150 200 250 300 350 400 450 500])
xticklabels({'-100','-50', '0','50', '100','150', '200','250','300'})
set(gca,'fontsize',12,'FontWeight','bold')
xlabel('Time [ms]','FontSize',16,'FontWeight','bold')
ylabel('RelResCSD [%]','FontSize',16,'FontWeight','bold')
title('RelResCSD Control n = 7','FontSize',20,'FontWeight','bold')
legend('Early','Late','Pre','Opto','Post','Location','best')


AVREC_1 = nan(800,12,3); RelRes_1 = nan(800,12,3);
AVREC_2 = nan(800,12); RelRes_2 = nan(800,12);
AVREC_3 = nan(800,12,3); RelRes_3 = nan(800,12,3);
AVREC_4 = nan(800,12,2); RelRes_4 = nan(800,12,2);
 AVREC_5 = nan(800,12,2); RelRes_5 = nan(800,12,2);
for i1 = 1:length(ID)
%     BF = find(Data(3).(ID{i1}).Frqz == Data(3).(ID{i1}).GS_BF);
    BF =find(Data(1).(ID{i1}).Early_RMS_AVREC == max(Data(1).(ID{i1}).Early_RMS_AVREC));
    AVREC_1(1:length(Data(1).(ID{i1}).AVREC_raw {BF}),i1,1)=Data(M(1)).(ID{i1}).AVREC_raw {BF};
    RelRes_1 (1:length(Data(1).(ID{i1}).RELRES_raw {BF}),i1,1)=Data(M(1)).(ID{i1}).RELRES_raw {BF};
    
    BF =find(Data(2).(ID{i1}).Early_RMS_AVREC == max(Data(2).(ID{i1}).Early_RMS_AVREC));
    AVREC_1(1:length(Data(2).(ID{i1}).AVREC_raw {BF}),i1,2)=Data(2).(ID{i1}).AVREC_raw {BF};
    RelRes_1 (1:length(Data(2).(ID{i1}).RELRES_raw {BF}),i1,2)=Data(2).(ID{i1}).RELRES_raw {BF};
    
    BF =find(Data(3).(ID{i1}).Early_RMS_AVREC == max(Data(3).(ID{i1}).Early_RMS_AVREC));
    AVREC_1(1:length(Data(3).(ID{i1}).AVREC_raw {BF}),i1,3)=Data(3).(ID{i1}).AVREC_raw {BF};
    RelRes_1 (1:length(Data(3).(ID{i1}).RELRES_raw {BF}),i1,3)=Data(3).(ID{i1}).RELRES_raw {BF};
    
    BF =find(Data(4).(ID{i1}).Early_RMS_AVREC == max(Data(4).(ID{i1}).Early_RMS_AVREC));
    AVREC_2(1:length(Data(4).(ID{i1}).AVREC_raw {BF}),i1)=Data(4).(ID{i1}).AVREC_raw {BF};
    RelRes_2 (1:length(Data(4).(ID{i1}).RELRES_raw {BF}),i1)=Data(4).(ID{i1}).RELRES_raw {BF};
    
    BF =find(Data(5).(ID{i1}).Early_RMS_AVREC == max(Data(5).(ID{i1}).Early_RMS_AVREC));
    AVREC_3(1:length(Data(5).(ID{i1}).AVREC_raw {BF}),i1,1)=Data(5).(ID{i1}).AVREC_raw {BF};
    RelRes_3 (1:length(Data(5).(ID{i1}).RELRES_raw {BF}),i1,1)=Data(5).(ID{i1}).RELRES_raw {BF};
    
    BF =find(Data(6).(ID{i1}).Early_RMS_AVREC == max(Data(6).(ID{i1}).Early_RMS_AVREC));
    AVREC_3(1:length(Data(6).(ID{i1}).AVREC_raw {BF}),i1,2)=Data(6).(ID{i1}).AVREC_raw {BF};
    RelRes_3 (1:length(Data(6).(ID{i1}).RELRES_raw {BF}),i1,2)=Data(6).(ID{i1}).RELRES_raw {BF};
    
    BF =find(Data(7).(ID{i1}).Early_RMS_AVREC == max(Data(7).(ID{i1}).Early_RMS_AVREC));
    AVREC_3(1:length(Data(7).(ID{i1}).AVREC_raw {BF}),i1,3)=Data(7).(ID{i1}).AVREC_raw {BF};
    RelRes_3 (1:length(Data(7).(ID{i1}).RELRES_raw {BF}),i1,3)=Data(7).(ID{i1}).RELRES_raw {BF};
    
    BF =find(Data(8).(ID{i1}).Early_RMS_AVREC == max(Data(8).(ID{i1}).Early_RMS_AVREC));
    AVREC_4(1:length(Data(8).(ID{i1}).AVREC_raw {BF}),i1,1)=Data(8).(ID{i1}).AVREC_raw {BF};
    RelRes_4 (1:length(Data(8).(ID{i1}).RELRES_raw {BF}),i1,1)=Data(8).(ID{i1}).RELRES_raw {BF};
    
    BF =find(Data(9).(ID{i1}).Early_RMS_AVREC == max(Data(9).(ID{i1}).Early_RMS_AVREC));
    AVREC_5(1:length(Data(9).(ID{i1}).AVREC_raw {BF}),i1,2)=Data(9).(ID{i1}).AVREC_raw {BF};
    RelRes_5 (1:length(Data(9).(ID{i1}).RELRES_raw {BF}),i1,2)=Data(9).(ID{i1}).RELRES_raw {BF};
    try
        BF =find(Data(10).(ID{i1}).Early_RMS_AVREC == max(Data(10).(ID{i1}).Early_RMS_AVREC));
    AVREC_5(1:length(Data(10).(ID{i1}).AVREC_raw {BF}),i1,1)=Data(10).(ID{i1}).AVREC_raw {BF};
    RelRes_5 (1:length(Data(10).(ID{i1}).RELRES_raw {BF}),i1,1)=Data(10).(ID{i1}).RELRES_raw {BF};
    catch,end
    try
        BF =find(Data(11).(ID{i1}).Early_RMS_AVREC == max(Data(11).(ID{i1}).Early_RMS_AVREC));
    AVREC_5(1:length(Data(11).(ID{i1}).AVREC_raw {BF}),i1,2)=Data(11).(ID{i1}).AVREC_raw {BF};
    RelRes_5 (1:length(Data(11).(ID{i1}).RELRES_raw {BF}),i1,2)=Data(11).(ID{i1}).RELRES_raw {BF};
    catch,end
end
AVREC_1 = nanmean(AVREC_1,3);
AVREC_3 = nanmean(AVREC_3,3);
AVREC_4 = nanmean(AVREC_4,3);
AVREC_5 = nanmean(AVREC_5,3);
RelRes_1 = nanmean(RelRes_1,3);
RelRes_3 = nanmean(RelRes_3,3);
RelRes_4 = nanmean(RelRes_4,3);
RelRes_5 = nanmean(RelRes_5,3);

figure
subplot(2,1,1)
hold on
area([200 250 ],[0.002 0.002 ],'FaceColor','y','FaceAlpha',0.05)
area([280 500 ],[0.002 0.002 ],'FaceColor','r','FaceAlpha',0.05)
plot(nanmean(AVREC_1'),'LineWidth',2)
plot(nanmean(AVREC_2'),'LineWidth',2)
plot(nanmean(AVREC_3'),'LineWidth',2)
plot(nanmean(AVREC_4'),'LineWidth',2)
plot(nanmean(AVREC_5'),'LineWidth',2)
plot([200 200],[0 0.0025],'--w')
plot([400 400],[0 0.0025],'--w')
xlim([100 500])
xticks([100 150 200 250 300 350 400 450 500])
xticklabels({'-100','-50', '0','50', '100','150', '200','250','300'})
set(gca,'fontsize',12,'FontWeight','bold')
xlabel('Time [ms]','FontSize',16,'FontWeight','bold')
ylabel('AVREC [mV/mm²]','FontSize',16,'FontWeight','bold')
title('AVREC Control n = 7','FontSize',20,'FontWeight','bold')
legend('Early','Late','Pre 1-3','Opto','Post 1-3','Post 4-5','Post 6-7','Location','best')

subplot(2,1,2)
hold on
area([200  200 250 250],[-30 15 15 -30],'FaceColor','y','FaceAlpha',0.05)
area([280 280 500 500],[-30 15 15 -30],'FaceColor','r','FaceAlpha',0.05)
plot(nanmean(RelRes_1'*100),'LineWidth',2)
plot(nanmean(RelRes_2'*100),'LineWidth',2)
plot(nanmean(RelRes_3'*100),'LineWidth',2)
plot(nanmean(RelRes_4'*100),'LineWidth',2)
plot(nanmean(RelRes_5'*100),'LineWidth',2)
plot([200 200],[-30 15],'--w')
plot([400 400],[-30 15],'--w')
xlim([100 500])
xticks([100 150 200 250 300 350 400 450 500])
xticklabels({'-100','-50', '0','50', '100','150', '200','250','300'})
set(gca,'fontsize',12,'FontWeight','bold')
xlabel('Time [ms]','FontSize',16,'FontWeight','bold')
ylabel('RelResCSD [%]','FontSize',16,'FontWeight','bold')
title('RelResCSD Control n = 7','FontSize',20,'FontWeight','bold')
legend('Early','Late','Pre 1-3','Opto','Post 1-3','Post 4-5','Post 6-7','Location','best')
















% %%%keyboard
clear


% Fig3
SORT = 'ST_based'; % GS_based

load('K:\CSD_dynamic_analysis\DATA\Output\Output_Input_Control_7post_n7_Data_Threshold_25_Zscore_0_binned_1_mirror_1.mat');
BF =  Data.BF_Pos;

Control_M= nan(length(Data.names),length(Data.GS_based),3);

CNorm = [];

for i1 = 1:3
    CNorm(:,i1,1) = Data.(SORT)(i1).Full_RMS_AVREC(:,BF); 
    CNorm(:,i1,2) = Data.(SORT)(i1).Full_RMS_AVREC(:,BF+1);
    CNorm(:,i1,3) = Data.(SORT)(i1).Full_RMS_AVREC(:,BF+2);
end

for i1 = 1:length(Data.names)
    Norm = nanmean(CNorm(i1,:,:));
    for i2 = 1:length(Data.GS_based)
        Control_M(i1,i2,1) = Data.(SORT)(i2).Full_RMS_AVREC(i1,BF)/Norm(1);
        Control_M(i1,i2,2) = Data.(SORT)(i2).Full_RMS_AVREC(i1,BF+1)/Norm(2);
        Control_M(i1,i2,3) = Data.(SORT)(i2).Full_RMS_AVREC(i1,BF+2)/Norm(3);
    end
end

C_SEM = nanstd(Control_M)./sqrt(sum(~isnan(Control_M)));

load('K:\CSD_dynamic_analysis\DATA\Output\Output_Input_HighP_7post_Data_Threshold_25_Zscore_0_binned_1_mirror_1.mat');

BF =  Data.BF_Pos;

Laser_M= nan(length(Data.names),length(Data.GS_based),3); 

lNorm = [];

for i1 = 1:3
    lNorm(:,i1,1) = Data.(SORT)(i1).Full_RMS_AVREC(:,BF); 
    lNorm(:,i1,2) = Data.(SORT)(i1).Full_RMS_AVREC(:,BF+1);
    lNorm(:,i1,3) = Data.(SORT)(i1).Full_RMS_AVREC(:,BF+2);
end

for i1 = 1:length(Data.names)
    Norm = nanmean(lNorm(i1,:,:));
    for i2 = 1:length(Data.GS_based)
        Laser_M(i1,i2,1) = Data.(SORT)(i2).Full_RMS_AVREC(i1,BF)/Norm(1);
        Laser_M(i1,i2,2) = Data.(SORT)(i2).Full_RMS_AVREC(i1,BF+1)/Norm(2);
        Laser_M(i1,i2,3) = Data.(SORT)(i2).Full_RMS_AVREC(i1,BF+2)/Norm(3);
    end
end

All =vertcat(CNorm,lNorm);
All2 = nanmean(All);
All3 = All./All2;%normalized
AllMean = nanmean(nanmean(All3,2),1);
AllSTD_1 = mean(std(All3(:,:,1)'))./sqrt(19);
AllSTD_2 = mean(std(All3(:,:,2)'))./sqrt(19);
AllSTD_3 = mean(std(All3(:,:,3)'))./sqrt(19);
% %%%keyboard
% ALLSTD = nanmean(nanstd(common'))

L_SEM = nanstd(Laser_M)./sqrt(sum(~isnan(Laser_M)));
STD_prime =vertcat(Control_M(:,1:3,:),Laser_M(:,1:3,:));
% %%%keyboard

X1 = Laser_M(:,:,1); 
X2 = nan(size(X1));
X2(1:size(Control_M(:,:,1),1),:)= Control_M(:,:,1);
PRE_L = mean(X1(:,1:3),2);
PRE_C = mean(X2(:,1:3),2);
pVal_L =[];
pVal_C =[];
pVal_Both = [];

for i1 =1:11
    
    %Laser with Pre
%      X3 =[PRE_L,X1(:,i1)];
%     L = [1,2];
%     O = teg_repeated_measures_ANOVA(X3,L,{'Time','Measurement'})
%      S=anovan([1:12]',{PRE_L,X1(:,i1)});
%     S=anovan(X3,{'Pre','M'});
    [h,p] = ttest2(PRE_L,X1(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_L =[pVal_L Y];
    
    %Control with Pre
%     X3 =[PRE_C,X2(:,i1)];
    [h,p] = ttest2(PRE_C,X2(:,i1));
%     S=anovan(X3,{'Pre','M'});
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_C =[pVal_C Y];
    
    % Both
%     X3 =[X1(:,i1),X2(:,i1)];   
%     S=anovan(X3,{'Opto','Control'});
    [h,p] = ttest2(X1(:,i1),X2(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_Both =[pVal_Both Y]; 
end
 
stats_BF.pVal_L = pVal_L; stats_BF.pVal_C = pVal_C; stats_BF.pVal_Both = pVal_Both; 

X1 = Laser_M(:,:,2); 
X2 = nan(size(X1));
X2(1:size(Control_M(:,:,1),1),:)= Control_M(:,:,2);
PRE_L = mean(X1(:,1:3),2);
PRE_C = mean(X2(:,1:3),2);
pVal_L =[];
pVal_C =[];
pVal_Both = [];

for i1 =1:11
    
    %Laser with Pre
%     X3 =[PRE_L,X1(:,i1)];
    %L = [1,2];
    %O = teg_repeated_measures_ANOVA(X3,L,{'Time','Measurement'})
%     S=anovan([1:12]',{PRE_L,X1(:,i1)});
%     S=anovan(X3,{'Pre','M'});
    [h,p] = ttest2(PRE_L,X1(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_L =[pVal_L Y];
    
    %Control with Pre
%     X3 =[PRE_C,X2(:,i1)];
    [h,p] = ttest2(PRE_C,X2(:,i1));
%     S=anovan(X3,{'Pre','M'});
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_C =[pVal_C Y];
    
    % Both
%     X3 =[X1(:,i1),X2(:,i1)];   
%     S=anovan(X3,{'Opto','Control'});
    [h,p] = ttest2(X1(:,i1),X2(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_Both =[pVal_Both Y]; 
end
    
stats_NearBF.pVal_L = pVal_L; stats_NearBF.pVal_C = pVal_C; stats_NearBF.pVal_Both = pVal_Both; 


X1 = Laser_M(:,:,3); 
X2 = nan(size(X1));
X2(1:size(Control_M(:,:,1),1),:)= Control_M(:,:,3);
PRE_L = mean(X1(:,1:3),2);
PRE_C = mean(X2(:,1:3),2);
pVal_L =[];
pVal_C =[];
pVal_Both = [];

for i1 =1:11
    
    %Laser with Pre
%     X3 =[PRE_L,X1(:,i1)];
    %L = [1,2];
    %O = teg_repeated_measures_ANOVA(X3,L,{'Time','Measurement'})
%     S=anovan([1:12]',{PRE_L,X1(:,i1)});
%     S=anovan(X3,{'Pre','M'});
    [h,p] = ttest2(PRE_L,X1(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_L =[pVal_L Y];
    
    %Control with Pre
%     X3 =[PRE_C,X2(:,i1)];
    [h,p] = ttest2(PRE_C,X2(:,i1));
%     S=anovan(X3,{'Pre','M'});
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_C =[pVal_C Y];
    
    % Both
%     X3 =[X1(:,i1),X2(:,i1)];   
%     S=anovan(X3,{'Opto','Control'});
    [h,p] = ttest2(X1(:,i1),X2(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_Both =[pVal_Both Y]; 
end
    
stats_NonBF.pVal_L = pVal_L; stats_NonBF.pVal_C = pVal_C; stats_NonBF.pVal_Both = pVal_Both; 

% %%keyboard

figure
subplot(2,3,1)
S1 =Laser_M(:,:,1);
S2 = nan(size(S1));
S2(1:7,:) = Control_M(:,:,1);
[stat,pperm,hperm,statperm]= elbe_permtest_005(S1,S2,1,'raw',1000,0,0.05,0);
% keyboard
SigPos=find(pperm <= 0.025 | pperm >= 1.025);;
shadedErrorBar([],nanmean(Control_M(:,:,1),1),C_SEM(:,:,1),'lineprops',{'b','markerfacecolor','b'})
hold on
shadedErrorBar([],nanmean(Laser_M(:,:,1),1),L_SEM(:,:,1),'lineprops',{'r','markerfacecolor','r'})
prime_M = mean(nanmean(STD_prime(:,:,1))');
prime_STD = mean(nanstd(STD_prime(:,:,1)'));
plot([SigPos],[ones(size(SigPos))*1.12],'k*')
% prime_STD = std(nanmean(STD_prime(:,:,1)))./sqrt(sum(~isnan(STD_prime(:,:,1))));
% prime_STD = mean(prime_STD);
% plot([1 length(Data.GS_based)],[prime_M+3*AllSTD_1 prime_M+3*AllSTD_1 ],':k')
% plot([1 length(Data.GS_based)],[prime_M-3*AllSTD_1 prime_M-3*AllSTD_1 ],':k')
% y1 = ones(1,11)*1.15; y2 = ones(1,11)*1.17;y3 = ones(1,11)*1.19;
% 
% plot([1:11].*(~isnan(stats_BF.pVal_L)),y1.*(~isnan(stats_BF.pVal_L)),'y*',...
%     [1:11].*(~isnan(stats_BF.pVal_C)),y2.*(~isnan(stats_BF.pVal_C)),'b*',...
%     [1:11].*(~isnan(stats_BF.pVal_Both)).*(~isnan(stats_BF.pVal_Both)),y3,'k*')
ylim([0.85 1.2])
xlim ([1 length(Data.GS_based)])
xticks([1 2 3 4 5 6 7 8 9 10 11])
xticklabels({'Pre(n-2)','Pre(n-1)', 'Pre(n)','Laser', 'Post 1','Post 2', 'Post 3','Post 4','Post 5','Post 6','Post 7'})
xtickangle(270)
set(gca,'fontsize',12,'FontWeight','bold')

title('FULL RMS AVREC BF')
legend ('Control','','','','Opto','','','','Location','best')


subplot(2,3,2)
S1 =Laser_M(:,:,2);
S2 = nan(size(S1));
S2(1:7,:) = Control_M(:,:,2);
[stat,pperm,hperm,statperm]= elbe_permtest_005(S1,S2,1,'raw',1000,0,0.05,0);
%keyboard
SigPos=find(pperm <= 0.025 | pperm >= 1.025);
shadedErrorBar([],nanmean(Control_M(:,:,2),1),C_SEM(:,:,2),'lineprops',{'b','markerfacecolor','b'})
hold on
shadedErrorBar([],nanmean(Laser_M(:,:,2),1),L_SEM(:,:,2),'lineprops',{'r','markerfacecolor','r'})
prime_M = mean(nanmean(STD_prime(:,:,2))');
prime_STD = mean(nanstd(STD_prime(:,:,2)'));
plot([SigPos],[ones(size(SigPos))*1.12],'k*')
% prime_STD = std(nanmean(STD_prime(:,:,2)))./sqrt(sum(~isnan(STD_prime(:,:,2))));
% prime_STD = mean(prime_STD);
% plot([1 length(Data.GS_based)],[prime_M+3*AllSTD_2 prime_M+3*AllSTD_2 ],':k')
% plot([1 length(Data.GS_based)],[prime_M-3*AllSTD_2 prime_M-3*AllSTD_2 ],':k')
% y1 = ones(1,11)*1.15; y2 = ones(1,11)*1.17;y3 = ones(1,11)*1.19;
% 
% plot([1:11].*(~isnan(stats_NearBF.pVal_L)),y1.*(~isnan(stats_NearBF.pVal_L)),'y*',...
%     [1:11].*(~isnan(stats_NearBF.pVal_C)),y2.*(~isnan(stats_NearBF.pVal_C)),'b*',...
%     [1:11].*(~isnan(stats_NearBF.pVal_Both)).*(~isnan(stats_NearBF.pVal_Both)),y3,'k*')
ylim([0.85 1.2])
xlim ([1 length(Data.GS_based)])
xticks([1 2 3 4 5 6 7 8 9 10 11])
xticklabels({'Pre(n-2)','Pre(n-1)', 'Pre(n)','Laser', 'Post 1','Post 2', 'Post 3','Post 4','Post 5','Post 6','Post 7'})
xtickangle(270)
set(gca,'fontsize',12,'FontWeight','bold')

title('± near BF (1-2 Oct)')
legend ('Control','','','','Opto','','','','Location','best')

subplot(2,3,3)
S1 =Laser_M(:,:,3);
S2 = nan(size(S1));
S2(1:7,:) = Control_M(:,:,3);
[stat,pperm,hperm,statperm]= elbe_permtest_005(S1,S2,1,'raw',1000,0,0.05,0);
SigPos=find(pperm <= 0.025 | pperm >= 1.025);
shadedErrorBar([],nanmean(Control_M(:,:,3),1),C_SEM(:,:,3),'lineprops',{'b','markerfacecolor','b'})
hold on
shadedErrorBar([],nanmean(Laser_M(:,:,3),1),L_SEM(:,:,3),'lineprops',{'r','markerfacecolor','r'})
prime_M = mean(nanmean(STD_prime(:,:,3)'));
prime_STD = mean(nanstd(STD_prime(:,:,3)'));
plot([SigPos],[ones(size(SigPos))*1.12],'k*')
% prime_STD = std(nanmean(STD_prime(:,:,3)))./sqrt(sum(~isnan(STD_prime(:,:,3))));
% prime_STD = mean(prime_STD);
% plot([1 length(Data.GS_based)],[prime_M+3*AllSTD_3 prime_M+3*AllSTD_3 ],':k')
% plot([1 length(Data.GS_based)],[prime_M-3*AllSTD_3 prime_M-3*AllSTD_3 ],':k')
% y1 = ones(1,11)*1.15; y2 = ones(1,11)*1.17;y3 = ones(1,11)*1.19;
% 
% plot([1:11].*(~isnan(stats_NonBF.pVal_L)),y1.*(~isnan(stats_NonBF.pVal_L)),'y*',...
%     [1:11].*(~isnan(stats_NonBF.pVal_C)),y2.*(~isnan(stats_NonBF.pVal_C)),'b*',...
%     [1:11].*(~isnan(stats_NonBF.pVal_Both)).*(~isnan(stats_NonBF.pVal_Both)),y3,'k*')
ylim([0.85 1.2])
xlim ([1 length(Data.GS_based)])
xticks([1 2 3 4 5 6 7 8 9 10 11])
xticklabels({'Pre(n-2)','Pre(n-1)', 'Pre(n)','Laser', 'Post 1','Post 2', 'Post 3','Post 4','Post 5','Post 6','Post 7'})
xtickangle(270)
set(gca,'fontsize',12,'FontWeight','bold')

title('± non BF (3-4 Oct)')
legend ('Control','','','','Opto','','','','Location','best')
%  %%%keyboard


%%% RELRES
SORT = 'ST_based'; % GS_based

load('K:\CSD_dynamic_analysis\DATA\Output\Output_Input_Control_7post_n7_Data_Threshold_25_Zscore_0_binned_1_mirror_1.mat');
BF =  Data.BF_Pos;

Control_M= nan(length(Data.names),length(Data.GS_based),3);

CNorm = [];

for i1 = 1:3
    CNorm(:,i1,1) = Data.(SORT)(i1).Full_RMS_RELRES(:,BF); 
    CNorm(:,i1,2) = Data.(SORT)(i1).Full_RMS_RELRES(:,BF+1);
    CNorm(:,i1,3) = Data.(SORT)(i1).Full_RMS_RELRES(:,BF+2);
end

for i1 = 1:length(Data.names)
    Norm = nanmean(CNorm(i1,:,:));
    for i2 = 1:length(Data.GS_based)
        Control_M(i1,i2,1) = Data.(SORT)(i2).Full_RMS_RELRES(i1,BF)/Norm(1);
        Control_M(i1,i2,2) = Data.(SORT)(i2).Full_RMS_RELRES(i1,BF+1)/Norm(2);
        Control_M(i1,i2,3) = Data.(SORT)(i2).Full_RMS_RELRES(i1,BF+2)/Norm(3);
    end
end

C_SEM = nanstd(Control_M)./sqrt(sum(~isnan(Control_M)));

load('K:\CSD_dynamic_analysis\DATA\Output\Output_Input_HighP_7post_Data_Threshold_25_Zscore_0_binned_1_mirror_1.mat');

BF =  Data.BF_Pos;

Laser_M= nan(length(Data.names),length(Data.GS_based),3); 

lNorm = [];

for i1 = 1:3
    lNorm(:,i1,1) = Data.(SORT)(i1).Full_RMS_RELRES(:,BF); 
    lNorm(:,i1,2) = Data.(SORT)(i1).Full_RMS_RELRES(:,BF+1);
    lNorm(:,i1,3) = Data.(SORT)(i1).Full_RMS_RELRES(:,BF+2);
end

for i1 = 1:length(Data.names)
    Norm = nanmean(lNorm(i1,:,:));
    for i2 = 1:length(Data.GS_based)
        Laser_M(i1,i2,1) = Data.(SORT)(i2).Full_RMS_RELRES(i1,BF)/Norm(1);
        Laser_M(i1,i2,2) = Data.(SORT)(i2).Full_RMS_RELRES(i1,BF+1)/Norm(2);
        Laser_M(i1,i2,3) = Data.(SORT)(i2).Full_RMS_RELRES(i1,BF+2)/Norm(3);
    end
end

All =vertcat(CNorm,lNorm);
All2 = nanmean(All);
All3 = All./All2;%normalized
AllMean = nanmean(nanmean(All3,2),1);
AllSTD_1 = mean(std(All3(:,:,1)'))./sqrt(19);
AllSTD_2 = mean(std(All3(:,:,2)'))./sqrt(19);
AllSTD_3 = mean(std(All3(:,:,3)'))./sqrt(19);
L_SEM = nanstd(Laser_M)./sqrt(sum(~isnan(Laser_M)));
STD_prime =vertcat(Control_M(:,1:3,:),Laser_M(:,1:3,:));


X1 = Laser_M(:,:,1); 
X2 = nan(size(X1));
X2(1:size(Control_M(:,:,1),1),:)= Control_M(:,:,1);
PRE_L = mean(X1(:,1:3),2);
PRE_C = mean(X2(:,1:3),2);
pVal_L =[];
pVal_C =[];
pVal_Both = [];

for i1 =1:11
    
    %Laser with Pre
%     X3 =[PRE_L,X1(:,i1)];
    %L = [1,2];
    %O = teg_repeated_measures_ANOVA(X3,L,{'Time','Measurement'})
%     S=anovan([1:12]',{PRE_L,X1(:,i1)});
%     S=anovan(X3,{'Pre','M'});
    [h,p] = ttest2(PRE_L,X1(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_L =[pVal_L Y];
    
    %Control with Pre
%     X3 =[PRE_C,X2(:,i1)];
    [h,p] = ttest2(PRE_C,X2(:,i1));
%     S=anovan(X3,{'Pre','M'});
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_C =[pVal_C Y];
    
    % Both
%     X3 =[X1(:,i1),X2(:,i1)];   
%     S=anovan(X3,{'Opto','Control'});
    [h,p] = ttest2(X1(:,i1),X2(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_Both =[pVal_Both Y]; 
end
    
stats_BF.pVal_L = pVal_L; stats_BF.pVal_C = pVal_C; stats_BF.pVal_Both = pVal_Both; 

X1 = Laser_M(:,:,2); 
X2 = nan(size(X1));
X2(1:size(Control_M(:,:,1),1),:)= Control_M(:,:,2);
PRE_L = mean(X1(:,1:3),2);
PRE_C = mean(X2(:,1:3),2);
pVal_L =[];
pVal_C =[];
pVal_Both = [];

for i1 =1:11
    
    %Laser with Pre
%     X3 =[PRE_L,X1(:,i1)];
    %L = [1,2];
    %O = teg_repeated_measures_ANOVA(X3,L,{'Time','Measurement'})
%     S=anovan([1:12]',{PRE_L,X1(:,i1)});
%     S=anovan(X3,{'Pre','M'});
    [h,p] = ttest2(PRE_L,X1(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_L =[pVal_L Y];
    
    %Control with Pre
%     X3 =[PRE_C,X2(:,i1)];
    [h,p] = ttest2(PRE_C,X2(:,i1));
%     S=anovan(X3,{'Pre','M'});
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_C =[pVal_C Y];
    
    % Both
%     X3 =[X1(:,i1),X2(:,i1)];   
%     S=anovan(X3,{'Opto','Control'});
    [h,p] = ttest2(X1(:,i1),X2(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_Both =[pVal_Both Y]; 
end
    
stats_NearBF.pVal_L = pVal_L; stats_NearBF.pVal_C = pVal_C; stats_NearBF.pVal_Both = pVal_Both; 


X1 = Laser_M(:,:,3); 
X2 = nan(size(X1));
X2(1:size(Control_M(:,:,1),1),:)= Control_M(:,:,3);
PRE_L = mean(X1(:,1:3),2);
PRE_C = mean(X2(:,1:3),2);
pVal_L =[];
pVal_C =[];
pVal_Both = [];

for i1 =1:11
    
    %Laser with Pre
%     X3 =[PRE_L,X1(:,i1)];
    %L = [1,2];
    %O = teg_repeated_measures_ANOVA(X3,L,{'Time','Measurement'})
%     S=anovan([1:12]',{PRE_L,X1(:,i1)});
%     S=anovan(X3,{'Pre','M'});
    [h,p] = ttest2(PRE_L,X1(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_L =[pVal_L Y];
    
    %Control with Pre
%     X3 =[PRE_C,X2(:,i1)];
    [h,p] = ttest2(PRE_C,X2(:,i1));
%     S=anovan(X3,{'Pre','M'});
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_C =[pVal_C Y];
    
    % Both
%     X3 =[X1(:,i1),X2(:,i1)];   
%     S=anovan(X3,{'Opto','Control'});
    [h,p] = ttest2(X1(:,i1),X2(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_Both =[pVal_Both Y]; 
end
    
stats_NonBF.pVal_L = pVal_L; stats_NonBF.pVal_C = pVal_C; stats_NonBF.pVal_Both = pVal_Both; 
% %%keyboard
subplot(2,3,4)
S1 =Laser_M(:,:,1);
S2 = nan(size(S1));
S2(1:7,:) = Control_M(:,:,1);
[stat,pperm,hperm,statperm]= elbe_permtest_005(S1,S2,1,'raw',1000,0,0.05,0);
SigPos=find(pperm <= 0.025 | pperm >= 1.025);
shadedErrorBar([],nanmean(Control_M(:,:,1),1),C_SEM(:,:,1),'lineprops',{'b','markerfacecolor','b'})
hold on
shadedErrorBar([],nanmean(Laser_M(:,:,1),1),L_SEM(:,:,1),'lineprops',{'r','markerfacecolor','r'})
prime_M = mean(nanmean(STD_prime(:,:,1)'));
prime_STD = mean(nanstd(STD_prime(:,:,1)'));
plot([SigPos],[ones(size(SigPos))*1.4],'k*')
% prime_STD = std(nanmean(STD_prime(:,:,1)))./sqrt(sum(~isnan(STD_prime(:,:,1))));
% prime_STD = mean(prime_STD);
% plot([1 length(Data.GS_based)],[prime_M+3*AllSTD_1 prime_M+3*AllSTD_1],':k')
% plot([1 length(Data.GS_based)],[prime_M-3*AllSTD_1 prime_M-3*AllSTD_1],':k')
% y1 = ones(1,11)*1.5; y2 = ones(1,11)*1.7;y3 = ones(1,11)*1.9;
% 
% plot([1:11].*(~isnan(stats_BF.pVal_L)),y1.*(~isnan(stats_BF.pVal_L)),'r*',...
%     [1:11].*(~isnan(stats_BF.pVal_C)),y2.*(~isnan(stats_BF.pVal_C)),'b*',...
%     [1:11].*(~isnan(stats_BF.pVal_Both)).*(~isnan(stats_BF.pVal_Both)),y3,'k*')
ylim([0.8 1.6])
xlim ([1 length(Data.GS_based)])
xticks([1 2 3 4 5 6 7 8 9 10 11])
xticklabels({'Pre(n-2)','Pre(n-1)', 'Pre(n)','Laser', 'Post 1','Post 2', 'Post 3','Post 4','Post 5','Post 6','Post 7'})
xtickangle(270)
set(gca,'fontsize',12,'FontWeight','bold')

title('FULL RMS RELRESCSD BF')
legend ('Control','','','','Opto','','','','Location','best')


subplot(2,3,5)
S1 =Laser_M(:,:,2);
S2 = nan(size(S1));
S2(1:7,:) = Control_M(:,:,2);
[stat,pperm,hperm,statperm]= elbe_permtest_005(S1,S2,1,'raw',1000,0,0.05,0);
SigPos=find(pperm <= 0.025 | pperm >= 1.025);
shadedErrorBar([],nanmean(Control_M(:,:,2),1),C_SEM(:,:,2),'lineprops',{'b','markerfacecolor','b'})
hold on
shadedErrorBar([],nanmean(Laser_M(:,:,2),1),L_SEM(:,:,2),'lineprops',{'r','markerfacecolor','r'})
prime_M = mean(nanmean(STD_prime(:,:,2)'));
prime_STD = mean(nanstd(STD_prime(:,:,2)'));
plot([SigPos],[ones(size(SigPos))*1.4],'k*')
% prime_STD = std(nanmean(STD_prime(:,:,2)))./sqrt(sum(~isnan(STD_prime(:,:,2))));
% prime_STD = mean(prime_STD);
% plot([1 length(Data.GS_based)],[prime_M+3*AllSTD_2 prime_M+3*AllSTD_2],':k')
% plot([1 length(Data.GS_based)],[prime_M-3*AllSTD_2 prime_M-3*AllSTD_2],':k')
% y1 = ones(1,11)*1.5; y2 = ones(1,11)*1.7;y3 = ones(1,11)*1.9;
% 
% plot([1:11].*(~isnan(stats_NearBF.pVal_L)),y1.*(~isnan(stats_NearBF.pVal_L)),'r*',...
%     [1:11].*(~isnan(stats_NearBF.pVal_C)),y2.*(~isnan(stats_NearBF.pVal_C)),'b*',...
%     [1:11].*(~isnan(stats_NearBF.pVal_Both)).*(~isnan(stats_NearBF.pVal_Both)),y3,'k*')
ylim([0.8 1.6])
xlim ([1 length(Data.GS_based)])
xticks([1 2 3 4 5 6 7 8 9 10 11])
xticklabels({'Pre(n-2)','Pre(n-1)', 'Pre(n)','Laser', 'Post 1','Post 2', 'Post 3','Post 4','Post 5','Post 6','Post 7'})
xtickangle(270)
set(gca,'fontsize',12,'FontWeight','bold')
title('± near BF (1-2 Oct)')
legend ('Control','','','','Opto','','','','Location','best')

subplot(2,3,6)
S1 =Laser_M(:,:,3);
S2 = nan(size(S1));
S2(1:7,:) = Control_M(:,:,3);
[stat,pperm,hperm,statperm]= elbe_permtest_005(S1,S2,1,'raw',1000,0,0.05,0);
SigPos=find(pperm <= 0.025 | pperm >= 1.025);
shadedErrorBar([],nanmean(Control_M(:,:,3),1),C_SEM(:,:,3),'lineprops',{'b','markerfacecolor','b'})
hold on
shadedErrorBar([],nanmean(Laser_M(:,:,3),1),L_SEM(:,:,3),'lineprops',{'r','markerfacecolor','r'})
prime_M = mean(nanmean(STD_prime(:,:,3)'));
prime_STD = mean(nanstd(STD_prime(:,:,3)'));
plot([SigPos],[ones(size(SigPos))*1.4],'k*');
% prime_STD = std(nanmean(STD_prime(:,:,3)))./sqrt(sum(~isnan(STD_prime(:,:,3))));
% prime_STD = mean(prime_STD);
% plot([1 length(Data.GS_based)],[prime_M+3*AllSTD_3 prime_M+3*AllSTD_3],':k')
% plot([1 length(Data.GS_based)],[prime_M-3*AllSTD_3 prime_M-3*AllSTD_3],':k')
% y1 = ones(1,11)*1.5; y2 = ones(1,11)*1.7;y3 = ones(1,11)*1.9;
% 
% plot([1:11].*(~isnan(stats_NonBF.pVal_L)),y1.*(~isnan(stats_NonBF.pVal_L)),'r*',...
%     [1:11].*(~isnan(stats_NonBF.pVal_C)),y2.*(~isnan(stats_NonBF.pVal_C)),'b*',...
%     [1:11].*(~isnan(stats_NonBF.pVal_Both)).*(~isnan(stats_NonBF.pVal_Both)),y3,'k*')
ylim([0.8 1.6])
xlim ([1 length(Data.GS_based)])
xticks([1 2 3 4 5 6 7 8 9 10 11])
xticklabels({'Pre(n-2)','Pre(n-1)', 'Pre(n)','Laser', 'Post 1','Post 2', 'Post 3','Post 4','Post 5','Post 6','Post 7'})
xtickangle(270)
set(gca,'fontsize',12,'FontWeight','bold')

title('± non BF (3-4 Oct)')
legend ('Control','','','','Opto','','','','Location','best')
h=gcf;

set(h,'PaperOrientation','landscape');

set(h,'PaperPosition', [1 1 28 19]);



SORT = 'ST_based'; % GS_based

load('K:\CSD_dynamic_analysis\DATA\Output\Output_Input_Control_7post_n7_Data_Threshold_25_Zscore_0_binned_1_mirror_1.mat');
BF =  Data.BF_Pos;

Control_M= nan(length(Data.names),length(Data.GS_based),3);

CNorm = [];

for i1 = 1:3
    CNorm(:,i1,1) = Data.(SORT)(i1).Early_RMS_AVREC(:,BF); 
    CNorm(:,i1,2) = Data.(SORT)(i1).Early_RMS_AVREC(:,BF+1);
    CNorm(:,i1,3) = Data.(SORT)(i1).Early_RMS_AVREC(:,BF+2);
end

for i1 = 1:length(Data.names)
    Norm = nanmean(CNorm(i1,:,:));
    for i2 = 1:length(Data.GS_based)
        Control_M(i1,i2,1) = Data.(SORT)(i2).Early_RMS_AVREC(i1,BF)/Norm(1);
        Control_M(i1,i2,2) = Data.(SORT)(i2).Early_RMS_AVREC(i1,BF+1)/Norm(2);
        Control_M(i1,i2,3) = Data.(SORT)(i2).Early_RMS_AVREC(i1,BF+2)/Norm(3);
    end
end

C_SEM = nanstd(Control_M)./sqrt(sum(~isnan(Control_M)));

load('K:\CSD_dynamic_analysis\DATA\Output\Output_Input_HighP_7post_Data_Threshold_25_Zscore_0_binned_1_mirror_1.mat');

BF =  Data.BF_Pos;

Laser_M= nan(length(Data.names),length(Data.GS_based),3); 

lNorm = [];

for i1 = 1:3
    lNorm(:,i1,1) = Data.(SORT)(i1).Early_RMS_AVREC(:,BF); 
    lNorm(:,i1,2) = Data.(SORT)(i1).Early_RMS_AVREC(:,BF+1);
    lNorm(:,i1,3) = Data.(SORT)(i1).Early_RMS_AVREC(:,BF+2);
end

for i1 = 1:length(Data.names)
    Norm = nanmean(lNorm(i1,:,:));
    for i2 = 1:length(Data.GS_based)
        Laser_M(i1,i2,1) = Data.(SORT)(i2).Early_RMS_AVREC(i1,BF)/Norm(1);
        Laser_M(i1,i2,2) = Data.(SORT)(i2).Early_RMS_AVREC(i1,BF+1)/Norm(2);
        Laser_M(i1,i2,3) = Data.(SORT)(i2).Early_RMS_AVREC(i1,BF+2)/Norm(3);
    end
end

All =vertcat(CNorm,lNorm);
All2 = nanmean(All);
All3 = All./All2;%normalized
AllMean = nanmean(nanmean(All3,2),1);
AllSTD_1 = mean(std(All3(:,:,1)'))./sqrt(19);
AllSTD_2 = mean(std(All3(:,:,2)'))./sqrt(19);
AllSTD_3 = mean(std(All3(:,:,3)'))./sqrt(19);
% %%%keyboard
% ALLSTD = nanmean(nanstd(common'))

L_SEM = nanstd(Laser_M)./sqrt(sum(~isnan(Laser_M)));
STD_prime =vertcat(Control_M(:,1:3,:),Laser_M(:,1:3,:));
% %%%keyboard

X1 = Laser_M(:,:,1); 
X2 = nan(size(X1));
X2(1:size(Control_M(:,:,1),1),:)= Control_M(:,:,1);
PRE_L = mean(X1(:,1:3),2);
PRE_C = mean(X2(:,1:3),2);
pVal_L =[];
pVal_C =[];
pVal_Both = [];

for i1 =1:11
    
    %Laser with Pre
%      X3 =[PRE_L,X1(:,i1)];
%     L = [1,2];
%     O = teg_repeated_measures_ANOVA(X3,L,{'Time','Measurement'})
%      S=anovan([1:12]',{PRE_L,X1(:,i1)});
%     S=anovan(X3,{'Pre','M'});
    [h,p] = ttest2(PRE_L,X1(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_L =[pVal_L Y];
    
    %Control with Pre
%     X3 =[PRE_C,X2(:,i1)];
    [h,p] = ttest2(PRE_C,X2(:,i1));
%     S=anovan(X3,{'Pre','M'});
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_C =[pVal_C Y];
    
    % Both
%     X3 =[X1(:,i1),X2(:,i1)];   
%     S=anovan(X3,{'Opto','Control'});
    [h,p] = ttest2(X1(:,i1),X2(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_Both =[pVal_Both Y]; 
end
 
stats_BF.pVal_L = pVal_L; stats_BF.pVal_C = pVal_C; stats_BF.pVal_Both = pVal_Both; 

X1 = Laser_M(:,:,2); 
X2 = nan(size(X1));
X2(1:size(Control_M(:,:,1),1),:)= Control_M(:,:,2);
PRE_L = mean(X1(:,1:3),2);
PRE_C = mean(X2(:,1:3),2);
pVal_L =[];
pVal_C =[];
pVal_Both = [];

for i1 =1:11
    
    %Laser with Pre
%     X3 =[PRE_L,X1(:,i1)];
    %L = [1,2];
    %O = teg_repeated_measures_ANOVA(X3,L,{'Time','Measurement'})
%     S=anovan([1:12]',{PRE_L,X1(:,i1)});
%     S=anovan(X3,{'Pre','M'});
    [h,p] = ttest2(PRE_L,X1(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_L =[pVal_L Y];
    
    %Control with Pre
%     X3 =[PRE_C,X2(:,i1)];
    [h,p] = ttest2(PRE_C,X2(:,i1));
%     S=anovan(X3,{'Pre','M'});
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_C =[pVal_C Y];
    
    % Both
%     X3 =[X1(:,i1),X2(:,i1)];   
%     S=anovan(X3,{'Opto','Control'});
    [h,p] = ttest2(X1(:,i1),X2(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_Both =[pVal_Both Y]; 
end
    
stats_NearBF.pVal_L = pVal_L; stats_NearBF.pVal_C = pVal_C; stats_NearBF.pVal_Both = pVal_Both; 


X1 = Laser_M(:,:,3); 
X2 = nan(size(X1));
X2(1:size(Control_M(:,:,1),1),:)= Control_M(:,:,3);
PRE_L = mean(X1(:,1:3),2);
PRE_C = mean(X2(:,1:3),2);
pVal_L =[];
pVal_C =[];
pVal_Both = [];

for i1 =1:11
    
    %Laser with Pre
%     X3 =[PRE_L,X1(:,i1)];
    %L = [1,2];
    %O = teg_repeated_measures_ANOVA(X3,L,{'Time','Measurement'})
%     S=anovan([1:12]',{PRE_L,X1(:,i1)});
%     S=anovan(X3,{'Pre','M'});
    [h,p] = ttest2(PRE_L,X1(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_L =[pVal_L Y];
    
    %Control with Pre
%     X3 =[PRE_C,X2(:,i1)];
    [h,p] = ttest2(PRE_C,X2(:,i1));
%     S=anovan(X3,{'Pre','M'});
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_C =[pVal_C Y];
    
    % Both
%     X3 =[X1(:,i1),X2(:,i1)];   
%     S=anovan(X3,{'Opto','Control'});
    [h,p] = ttest2(X1(:,i1),X2(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_Both =[pVal_Both Y]; 
end
    
stats_NonBF.pVal_L = pVal_L; stats_NonBF.pVal_C = pVal_C; stats_NonBF.pVal_Both = pVal_Both; 

% keyboard %Baustelle
figure
clear V
subplot(2,3,1)
S1 =Laser_M(:,:,1);
S2 = nan(size(S1));
S2(1:7,:) = Control_M(:,:,1);
% [stat,pperm,hperm,statperm]= elbe_permtest_005(S1,S2,1,'raw',1000,0,0.05,0);
% SigPos=find(pperm <= 0.025 | pperm >= 1.025);
L = [2 11];
S3 = [S1, S2];
O = teg_repeated_measures_ANOVA(S3,L,{'Group','Measurement'});
[p,n]=ttest(S1,S2);
try 
V =reshape(O.groups,2,[])'; Vlog =(V <= 11);
Vlog =Vlog(:,1)+Vlog(:,2);
Vlog = Vlog ~= 2;V = V.*Vlog; V(V == 0) = nan;V(:,1) =V(:,1)-11;
V =V(V(:,1) == V(:,2));
catch,
end

shadedErrorBar([],nanmean(Control_M(:,:,1),1),C_SEM(:,:,1),'lineprops',{'b','markerfacecolor','b'})
hold on
shadedErrorBar([],nanmean(Laser_M(:,:,1),1),L_SEM(:,:,1),'lineprops',{'r','markerfacecolor','r'})
prime_M = mean(nanmean(STD_prime(:,:,1))');
prime_STD = mean(nanstd(STD_prime(:,:,1)'));
plot(1:length(p),p*1.2,'k*')
plot(V,ones(size(V))*1.3,'r*')

T1 = [O.labels{1} ' F(' num2str(O.R(1,2)) ','  num2str(O.R(1,3)) ') = ' num2str(O.R(1,1)) ' p = ' num2str(O.R(1,4))];
T2 = [O.labels{2} ' F(' num2str(O.R(2,2)) ','  num2str(O.R(2,3)) ') = ' num2str(O.R(2,1)) ' p = ' num2str(O.R(2,4))];
T3 = [O.labels{3} ' F(' num2str(O.R(3,2)) ','  num2str(O.R(3,3)) ') = ' num2str(O.R(3,1)) ' p = ' num2str(O.R(3,4))];
T =strvcat(T1, T2, T3);
text(1.2,0.9, T);
 
% prime_STD = std(nanmean(STD_prime(:,:,1)))./sqrt(sum(~isnan(STD_prime(:,:,1))));
% prime_STD = mean(prime_STD);
% plot([1 length(Data.GS_based)],[prime_M+3*AllSTD_1 prime_M+3*AllSTD_1 ],':k')
% plot([1 length(Data.GS_based)],[prime_M-3*AllSTD_1 prime_M-3*AllSTD_1 ],':k')
% y1 = ones(1,11)*1.15; y2 = ones(1,11)*1.17;y3 = ones(1,11)*1.19;
% 
% plot([1:11].*(~isnan(stats_BF.pVal_L)),y1.*(~isnan(stats_BF.pVal_L)),'y*',...
%     [1:11].*(~isnan(stats_BF.pVal_C)),y2.*(~isnan(stats_BF.pVal_C)),'b*',...
%     [1:11].*(~isnan(stats_BF.pVal_Both)).*(~isnan(stats_BF.pVal_Both)),y3,'k*')
ylim([0.85 1.3])
xlim ([1 length(Data.GS_based)])
xticks([1 2 3 4 5 6 7 8 9 10 11])
xticklabels({'Pre(n-2)','Pre(n-1)', 'Pre(n)','Laser', 'Post 1','Post 2', 'Post 3','Post 4','Post 5','Post 6','Post 7'})
xtickangle(270)
set(gca,'fontsize',12,'FontWeight','bold')

title('Early RMS AVREC BF')
legend ('Control','','','','Opto','','','','Location','best')


clear V
subplot(2,3,2)
S1 =Laser_M(:,:,2);
S2 = nan(size(S1));
S2(1:7,:) = Control_M(:,:,2);
% [stat,pperm,hperm,statperm]= elbe_permtest_005(S1,S2,1,'raw',1000,0,0.05,0);
% SigPos=find(pperm <= 0.025 | pperm >= 1.025);
L = [2 11];
S3 = [S1, S2];
O = teg_repeated_measures_ANOVA(S3,L,{'Group','Measurement'});
[p,n]=ttest(S1,S2);
V =reshape(O.groups,2,[])'; Vlog =(V <= 11);
Vlog =Vlog(:,1)+Vlog(:,2);
Vlog = Vlog ~= 2;V = V.*Vlog; V(V == 0) = nan;V(:,1) =V(:,1)-11;
V =V(V(:,1) == V(:,2));

shadedErrorBar([],nanmean(Control_M(:,:,2),1),C_SEM(:,:,2),'lineprops',{'b','markerfacecolor','b'})
hold on
shadedErrorBar([],nanmean(Laser_M(:,:,2),1),L_SEM(:,:,2),'lineprops',{'r','markerfacecolor','r'})
prime_M = mean(nanmean(STD_prime(:,:,2))');
prime_STD = mean(nanstd(STD_prime(:,:,2)'));
plot(1:length(p),p*1.2,'k*')
plot(V,ones(size(V))*1.3,'r*')
T1 = [O.labels{1} ' F(' num2str(O.R(1,2)) ','  num2str(O.R(1,3)) ') = ' num2str(O.R(1,1)) ' p = ' num2str(O.R(1,4))];
T2 = [O.labels{2} ' F(' num2str(O.R(2,2)) ','  num2str(O.R(2,3)) ') = ' num2str(O.R(2,1)) ' p = ' num2str(O.R(2,4))];
T3 = [O.labels{3} ' F(' num2str(O.R(3,2)) ','  num2str(O.R(3,3)) ') = ' num2str(O.R(3,1)) ' p = ' num2str(O.R(3,4))];
T =strvcat(T1, T2, T3);
text(1.2,0.9, T);
% prime_STD = std(nanmean(STD_prime(:,:,2)))./sqrt(sum(~isnan(STD_prime(:,:,2))));
% prime_STD = mean(prime_STD);
% plot([1 length(Data.GS_based)],[prime_M+3*AllSTD_2 prime_M+3*AllSTD_2 ],':k')
% plot([1 length(Data.GS_based)],[prime_M-3*AllSTD_2 prime_M-3*AllSTD_2 ],':k')
% y1 = ones(1,11)*1.15; y2 = ones(1,11)*1.17;y3 = ones(1,11)*1.19;
% 
% plot([1:11].*(~isnan(stats_NearBF.pVal_L)),y1.*(~isnan(stats_NearBF.pVal_L)),'y*',...
%     [1:11].*(~isnan(stats_NearBF.pVal_C)),y2.*(~isnan(stats_NearBF.pVal_C)),'b*',...
%     [1:11].*(~isnan(stats_NearBF.pVal_Both)).*(~isnan(stats_NearBF.pVal_Both)),y3,'k*')
ylim([0.85 1.3])
xlim ([1 length(Data.GS_based)])
xticks([1 2 3 4 5 6 7 8 9 10 11])
xticklabels({'Pre(n-2)','Pre(n-1)', 'Pre(n)','Laser', 'Post 1','Post 2', 'Post 3','Post 4','Post 5','Post 6','Post 7'})
xtickangle(270)
set(gca,'fontsize',12,'FontWeight','bold')

title('± near BF (1-2 Oct)')
% legend ('Control','','','','Opto','','','','Location','best')

clear V
subplot(2,3,3)
S1 =Laser_M(:,:,3);
S2 = nan(size(S1));
S2(1:7,:) = Control_M(:,:,3);
L = [2 11];
S3 = [S1, S2];
O = teg_repeated_measures_ANOVA(S3,L,{'Group','Measurement'});
[p,n]=ttest(S1,S2);
try
V =reshape(O.groups,2,[])'; Vlog =(V <= 11);
Vlog =Vlog(:,1)+Vlog(:,2);
Vlog = Vlog ~= 2;V = V.*Vlog; V(V == 0) = nan;V(:,1) =V(:,1)-11;
V =V(V(:,1) == V(:,2));
end
% [stat,pperm,hperm,statperm]= elbe_permtest_005(S1,S2,1,'raw',1000,0,0.05,0);
% SigPos=find(pperm <= 0.025 | pperm >= 1.025);
shadedErrorBar([],nanmean(Control_M(:,:,3),1),C_SEM(:,:,3),'lineprops',{'b','markerfacecolor','b'})
hold on
shadedErrorBar([],nanmean(Laser_M(:,:,3),1),L_SEM(:,:,3),'lineprops',{'r','markerfacecolor','r'})
prime_M = mean(nanmean(STD_prime(:,:,3)'));
prime_STD = mean(nanstd(STD_prime(:,:,3)'));
plot(1:length(p),p*1.2,'k*')
try
 plot(V,ones(size(V))*1.3,'r*')
end
T1 = [O.labels{1} ' F(' num2str(O.R(1,2)) ','  num2str(O.R(1,3)) ') = ' num2str(O.R(1,1)) ' p = ' num2str(O.R(1,4))];
T2 = [O.labels{2} ' F(' num2str(O.R(2,2)) ','  num2str(O.R(2,3)) ') = ' num2str(O.R(2,1)) ' p = ' num2str(O.R(2,4))];
T3 = [O.labels{3} ' F(' num2str(O.R(3,2)) ','  num2str(O.R(3,3)) ') = ' num2str(O.R(3,1)) ' p = ' num2str(O.R(3,4))];
T =strvcat(T1, T2, T3);
text(1.2,0.9, T);

% prime_STD = std(nanmean(STD_prime(:,:,3)))./sqrt(sum(~isnan(STD_prime(:,:,3))));
% prime_STD = mean(prime_STD);
% plot([1 length(Data.GS_based)],[prime_M+3*AllSTD_3 prime_M+3*AllSTD_3 ],':k')
% plot([1 length(Data.GS_based)],[prime_M-3*AllSTD_3 prime_M-3*AllSTD_3 ],':k')
% y1 = ones(1,11)*1.15; y2 = ones(1,11)*1.17;y3 = ones(1,11)*1.19;
% 
% plot([1:11].*(~isnan(stats_NonBF.pVal_L)),y1.*(~isnan(stats_NonBF.pVal_L)),'y*',...
%     [1:11].*(~isnan(stats_NonBF.pVal_C)),y2.*(~isnan(stats_NonBF.pVal_C)),'b*',...
%     [1:11].*(~isnan(stats_NonBF.pVal_Both)).*(~isnan(stats_NonBF.pVal_Both)),y3,'k*')
ylim([0.85 1.3])
xlim ([1 length(Data.GS_based)])
xticks([1 2 3 4 5 6 7 8 9 10 11])
xticklabels({'Pre(n-2)','Pre(n-1)', 'Pre(n)','Laser', 'Post 1','Post 2', 'Post 3','Post 4','Post 5','Post 6','Post 7'})
xtickangle(270)
set(gca,'fontsize',12,'FontWeight','bold')

title('± non BF (3-4 Oct)')
% legend ('Control','','','','Opto','','','','Location','best')
%  %%%keyboard


%%% RELRES
SORT = 'ST_based'; % GS_based

load('K:\CSD_dynamic_analysis\DATA\Output\Output_Input_Control_7post_n7_Data_Threshold_25_Zscore_0_binned_1_mirror_1.mat');
BF =  Data.BF_Pos;

Control_M= nan(length(Data.names),length(Data.GS_based),3);

CNorm = [];

for i1 = 1:3
    CNorm(:,i1,1) = Data.(SORT)(i1).Early_RMS_RELRES(:,BF); 
    CNorm(:,i1,2) = Data.(SORT)(i1).Early_RMS_RELRES(:,BF+1);
    CNorm(:,i1,3) = Data.(SORT)(i1).Early_RMS_RELRES(:,BF+2);
end

for i1 = 1:length(Data.names)
    Norm = nanmean(CNorm(i1,:,:));
    for i2 = 1:length(Data.GS_based)
        Control_M(i1,i2,1) = Data.(SORT)(i2).Early_RMS_RELRES(i1,BF)/Norm(1);
        Control_M(i1,i2,2) = Data.(SORT)(i2).Early_RMS_RELRES(i1,BF+1)/Norm(2);
        Control_M(i1,i2,3) = Data.(SORT)(i2).Early_RMS_RELRES(i1,BF+2)/Norm(3);
    end
end

C_SEM = nanstd(Control_M)./sqrt(sum(~isnan(Control_M)));

load('K:\CSD_dynamic_analysis\DATA\Output\Output_Input_HighP_7post_Data_Threshold_25_Zscore_0_binned_1_mirror_1.mat');

BF =  Data.BF_Pos;

Laser_M= nan(length(Data.names),length(Data.GS_based),3); 

lNorm = [];

for i1 = 1:3
    lNorm(:,i1,1) = Data.(SORT)(i1).Early_RMS_RELRES(:,BF); 
    lNorm(:,i1,2) = Data.(SORT)(i1).Early_RMS_RELRES(:,BF+1);
    lNorm(:,i1,3) = Data.(SORT)(i1).Early_RMS_RELRES(:,BF+2);
end

for i1 = 1:length(Data.names)
    Norm = nanmean(lNorm(i1,:,:));
    for i2 = 1:length(Data.GS_based)
        Laser_M(i1,i2,1) = Data.(SORT)(i2).Early_RMS_RELRES(i1,BF)/Norm(1);
        Laser_M(i1,i2,2) = Data.(SORT)(i2).Early_RMS_RELRES(i1,BF+1)/Norm(2);
        Laser_M(i1,i2,3) = Data.(SORT)(i2).Early_RMS_RELRES(i1,BF+2)/Norm(3);
    end
end

All =vertcat(CNorm,lNorm);
All2 = nanmean(All);
All3 = All./All2;%normalized
AllMean = nanmean(nanmean(All3,2),1);
AllSTD_1 = mean(std(All3(:,:,1)'))./sqrt(19);
AllSTD_2 = mean(std(All3(:,:,2)'))./sqrt(19);
AllSTD_3 = mean(std(All3(:,:,3)'))./sqrt(19);
L_SEM = nanstd(Laser_M)./sqrt(sum(~isnan(Laser_M)));
STD_prime =vertcat(Control_M(:,1:3,:),Laser_M(:,1:3,:));


X1 = Laser_M(:,:,1); 
X2 = nan(size(X1));
X2(1:size(Control_M(:,:,1),1),:)= Control_M(:,:,1);
PRE_L = mean(X1(:,1:3),2);
PRE_C = mean(X2(:,1:3),2);
pVal_L =[];
pVal_C =[];
pVal_Both = [];

for i1 =1:11
    
    %Laser with Pre
%     X3 =[PRE_L,X1(:,i1)];
    %L = [1,2];
    %O = teg_repeated_measures_ANOVA(X3,L,{'Time','Measurement'})
%     S=anovan([1:12]',{PRE_L,X1(:,i1)});
%     S=anovan(X3,{'Pre','M'});
    [h,p] = ttest2(PRE_L,X1(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_L =[pVal_L Y];
    
    %Control with Pre
%     X3 =[PRE_C,X2(:,i1)];
    [h,p] = ttest2(PRE_C,X2(:,i1));
%     S=anovan(X3,{'Pre','M'});
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_C =[pVal_C Y];
    
    % Both
%     X3 =[X1(:,i1),X2(:,i1)];   
%     S=anovan(X3,{'Opto','Control'});
    [h,p] = ttest2(X1(:,i1),X2(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_Both =[pVal_Both Y]; 
end
    
stats_BF.pVal_L = pVal_L; stats_BF.pVal_C = pVal_C; stats_BF.pVal_Both = pVal_Both; 

X1 = Laser_M(:,:,2); 
X2 = nan(size(X1));
X2(1:size(Control_M(:,:,1),1),:)= Control_M(:,:,2);
PRE_L = mean(X1(:,1:3),2);
PRE_C = mean(X2(:,1:3),2);
pVal_L =[];
pVal_C =[];
pVal_Both = [];

for i1 =1:11
    
    %Laser with Pre
%     X3 =[PRE_L,X1(:,i1)];
    %L = [1,2];
    %O = teg_repeated_measures_ANOVA(X3,L,{'Time','Measurement'})
%     S=anovan([1:12]',{PRE_L,X1(:,i1)});
%     S=anovan(X3,{'Pre','M'});
    [h,p] = ttest2(PRE_L,X1(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_L =[pVal_L Y];
    
    %Control with Pre
%     X3 =[PRE_C,X2(:,i1)];
    [h,p] = ttest2(PRE_C,X2(:,i1));
%     S=anovan(X3,{'Pre','M'});
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_C =[pVal_C Y];
    
    % Both
%     X3 =[X1(:,i1),X2(:,i1)];   
%     S=anovan(X3,{'Opto','Control'});
    [h,p] = ttest2(X1(:,i1),X2(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_Both =[pVal_Both Y]; 
end
    
stats_NearBF.pVal_L = pVal_L; stats_NearBF.pVal_C = pVal_C; stats_NearBF.pVal_Both = pVal_Both; 


X1 = Laser_M(:,:,3); 
X2 = nan(size(X1));
X2(1:size(Control_M(:,:,1),1),:)= Control_M(:,:,3);
PRE_L = mean(X1(:,1:3),2);
PRE_C = mean(X2(:,1:3),2);
pVal_L =[];
pVal_C =[];
pVal_Both = [];

for i1 =1:11
    
    %Laser with Pre
%     X3 =[PRE_L,X1(:,i1)];
    %L = [1,2];
    %O = teg_repeated_measures_ANOVA(X3,L,{'Time','Measurement'})
%     S=anovan([1:12]',{PRE_L,X1(:,i1)});
%     S=anovan(X3,{'Pre','M'});
    [h,p] = ttest2(PRE_L,X1(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_L =[pVal_L Y];
    
    %Control with Pre
%     X3 =[PRE_C,X2(:,i1)];
    [h,p] = ttest2(PRE_C,X2(:,i1));
%     S=anovan(X3,{'Pre','M'});
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_C =[pVal_C Y];
    
    % Both
%     X3 =[X1(:,i1),X2(:,i1)];   
%     S=anovan(X3,{'Opto','Control'});
    [h,p] = ttest2(X1(:,i1),X2(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_Both =[pVal_Both Y]; 
end
    
stats_NonBF.pVal_L = pVal_L; stats_NonBF.pVal_C = pVal_C; stats_NonBF.pVal_Both = pVal_Both; 
% %%keyboard%early relsres

clear V
subplot(2,3,4)
S1 =Laser_M(:,:,1);
S2 = nan(size(S1));
S2(1:7,:) = Control_M(:,:,1);
% [stat,pperm,hperm,statperm]= elbe_permtest_005(S1,S2,1,'raw',1000,0,0.05,0);
% SigPos=find(pperm <= 0.025 | pperm >= 1.025);
L = [2 11];
S3 = [S1, S2];
O = teg_repeated_measures_ANOVA(S3,L,{'Group','Measurement'});
[p,n]=ttest(S1,S2);
try
V =reshape(O.groups,2,[])'; Vlog =(V <= 11);
Vlog =Vlog(:,1)+Vlog(:,2);
Vlog = Vlog ~= 2;V = V.*Vlog; V(V == 0) = nan;V(:,1) =V(:,1)-11;
V =V(V(:,1) == V(:,2));
end
shadedErrorBar([],nanmean(Control_M(:,:,1),1),C_SEM(:,:,1),'lineprops',{'b','markerfacecolor','b'})
hold on
shadedErrorBar([],nanmean(Laser_M(:,:,1),1),L_SEM(:,:,1),'lineprops',{'r','markerfacecolor','r'})
prime_M = mean(nanmean(STD_prime(:,:,1)'));
prime_STD = mean(nanstd(STD_prime(:,:,1)'));
plot(1:length(p),p*1.4,'k*')
try
 plot(V,ones(size(V))*1.5,'r*')
end
T1 = [O.labels{1} ' F(' num2str(O.R(1,2)) ','  num2str(O.R(1,3)) ') = ' num2str(O.R(1,1)) ' p = ' num2str(O.R(1,4))];
T2 = [O.labels{2} ' F(' num2str(O.R(2,2)) ','  num2str(O.R(2,3)) ') = ' num2str(O.R(2,1)) ' p = ' num2str(O.R(2,4))];
T3 = [O.labels{3} ' F(' num2str(O.R(3,2)) ','  num2str(O.R(3,3)) ') = ' num2str(O.R(3,1)) ' p = ' num2str(O.R(3,4))];
T =strvcat(T1, T2, T3);
text(1.2,0.9, T);

% plot([SigPos],[ones(size(SigPos))*1.4],'k*')
% prime_STD = std(nanmean(STD_prime(:,:,1)))./sqrt(sum(~isnan(STD_prime(:,:,1))));
% prime_STD = mean(prime_STD);
% plot([1 length(Data.GS_based)],[prime_M+3*AllSTD_1 prime_M+3*AllSTD_1],':k')
% plot([1 length(Data.GS_based)],[prime_M-3*AllSTD_1 prime_M-3*AllSTD_1],':k')
% y1 = ones(1,11)*1.5; y2 = ones(1,11)*1.7;y3 = ones(1,11)*1.9;
% 
% plot([1:11].*(~isnan(stats_BF.pVal_L)),y1.*(~isnan(stats_BF.pVal_L)),'r*',...
%     [1:11].*(~isnan(stats_BF.pVal_C)),y2.*(~isnan(stats_BF.pVal_C)),'b*',...
%     [1:11].*(~isnan(stats_BF.pVal_Both)).*(~isnan(stats_BF.pVal_Both)),y3,'k*')
ylim([0.8 1.5])
xlim ([1 length(Data.GS_based)])
xticks([1 2 3 4 5 6 7 8 9 10 11])
xticklabels({'Pre(n-2)','Pre(n-1)', 'Pre(n)','Laser', 'Post 1','Post 2', 'Post 3','Post 4','Post 5','Post 6','Post 7'})
xtickangle(270)
set(gca,'fontsize',12,'FontWeight','bold')

title('Early RMS RELRESCSD BF')
legend ('Control','','','','Opto','','','','Location','best')

clear V
subplot(2,3,5)
S1 =Laser_M(:,:,2);
S2 = nan(size(S1));
S2(1:7,:) = Control_M(:,:,2);
L = [2 11];
S3 = [S1, S2];
O = teg_repeated_measures_ANOVA(S3,L,{'Group','Measurement'});
[p,n]=ttest(S1,S2);
try
V =reshape(O.groups,2,[])'; Vlog =(V <= 11);
Vlog =Vlog(:,1)+Vlog(:,2);
Vlog = Vlog ~= 2;V = V.*Vlog; V(V == 0) = nan;V(:,1) =V(:,1)-11;
V =V(V(:,1) == V(:,2));
end

% [stat,pperm,hperm,statperm]= elbe_permtest_005(S1,S2,1,'raw',1000,0,0.05,0);
% SigPos=find(pperm <= 0.025 | pperm >= 1.025);
shadedErrorBar([],nanmean(Control_M(:,:,2),1),C_SEM(:,:,2),'lineprops',{'b','markerfacecolor','b'})
hold on
shadedErrorBar([],nanmean(Laser_M(:,:,2),1),L_SEM(:,:,2),'lineprops',{'r','markerfacecolor','r'})
prime_M = mean(nanmean(STD_prime(:,:,2)'));
prime_STD = mean(nanstd(STD_prime(:,:,2)'));
plot(1:length(p),p*1.4,'k*')
try
 plot(V,ones(size(V))*1.5,'r*')
end
T1 = [O.labels{1} ' F(' num2str(O.R(1,2)) ','  num2str(O.R(1,3)) ') = ' num2str(O.R(1,1)) ' p = ' num2str(O.R(1,4))];
T2 = [O.labels{2} ' F(' num2str(O.R(2,2)) ','  num2str(O.R(2,3)) ') = ' num2str(O.R(2,1)) ' p = ' num2str(O.R(2,4))];
T3 = [O.labels{3} ' F(' num2str(O.R(3,2)) ','  num2str(O.R(3,3)) ') = ' num2str(O.R(3,1)) ' p = ' num2str(O.R(3,4))];
T =strvcat(T1, T2, T3);
text(1.2,0.9, T);

% plot([SigPos],[ones(size(SigPos))*1.4],'k*')
% prime_STD = std(nanmean(STD_prime(:,:,2)))./sqrt(sum(~isnan(STD_prime(:,:,2))));
% prime_STD = mean(prime_STD);
% plot([1 length(Data.GS_based)],[prime_M+3*AllSTD_2 prime_M+3*AllSTD_2],':k')
% plot([1 length(Data.GS_based)],[prime_M-3*AllSTD_2 prime_M-3*AllSTD_2],':k')
% y1 = ones(1,11)*1.5; y2 = ones(1,11)*1.7;y3 = ones(1,11)*1.9;
% 
% plot([1:11].*(~isnan(stats_NearBF.pVal_L)),y1.*(~isnan(stats_NearBF.pVal_L)),'r*',...
%     [1:11].*(~isnan(stats_NearBF.pVal_C)),y2.*(~isnan(stats_NearBF.pVal_C)),'b*',...
%     [1:11].*(~isnan(stats_NearBF.pVal_Both)).*(~isnan(stats_NearBF.pVal_Both)),y3,'k*')
ylim([0.8 1.5])
xlim ([1 length(Data.GS_based)])
xticks([1 2 3 4 5 6 7 8 9 10 11])
xticklabels({'Pre(n-2)','Pre(n-1)', 'Pre(n)','Laser', 'Post 1','Post 2', 'Post 3','Post 4','Post 5','Post 6','Post 7'})
xtickangle(270)
set(gca,'fontsize',12,'FontWeight','bold')
title('± near BF (1-2 Oct)')
% legend ('Control','','','','Opto','','','','Location','best')

clear V
subplot(2,3,6)
S1 =Laser_M(:,:,3);
S2 = nan(size(S1));
S2(1:7,:) = Control_M(:,:,3);
L = [2 11];
S3 = [S1, S2];
O = teg_repeated_measures_ANOVA(S3,L,{'Group','Measurement'});
[p,n]=ttest(S1,S2);
try
V =reshape(O.groups,2,[])'; Vlog =(V <= 11);
Vlog =Vlog(:,1)+Vlog(:,2);
Vlog = Vlog ~= 2;V = V.*Vlog; V(V == 0) = nan;V(:,1) =V(:,1)-11;
V =V(V(:,1) == V(:,2));
end
% [stat,pperm,hperm,statperm]= elbe_permtest_005(S1,S2,1,'raw',1000,0,0.05,0);
% SigPos=find(pperm <= 0.025 | pperm >= 1.025);
shadedErrorBar([],nanmean(Control_M(:,:,3),1),C_SEM(:,:,3),'lineprops',{'b','markerfacecolor','b'})
hold on
shadedErrorBar([],nanmean(Laser_M(:,:,3),1),L_SEM(:,:,3),'lineprops',{'r','markerfacecolor','r'})
prime_M = mean(nanmean(STD_prime(:,:,3)'));
prime_STD = mean(nanstd(STD_prime(:,:,3)'));
plot(1:length(p),p*1.4,'k*')
try
 plot(V,ones(size(V))*1.5,'r*')
end
T1 = [O.labels{1} ' F(' num2str(O.R(1,2)) ','  num2str(O.R(1,3)) ') = ' num2str(O.R(1,1)) ' p = ' num2str(O.R(1,4))];
T2 = [O.labels{2} ' F(' num2str(O.R(2,2)) ','  num2str(O.R(2,3)) ') = ' num2str(O.R(2,1)) ' p = ' num2str(O.R(2,4))];
T3 = [O.labels{3} ' F(' num2str(O.R(3,2)) ','  num2str(O.R(3,3)) ') = ' num2str(O.R(3,1)) ' p = ' num2str(O.R(3,4))];
T =strvcat(T1, T2, T3);
text(1.2,0.9, T);

% plot([SigPos],[ones(size(SigPos))*1.4],'k*')
% prime_STD = std(nanmean(STD_prime(:,:,3)))./sqrt(sum(~isnan(STD_prime(:,:,3))));
% prime_STD = mean(prime_STD);
% plot([1 length(Data.GS_based)],[prime_M+3*AllSTD_3 prime_M+3*AllSTD_3],':k')
% plot([1 length(Data.GS_based)],[prime_M-3*AllSTD_3 prime_M-3*AllSTD_3],':k')
% y1 = ones(1,11)*1.5; y2 = ones(1,11)*1.7;y3 = ones(1,11)*1.9;
% 
% plot([1:11].*(~isnan(stats_NonBF.pVal_L)),y1.*(~isnan(stats_NonBF.pVal_L)),'r*',...
%     [1:11].*(~isnan(stats_NonBF.pVal_C)),y2.*(~isnan(stats_NonBF.pVal_C)),'b*',...
%     [1:11].*(~isnan(stats_NonBF.pVal_Both)).*(~isnan(stats_NonBF.pVal_Both)),y3,'k*')
ylim([0.8 1.5])
xlim ([1 length(Data.GS_based)])
xticks([1 2 3 4 5 6 7 8 9 10 11])
xticklabels({'Pre(n-2)','Pre(n-1)', 'Pre(n)','Laser', 'Post 1','Post 2', 'Post 3','Post 4','Post 5','Post 6','Post 7'})
xtickangle(270)
set(gca,'fontsize',12,'FontWeight','bold')

title('± non BF (3-4 Oct)')
% legend ('Control','','','','Opto','','','','Location','best')
h=gcf;

set(h,'PaperOrientation','landscape');

set(h,'PaperPosition', [1 1 28 19]);

SORT = 'ST_based'; % GS_based

load('K:\CSD_dynamic_analysis\DATA\Output\Output_Input_Control_7post_n7_Data_Threshold_25_Zscore_0_binned_1_mirror_1.mat');
BF =  Data.BF_Pos;

Control_M= nan(length(Data.names),length(Data.GS_based),3);

CNorm = [];

for i1 = 1:3
    CNorm(:,i1,1) = Data.(SORT)(i1).Late_RMS_AVREC(:,BF); 
    CNorm(:,i1,2) = Data.(SORT)(i1).Late_RMS_AVREC(:,BF+1);
    CNorm(:,i1,3) = Data.(SORT)(i1).Late_RMS_AVREC(:,BF+2);
end

for i1 = 1:length(Data.names)
    Norm = nanmean(CNorm(i1,:,:));
    for i2 = 1:length(Data.GS_based)
        Control_M(i1,i2,1) = Data.(SORT)(i2).Late_RMS_AVREC(i1,BF)/Norm(1);
        Control_M(i1,i2,2) = Data.(SORT)(i2).Late_RMS_AVREC(i1,BF+1)/Norm(2);
        Control_M(i1,i2,3) = Data.(SORT)(i2).Late_RMS_AVREC(i1,BF+2)/Norm(3);
    end
end

C_SEM = nanstd(Control_M)./sqrt(sum(~isnan(Control_M)));

load('K:\CSD_dynamic_analysis\DATA\Output\Output_Input_HighP_7post_Data_Threshold_25_Zscore_0_binned_1_mirror_1.mat');

BF =  Data.BF_Pos;

Laser_M= nan(length(Data.names),length(Data.GS_based),3); 

lNorm = [];

for i1 = 1:3
    lNorm(:,i1,1) = Data.(SORT)(i1).Late_RMS_AVREC(:,BF); 
    lNorm(:,i1,2) = Data.(SORT)(i1).Late_RMS_AVREC(:,BF+1);
    lNorm(:,i1,3) = Data.(SORT)(i1).Late_RMS_AVREC(:,BF+2);
end

for i1 = 1:length(Data.names)
    Norm = nanmean(lNorm(i1,:,:));
    for i2 = 1:length(Data.GS_based)
        Laser_M(i1,i2,1) = Data.(SORT)(i2).Late_RMS_AVREC(i1,BF)/Norm(1);
        Laser_M(i1,i2,2) = Data.(SORT)(i2).Late_RMS_AVREC(i1,BF+1)/Norm(2);
        Laser_M(i1,i2,3) = Data.(SORT)(i2).Late_RMS_AVREC(i1,BF+2)/Norm(3);
    end
end

All =vertcat(CNorm,lNorm);
All2 = nanmean(All);
All3 = All./All2;%normalized
AllMean = nanmean(nanmean(All3,2),1);
AllSTD_1 = mean(std(All3(:,:,1)'))./sqrt(19);
AllSTD_2 = mean(std(All3(:,:,2)'))./sqrt(19);
AllSTD_3 = mean(std(All3(:,:,3)'))./sqrt(19);
% %%%keyboard
% ALLSTD = nanmean(nanstd(common'))

L_SEM = nanstd(Laser_M)./sqrt(sum(~isnan(Laser_M)));
STD_prime =vertcat(Control_M(:,1:3,:),Laser_M(:,1:3,:));
% %%%keyboard

X1 = Laser_M(:,:,1); 
X2 = nan(size(X1));
X2(1:size(Control_M(:,:,1),1),:)= Control_M(:,:,1);
PRE_L = mean(X1(:,1:3),2);
PRE_C = mean(X2(:,1:3),2);
pVal_L =[];
pVal_C =[];
pVal_Both = [];

for i1 =1:11
    
    %Laser with Pre
%      X3 =[PRE_L,X1(:,i1)];
%     L = [1,2];
%     O = teg_repeated_measures_ANOVA(X3,L,{'Time','Measurement'})
%      S=anovan([1:12]',{PRE_L,X1(:,i1)});
%     S=anovan(X3,{'Pre','M'});
    [h,p] = ttest2(PRE_L,X1(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_L =[pVal_L Y];
    
    %Control with Pre
%     X3 =[PRE_C,X2(:,i1)];
    [h,p] = ttest2(PRE_C,X2(:,i1));
%     S=anovan(X3,{'Pre','M'});
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_C =[pVal_C Y];
    
    % Both
%     X3 =[X1(:,i1),X2(:,i1)];   
%     S=anovan(X3,{'Opto','Control'});
    [h,p] = ttest2(X1(:,i1),X2(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_Both =[pVal_Both Y]; 
end
 
stats_BF.pVal_L = pVal_L; stats_BF.pVal_C = pVal_C; stats_BF.pVal_Both = pVal_Both; 

X1 = Laser_M(:,:,2); 
X2 = nan(size(X1));
X2(1:size(Control_M(:,:,1),1),:)= Control_M(:,:,2);
PRE_L = mean(X1(:,1:3),2);
PRE_C = mean(X2(:,1:3),2);
pVal_L =[];
pVal_C =[];
pVal_Both = [];

for i1 =1:11
    
    %Laser with Pre
%     X3 =[PRE_L,X1(:,i1)];
    %L = [1,2];
    %O = teg_repeated_measures_ANOVA(X3,L,{'Time','Measurement'})
%     S=anovan([1:12]',{PRE_L,X1(:,i1)});
%     S=anovan(X3,{'Pre','M'});
    [h,p] = ttest2(PRE_L,X1(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_L =[pVal_L Y];
    
    %Control with Pre
%     X3 =[PRE_C,X2(:,i1)];
    [h,p] = ttest2(PRE_C,X2(:,i1));
%     S=anovan(X3,{'Pre','M'});
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_C =[pVal_C Y];
    
    % Both
%     X3 =[X1(:,i1),X2(:,i1)];   
%     S=anovan(X3,{'Opto','Control'});
    [h,p] = ttest2(X1(:,i1),X2(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_Both =[pVal_Both Y]; 
end
    
stats_NearBF.pVal_L = pVal_L; stats_NearBF.pVal_C = pVal_C; stats_NearBF.pVal_Both = pVal_Both; 


X1 = Laser_M(:,:,3); 
X2 = nan(size(X1));
X2(1:size(Control_M(:,:,1),1),:)= Control_M(:,:,3);
PRE_L = mean(X1(:,1:3),2);
PRE_C = mean(X2(:,1:3),2);
pVal_L =[];
pVal_C =[];
pVal_Both = [];

for i1 =1:11
    
    %Laser with Pre
%     X3 =[PRE_L,X1(:,i1)];
    %L = [1,2];
    %O = teg_repeated_measures_ANOVA(X3,L,{'Time','Measurement'})
%     S=anovan([1:12]',{PRE_L,X1(:,i1)});
%     S=anovan(X3,{'Pre','M'});
    [h,p] = ttest2(PRE_L,X1(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_L =[pVal_L Y];
    
    %Control with Pre
%     X3 =[PRE_C,X2(:,i1)];
    [h,p] = ttest2(PRE_C,X2(:,i1));
%     S=anovan(X3,{'Pre','M'});
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_C =[pVal_C Y];
    
    % Both
%     X3 =[X1(:,i1),X2(:,i1)];   
%     S=anovan(X3,{'Opto','Control'});
    [h,p] = ttest2(X1(:,i1),X2(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_Both =[pVal_Both Y]; 
end
    
stats_NonBF.pVal_L = pVal_L; stats_NonBF.pVal_C = pVal_C; stats_NonBF.pVal_Both = pVal_Both; 

% keyboard % late
figure
subplot(2,3,1)
S1 =Laser_M(:,:,1);
S2 = nan(size(S1));
S2(1:7,:) = Control_M(:,:,1);
L = [2 11];
S3 = [S1, S2];
O = teg_repeated_measures_ANOVA(S3,L,{'Group','Measurement'});
[p,n]=ttest(S1,S2);
try
V =reshape(O.groups,2,[])'; Vlog =(V <= 11);
Vlog =Vlog(:,1)+Vlog(:,2);
Vlog = Vlog ~= 2;V = V.*Vlog; V(V == 0) = nan;V(:,1) =V(:,1)-11;
V =V(V(:,1) == V(:,2));
end
% [stat,pperm,hperm,statperm]= elbe_permtest_005(S1,S2,1,'raw',1000,0,0.05,0);
% SigPos=find(pperm <= 0.025 | pperm >= 1.025);
shadedErrorBar([],nanmean(Control_M(:,:,1),1),C_SEM(:,:,1),'lineprops',{'b','markerfacecolor','b'})
hold on
shadedErrorBar([],nanmean(Laser_M(:,:,1),1),L_SEM(:,:,1),'lineprops',{'r','markerfacecolor','r'})
prime_M = mean(nanmean(STD_prime(:,:,1))');
prime_STD = mean(nanstd(STD_prime(:,:,1)'));
plot(1:length(p),p*1.15,'k*')
try
 plot(V,ones(size(V))*1.25,'r*')
end
T1 = [O.labels{1} ' F(' num2str(O.R(1,2)) ','  num2str(O.R(1,3)) ') = ' num2str(O.R(1,1)) ' p = ' num2str(O.R(1,4))];
T2 = [O.labels{2} ' F(' num2str(O.R(2,2)) ','  num2str(O.R(2,3)) ') = ' num2str(O.R(2,1)) ' p = ' num2str(O.R(2,4))];
T3 = [O.labels{3} ' F(' num2str(O.R(3,2)) ','  num2str(O.R(3,3)) ') = ' num2str(O.R(3,1)) ' p = ' num2str(O.R(3,4))];
T =strvcat(T1, T2, T3);
text(1.2,0.9, T);
% plot([SigPos],[ones(size(SigPos))*1.15],'k*')
% prime_STD = std(nanmean(STD_prime(:,:,1)))./sqrt(sum(~isnan(STD_prime(:,:,1))));
% prime_STD = mean(prime_STD);
% plot([1 length(Data.GS_based)],[prime_M+3*AllSTD_1 prime_M+3*AllSTD_1 ],':k')
% plot([1 length(Data.GS_based)],[prime_M-3*AllSTD_1 prime_M-3*AllSTD_1 ],':k')
% y1 = ones(1,11)*1.15; y2 = ones(1,11)*1.17;y3 = ones(1,11)*1.19;
% 
% plot([1:11].*(~isnan(stats_BF.pVal_L)),y1.*(~isnan(stats_BF.pVal_L)),'y*',...
%     [1:11].*(~isnan(stats_BF.pVal_C)),y2.*(~isnan(stats_BF.pVal_C)),'b*',...
%     [1:11].*(~isnan(stats_BF.pVal_Both)).*(~isnan(stats_BF.pVal_Both)),y3,'k*')
ylim([0.8 1.3])
xlim ([1 length(Data.GS_based)])
xticks([1 2 3 4 5 6 7 8 9 10 11])
xticklabels({'Pre(n-2)','Pre(n-1)', 'Pre(n)','Laser', 'Post 1','Post 2', 'Post 3','Post 4','Post 5','Post 6','Post 7'})
xtickangle(270)
set(gca,'fontsize',12,'FontWeight','bold')

title('Late RMS AVREC BF')
legend ('Control','','','','Opto','','','','Location','best')


subplot(2,3,2)
S1 =Laser_M(:,:,2);
S2 = nan(size(S1));
S2(1:7,:) = Control_M(:,:,2);
L = [2 11];
S3 = [S1, S2];
O = teg_repeated_measures_ANOVA(S3,L,{'Group','Measurement'});
[p,n]=ttest(S1,S2);
try
V =reshape(O.groups,2,[])'; Vlog =(V <= 11);
Vlog =Vlog(:,1)+Vlog(:,2);
Vlog = Vlog ~= 2;V = V.*Vlog; V(V == 0) = nan;V(:,1) =V(:,1)-11;
V =V(V(:,1) == V(:,2));
end
% [stat,pperm,hperm,statperm]= elbe_permtest_005(S1,S2,1,'raw',1000,0,0.05,0);
% SigPos=find(pperm <= 0.025 | pperm >= 1.025);
shadedErrorBar([],nanmean(Control_M(:,:,2),1),C_SEM(:,:,2),'lineprops',{'b','markerfacecolor','b'})
hold on
shadedErrorBar([],nanmean(Laser_M(:,:,2),1),L_SEM(:,:,2),'lineprops',{'r','markerfacecolor','r'})
prime_M = mean(nanmean(STD_prime(:,:,2))');
prime_STD = mean(nanstd(STD_prime(:,:,2)'));
plot(1:length(p),p*1.15,'k*')
try
 plot(V,ones(size(V))*1.25,'r*')
end
T1 = [O.labels{1} ' F(' num2str(O.R(1,2)) ','  num2str(O.R(1,3)) ') = ' num2str(O.R(1,1)) ' p = ' num2str(O.R(1,4))];
T2 = [O.labels{2} ' F(' num2str(O.R(2,2)) ','  num2str(O.R(2,3)) ') = ' num2str(O.R(2,1)) ' p = ' num2str(O.R(2,4))];
T3 = [O.labels{3} ' F(' num2str(O.R(3,2)) ','  num2str(O.R(3,3)) ') = ' num2str(O.R(3,1)) ' p = ' num2str(O.R(3,4))];
T =strvcat(T1, T2, T3);
text(1.2,0.9, T);
% plot([SigPos],[ones(size(SigPos))*1.15],'k*')
% prime_STD = std(nanmean(STD_prime(:,:,2)))./sqrt(sum(~isnan(STD_prime(:,:,2))));
% prime_STD = mean(prime_STD);
% plot([1 length(Data.GS_based)],[prime_M+3*AllSTD_2 prime_M+3*AllSTD_2 ],':k')
% plot([1 length(Data.GS_based)],[prime_M-3*AllSTD_2 prime_M-3*AllSTD_2 ],':k')
% y1 = ones(1,11)*1.15; y2 = ones(1,11)*1.17;y3 = ones(1,11)*1.19;
% 
% plot([1:11].*(~isnan(stats_NearBF.pVal_L)),y1.*(~isnan(stats_NearBF.pVal_L)),'y*',...
%     [1:11].*(~isnan(stats_NearBF.pVal_C)),y2.*(~isnan(stats_NearBF.pVal_C)),'b*',...
%     [1:11].*(~isnan(stats_NearBF.pVal_Both)).*(~isnan(stats_NearBF.pVal_Both)),y3,'k*')
ylim([0.8 1.3])
xlim ([1 length(Data.GS_based)])
xticks([1 2 3 4 5 6 7 8 9 10 11])
xticklabels({'Pre(n-2)','Pre(n-1)', 'Pre(n)','Laser', 'Post 1','Post 2', 'Post 3','Post 4','Post 5','Post 6','Post 7'})
xtickangle(270)
set(gca,'fontsize',12,'FontWeight','bold')

title('± near BF (1-2 Oct)')
% legend ('Control','','','','Opto','','','','Location','best')

subplot(2,3,3)
S1 =Laser_M(:,:,3);
S2 = nan(size(S1));
S2(1:7,:) = Control_M(:,:,3);
S3 = [S1, S2];
O = teg_repeated_measures_ANOVA(S3,L,{'Group','Measurement'});
[p,n]=ttest(S1,S2);
try
V =reshape(O.groups,2,[])'; Vlog =(V <= 11);
Vlog =Vlog(:,1)+Vlog(:,2);
Vlog = Vlog ~= 2;V = V.*Vlog; V(V == 0) = nan;V(:,1) =V(:,1)-11;
V =V(V(:,1) == V(:,2));
end
% [stat,pperm,hperm,statperm]= elbe_permtest_005(S1,S2,1,'raw',1000,0,0.05,0);
% SigPos=find(pperm <= 0.025 | pperm >= 1.025);
shadedErrorBar([],nanmean(Control_M(:,:,3),1),C_SEM(:,:,3),'lineprops',{'b','markerfacecolor','b'})
hold on
shadedErrorBar([],nanmean(Laser_M(:,:,3),1),L_SEM(:,:,3),'lineprops',{'r','markerfacecolor','r'})
prime_M = mean(nanmean(STD_prime(:,:,3)'));
prime_STD = mean(nanstd(STD_prime(:,:,3)'));
plot(1:length(p),p*1.15,'k*')
try
 plot(V,ones(size(V))*1.25,'r*')
end
T1 = [O.labels{1} ' F(' num2str(O.R(1,2)) ','  num2str(O.R(1,3)) ') = ' num2str(O.R(1,1)) ' p = ' num2str(O.R(1,4))];
T2 = [O.labels{2} ' F(' num2str(O.R(2,2)) ','  num2str(O.R(2,3)) ') = ' num2str(O.R(2,1)) ' p = ' num2str(O.R(2,4))];
T3 = [O.labels{3} ' F(' num2str(O.R(3,2)) ','  num2str(O.R(3,3)) ') = ' num2str(O.R(3,1)) ' p = ' num2str(O.R(3,4))];
T =strvcat(T1, T2, T3);
text(1.2,0.9, T);
% plot([SigPos],[ones(size(SigPos))*1.15],'k*')
% prime_STD = std(nanmean(STD_prime(:,:,3)))./sqrt(sum(~isnan(STD_prime(:,:,3))));
% prime_STD = mean(prime_STD);
% plot([1 length(Data.GS_based)],[prime_M+3*AllSTD_3 prime_M+3*AllSTD_3 ],':k')
% plot([1 length(Data.GS_based)],[prime_M-3*AllSTD_3 prime_M-3*AllSTD_3 ],':k')
% y1 = ones(1,11)*1.15; y2 = ones(1,11)*1.17;y3 = ones(1,11)*1.19;
% 
% plot([1:11].*(~isnan(stats_NonBF.pVal_L)),y1.*(~isnan(stats_NonBF.pVal_L)),'y*',...
%     [1:11].*(~isnan(stats_NonBF.pVal_C)),y2.*(~isnan(stats_NonBF.pVal_C)),'b*',...
%     [1:11].*(~isnan(stats_NonBF.pVal_Both)).*(~isnan(stats_NonBF.pVal_Both)),y3,'k*')
ylim([0.8 1.3])
xlim ([1 length(Data.GS_based)])
xticks([1 2 3 4 5 6 7 8 9 10 11])
xticklabels({'Pre(n-2)','Pre(n-1)', 'Pre(n)','Laser', 'Post 1','Post 2', 'Post 3','Post 4','Post 5','Post 6','Post 7'})
xtickangle(270)
set(gca,'fontsize',12,'FontWeight','bold')

title('± non BF (3-4 Oct)')
% legend ('Control','','','','Opto','','','','Location','best')
%  %%%keyboard


%%% RELRES
SORT = 'ST_based'; % GS_based

load('K:\CSD_dynamic_analysis\DATA\Output\Output_Input_Control_7post_n7_Data_Threshold_25_Zscore_0_binned_1_mirror_1.mat');
BF =  Data.BF_Pos;

Control_M= nan(length(Data.names),length(Data.GS_based),3);

CNorm = [];

for i1 = 1:3
    CNorm(:,i1,1) = Data.(SORT)(i1).Late_RMS_RELRES(:,BF); 
    CNorm(:,i1,2) = Data.(SORT)(i1).Late_RMS_RELRES(:,BF+1);
    CNorm(:,i1,3) = Data.(SORT)(i1).Late_RMS_RELRES(:,BF+2);
end

for i1 = 1:length(Data.names)
    Norm = nanmean(CNorm(i1,:,:));
    for i2 = 1:length(Data.GS_based)
        Control_M(i1,i2,1) = Data.(SORT)(i2).Late_RMS_RELRES(i1,BF)/Norm(1);
        Control_M(i1,i2,2) = Data.(SORT)(i2).Late_RMS_RELRES(i1,BF+1)/Norm(2);
        Control_M(i1,i2,3) = Data.(SORT)(i2).Late_RMS_RELRES(i1,BF+2)/Norm(3);
    end
end

C_SEM = nanstd(Control_M)./sqrt(sum(~isnan(Control_M)));

load('K:\CSD_dynamic_analysis\DATA\Output\Output_Input_HighP_7post_Data_Threshold_25_Zscore_0_binned_1_mirror_1.mat');

BF =  Data.BF_Pos;

Laser_M= nan(length(Data.names),length(Data.GS_based),3); 

lNorm = [];

for i1 = 1:3
    lNorm(:,i1,1) = Data.(SORT)(i1).Late_RMS_RELRES(:,BF); 
    lNorm(:,i1,2) = Data.(SORT)(i1).Late_RMS_RELRES(:,BF+1);
    lNorm(:,i1,3) = Data.(SORT)(i1).Late_RMS_RELRES(:,BF+2);
end

for i1 = 1:length(Data.names)
    Norm = nanmean(lNorm(i1,:,:));
    for i2 = 1:length(Data.GS_based)
        Laser_M(i1,i2,1) = Data.(SORT)(i2).Late_RMS_RELRES(i1,BF)/Norm(1);
        Laser_M(i1,i2,2) = Data.(SORT)(i2).Late_RMS_RELRES(i1,BF+1)/Norm(2);
        Laser_M(i1,i2,3) = Data.(SORT)(i2).Late_RMS_RELRES(i1,BF+2)/Norm(3);
    end
end

All =vertcat(CNorm,lNorm);
All2 = nanmean(All);
All3 = All./All2;%normalized
AllMean = nanmean(nanmean(All3,2),1);
AllSTD_1 = mean(std(All3(:,:,1)'))./sqrt(19);
AllSTD_2 = mean(std(All3(:,:,2)'))./sqrt(19);
AllSTD_3 = mean(std(All3(:,:,3)'))./sqrt(19);
L_SEM = nanstd(Laser_M)./sqrt(sum(~isnan(Laser_M)));
STD_prime =vertcat(Control_M(:,1:3,:),Laser_M(:,1:3,:));


X1 = Laser_M(:,:,1); 
X2 = nan(size(X1));
X2(1:size(Control_M(:,:,1),1),:)= Control_M(:,:,1);
PRE_L = mean(X1(:,1:3),2);
PRE_C = mean(X2(:,1:3),2);
pVal_L =[];
pVal_C =[];
pVal_Both = [];

for i1 =1:11
    
    %Laser with Pre
%     X3 =[PRE_L,X1(:,i1)];
    %L = [1,2];
    %O = teg_repeated_measures_ANOVA(X3,L,{'Time','Measurement'})
%     S=anovan([1:12]',{PRE_L,X1(:,i1)});
%     S=anovan(X3,{'Pre','M'});
    [h,p] = ttest2(PRE_L,X1(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_L =[pVal_L Y];
    
    %Control with Pre
%     X3 =[PRE_C,X2(:,i1)];
    [h,p] = ttest2(PRE_C,X2(:,i1));
%     S=anovan(X3,{'Pre','M'});
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_C =[pVal_C Y];
    
    % Both
%     X3 =[X1(:,i1),X2(:,i1)];   
%     S=anovan(X3,{'Opto','Control'});
    [h,p] = ttest2(X1(:,i1),X2(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_Both =[pVal_Both Y]; 
end
    
stats_BF.pVal_L = pVal_L; stats_BF.pVal_C = pVal_C; stats_BF.pVal_Both = pVal_Both; 

X1 = Laser_M(:,:,2); 
X2 = nan(size(X1));
X2(1:size(Control_M(:,:,1),1),:)= Control_M(:,:,2);
PRE_L = mean(X1(:,1:3),2);
PRE_C = mean(X2(:,1:3),2);
pVal_L =[];
pVal_C =[];
pVal_Both = [];

for i1 =1:11
    
    %Laser with Pre
%     X3 =[PRE_L,X1(:,i1)];
    %L = [1,2];
    %O = teg_repeated_measures_ANOVA(X3,L,{'Time','Measurement'})
%     S=anovan([1:12]',{PRE_L,X1(:,i1)});
%     S=anovan(X3,{'Pre','M'});
    [h,p] = ttest2(PRE_L,X1(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_L =[pVal_L Y];
    
    %Control with Pre
%     X3 =[PRE_C,X2(:,i1)];
    [h,p] = ttest2(PRE_C,X2(:,i1));
%     S=anovan(X3,{'Pre','M'});
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_C =[pVal_C Y];
    
    % Both
%     X3 =[X1(:,i1),X2(:,i1)];   
%     S=anovan(X3,{'Opto','Control'});
    [h,p] = ttest2(X1(:,i1),X2(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_Both =[pVal_Both Y]; 
end
    
stats_NearBF.pVal_L = pVal_L; stats_NearBF.pVal_C = pVal_C; stats_NearBF.pVal_Both = pVal_Both; 


X1 = Laser_M(:,:,3); 
X2 = nan(size(X1));
X2(1:size(Control_M(:,:,1),1),:)= Control_M(:,:,3);
PRE_L = mean(X1(:,1:3),2);
PRE_C = mean(X2(:,1:3),2);
pVal_L =[];
pVal_C =[];
pVal_Both = [];

for i1 =1:11
    
    %Laser with Pre
%     X3 =[PRE_L,X1(:,i1)];
    %L = [1,2];
    %O = teg_repeated_measures_ANOVA(X3,L,{'Time','Measurement'})
%     S=anovan([1:12]',{PRE_L,X1(:,i1)});
%     S=anovan(X3,{'Pre','M'});
    [h,p] = ttest2(PRE_L,X1(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_L =[pVal_L Y];
    
    %Control with Pre
%     X3 =[PRE_C,X2(:,i1)];
    [h,p] = ttest2(PRE_C,X2(:,i1));
%     S=anovan(X3,{'Pre','M'});
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_C =[pVal_C Y];
    
    % Both
%     X3 =[X1(:,i1),X2(:,i1)];   
%     S=anovan(X3,{'Opto','Control'});
    [h,p] = ttest2(X1(:,i1),X2(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_Both =[pVal_Both Y]; 
end
    
stats_NonBF.pVal_L = pVal_L; stats_NonBF.pVal_C = pVal_C; stats_NonBF.pVal_Both = pVal_Both; 
% keyboard
subplot(2,3,4)
S1 =Laser_M(:,:,1);
S2 = nan(size(S1));
S2(1:7,:) = Control_M(:,:,1);
S3 = [S1, S2];
O = teg_repeated_measures_ANOVA(S3,L,{'Group','Measurement'});
[p,n]=ttest(S1,S2);
try
V =reshape(O.groups,2,[])'; Vlog =(V <= 11);
Vlog =Vlog(:,1)+Vlog(:,2);
Vlog = Vlog ~= 2;V = V.*Vlog; V(V == 0) = nan;V(:,1) =V(:,1)-11;
V =V(V(:,1) == V(:,2));
end
% [stat,pperm,hperm,statperm]= elbe_permtest_005(S1,S2,1,'raw',1000,0,0.05,0);
% SigPos=find(pperm <= 0.025 | pperm >= 1.025);
shadedErrorBar([],nanmean(Control_M(:,:,1),1),C_SEM(:,:,1),'lineprops',{'b','markerfacecolor','b'})
hold on
shadedErrorBar([],nanmean(Laser_M(:,:,1),1),L_SEM(:,:,1),'lineprops',{'r','markerfacecolor','r'})
prime_M = mean(nanmean(STD_prime(:,:,1)'));
prime_STD = mean(nanstd(STD_prime(:,:,1)'));
plot(1:length(p),p*1.5,'k*')
try
 plot(V,ones(size(V))*1.6,'r*')
end
T1 = [O.labels{1} ' F(' num2str(O.R(1,2)) ','  num2str(O.R(1,3)) ') = ' num2str(O.R(1,1)) ' p = ' num2str(O.R(1,4))];
T2 = [O.labels{2} ' F(' num2str(O.R(2,2)) ','  num2str(O.R(2,3)) ') = ' num2str(O.R(2,1)) ' p = ' num2str(O.R(2,4))];
T3 = [O.labels{3} ' F(' num2str(O.R(3,2)) ','  num2str(O.R(3,3)) ') = ' num2str(O.R(3,1)) ' p = ' num2str(O.R(3,4))];
T =strvcat(T1, T2, T3);
text(1.2,0.9, T);
% plot([SigPos],[ones(size(SigPos))*1.5],'k*')
% prime_STD = std(nanmean(STD_prime(:,:,1)))./sqrt(sum(~isnan(STD_prime(:,:,1))));
% prime_STD = mean(prime_STD);
% plot([1 length(Data.GS_based)],[prime_M+3*AllSTD_1 prime_M+3*AllSTD_1],':k')
% plot([1 length(Data.GS_based)],[prime_M-3*AllSTD_1 prime_M-3*AllSTD_1],':k')
% y1 = ones(1,11)*1.5; y2 = ones(1,11)*1.7;y3 = ones(1,11)*1.9;
% 
% plot([1:11].*(~isnan(stats_BF.pVal_L)),y1.*(~isnan(stats_BF.pVal_L)),'r*',...
%     [1:11].*(~isnan(stats_BF.pVal_C)),y2.*(~isnan(stats_BF.pVal_C)),'b*',...
%     [1:11].*(~isnan(stats_BF.pVal_Both)).*(~isnan(stats_BF.pVal_Both)),y3,'k*')
ylim([0.8 2])
xlim ([1 length(Data.GS_based)])
xticks([1 2 3 4 5 6 7 8 9 10 11])
xticklabels({'Pre(n-2)','Pre(n-1)', 'Pre(n)','Laser', 'Post 1','Post 2', 'Post 3','Post 4','Post 5','Post 6','Post 7'})
xtickangle(270)
set(gca,'fontsize',12,'FontWeight','bold')

title('Late RMS RELRESCSD BF')
% legend ('Control','','','','Opto','','','','Location','best')


subplot(2,3,5)
S1 =Laser_M(:,:,2);
S2 = nan(size(S1));
S2(1:7,:) = Control_M(:,:,2);
S3 = [S1, S2];
O = teg_repeated_measures_ANOVA(S3,L,{'Group','Measurement'});
[p,n]=ttest(S1,S2);
try
V =reshape(O.groups,2,[])'; Vlog =(V <= 11);
Vlog =Vlog(:,1)+Vlog(:,2);
Vlog = Vlog ~= 2;V = V.*Vlog; V(V == 0) = nan;V(:,1) =V(:,1)-11;
V =V(V(:,1) == V(:,2));
end

% [stat,pperm,hperm,statperm]= elbe_permtest_005(S1,S2,1,'raw',1000,0,0.05,0);
% SigPos=find(pperm <= 0.025 | pperm >= 1.025);
shadedErrorBar([],nanmean(Control_M(:,:,2),1),C_SEM(:,:,2),'lineprops',{'b','markerfacecolor','b'})
hold on
shadedErrorBar([],nanmean(Laser_M(:,:,2),1),L_SEM(:,:,2),'lineprops',{'r','markerfacecolor','r'})
prime_M = mean(nanmean(STD_prime(:,:,2)'));
prime_STD = mean(nanstd(STD_prime(:,:,2)'));
plot(1:length(p),p*1.5,'k*')
try
 plot(V,ones(size(V))*1.6,'r*')
end
T1 = [O.labels{1} ' F(' num2str(O.R(1,2)) ','  num2str(O.R(1,3)) ') = ' num2str(O.R(1,1)) ' p = ' num2str(O.R(1,4))];
T2 = [O.labels{2} ' F(' num2str(O.R(2,2)) ','  num2str(O.R(2,3)) ') = ' num2str(O.R(2,1)) ' p = ' num2str(O.R(2,4))];
T3 = [O.labels{3} ' F(' num2str(O.R(3,2)) ','  num2str(O.R(3,3)) ') = ' num2str(O.R(3,1)) ' p = ' num2str(O.R(3,4))];
T =strvcat(T1, T2, T3);
text(1.2,0.9, T);
% plot([SigPos],[ones(size(SigPos))*1.5],'k*')
% prime_STD = std(nanmean(STD_prime(:,:,2)))./sqrt(sum(~isnan(STD_prime(:,:,2))));
% prime_STD = mean(prime_STD);
% plot([1 length(Data.GS_based)],[prime_M+3*AllSTD_2 prime_M+3*AllSTD_2],':k')
% plot([1 length(Data.GS_based)],[prime_M-3*AllSTD_2 prime_M-3*AllSTD_2],':k')
% y1 = ones(1,11)*1.5; y2 = ones(1,11)*1.7;y3 = ones(1,11)*1.9;
% 
% plot([1:11].*(~isnan(stats_NearBF.pVal_L)),y1.*(~isnan(stats_NearBF.pVal_L)),'r*',...
%     [1:11].*(~isnan(stats_NearBF.pVal_C)),y2.*(~isnan(stats_NearBF.pVal_C)),'b*',...
%     [1:11].*(~isnan(stats_NearBF.pVal_Both)).*(~isnan(stats_NearBF.pVal_Both)),y3,'k*')
ylim([0.8 2])
xlim ([1 length(Data.GS_based)])
xticks([1 2 3 4 5 6 7 8 9 10 11])
xticklabels({'Pre(n-2)','Pre(n-1)', 'Pre(n)','Laser', 'Post 1','Post 2', 'Post 3','Post 4','Post 5','Post 6','Post 7'})
xtickangle(270)
set(gca,'fontsize',12,'FontWeight','bold')
title('± near BF (1-2 Oct)')
legend ('Control','','','','Opto','','','','Location','best')

subplot(2,3,6)
S1 =Laser_M(:,:,3);
S2 = nan(size(S1));
S2(1:7,:) = Control_M(:,:,3);
S3 = [S1, S2];
O = teg_repeated_measures_ANOVA(S3,L,{'Group','Measurement'});
[p,n]=ttest(S1,S2);
try
V =reshape(O.groups,2,[])'; Vlog =(V <= 11);
Vlog =Vlog(:,1)+Vlog(:,2);
Vlog = Vlog ~= 2;V = V.*Vlog; V(V == 0) = nan;V(:,1) =V(:,1)-11;
V =V(V(:,1) == V(:,2));
end

% [stat,pperm,hperm,statperm]= elbe_permtest_005(S1,S2,1,'raw',1000,0,0.05,0);
% SigPos=find(pperm <= 0.025 | pperm >= 1.025);
shadedErrorBar([],nanmean(Control_M(:,:,3),1),C_SEM(:,:,3),'lineprops',{'b','markerfacecolor','b'})
hold on
shadedErrorBar([],nanmean(Laser_M(:,:,3),1),L_SEM(:,:,3),'lineprops',{'r','markerfacecolor','r'})
prime_M = mean(nanmean(STD_prime(:,:,3)'));
prime_STD = mean(nanstd(STD_prime(:,:,3)'));
plot(1:length(p),p*1.5,'k*')
try
 plot(V,ones(size(V))*1.6,'r*')
end
T1 = [O.labels{1} ' F(' num2str(O.R(1,2)) ','  num2str(O.R(1,3)) ') = ' num2str(O.R(1,1)) ' p = ' num2str(O.R(1,4))];
T2 = [O.labels{2} ' F(' num2str(O.R(2,2)) ','  num2str(O.R(2,3)) ') = ' num2str(O.R(2,1)) ' p = ' num2str(O.R(2,4))];
T3 = [O.labels{3} ' F(' num2str(O.R(3,2)) ','  num2str(O.R(3,3)) ') = ' num2str(O.R(3,1)) ' p = ' num2str(O.R(3,4))];
T =strvcat(T1, T2, T3);
text(1.2,0.9, T);
% plot([SigPos],[ones(size(SigPos))*1.5],'k*')
% prime_STD = std(nanmean(STD_prime(:,:,3)))./sqrt(sum(~isnan(STD_prime(:,:,3))));
% prime_STD = mean(prime_STD);
% plot([1 length(Data.GS_based)],[prime_M+3*AllSTD_3 prime_M+3*AllSTD_3],':k')
% plot([1 length(Data.GS_based)],[prime_M-3*AllSTD_3 prime_M-3*AllSTD_3],':k')
% y1 = ones(1,11)*1.5; y2 = ones(1,11)*1.7;y3 = ones(1,11)*1.9;
% 
% plot([1:11].*(~isnan(stats_NonBF.pVal_L)),y1.*(~isnan(stats_NonBF.pVal_L)),'r*',...
%     [1:11].*(~isnan(stats_NonBF.pVal_C)),y2.*(~isnan(stats_NonBF.pVal_C)),'b*',...
%     [1:11].*(~isnan(stats_NonBF.pVal_Both)).*(~isnan(stats_NonBF.pVal_Both)),y3,'k*')
ylim([0.8 2])
xlim ([1 length(Data.GS_based)])
xticks([1 2 3 4 5 6 7 8 9 10 11])
xticklabels({'Pre(n-2)','Pre(n-1)', 'Pre(n)','Laser', 'Post 1','Post 2', 'Post 3','Post 4','Post 5','Post 6','Post 7'})
xtickangle(270)
set(gca,'fontsize',12,'FontWeight','bold')

title('± non BF (3-4 Oct)')
% legend ('Control','','','','Opto','','','','Location','best')
h=gcf;

set(h,'PaperOrientation','landscape');

set(h,'PaperPosition', [1 1 28 19]);




























% %%keyboard

%%% RELRES
SORT = 'ST_based'; % GS_based

load('K:\CSD_dynamic_analysis\DATA\Output\Output_Input_Control_7post_n7_Data_Threshold_25_Zscore_0_binned_1_mirror_1.mat');
BF =  Data.BF_Pos;

Control_M= nan(length(Data.names),length(Data.GS_based),3);

CNorm = [];

for i1 = 1:3
    CNorm(:,i1,1) = Data.(SORT)(i1).Full_RMS_RELRES(:,BF); 
    CNorm(:,i1,2) = Data.(SORT)(i1).Full_RMS_RELRES(:,BF+1);
    CNorm(:,i1,3) = Data.(SORT)(i1).Full_RMS_RELRES(:,BF+2);
end

for i1 = 1:length(Data.names)
    Norm = nanmean(CNorm(i1,:,:));
    for i2 = 1:length(Data.GS_based)
        Control_M(i1,i2,1) = Data.(SORT)(i2).Full_RMS_RELRES(i1,BF)/Norm(1);
        Control_M(i1,i2,2) = Data.(SORT)(i2).Full_RMS_RELRES(i1,BF+1)/Norm(2);
        Control_M(i1,i2,3) = Data.(SORT)(i2).Full_RMS_RELRES(i1,BF+2)/Norm(3);
    end
end

C_SEM = nanstd(Control_M)./sqrt(sum(~isnan(Control_M)));

load('K:\CSD_dynamic_analysis\DATA\Output\Output_Input_HighP_7post_Data_Threshold_25_Zscore_0_binned_1_mirror_1.mat');

BF =  Data.BF_Pos;

Laser_M= nan(length(Data.names),length(Data.GS_based),3); 

lNorm = [];

for i1 = 1:3
    lNorm(:,i1,1) = Data.(SORT)(i1).Full_RMS_RELRES(:,BF); 
    lNorm(:,i1,2) = Data.(SORT)(i1).Full_RMS_RELRES(:,BF+1);
    lNorm(:,i1,3) = Data.(SORT)(i1).Full_RMS_RELRES(:,BF+2);
end

for i1 = 1:length(Data.names)
    Norm = nanmean(lNorm(i1,:,:));
    for i2 = 1:length(Data.GS_based)
        Laser_M(i1,i2,1) = Data.(SORT)(i2).Full_RMS_RELRES(i1,BF)/Norm(1);
        Laser_M(i1,i2,2) = Data.(SORT)(i2).Full_RMS_RELRES(i1,BF+1)/Norm(2);
        Laser_M(i1,i2,3) = Data.(SORT)(i2).Full_RMS_RELRES(i1,BF+2)/Norm(3);
    end
end
All =vertcat(CNorm,lNorm);
All2 = nanmean(All);
All3 = All./All2;%normalized
AllMean = nanmean(nanmean(All3,2),1);
AllSTD_1 = mean(std(All3(:,:,1)'))./sqrt(19);
AllSTD_2 = mean(std(All3(:,:,2)'))./sqrt(19);
AllSTD_3 = mean(std(All3(:,:,3)'))./sqrt(19);
L_SEM = nanstd(Laser_M)./sqrt(sum(~isnan(Laser_M)));
STD_prime =vertcat(Control_M(:,1:3,:),Laser_M(:,1:3,:));

X1 = Laser_M(:,:,1); 
X2 = nan(size(X1));
X2(1:size(Control_M(:,:,1),1),:)= Control_M(:,:,1);
PRE_L = mean(X1(:,1:3),2);
PRE_C = mean(X2(:,1:3),2);
pVal_L =[];
pVal_C =[];
pVal_Both = [];

for i1 =1:11
    
    %Laser with Pre
%     X3 =[PRE_L,X1(:,i1)];
    %L = [1,2];
    %O = teg_repeated_measures_ANOVA(X3,L,{'Time','Measurement'})
%     S=anovan([1:12]',{PRE_L,X1(:,i1)});
%     S=anovan(X3,{'Pre','M'});
    [h,p] = ttest2(PRE_L,X1(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_L =[pVal_L Y];
    
    %Control with Pre
%     X3 =[PRE_C,X2(:,i1)];
    [h,p] = ttest2(PRE_C,X2(:,i1));
%     S=anovan(X3,{'Pre','M'});
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_C =[pVal_C Y];
    
    % Both
%     X3 =[X1(:,i1),X2(:,i1)];   
%     S=anovan(X3,{'Opto','Control'});
    [h,p] = ttest2(X1(:,i1),X2(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_Both =[pVal_Both Y]; 
end
    
stats_BF.pVal_L = pVal_L; stats_BF.pVal_C = pVal_C; stats_BF.pVal_Both = pVal_Both; 

X1 = Laser_M(:,:,2); 
X2 = nan(size(X1));
X2(1:size(Control_M(:,:,1),1),:)= Control_M(:,:,2);
PRE_L = mean(X1(:,1:3),2);
PRE_C = mean(X2(:,1:3),2);
pVal_L =[];
pVal_C =[];
pVal_Both = [];

for i1 =1:11
    
    %Laser with Pre
%     X3 =[PRE_L,X1(:,i1)];
    %L = [1,2];
    %O = teg_repeated_measures_ANOVA(X3,L,{'Time','Measurement'})
%     S=anovan([1:12]',{PRE_L,X1(:,i1)});
%     S=anovan(X3,{'Pre','M'});
    [h,p] = ttest2(PRE_L,X1(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_L =[pVal_L Y];
    
    %Control with Pre
%     X3 =[PRE_C,X2(:,i1)];
    [h,p] = ttest2(PRE_C,X2(:,i1));
%     S=anovan(X3,{'Pre','M'});
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_C =[pVal_C Y];
    
    % Both
%     X3 =[X1(:,i1),X2(:,i1)];   
%     S=anovan(X3,{'Opto','Control'});
    [h,p] = ttest2(X1(:,i1),X2(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_Both =[pVal_Both Y]; 
end
    
stats_NearBF.pVal_L = pVal_L; stats_NearBF.pVal_C = pVal_C; stats_NearBF.pVal_Both = pVal_Both; 


X1 = Laser_M(:,:,3); 
X2 = nan(size(X1));
X2(1:size(Control_M(:,:,1),1),:)= Control_M(:,:,3);
PRE_L = mean(X1(:,1:3),2);
PRE_C = mean(X2(:,1:3),2);
pVal_L =[];
pVal_C =[];
pVal_Both = [];

for i1 =1:11
    
    %Laser with Pre
%     X3 =[PRE_L,X1(:,i1)];
    %L = [1,2];
    %O = teg_repeated_measures_ANOVA(X3,L,{'Time','Measurement'})
%     S=anovan([1:12]',{PRE_L,X1(:,i1)});
%     S=anovan(X3,{'Pre','M'});
    [h,p] = ttest2(PRE_L,X1(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_L =[pVal_L Y];
    
    %Control with Pre
%     X3 =[PRE_C,X2(:,i1)];
    [h,p] = ttest2(PRE_C,X2(:,i1));
%     S=anovan(X3,{'Pre','M'});
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_C =[pVal_C Y];
    
    % Both
%     X3 =[X1(:,i1),X2(:,i1)];   
%     S=anovan(X3,{'Opto','Control'});
    [h,p] = ttest2(X1(:,i1),X2(:,i1));
    Y = h*p;
    if Y == 0, Y = nan; end
    pVal_Both =[pVal_Both Y]; 
end
    
stats_NonBF.pVal_L = pVal_L; stats_NonBF.pVal_C = pVal_C; stats_NonBF.pVal_Both = pVal_Both; 

figure
subplot(1,3,1)
S1 =Laser_M(:,:,1);
S2 = nan(size(S1));
S2(1:7,:) = Control_M(:,:,1);
[stat,pperm,hperm,statperm]= elbe_permtest_005(S1,S2,1,'raw',1000,0,0.05,0);
SigPos=find(pperm <= 0.025 | pperm >= 1.025);
shadedErrorBar([],nanmean(Control_M(:,:,1),1),C_SEM(:,:,1),'lineprops',{'b','markerfacecolor','b'})
hold on
shadedErrorBar([],nanmean(Laser_M(:,:,1),1),L_SEM(:,:,1),'lineprops',{'r','markerfacecolor','r'})
prime_M = mean(nanmean(STD_prime(:,:,1)'));
prime_STD = mean(nanstd(STD_prime(:,:,1)'));
plot([SigPos],[ones(size(SigPos))*1.5],'k*')
% prime_STD = std(nanmean(STD_prime(:,:,1)))./sqrt(sum(~isnan(STD_prime(:,:,1))));
% prime_STD = mean(prime_STD);
% plot([1 length(Data.GS_based)],[prime_M+3*AllSTD_1 prime_M+3*AllSTD_1],':k')
% plot([1 length(Data.GS_based)],[prime_M-3*AllSTD_1 prime_M-3*AllSTD_1],':k')
% y1 = ones(1,11)*1.6; y2 = ones(1,11)*1.65;y3 = ones(1,11)*1.7;
% 
% plot([1:11].*(~isnan(stats_BF.pVal_L)),y1.*(~isnan(stats_BF.pVal_L)),'r*',...
%     [1:11].*(~isnan(stats_BF.pVal_C)),y2.*(~isnan(stats_BF.pVal_C)),'b*',...
%     [1:11].*(~isnan(stats_BF.pVal_Both)).*(~isnan(stats_BF.pVal_Both)),y3,'k*')
ylim([0.8 1.8])
xlim ([1 length(Data.GS_based)])
xticks([1 2 3 4 5 6 7 8 9 10 11])
xticklabels({'Pre(n-2)','Pre(n-1)', 'Pre(n)','Laser', 'Post 1','Post 2', 'Post 3','Post 4','Post 5','Post 6','Post 7'})
xtickangle(270)
set(gca,'fontsize',12,'FontWeight','bold')

title('BF')
legend ('Control','','','','Opto','','','','Location','best')


subplot(1,3,2)
S1 =Laser_M(:,:,2);
S2 = nan(size(S1));
S2(1:7,:) = Control_M(:,:,2);
[stat,pperm,hperm,statperm]= elbe_permtest_005(S1,S2,1,'raw',1000,0,0.05,0);
SigPos=find(pperm <= 0.025 | pperm >= 1.025);
shadedErrorBar([],nanmean(Control_M(:,:,2),1),C_SEM(:,:,2),'lineprops',{'b','markerfacecolor','b'})
hold on
shadedErrorBar([],nanmean(Laser_M(:,:,2),1),L_SEM(:,:,2),'lineprops',{'r','markerfacecolor','r'})
prime_M = mean(nanmean(STD_prime(:,:,2)'));
prime_STD = mean(nanstd(STD_prime(:,:,2)'));
plot([SigPos],[ones(size(SigPos))*1.5],'k*')
% prime_STD = std(nanmean(STD_prime(:,:,2)))./sqrt(sum(~isnan(STD_prime(:,:,2))));
% prime_STD = mean(prime_STD);
% plot([1 length(Data.GS_based)],[prime_M+3*AllSTD_2 prime_M+3*AllSTD_2],':k')
% plot([1 length(Data.GS_based)],[prime_M-3*AllSTD_2 prime_M-3*AllSTD_2],':k')
% y1 = ones(1,11)*1.6; y2 = ones(1,11)*1.65;y3 = ones(1,11)*1.7;
% 
% plot([1:11].*(~isnan(stats_NearBF.pVal_L)),y1.*(~isnan(stats_NearBF.pVal_L)),'r*',...
%     [1:11].*(~isnan(stats_NearBF.pVal_C)),y2.*(~isnan(stats_NearBF.pVal_C)),'b*',...
%     [1:11].*(~isnan(stats_NearBF.pVal_Both)).*(~isnan(stats_NearBF.pVal_Both)),y3,'k*')
ylim([0.8 1.8])
xlim ([1 length(Data.GS_based)])
xticks([1 2 3 4 5 6 7 8 9 10 11])
xticklabels({'Pre(n-2)','Pre(n-1)', 'Pre(n)','Laser', 'Post 1','Post 2', 'Post 3','Post 4','Post 5','Post 6','Post 7'})
xtickangle(270)
title('± near BF (1-2 Oct)')
set(gca,'fontsize',12,'FontWeight','bold')

legend ('Control','','','','Opto','','','','Location','best')

subplot(1,3,3)
S1 =Laser_M(:,:,3);
S2 = nan(size(S1));
S2(1:7,:) = Control_M(:,:,3);
[stat,pperm,hperm,statperm]= elbe_permtest_005(S1,S2,1,'raw',1000,0,0.05,0);
SigPos=find(pperm <= 0.025 | pperm >= 1.025);
shadedErrorBar([],nanmean(Control_M(:,:,3),1),C_SEM(:,:,3),'lineprops',{'b','markerfacecolor','b'})
hold on
shadedErrorBar([],nanmean(Laser_M(:,:,3),1),L_SEM(:,:,3),'lineprops',{'r','markerfacecolor','r'})
prime_M = mean(nanmean(STD_prime(:,:,3)'));
prime_STD = mean(nanstd(STD_prime(:,:,3)'));
plot([SigPos],[ones(size(SigPos))*1.5],'k*')
% prime_STD = std(nanmean(STD_prime(:,:,3)))./sqrt(sum(~isnan(STD_prime(:,:,3))));
% prime_STD = mean(prime_STD);
% plot([1 length(Data.GS_based)],[prime_M+3*AllSTD_3 prime_M+3*AllSTD_3],':k')
% plot([1 length(Data.GS_based)],[prime_M-3*AllSTD_3 prime_M-3*AllSTD_3],':k')
% y1 = ones(1,11)*1.6; y2 = ones(1,11)*1.65;y3 = ones(1,11)*1.7;
% plot([1:11].*(~isnan(stats_NonBF.pVal_L)),y1.*(~isnan(stats_NonBF.pVal_L)),'r*',...
%     [1:11].*(~isnan(stats_NonBF.pVal_C)),y2.*(~isnan(stats_NonBF.pVal_C)),'b*',...
%     [1:11].*(~isnan(stats_NonBF.pVal_Both)).*(~isnan(stats_NonBF.pVal_Both)),y3,'k*')
ylim([0.8 1.8])
xlim ([1 length(Data.GS_based)])
xticks([1 2 3 4 5 6 7 8 9 10 11])
xticklabels({'Pre(n-2)','Pre(n-1)', 'Pre(n)','Laser', 'Post 1','Post 2', 'Post 3','Post 4','Post 5','Post 6','Post 7'})
xtickangle(270)
title('± non BF (3-4 Oct)')
set(gca,'fontsize',12,'FontWeight','bold')

legend ('Control','','','','Opto','','','','Location','best')

%  %%%keyboard
clear


load('K:\CSD_dynamic_analysis\DATA\Output\Output_Input_HighP_7post_Data_Threshold_25_Zscore_0_binned_1_mirror_0.mat')
% measurements pre, laser, post
M = [1, 4, 9]
SORT = 'ST_based'; 
BF = Data.BF_Pos;
BFlim= 2;

LNorm = nan(length(Data.names),3);
for i1 = 1:3
    LNorm(:,i1,1) =  Data.(SORT)(i1).Full_RMS_AVREC(:,BF); %Full AVREC RMS
    LNorm(:,i1,2) =  Data.(SORT)(i1).Early_RMS_RELRES(:,BF);% Full RELRES RMS   
end
LNorm = nanmean(LNorm,2);

figure
subplot(2,1,1)
Pre(1,:)= nanmean(Data.(SORT)(M(1)).Full_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1));
Pre(2,:)= nanstd(Data.(SORT)(M(1)).Full_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1))./sqrt(sum(~isnan(Data.(SORT)(M(1)).Full_RMS_AVREC(:,BF-BFlim:BF+BFlim))));

L(1,:)= nanmean(Data.(SORT)(M(2)).Full_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1));
L(2,:)= nanstd(Data.(SORT)(M(2)).Full_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1))./sqrt(sum(~isnan(Data.(SORT)(M(2)).Full_RMS_AVREC(:,BF-BFlim:BF+BFlim))));

Post(1,:)= nanmean(Data.(SORT)(M(3)).Full_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1));
Post(2,:)= nanstd(Data.(SORT)(M(3)).Full_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1))./sqrt(sum(~isnan(Data.(SORT)(M(3)).Full_RMS_AVREC(:,BF-BFlim:BF+BFlim))));

shadedErrorBar([],Pre(1,:),Pre(2,:),'lineprops',{'b','markerfacecolor','b'})
hold on
shadedErrorBar([],L(1,:),L(2,:),'lineprops',{'g','markerfacecolor','g'})
shadedErrorBar([],Post(1,:),Post(2,:),'lineprops',{'r','markerfacecolor','r'})
xticks([1 2 3 4 5 ])
xticklabels({'- nonBF','- nearBF', 'BF','+ nearBF', '+ nonBF'})
ylabel('normalized Full RMS AVREC','FontSize',16,'FontWeight','bold')
set(gca,'fontsize',12,'FontWeight','bold')
xlabel('Frequency bins','FontSize',16,'FontWeight','bold')
legend('Pre','','','','Opto','','','','Post','','','','Location','best')
title('Full RMS AVREC [3 Pres BF normalized]')

subplot(2,1,2)
Pre(1,:)= nanmean(Data.(SORT)(M(1)).Full_RMS_RELRES(:,BF-BFlim:BF+BFlim));
Pre(2,:)= nanstd(Data.(SORT)(M(1)).Full_RMS_RELRES(:,BF-BFlim:BF+BFlim))./sqrt(sum(~isnan(Data.(SORT)(M(1)).Full_RMS_RELRES(:,BF-BFlim:BF+BFlim))));
L(1,:)= nanmean(Data.(SORT)(M(2)).Full_RMS_RELRES(:,BF-BFlim:BF+BFlim));
L(2,:)= nanstd(Data.(SORT)(M(2)).Full_RMS_RELRES(:,BF-BFlim:BF+BFlim))./sqrt(sum(~isnan(Data.(SORT)(M(2)).Full_RMS_RELRES(:,BF-BFlim:BF+BFlim))));
Post(1,:)= nanmean(Data.(SORT)(M(3)).Full_RMS_RELRES(:,BF-BFlim:BF+BFlim));
Post(2,:)= nanstd(Data.(SORT)(M(3)).Full_RMS_RELRES(:,BF-BFlim:BF+BFlim))./sqrt(sum(~isnan(Data.(SORT)(M(3)).Full_RMS_RELRES(:,BF-BFlim:BF+BFlim))));

shadedErrorBar([],Pre(1,:),Pre(2,:),'lineprops',{'b','markerfacecolor','b'})
hold on
shadedErrorBar([],L(1,:),L(2,:),'lineprops',{'g','markerfacecolor','g'})
shadedErrorBar([],Post(1,:),Post(2,:),'lineprops',{'r','markerfacecolor','r'})
xticks([1 2 3 4 5 ])
xticklabels({'- nonBF','- nearBF', 'BF','+ nearBF', '+ nonBF'})
ylabel('Full RMS RelResCSD [%]','FontSize',16,'FontWeight','bold')
set(gca,'fontsize',12,'FontWeight','bold')
xlabel('Frequency bins','FontSize',16,'FontWeight','bold')
legend('Pre','','','','Opto','','','','Post','','','','Location','best')
title('Full RMS RELRES')
% keyboard


%%% Lippi
LNorm = nan(length(Data.names),3);
for i1 = 1:3
    LNorm(:,i1,1) =  Data.(SORT)(i1).Full_RMS_AVREC(:,BF); %Full AVREC RMS
    LNorm(:,i1,2) =  Data.(SORT)(i1).Full_RMS_RELRES(:,BF);% Full RELRES RMS   
end
LNorm = nanmean(LNorm,2);

Pre_1 = nan(12,5,3);
Pre_1(:,1:5,1) = Data.(SORT)(1).Full_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);
Pre_1(:,1:5,2) = Data.(SORT)(2).Full_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);
Pre_1(:,1:5,3) = Data.(SORT)(3).Full_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);

PRE_1_STD = vertcat(Pre_1(:,:,1),Pre_1(:,:,2),Pre_1(:,:,3));
PRE_1_STD = nanstd(PRE_1_STD)./sum(sum(~isnan(Pre_1),3));

Post_1 = nan(12,5,3);
Post_1(:,1:5,1) = Data.(SORT)(5).Full_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);
Post_1(:,1:5,2) = Data.(SORT)(6).Full_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);
Post_1(:,1:5,3) = Data.(SORT)(7).Full_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);

Post_1_STD = vertcat(Post_1(:,:,1),Post_1(:,:,2),Post_1(:,:,3));
Post_1_STD = nanstd(Post_1_STD)./sum(sum(~isnan(Post_1),3));


Post_2 = nan(12,5,2);
Post_2(:,1:5,1) = Data.(SORT)(8).Full_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);
Post_2(:,1:5,2) = Data.(SORT)(9).Full_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);
% Post_2(:,1:5,3) = Data.(SORT)(10).Full_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);
% Post_2(:,1:5,3) = Data.(SORT)(11).Full_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);

Post_2_STD = vertcat(Post_2(:,:,1),Post_2(:,:,2));
Post_2_STD = nanstd(Post_2_STD)./sum(sum(~isnan(Post_2),3));

Post_3 = nan(12,5,2);
Post_3(:,1:5,1) = Data.(SORT)(10).Full_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);
Post_3(:,1:5,2) = Data.(SORT)(11).Full_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);
% Post_3(:,1:5,3) = Data.(SORT)(11).Full_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);

Post_3_STD = vertcat(Post_3(:,:,1),Post_3(:,:,2));
Post_3_STD = nanstd(Post_3_STD)./sum(sum(~isnan(Post_3),3));


figure
subplot(2,1,1)


L(1,:)= nanmean(Data.(SORT)(M(2)).Full_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1));
L(2,:)= nanstd(Data.(SORT)(M(2)).Full_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1))./sqrt(sum(~isnan(Data.(SORT)(M(2)).Full_RMS_AVREC(:,BF-BFlim:BF+BFlim))));

plot(nanmean(nanmean(Pre_1,3)),'LineWidth',2)% shadedErrorBar([],nanmean(nanmean(Pre_1,3)),PRE_1_STD,'lineprops',{'b','markerfacecolor','b'})
hold on
plot(L(1,:),'LineWidth',2)% shadedErrorBar([],L(1,:),L(2,:),'lineprops',{'g','markerfacecolor','g'})
plot(nanmean(nanmean(Post_1,3)),'LineWidth',2)% shadedErrorBar([],nanmean(nanmean(Post_1,3)),Post_1_STD,'lineprops',{'r','markerfacecolor','r'})
plot(nanmean(nanmean(Post_2,3)),'LineWidth',2)% shadedErrorBar([],nanmean(nanmean(Post_2,3)),Post_2_STD,'lineprops',{'y','markerfacecolor','y'})
plot(nanmean(nanmean(Post_3,3)),'LineWidth',2)% shadedErrorBar([],nanmean(nanmean(Post_3,3)),Post_3_STD,'lineprops',{'c','markerfacecolor','c'})
xticks([1 2 3 4 5 ])
xticklabels({'- nonBF','- nearBF', 'BF','+ nearBF', '+ nonBF'})
ylabel('normalized Full RMS AVREC','FontSize',16,'FontWeight','bold')
set(gca,'fontsize',12,'FontWeight','bold')
xlabel('Frequency bins','FontSize',16,'FontWeight','bold')
legend('Pre 1-3','Opto','Post 1-3','Post 4-5','Post 6-7','Location','best')
title('Full RMS AVREC [3 Pres BF normalized]')


subplot(2,1,2)
Pre_1 = nan(12,5,3);
Pre_1(:,1:5,1) = Data.(SORT)(1).Full_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);
Pre_1(:,1:5,2) = Data.(SORT)(2).Full_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);
Pre_1(:,1:5,3) = Data.(SORT)(3).Full_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);

PRE_1_STD = vertcat(Pre_1(:,:,1),Pre_1(:,:,2),Pre_1(:,:,3));
PRE_1_STD = nanstd(PRE_1_STD)./sum(sum(~isnan(Pre_1),3));

Post_1 = nan(12,5,3);
Post_1(:,1:5,1) = Data.(SORT)(5).Full_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);
Post_1(:,1:5,2) = Data.(SORT)(6).Full_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);
Post_1(:,1:5,3) = Data.(SORT)(7).Full_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);

Post_1_STD = vertcat(Post_1(:,:,1),Post_1(:,:,2),Post_1(:,:,3));
Post_1_STD = nanstd(Post_1_STD)./sum(sum(~isnan(Post_1),3));


Post_2 = nan(12,5,2);
Post_2(:,1:5,1) = Data.(SORT)(8).Full_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);
Post_2(:,1:5,2) = Data.(SORT)(9).Full_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);
% Post_2(:,1:5,3) = Data.(SORT)(10).Full_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);
% Post_2(:,1:5,3) = Data.(SORT)(11).Full_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);

Post_2_STD = vertcat(Post_2(:,:,1),Post_2(:,:,2));
Post_2_STD = nanstd(Post_2_STD)./sum(sum(~isnan(Post_2),3));

Post_3 = nan(12,5,2);
Post_3(:,1:5,1) = Data.(SORT)(10).Full_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);
Post_3(:,1:5,2) = Data.(SORT)(11).Full_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);
% Post_3(:,1:5,3) = Data.(SORT)(11).Full_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);

Post_3_STD = vertcat(Post_3(:,:,1),Post_3(:,:,2));
Post_3_STD = nanstd(Post_3_STD)./sum(sum(~isnan(Post_3),3));

clear L
L(1,:)= nanmean(Data.(SORT)(M(2)).Full_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2));
L(2,:)= nanstd(Data.(SORT)(M(2)).Full_RMS_RELRES(:,BF-BFlim:BF+BFlim))./sqrt(sum(~isnan(Data.(SORT)(M(2)).Full_RMS_RELRES(:,BF-BFlim:BF+BFlim))));

plot(nanmean(nanmean(Pre_1,3)),'LineWidth',2)% shadedErrorBar([],nanmean(nanmean(Pre_1,3)),PRE_1_STD,'lineprops',{'b','markerfacecolor','b'})
hold on
plot(L(1,:),'LineWidth',2)% shadedErrorBar([],L(1,:),L(2,:),'lineprops',{'g','markerfacecolor','g'})
plot(nanmean(nanmean(Post_1,3)),'LineWidth',2)% shadedErrorBar([],nanmean(nanmean(Post_1,3)),Post_1_STD,'lineprops',{'r','markerfacecolor','r'})
plot(nanmean(nanmean(Post_2,3)),'LineWidth',2)% shadedErrorBar([],nanmean(nanmean(Post_2,3)),Post_2_STD,'lineprops',{'y','markerfacecolor','y'})
plot(nanmean(nanmean(Post_3,3)),'LineWidth',2)% shadedErrorBar([],nanmean(nanmean(Post_3,3)),Post_3_STD,'lineprops',{'c','markerfacecolor','c'})
% ylim([0.35 0.8])
xticks([1 2 3 4 5 ])
xticklabels({'- nonBF','- nearBF', 'BF','+ nearBF', '+ nonBF'})
ylabel('normalized Full RMS AVREC','FontSize',16,'FontWeight','bold')
set(gca,'fontsize',12,'FontWeight','bold')
xlabel('Frequency bins','FontSize',16,'FontWeight','bold')
legend('Pre 1-3','Opto','Post 1-3','Post 4-5','Post 6-7','Location','best')
title('Full RMS RELRES [3 Pres BF normalized]')



LNorm = nan(length(Data.names),3);
for i1 = 1:3
    LNorm(:,i1,1) =  Data.(SORT)(i1).Early_RMS_AVREC(:,BF); %Full AVREC RMS
    LNorm(:,i1,2) =  Data.(SORT)(i1).Early_RMS_RELRES(:,BF);% Full RELRES RMS   
end
LNorm = nanmean(LNorm,2);

Pre_1 = nan(12,5,3);
Pre_1(:,1:5,1) = Data.(SORT)(1).Early_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);
Pre_1(:,1:5,2) = Data.(SORT)(2).Early_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);
Pre_1(:,1:5,3) = Data.(SORT)(3).Early_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);

PRE_1_STD = vertcat(Pre_1(:,:,1),Pre_1(:,:,2),Pre_1(:,:,3));
PRE_1_STD = nanstd(PRE_1_STD)./sum(sum(~isnan(Pre_1),3));

Post_1 = nan(12,5,3);
Post_1(:,1:5,1) = Data.(SORT)(5).Early_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);
Post_1(:,1:5,2) = Data.(SORT)(6).Early_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);
Post_1(:,1:5,3) = Data.(SORT)(7).Early_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);

Post_1_STD = vertcat(Post_1(:,:,1),Post_1(:,:,2),Post_1(:,:,3));
Post_1_STD = nanstd(Post_1_STD)./sum(sum(~isnan(Post_1),3));


Post_2 = nan(12,5,2);
Post_2(:,1:5,1) = Data.(SORT)(8).Early_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);
Post_2(:,1:5,2) = Data.(SORT)(9).Early_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);
% Post_2(:,1:5,3) = Data.(SORT)(10).Early_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);
% Post_2(:,1:5,3) = Data.(SORT)(11).Early_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);

Post_2_STD = vertcat(Post_2(:,:,1),Post_2(:,:,2));
Post_2_STD = nanstd(Post_2_STD)./sum(sum(~isnan(Post_2),3));

Post_3 = nan(12,5,2);
Post_3(:,1:5,1) = Data.(SORT)(10).Early_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);
Post_3(:,1:5,2) = Data.(SORT)(11).Early_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);
% Post_3(:,1:5,3) = Data.(SORT)(11).Early_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);

Post_3_STD = vertcat(Post_3(:,:,1),Post_3(:,:,2));
Post_3_STD = nanstd(Post_3_STD)./sum(sum(~isnan(Post_3),3));


figure
subplot(2,1,1)


L(1,:)= nanmean(Data.(SORT)(M(2)).Early_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1));
L(2,:)= nanstd(Data.(SORT)(M(2)).Early_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1))./sqrt(sum(~isnan(Data.(SORT)(M(2)).Early_RMS_AVREC(:,BF-BFlim:BF+BFlim))));

plot(nanmean(nanmean(Pre_1,3)),'LineWidth',2)% shadedErrorBar([],nanmean(nanmean(Pre_1,3)),PRE_1_STD,'lineprops',{'b','markerfacecolor','b'})
hold on
plot(L(1,:),'LineWidth',2)% shadedErrorBar([],L(1,:),L(2,:),'lineprops',{'g','markerfacecolor','g'})
plot(nanmean(nanmean(Post_1,3)),'LineWidth',2)% shadedErrorBar([],nanmean(nanmean(Post_1,3)),Post_1_STD,'lineprops',{'r','markerfacecolor','r'})
plot(nanmean(nanmean(Post_2,3)),'LineWidth',2)% shadedErrorBar([],nanmean(nanmean(Post_2,3)),Post_2_STD,'lineprops',{'y','markerfacecolor','y'})
plot(nanmean(nanmean(Post_3,3)),'LineWidth',2)% shadedErrorBar([],nanmean(nanmean(Post_3,3)),Post_3_STD,'lineprops',{'c','markerfacecolor','c'})
xticks([1 2 3 4 5 ])
xticklabels({'- nonBF','- nearBF', 'BF','+ nearBF', '+ nonBF'})
ylabel('normalized RMS AVREC','FontSize',16,'FontWeight','bold')
set(gca,'fontsize',12,'FontWeight','bold')
xlabel('Frequency bins','FontSize',16,'FontWeight','bold')
legend('Pre 1-3','Opto','Post 1-3','Post 4-5','Post 6-7','Location','best')
title('Early RMS AVREC [3 Pres BF normalized]')


subplot(2,1,2)
Pre_1 = nan(12,5,3);
Pre_1(:,1:5,1) = Data.(SORT)(1).Early_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);
Pre_1(:,1:5,2) = Data.(SORT)(2).Early_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);
Pre_1(:,1:5,3) = Data.(SORT)(3).Early_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);

PRE_1_STD = vertcat(Pre_1(:,:,1),Pre_1(:,:,2),Pre_1(:,:,3));
PRE_1_STD = nanstd(PRE_1_STD)./sum(sum(~isnan(Pre_1),3));

Post_1 = nan(12,5,3);
Post_1(:,1:5,1) = Data.(SORT)(5).Early_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);
Post_1(:,1:5,2) = Data.(SORT)(6).Early_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);
Post_1(:,1:5,3) = Data.(SORT)(7).Early_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);

Post_1_STD = vertcat(Post_1(:,:,1),Post_1(:,:,2),Post_1(:,:,3));
Post_1_STD = nanstd(Post_1_STD)./sum(sum(~isnan(Post_1),3));


Post_2 = nan(12,5,2);
Post_2(:,1:5,1) = Data.(SORT)(8).Early_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);
Post_2(:,1:5,2) = Data.(SORT)(9).Early_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);
% Post_2(:,1:5,3) = Data.(SORT)(10).Early_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);
% Post_2(:,1:5,3) = Data.(SORT)(11).Early_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);

Post_2_STD = vertcat(Post_2(:,:,1),Post_2(:,:,2));
Post_2_STD = nanstd(Post_2_STD)./sum(sum(~isnan(Post_2),3));

Post_3 = nan(12,5,2);
Post_3(:,1:5,1) = Data.(SORT)(10).Early_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);
Post_3(:,1:5,2) = Data.(SORT)(11).Early_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);
% Post_3(:,1:5,3) = Data.(SORT)(11).Early_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);

Post_3_STD = vertcat(Post_3(:,:,1),Post_3(:,:,2));
Post_3_STD = nanstd(Post_3_STD)./sum(sum(~isnan(Post_3),3));

clear L
L(1,:)= nanmean(Data.(SORT)(M(2)).Early_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2));
L(2,:)= nanstd(Data.(SORT)(M(2)).Early_RMS_RELRES(:,BF-BFlim:BF+BFlim))./sqrt(sum(~isnan(Data.(SORT)(M(2)).Early_RMS_RELRES(:,BF-BFlim:BF+BFlim))));

plot(nanmean(nanmean(Pre_1,3)),'LineWidth',2)% shadedErrorBar([],nanmean(nanmean(Pre_1,3)),PRE_1_STD,'lineprops',{'b','markerfacecolor','b'})
hold on
plot(L(1,:),'LineWidth',2)% shadedErrorBar([],L(1,:),L(2,:),'lineprops',{'g','markerfacecolor','g'})
plot(nanmean(nanmean(Post_1,3)),'LineWidth',2)% shadedErrorBar([],nanmean(nanmean(Post_1,3)),Post_1_STD,'lineprops',{'r','markerfacecolor','r'})
plot(nanmean(nanmean(Post_2,3)),'LineWidth',2)% shadedErrorBar([],nanmean(nanmean(Post_2,3)),Post_2_STD,'lineprops',{'y','markerfacecolor','y'})
plot(nanmean(nanmean(Post_3,3)),'LineWidth',2)% shadedErrorBar([],nanmean(nanmean(Post_3,3)),Post_3_STD,'lineprops',{'c','markerfacecolor','c'})
% ylim([0.35 0.8])
xticks([1 2 3 4 5 ])
xticklabels({'- nonBF','- nearBF', 'BF','+ nearBF', '+ nonBF'})
ylabel('normalized RMS RelRes','FontSize',16,'FontWeight','bold')
set(gca,'fontsize',12,'FontWeight','bold')
xlabel('Frequency bins','FontSize',16,'FontWeight','bold')
legend('Pre 1-3','Opto','Post 1-3','Post 4-5','Post 6-7','Location','best')
title('Early RMS RELRES [3 Pres BF normalized]')



LNorm = nan(length(Data.names),3);
for i1 = 1:3
    LNorm(:,i1,1) =  Data.(SORT)(i1).Late_RMS_AVREC(:,BF); %Full AVREC RMS
    LNorm(:,i1,2) =  Data.(SORT)(i1).Late_RMS_RELRES(:,BF);% Full RELRES RMS   
end
LNorm = nanmean(LNorm,2);

Pre_1 = nan(12,5,3);
Pre_1(:,1:5,1) = Data.(SORT)(1).Late_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);
Pre_1(:,1:5,2) = Data.(SORT)(2).Late_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);
Pre_1(:,1:5,3) = Data.(SORT)(3).Late_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);

PRE_1_STD = vertcat(Pre_1(:,:,1),Pre_1(:,:,2),Pre_1(:,:,3));
PRE_1_STD = nanstd(PRE_1_STD)./sum(sum(~isnan(Pre_1),3));

Post_1 = nan(12,5,3);
Post_1(:,1:5,1) = Data.(SORT)(5).Late_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);
Post_1(:,1:5,2) = Data.(SORT)(6).Late_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);
Post_1(:,1:5,3) = Data.(SORT)(7).Late_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);

Post_1_STD = vertcat(Post_1(:,:,1),Post_1(:,:,2),Post_1(:,:,3));
Post_1_STD = nanstd(Post_1_STD)./sum(sum(~isnan(Post_1),3));


Post_2 = nan(12,5,2);
Post_2(:,1:5,1) = Data.(SORT)(8).Late_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);
Post_2(:,1:5,2) = Data.(SORT)(9).Late_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);
% Post_2(:,1:5,3) = Data.(SORT)(10).Late_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);
% Post_2(:,1:5,3) = Data.(SORT)(11).Late_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);

Post_2_STD = vertcat(Post_2(:,:,1),Post_2(:,:,2));
Post_2_STD = nanstd(Post_2_STD)./sum(sum(~isnan(Post_2),3));

Post_3 = nan(12,5,2);
Post_3(:,1:5,1) = Data.(SORT)(10).Late_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);
Post_3(:,1:5,2) = Data.(SORT)(11).Late_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);
% Post_3(:,1:5,3) = Data.(SORT)(11).Late_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1);

Post_3_STD = vertcat(Post_3(:,:,1),Post_3(:,:,2));
Post_3_STD = nanstd(Post_3_STD)./sum(sum(~isnan(Post_3),3));


figure
subplot(2,1,1)


L(1,:)= nanmean(Data.(SORT)(M(2)).Late_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1));
L(2,:)= nanstd(Data.(SORT)(M(2)).Late_RMS_AVREC(:,BF-BFlim:BF+BFlim)./LNorm(:,:,1))./sqrt(sum(~isnan(Data.(SORT)(M(2)).Late_RMS_AVREC(:,BF-BFlim:BF+BFlim))));

plot(nanmean(nanmean(Pre_1,3)),'LineWidth',2)% shadedErrorBar([],nanmean(nanmean(Pre_1,3)),PRE_1_STD,'lineprops',{'b','markerfacecolor','b'})
hold on
plot(L(1,:),'LineWidth',2)% shadedErrorBar([],L(1,:),L(2,:),'lineprops',{'g','markerfacecolor','g'})
plot(nanmean(nanmean(Post_1,3)),'LineWidth',2)% shadedErrorBar([],nanmean(nanmean(Post_1,3)),Post_1_STD,'lineprops',{'r','markerfacecolor','r'})
plot(nanmean(nanmean(Post_2,3)),'LineWidth',2)% shadedErrorBar([],nanmean(nanmean(Post_2,3)),Post_2_STD,'lineprops',{'y','markerfacecolor','y'})
plot(nanmean(nanmean(Post_3,3)),'LineWidth',2)% shadedErrorBar([],nanmean(nanmean(Post_3,3)),Post_3_STD,'lineprops',{'c','markerfacecolor','c'})
xticks([1 2 3 4 5 ])
xticklabels({'- nonBF','- nearBF', 'BF','+ nearBF', '+ nonBF'})
ylabel('normalized RMS AVREC','FontSize',16,'FontWeight','bold')
set(gca,'fontsize',12,'FontWeight','bold')
xlabel('Frequency bins','FontSize',16,'FontWeight','bold')
legend('Pre 1-3','Opto','Post 1-3','Post 4-5','Post 6-7','Location','best')
title('Late RMS AVREC [3 Pres BF normalized]')


subplot(2,1,2)
Pre_1 = nan(12,5,3);
Pre_1(:,1:5,1) = Data.(SORT)(1).Late_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);
Pre_1(:,1:5,2) = Data.(SORT)(2).Late_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);
Pre_1(:,1:5,3) = Data.(SORT)(3).Late_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);

PRE_1_STD = vertcat(Pre_1(:,:,1),Pre_1(:,:,2),Pre_1(:,:,3));
PRE_1_STD = nanstd(PRE_1_STD)./sum(sum(~isnan(Pre_1),3));

Post_1 = nan(12,5,3);
Post_1(:,1:5,1) = Data.(SORT)(5).Late_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);
Post_1(:,1:5,2) = Data.(SORT)(6).Late_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);
Post_1(:,1:5,3) = Data.(SORT)(7).Late_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);

Post_1_STD = vertcat(Post_1(:,:,1),Post_1(:,:,2),Post_1(:,:,3));
Post_1_STD = nanstd(Post_1_STD)./sum(sum(~isnan(Post_1),3));


Post_2 = nan(12,5,2);
Post_2(:,1:5,1) = Data.(SORT)(8).Late_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);
Post_2(:,1:5,2) = Data.(SORT)(9).Late_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);
% Post_2(:,1:5,3) = Data.(SORT)(10).Late_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);
% Post_2(:,1:5,3) = Data.(SORT)(11).Late_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);

Post_2_STD = vertcat(Post_2(:,:,1),Post_2(:,:,2));
Post_2_STD = nanstd(Post_2_STD)./sum(sum(~isnan(Post_2),3));

Post_3 = nan(12,5,2);
Post_3(:,1:5,1) = Data.(SORT)(10).Late_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);
Post_3(:,1:5,2) = Data.(SORT)(11).Late_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);
% Post_3(:,1:5,3) = Data.(SORT)(11).Late_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2);

Post_3_STD = vertcat(Post_3(:,:,1),Post_3(:,:,2));
Post_3_STD = nanstd(Post_3_STD)./sum(sum(~isnan(Post_3),3));

clear L
L(1,:)= nanmean(Data.(SORT)(M(2)).Late_RMS_RELRES(:,BF-BFlim:BF+BFlim)./LNorm(:,:,2));
L(2,:)= nanstd(Data.(SORT)(M(2)).Late_RMS_RELRES(:,BF-BFlim:BF+BFlim))./sqrt(sum(~isnan(Data.(SORT)(M(2)).Late_RMS_RELRES(:,BF-BFlim:BF+BFlim))));

plot(nanmean(nanmean(Pre_1,3)),'LineWidth',2)% shadedErrorBar([],nanmean(nanmean(Pre_1,3)),PRE_1_STD,'lineprops',{'b','markerfacecolor','b'})
hold on
plot(L(1,:),'LineWidth',2)% shadedErrorBar([],L(1,:),L(2,:),'lineprops',{'g','markerfacecolor','g'})
plot(nanmean(nanmean(Post_1,3)),'LineWidth',2)% shadedErrorBar([],nanmean(nanmean(Post_1,3)),Post_1_STD,'lineprops',{'r','markerfacecolor','r'})
plot(nanmean(nanmean(Post_2,3)),'LineWidth',2)% shadedErrorBar([],nanmean(nanmean(Post_2,3)),Post_2_STD,'lineprops',{'y','markerfacecolor','y'})
plot(nanmean(nanmean(Post_3,3)),'LineWidth',2)% shadedErrorBar([],nanmean(nanmean(Post_3,3)),Post_3_STD,'lineprops',{'c','markerfacecolor','c'})
% ylim([0.35 0.8])
xticks([1 2 3 4 5 ])
xticklabels({'- nonBF','- nearBF', 'BF','+ nearBF', '+ nonBF'})
ylabel('normalized RMS RelRes','FontSize',16,'FontWeight','bold')
set(gca,'fontsize',12,'FontWeight','bold')
xlabel('Frequency bins','FontSize',16,'FontWeight','bold')
legend('Pre 1-3','Opto','Post 1-3','Post 4-5','Post 6-7','Location','best')
title('Late RMS RELRES [3 Pres BF normalized]')





clear
































%%% Bargraph Fig 4
M = [1, 4, 9]
SORT = 'ST_based';
BFlim= 0;

c = load('K:\CSD_dynamic_analysis\DATA\Output\Output_Input_Control_7post_n7_Data_Threshold_25_Zscore_0_binned_1_mirror_1.mat');
BF_C = c.Data.BF_Pos;
l = load('K:\CSD_dynamic_analysis\DATA\Output\Output_Input_HighP_7post_Data_Threshold_25_Zscore_0_binned_1_mirror_1.mat');
BF_L = l.Data.BF_Pos;
y = load('K:\CSD_dynamic_analysis\DATA\Output\Output_Input_YFP_7post_Data_Threshold_25_Zscore_0_binned_1_mirror_1.mat');
BF_Y = y.Data.BF_Pos;

CNorm = nan(length(c.Data.names),3);
for i1 = 1:3
    CNorm(:,i1,1) =  c.Data.(SORT)(i1).Full_RMS_AVREC(:,BF_C+BFlim);%full
    CNorm(:,i1,2) =  c.Data.(SORT)(i1).Early_RMS_AVREC(:,BF_C+BFlim);%early
    CNorm(:,i1,3) =  c.Data.(SORT)(i1).Late_RMS_AVREC(:,BF_C+BFlim);%late
end
CNorm = nanmean(CNorm,2);

LNorm = nan(length(l.Data.names),3);
for i1 = 1:3
    LNorm(:,i1,1) =  l.Data.(SORT)(i1).Full_RMS_AVREC(:,BF_L+BFlim);%full
    LNorm(:,i1,2) =  l.Data.(SORT)(i1).Early_RMS_AVREC(:,BF_L+BFlim);%early
    LNorm(:,i1,3) =  l.Data.(SORT)(i1).Late_RMS_AVREC(:,BF_L+BFlim);%late
end
LNorm = nanmean(LNorm,2);

YNorm = nan(length(y.Data.names),3);
for i1 = 1:3
    YNorm(:,i1,1) =  y.Data.(SORT)(i1).Full_RMS_AVREC(:,BF_Y+BFlim);%full
    YNorm(:,i1,2) =  y.Data.(SORT)(i1).Early_RMS_AVREC(:,BF_Y+BFlim);%early
    YNorm(:,i1,3) =  y.Data.(SORT)(i1).Late_RMS_AVREC(:,BF_Y+BFlim);%late
end
YNorm = nanmean(YNorm,2);

%Pre
C_Pre_full = c.Data.(SORT)(M(1)).Full_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,1);%Full
C_Pre_early = c.Data.(SORT)(M(1)).Early_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,2);%early
C_Pre_late = c.Data.(SORT)(M(1)).Late_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,3);%late

L_Pre_full = l.Data.(SORT)(M(1)).Full_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,1);%Full
L_Pre_early = l.Data.(SORT)(M(1)).Early_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,2);%early
L_Pre_late = l.Data.(SORT)(M(1)).Late_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,3);%late

Y_Pre_full = y.Data.(SORT)(M(1)).Full_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,1);%Full
Y_Pre_early = y.Data.(SORT)(M(1)).Early_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,2);%early
Y_Pre_late = y.Data.(SORT)(M(1)).Late_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,3);%late

%Laser
C_L_full = c.Data.(SORT)(M(2)).Full_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,1);%Full
C_L_early = c.Data.(SORT)(M(2)).Early_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,2);%early
C_L_late = c.Data.(SORT)(M(2)).Late_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,3);%late

L_L_full = l.Data.(SORT)(M(2)).Full_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,1);%Full
L_L_early = l.Data.(SORT)(M(2)).Early_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,2);%early
L_L_late = l.Data.(SORT)(M(2)).Late_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,3);%late

Y_L_full = y.Data.(SORT)(M(2)).Full_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,1);%Full
Y_L_early = y.Data.(SORT)(M(2)).Early_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,2);%early
Y_L_late = y.Data.(SORT)(M(2)).Late_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,3);%late


%Post
C_Post_full = c.Data.(SORT)(M(3)).Full_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,1);%Full
C_Post_early = c.Data.(SORT)(M(3)).Early_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,2);%early
C_Post_late = c.Data.(SORT)(M(3)).Late_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,3);%late

L_Post_full = l.Data.(SORT)(M(3)).Full_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,1);%Full
L_Post_early = l.Data.(SORT)(M(3)).Early_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,2);%early
L_Post_late = l.Data.(SORT)(M(3)).Late_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,3);%late

Y_Post_full = y.Data.(SORT)(M(3)).Full_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,1);%Full
Y_Post_early = y.Data.(SORT)(M(3)).Early_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,2);%early
Y_Post_late = y.Data.(SORT)(M(3)).Late_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,3);%late
% keyboard
figure
subplot(2,3,1)
bar([1 4 7],[nanmean(C_Pre_full) nanmean(C_L_full) nanmean(C_Post_full)],0.3)
hold on
bar([2 5 8],[nanmean(L_Pre_full) nanmean(L_L_full) nanmean(L_Post_full)],0.3)
bar([3 6 9],[nanmean(Y_Pre_full) nanmean(Y_L_full) nanmean(Y_Post_full)],0.3)

errorbar([1 4 7],[nanmean(C_Pre_full) nanmean(C_L_full) nanmean(C_Post_full)],...
    [nanstd(C_Pre_full)/sqrt(sum(~isnan(C_Pre_full)))...
    nanstd(C_L_full)/sqrt(sum(~isnan(C_L_full)))...
    nanstd(C_Post_full)/sqrt(sum(~isnan(C_Post_full)))],'.', 'LineWidth', 2)

errorbar([2 5 8],[nanmean(L_Pre_full) nanmean(L_L_full) nanmean(L_Post_full)],...
    [nanstd(L_Pre_full)/sqrt(sum(~isnan(L_Pre_full)))...
    nanstd(L_L_full)/sqrt(sum(~isnan(L_L_full)))...
    nanstd(L_Post_full)/sqrt(sum(~isnan(L_Post_full)))],'.', 'LineWidth', 2)

errorbar([3 6 9],[nanmean(Y_Pre_full) nanmean(Y_L_full) nanmean(Y_Post_full)],...
    [nanstd(Y_Pre_full)/sqrt(sum(~isnan(Y_Pre_full)))...
    nanstd(Y_L_full)/sqrt(sum(~isnan(Y_L_full)))...
    nanstd(Y_Post_full)/sqrt(sum(~isnan(Y_Post_full)))],'.', 'LineWidth', 2)

xticks([2 5 8])
xticklabels({'Pre','Opto','Post'})
ylim([0.75 1.25])
set(gca,'fontsize',12,'FontWeight','bold')
title ('Full RMS AVREC')
legend('Control','Opto','YFP','Location','best')

subplot(2,3,2)
bar([1 4 7],[nanmean(C_Pre_early) nanmean(C_L_early) nanmean(C_Post_early)],0.3)
hold on
bar([2 5 8],[nanmean(L_Pre_early) nanmean(L_L_early) nanmean(L_Post_early)],0.3)
bar([3 6 9],[nanmean(Y_Pre_early) nanmean(Y_L_early) nanmean(Y_Post_early)],0.3)

errorbar([1 4 7],[nanmean(C_Pre_early) nanmean(C_L_early) nanmean(C_Post_early)],...
    [nanstd(C_Pre_early)/sqrt(sum(~isnan(C_Pre_early)))...
    nanstd(C_L_early)/sqrt(sum(~isnan(C_L_early)))...
    nanstd(C_Post_early)/sqrt(sum(~isnan(C_Post_early)))],'.', 'LineWidth', 2)

errorbar([2 5 8],[nanmean(L_Pre_early) nanmean(L_L_early) nanmean(L_Post_early)],...
    [nanstd(L_Pre_early)/sqrt(sum(~isnan(L_Pre_early)))...
    nanstd(L_L_early)/sqrt(sum(~isnan(L_L_early)))...
    nanstd(L_Post_early)/sqrt(sum(~isnan(L_Post_early)))],'.', 'LineWidth', 2)

errorbar([3 6 9],[nanmean(Y_Pre_early) nanmean(Y_L_early) nanmean(Y_Post_early)],...
    [nanstd(Y_Pre_early)/sqrt(sum(~isnan(Y_Pre_early)))...
    nanstd(Y_L_early)/sqrt(sum(~isnan(Y_L_early)))...
    nanstd(Y_Post_early)/sqrt(sum(~isnan(Y_Post_early)))],'.', 'LineWidth', 2)

xticks([2 5 8])
xticklabels({'Pre','Opto','Post'})
ylim([0.75 1.25])
set(gca,'fontsize',12,'FontWeight','bold')
title ('Early RMS AVREC')
legend('Control','Opto','YFP','Location','best')

subplot(2,3,3)
bar([1 4 7],[nanmean(C_Pre_late) nanmean(C_L_late) nanmean(C_Post_late)],0.3)
hold on
bar([2 5 8],[nanmean(L_Pre_late) nanmean(L_L_late) nanmean(L_Post_late)],0.3)
bar([3 6 9],[nanmean(Y_Pre_late) nanmean(Y_L_late) nanmean(Y_Post_late)],0.3)

errorbar([1 4 7],[nanmean(C_Pre_late) nanmean(C_L_late) nanmean(C_Post_late)],...
    [nanstd(C_Pre_late)/sqrt(sum(~isnan(C_Pre_late)))...
    nanstd(C_L_late)/sqrt(sum(~isnan(C_L_late)))...
    nanstd(C_Post_late)/sqrt(sum(~isnan(C_Post_late)))],'.', 'LineWidth', 2)

errorbar([2 5 8],[nanmean(L_Pre_late) nanmean(L_L_late) nanmean(L_Post_late)],...
    [nanstd(L_Pre_late)/sqrt(sum(~isnan(L_Pre_late)))...
    nanstd(L_L_late)/sqrt(sum(~isnan(L_L_late)))...
    nanstd(L_Post_late)/sqrt(sum(~isnan(L_Post_late)))],'.', 'LineWidth', 2)

errorbar([3 6 9],[nanmean(Y_Pre_late) nanmean(Y_L_late) nanmean(Y_Post_late)],...
    [nanstd(Y_Pre_late)/sqrt(sum(~isnan(Y_Pre_late)))...
    nanstd(Y_L_late)/sqrt(sum(~isnan(Y_L_late)))...
    nanstd(Y_Post_late)/sqrt(sum(~isnan(Y_Post_late)))],'.', 'LineWidth', 2)

xticks([2 5 8])
xticklabels({'Pre','Opto','Post'})
ylim([0.75 1.25])
set(gca,'fontsize',12,'FontWeight','bold')
title ('Late RMS AVREC')
legend('Control','Opto','YFP','Location','best')

%%%%%%%%%%%%RELRES
CNorm = nan(length(c.Data.names),3);
for i1 = 1:3
    CNorm(:,i1,1) =  c.Data.(SORT)(i1).Full_RMS_RELRES(:,BF_C+BFlim);%full
    CNorm(:,i1,2) =  c.Data.(SORT)(i1).Early_RMS_RELRES(:,BF_C+BFlim);%early
    CNorm(:,i1,3) =  c.Data.(SORT)(i1).Late_RMS_RELRES(:,BF_C+BFlim);%late
end
CNorm = nanmean(CNorm,2);

LNorm = nan(length(l.Data.names),3);
for i1 = 1:3
    LNorm(:,i1,1) =  l.Data.(SORT)(i1).Full_RMS_RELRES(:,BF_L+BFlim);%full
    LNorm(:,i1,2) =  l.Data.(SORT)(i1).Early_RMS_RELRES(:,BF_L+BFlim);%early
    LNorm(:,i1,3) =  l.Data.(SORT)(i1).Late_RMS_RELRES(:,BF_L+BFlim);%late
end
LNorm = nanmean(LNorm,2);

YNorm = nan(length(y.Data.names),3);
for i1 = 1:3
    YNorm(:,i1,1) =  y.Data.(SORT)(i1).Full_RMS_RELRES(:,BF_Y+BFlim);%full
    YNorm(:,i1,2) =  y.Data.(SORT)(i1).Early_RMS_RELRES(:,BF_Y+BFlim);%early
    YNorm(:,i1,3) =  y.Data.(SORT)(i1).Late_RMS_RELRES(:,BF_Y+BFlim);%late
end
YNorm = nanmean(YNorm,2);

%Pre
C_Pre_full = c.Data.(SORT)(M(1)).Full_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,1);%Full
C_Pre_early = c.Data.(SORT)(M(1)).Early_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,2);%early
C_Pre_late = c.Data.(SORT)(M(1)).Late_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,3);%late

L_Pre_full = l.Data.(SORT)(M(1)).Full_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,1);%Full
L_Pre_early = l.Data.(SORT)(M(1)).Early_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,2);%early
L_Pre_late = l.Data.(SORT)(M(1)).Late_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,3);%late

Y_Pre_full = y.Data.(SORT)(M(1)).Full_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,1);%Full
Y_Pre_early = y.Data.(SORT)(M(1)).Early_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,2);%early
Y_Pre_late = y.Data.(SORT)(M(1)).Late_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,3);%late

%Laser
C_L_full = c.Data.(SORT)(M(2)).Full_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,1);%Full
C_L_early = c.Data.(SORT)(M(2)).Early_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,2);%early
C_L_late = c.Data.(SORT)(M(2)).Late_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,3);%late

L_L_full = l.Data.(SORT)(M(2)).Full_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,1);%Full
L_L_early = l.Data.(SORT)(M(2)).Early_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,2);%early
L_L_late = l.Data.(SORT)(M(2)).Late_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,3);%late

Y_L_full = y.Data.(SORT)(M(2)).Full_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,1);%Full
Y_L_early = y.Data.(SORT)(M(2)).Early_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,2);%early
Y_L_late = y.Data.(SORT)(M(2)).Late_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,3);%late


%Post
C_Post_full = c.Data.(SORT)(M(3)).Full_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,1);%Full
C_Post_early = c.Data.(SORT)(M(3)).Early_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,2);%early
C_Post_late = c.Data.(SORT)(M(3)).Late_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,3);%late

L_Post_full = l.Data.(SORT)(M(3)).Full_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,1);%Full
L_Post_early = l.Data.(SORT)(M(3)).Early_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,2);%early
L_Post_late = l.Data.(SORT)(M(3)).Late_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,3);%late

Y_Post_full = y.Data.(SORT)(M(3)).Full_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,1);%Full
Y_Post_early = y.Data.(SORT)(M(3)).Early_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,2);%early
Y_Post_late = y.Data.(SORT)(M(3)).Late_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,3);%late


subplot(2,3,4)
bar([1 4 7],[nanmean(C_Pre_full) nanmean(C_L_full) nanmean(C_Post_full)],0.3)
hold on
bar([2 5 8],[nanmean(L_Pre_full) nanmean(L_L_full) nanmean(L_Post_full)],0.3)
bar([3 6 9],[nanmean(Y_Pre_full) nanmean(Y_L_full) nanmean(Y_Post_full)],0.3)

errorbar([1 4 7],[nanmean(C_Pre_full) nanmean(C_L_full) nanmean(C_Post_full)],...
    [nanstd(C_Pre_full)/sqrt(sum(~isnan(C_Pre_full)))...
    nanstd(C_L_full)/sqrt(sum(~isnan(C_L_full)))...
    nanstd(C_Post_full)/sqrt(sum(~isnan(C_Post_full)))],'.', 'LineWidth', 2)

errorbar([2 5 8],[nanmean(L_Pre_full) nanmean(L_L_full) nanmean(L_Post_full)],...
    [nanstd(L_Pre_full)/sqrt(sum(~isnan(L_Pre_full)))...
    nanstd(L_L_full)/sqrt(sum(~isnan(L_L_full)))...
    nanstd(L_Post_full)/sqrt(sum(~isnan(L_Post_full)))],'.', 'LineWidth', 2)

errorbar([3 6 9],[nanmean(Y_Pre_full) nanmean(Y_L_full) nanmean(Y_Post_full)],...
    [nanstd(Y_Pre_full)/sqrt(sum(~isnan(Y_Pre_full)))...
    nanstd(Y_L_full)/sqrt(sum(~isnan(Y_L_full)))...
    nanstd(Y_Post_full)/sqrt(sum(~isnan(Y_Post_full)))],'.', 'LineWidth', 2)

xticks([2 5 8])
xticklabels({'Pre','Opto','Post'})
ylim([0.75 1.5])
set(gca,'fontsize',12,'FontWeight','bold')
title ('Full RMS RelRes')
legend('Control','Opto','YFP','Location','best')

subplot(2,3,5)
bar([1 4 7],[nanmean(C_Pre_early) nanmean(C_L_early) nanmean(C_Post_early)],0.3)
hold on
bar([2 5 8],[nanmean(L_Pre_early) nanmean(L_L_early) nanmean(L_Post_early)],0.3)
bar([3 6 9],[nanmean(Y_Pre_early) nanmean(Y_L_early) nanmean(Y_Post_early)],0.3)

errorbar([1 4 7],[nanmean(C_Pre_early) nanmean(C_L_early) nanmean(C_Post_early)],...
    [nanstd(C_Pre_early)/sqrt(sum(~isnan(C_Pre_early)))...
    nanstd(C_L_early)/sqrt(sum(~isnan(C_L_early)))...
    nanstd(C_Post_early)/sqrt(sum(~isnan(C_Post_early)))],'.', 'LineWidth', 2)

errorbar([2 5 8],[nanmean(L_Pre_early) nanmean(L_L_early) nanmean(L_Post_early)],...
    [nanstd(L_Pre_early)/sqrt(sum(~isnan(L_Pre_early)))...
    nanstd(L_L_early)/sqrt(sum(~isnan(L_L_early)))...
    nanstd(L_Post_early)/sqrt(sum(~isnan(L_Post_early)))],'.', 'LineWidth', 2)

errorbar([3 6 9],[nanmean(Y_Pre_early) nanmean(Y_L_early) nanmean(Y_Post_early)],...
    [nanstd(Y_Pre_early)/sqrt(sum(~isnan(Y_Pre_early)))...
    nanstd(Y_L_early)/sqrt(sum(~isnan(Y_L_early)))...
    nanstd(Y_Post_early)/sqrt(sum(~isnan(Y_Post_early)))],'.', 'LineWidth', 2)

xticks([2 5 8])
xticklabels({'Pre','Opto','Post'})
ylim([0.75 1.5])
set(gca,'fontsize',12,'FontWeight','bold')
title ('Early RMS RelRes')
legend('Control','Opto','YFP','Location','best')

subplot(2,3,6)
bar([1 4 7],[nanmean(C_Pre_late) nanmean(C_L_late) nanmean(C_Post_late)],0.3)
hold on
bar([2 5 8],[nanmean(L_Pre_late) nanmean(L_L_late) nanmean(L_Post_late)],0.3)
bar([3 6 9],[nanmean(Y_Pre_late) nanmean(Y_L_late) nanmean(Y_Post_late)],0.3)

errorbar([1 4 7],[nanmean(C_Pre_late) nanmean(C_L_late) nanmean(C_Post_late)],...
    [nanstd(C_Pre_late)/sqrt(sum(~isnan(C_Pre_late)))...
    nanstd(C_L_late)/sqrt(sum(~isnan(C_L_late)))...
    nanstd(C_Post_late)/sqrt(sum(~isnan(C_Post_late)))],'.', 'LineWidth', 2)

errorbar([2 5 8],[nanmean(L_Pre_late) nanmean(L_L_late) nanmean(L_Post_late)],...
    [nanstd(L_Pre_late)/sqrt(sum(~isnan(L_Pre_late)))...
    nanstd(L_L_late)/sqrt(sum(~isnan(L_L_late)))...
    nanstd(L_Post_late)/sqrt(sum(~isnan(L_Post_late)))],'.', 'LineWidth', 2)

errorbar([3 6 9],[nanmean(Y_Pre_late) nanmean(Y_L_late) nanmean(Y_Post_late)],...
    [nanstd(Y_Pre_late)/sqrt(sum(~isnan(Y_Pre_late)))...
    nanstd(Y_L_late)/sqrt(sum(~isnan(Y_L_late)))...
    nanstd(Y_Post_late)/sqrt(sum(~isnan(Y_Post_late)))],'.', 'LineWidth', 2)

xticks([2 5 8])
xticklabels({'Pre','Opto','Post'})
ylim([0.75 1.5])
set(gca,'fontsize',12,'FontWeight','bold')
title ('Late RMS RelRes')
legend('Control','Opto','YFP','Location','best')

clear C_Pre_full  C_Pre_early C_Pre_late

%%% RELRES

CNorm = nan(length(c.Data.names),3);
for i1 = 1:3
    CNorm(:,i1,1) =  c.Data.(SORT)(i1).Full_RMS_RELRES(:,BF_C+BFlim);%full
    CNorm(:,i1,2) =  c.Data.(SORT)(i1).Early_RMS_RELRES(:,BF_C+BFlim);%early
    CNorm(:,i1,3) =  c.Data.(SORT)(i1).Late_RMS_RELRES(:,BF_C+BFlim);%late
end
CNorm = nanmean(CNorm,2);

LNorm = nan(length(l.Data.names),3);
for i1 = 1:3
    LNorm(:,i1,1) =  l.Data.(SORT)(i1).Full_RMS_RELRES(:,BF_L+BFlim);%full
    LNorm(:,i1,2) =  l.Data.(SORT)(i1).Early_RMS_RELRES(:,BF_L+BFlim);%early
    LNorm(:,i1,3) =  l.Data.(SORT)(i1).Late_RMS_RELRES(:,BF_L+BFlim);%late
end
LNorm = nanmean(LNorm,2);

YNorm = nan(length(y.Data.names),3);
for i1 = 1:3
    YNorm(:,i1,1) =  y.Data.(SORT)(i1).Full_RMS_RELRES(:,BF_Y+BFlim);%full
    YNorm(:,i1,2) =  y.Data.(SORT)(i1).Early_RMS_RELRES(:,BF_Y+BFlim);%early
    YNorm(:,i1,3) =  y.Data.(SORT)(i1).Late_RMS_RELRES(:,BF_Y+BFlim);%late
end
YNorm = nanmean(YNorm,2);

%Pre
C_Pre_full = c.Data.(SORT)(M(1)).Full_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,1);%Full
C_Pre_early = c.Data.(SORT)(M(1)).Early_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,2);%early
C_Pre_late = c.Data.(SORT)(M(1)).Late_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,3);%late

L_Pre_full = l.Data.(SORT)(M(1)).Full_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,1);%Full
L_Pre_early = l.Data.(SORT)(M(1)).Early_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,2);%early
L_Pre_late = l.Data.(SORT)(M(1)).Late_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,3);%late

Y_Pre_full = y.Data.(SORT)(M(1)).Full_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,1);%Full
Y_Pre_early = y.Data.(SORT)(M(1)).Early_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,2);%early
Y_Pre_late = y.Data.(SORT)(M(1)).Late_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,3);%late

%Laser
C_L_full = c.Data.(SORT)(M(2)).Full_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,1);%Full
C_L_early = c.Data.(SORT)(M(2)).Early_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,2);%early
C_L_late = c.Data.(SORT)(M(2)).Late_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,3);%late

L_L_full = l.Data.(SORT)(M(2)).Full_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,1);%Full
L_L_early = l.Data.(SORT)(M(2)).Early_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,2);%early
L_L_late = l.Data.(SORT)(M(2)).Late_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,3);%late

Y_L_full = y.Data.(SORT)(M(2)).Full_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,1);%Full
Y_L_early = y.Data.(SORT)(M(2)).Early_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,2);%early
Y_L_late = y.Data.(SORT)(M(2)).Late_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,3);%late


%Post
C_Post_full = c.Data.(SORT)(M(3)).Full_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,1);%Full
C_Post_early = c.Data.(SORT)(M(3)).Early_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,2);%early
C_Post_late = c.Data.(SORT)(M(3)).Late_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,3);%late

L_Post_full = l.Data.(SORT)(M(3)).Full_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,1);%Full
L_Post_early = l.Data.(SORT)(M(3)).Early_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,2);%early
L_Post_late = l.Data.(SORT)(M(3)).Late_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,3);%late

Y_Post_full = y.Data.(SORT)(M(3)).Full_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,1);%Full
Y_Post_early = y.Data.(SORT)(M(3)).Early_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,2);%early
Y_Post_late = y.Data.(SORT)(M(3)).Late_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,3);%late

figure
subplot(1,3,1)
bar([1 4 7],[nanmean(C_Pre_full) nanmean(C_L_full) nanmean(C_Post_full)],0.3)
hold on
bar([2 5 8],[nanmean(L_Pre_full) nanmean(L_L_full) nanmean(L_Post_full)],0.3)
bar([3 6 9],[nanmean(Y_Pre_full) nanmean(Y_L_full) nanmean(Y_Post_full)],0.3)

errorbar([1 4 7],[nanmean(C_Pre_full) nanmean(C_L_full) nanmean(C_Post_full)],...
    [nanstd(C_Pre_full)/sqrt(sum(~isnan(C_Pre_full)))...
    nanstd(C_L_full)/sqrt(sum(~isnan(C_L_full)))...
    nanstd(C_Post_full)/sqrt(sum(~isnan(C_Post_full)))],'.', 'LineWidth', 2)

errorbar([2 5 8],[nanmean(L_Pre_full) nanmean(L_L_full) nanmean(L_Post_full)],...
    [nanstd(L_Pre_full)/sqrt(sum(~isnan(L_Pre_full)))...
    nanstd(L_L_full)/sqrt(sum(~isnan(L_L_full)))...
    nanstd(L_Post_full)/sqrt(sum(~isnan(L_Post_full)))],'.', 'LineWidth', 2)

errorbar([3 6 9],[nanmean(Y_Pre_full) nanmean(Y_L_full) nanmean(Y_Post_full)],...
    [nanstd(Y_Pre_full)/sqrt(sum(~isnan(Y_Pre_full)))...
    nanstd(Y_L_full)/sqrt(sum(~isnan(Y_L_full)))...
    nanstd(Y_Post_full)/sqrt(sum(~isnan(Y_Post_full)))],'.', 'LineWidth', 2)

xticks([2 5 8])
xticklabels({'Pre','Opto','Post'})
ylim([0.75 2])
set(gca,'fontsize',12,'FontWeight','bold')
title ('Full RMS RELRES')
legend('Control','Opto','YFP','Location','best')

subplot(1,3,2)
bar([1 4 7],[nanmean(C_Pre_early) nanmean(C_L_early) nanmean(C_Post_early)],0.3)
hold on
bar([2 5 8],[nanmean(L_Pre_early) nanmean(L_L_early) nanmean(L_Post_early)],0.3)
bar([3 6 9],[nanmean(Y_Pre_early) nanmean(Y_L_early) nanmean(Y_Post_early)],0.3)

errorbar([1 4 7],[nanmean(C_Pre_early) nanmean(C_L_early) nanmean(C_Post_early)],...
    [nanstd(C_Pre_early)/sqrt(sum(~isnan(C_Pre_early)))...
    nanstd(C_L_early)/sqrt(sum(~isnan(C_L_early)))...
    nanstd(C_Post_early)/sqrt(sum(~isnan(C_Post_early)))],'.', 'LineWidth', 2)

errorbar([2 5 8],[nanmean(L_Pre_early) nanmean(L_L_early) nanmean(L_Post_early)],...
    [nanstd(L_Pre_early)/sqrt(sum(~isnan(L_Pre_early)))...
    nanstd(L_L_early)/sqrt(sum(~isnan(L_L_early)))...
    nanstd(L_Post_early)/sqrt(sum(~isnan(L_Post_early)))],'.', 'LineWidth', 2)

errorbar([3 6 9],[nanmean(Y_Pre_early) nanmean(Y_L_early) nanmean(Y_Post_early)],...
    [nanstd(Y_Pre_early)/sqrt(sum(~isnan(Y_Pre_early)))...
    nanstd(Y_L_early)/sqrt(sum(~isnan(Y_L_early)))...
    nanstd(Y_Post_early)/sqrt(sum(~isnan(Y_Post_early)))],'.', 'LineWidth', 2)

xticks([2 5 8])
xticklabels({'Pre','Opto','Post'})
ylim([0.75 2])
set(gca,'fontsize',12,'FontWeight','bold')
title ('Early RMS RELRES')
legend('Control','Opto','YFP','Location','best')

subplot(1,3,3)
bar([1 4 7],[nanmean(C_Pre_late) nanmean(C_L_late) nanmean(C_Post_late)],0.3)
hold on
bar([2 5 8],[nanmean(L_Pre_late) nanmean(L_L_late) nanmean(L_Post_late)],0.3)
bar([3 6 9],[nanmean(Y_Pre_late) nanmean(Y_L_late) nanmean(Y_Post_late)],0.3)

errorbar([1 4 7],[nanmean(C_Pre_late) nanmean(C_L_late) nanmean(C_Post_late)],...
    [nanstd(C_Pre_late)/sqrt(sum(~isnan(C_Pre_late)))...
    nanstd(C_L_late)/sqrt(sum(~isnan(C_L_late)))...
    nanstd(C_Post_late)/sqrt(sum(~isnan(C_Post_late)))],'.', 'LineWidth', 2)

errorbar([2 5 8],[nanmean(L_Pre_late) nanmean(L_L_late) nanmean(L_Post_late)],...
    [nanstd(L_Pre_late)/sqrt(sum(~isnan(L_Pre_late)))...
    nanstd(L_L_late)/sqrt(sum(~isnan(L_L_late)))...
    nanstd(L_Post_late)/sqrt(sum(~isnan(L_Post_late)))],'.', 'LineWidth', 2)

errorbar([3 6 9],[nanmean(Y_Pre_late) nanmean(Y_L_late) nanmean(Y_Post_late)],...
    [nanstd(Y_Pre_late)/sqrt(sum(~isnan(Y_Pre_late)))...
    nanstd(Y_L_late)/sqrt(sum(~isnan(Y_L_late)))...
    nanstd(Y_Post_late)/sqrt(sum(~isnan(Y_Post_late)))],'.', 'LineWidth', 2)

xticks([2 5 8])
xticklabels({'Pre','Opto','Post'})
ylim([0.75 2])
set(gca,'fontsize',12,'FontWeight','bold')
title ('Late RMS RELRES')
legend('Control','Opto','YFP','Location','best')

clear C_Pre_full  C_Pre_early C_Pre_late



%%% Lippi
CNorm = nan(length(c.Data.names),3);
for i1 = 1:3
    CNorm(:,i1,1) =  c.Data.(SORT)(i1).Full_RMS_AVREC(:,BF_C+BFlim);%full
    CNorm(:,i1,2) =  c.Data.(SORT)(i1).Early_RMS_AVREC(:,BF_C+BFlim);%early
    CNorm(:,i1,3) =  c.Data.(SORT)(i1).Late_RMS_AVREC(:,BF_C+BFlim);%late
end
CNorm = nanmean(CNorm,2);

LNorm = nan(length(l.Data.names),3);
for i1 = 1:3
    LNorm(:,i1,1) =  l.Data.(SORT)(i1).Full_RMS_AVREC(:,BF_L+BFlim);%full
    LNorm(:,i1,2) =  l.Data.(SORT)(i1).Early_RMS_AVREC(:,BF_L+BFlim);%early
    LNorm(:,i1,3) =  l.Data.(SORT)(i1).Late_RMS_AVREC(:,BF_L+BFlim);%late
end
LNorm = nanmean(LNorm,2);

YNorm = nan(length(y.Data.names),3);
for i1 = 1:3
    YNorm(:,i1,1) =  y.Data.(SORT)(i1).Full_RMS_AVREC(:,BF_Y+BFlim);%full
    YNorm(:,i1,2) =  y.Data.(SORT)(i1).Early_RMS_AVREC(:,BF_Y+BFlim);%early
    YNorm(:,i1,3) =  y.Data.(SORT)(i1).Late_RMS_AVREC(:,BF_Y+BFlim);%late
end

YNorm = nanmean(YNorm,2);
%Pre 1-3
C_Pre_full(:,:,1) = c.Data.(SORT)(1).Full_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,1);%Full
C_Pre_full(:,:,2) = c.Data.(SORT)(2).Full_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,1);
C_Pre_full(:,:,3) = c.Data.(SORT)(3).Full_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,1);
C_Pre_full = nanmean(nanmean(C_Pre_full,3),2);

C_Pre_early (:,:,1) = c.Data.(SORT)(1).Early_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,2);%Full
C_Pre_early (:,:,2) = c.Data.(SORT)(2).Early_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,2);
C_Pre_early (:,:,3) = c.Data.(SORT)(3).Early_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,2);
C_Pre_early  = nanmean(C_Pre_early ,3);

C_Pre_late (:,:,1) = c.Data.(SORT)(1).Late_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,3);%Full
C_Pre_late (:,:,2) = c.Data.(SORT)(2).Late_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,3);
C_Pre_late (:,:,3) = c.Data.(SORT)(3).Late_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,3);
C_Pre_late  = nanmean(C_Pre_late ,3);

L_Pre_full(:,:,1) = l.Data.(SORT)(1).Full_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,1);%Full
L_Pre_full(:,:,2) = l.Data.(SORT)(2).Full_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,1);
L_Pre_full(:,:,3) = l.Data.(SORT)(3).Full_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,1);
L_Pre_full = nanmean(L_Pre_full,3);

L_Pre_early (:,:,1) = l.Data.(SORT)(1).Early_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,2);%Full
L_Pre_early (:,:,2) = l.Data.(SORT)(2).Early_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,2);
L_Pre_early (:,:,3) = l.Data.(SORT)(3).Early_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,2);
L_Pre_early  = nanmean(L_Pre_early ,3);

L_Pre_late (:,:,1) = l.Data.(SORT)(1).Late_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,3);%Full
L_Pre_late (:,:,2) = l.Data.(SORT)(2).Late_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,3);
L_Pre_late (:,:,3) = l.Data.(SORT)(3).Late_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,3);
L_Pre_late  = nanmean(L_Pre_late ,3);

Y_Pre_full(:,:,1) = y.Data.(SORT)(1).Full_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,1);%Full
Y_Pre_full(:,:,2) = y.Data.(SORT)(2).Full_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,1);
Y_Pre_full(:,:,3) = y.Data.(SORT)(3).Full_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,1);
Y_Pre_full = nanmean(Y_Pre_full,3);

Y_Pre_early (:,:,1) = y.Data.(SORT)(1).Early_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,2);%Full
Y_Pre_early (:,:,2) = y.Data.(SORT)(2).Early_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,2);
Y_Pre_early (:,:,3) = y.Data.(SORT)(3).Early_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,2);
Y_Pre_early  = nanmean(Y_Pre_early ,3);

Y_Pre_late (:,:,1) = y.Data.(SORT)(1).Late_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,3);%Full
Y_Pre_late (:,:,2) = y.Data.(SORT)(2).Late_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,3);
Y_Pre_late (:,:,3) = y.Data.(SORT)(3).Late_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,3);
Y_Pre_late  = nanmean(Y_Pre_late ,3);

%Laser
C_L_full = c.Data.(SORT)(4).Full_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,1);%Full
C_L_early = c.Data.(SORT)(4).Early_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,2);%early
C_L_late = c.Data.(SORT)(4).Late_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,3);%late

L_L_full = l.Data.(SORT)(4).Full_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,1);%Full
L_L_early = l.Data.(SORT)(4).Early_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,2);%early
L_L_late = l.Data.(SORT)(4).Late_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,3);%late

Y_L_full = y.Data.(SORT)(4).Full_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,1);%Full
Y_L_early = y.Data.(SORT)(4).Early_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,2);%early
Y_L_late = y.Data.(SORT)(4).Late_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,3);%late


%Post 1-3
C_Post1_full(:,:,1) = c.Data.(SORT)(5).Full_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,1);%Full
C_Post1_full(:,:,2) = c.Data.(SORT)(6).Full_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,1);
C_Post1_full(:,:,3) = c.Data.(SORT)(7).Full_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,1);
C_Post1_full = nanmean(C_Post1_full,3);

C_Post1_early (:,:,1) = c.Data.(SORT)(5).Early_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,2);%Full
C_Post1_early (:,:,2) = c.Data.(SORT)(6).Early_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,2);
C_Post1_early (:,:,3) = c.Data.(SORT)(7).Early_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,2);
C_Post1_early  = nanmean(C_Post1_early ,3);

C_Post1_late (:,:,1) = c.Data.(SORT)(5).Late_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,3);%Full
C_Post1_late (:,:,2) = c.Data.(SORT)(6).Late_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,3);
C_Post1_late (:,:,3) = c.Data.(SORT)(7).Late_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,3);
C_Post1_late  = nanmean(C_Post1_late ,3);

L_Post1_full(:,:,1) = l.Data.(SORT)(5).Full_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,1);%Full
L_Post1_full(:,:,2) = l.Data.(SORT)(6).Full_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,1);
L_Post1_full(:,:,3) = l.Data.(SORT)(7).Full_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,1);
L_Post1_full = nanmean(L_Post1_full,3);

L_Post1_early (:,:,1) = l.Data.(SORT)(5).Early_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,2);%Full
L_Post1_early (:,:,2) = l.Data.(SORT)(6).Early_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,2);
L_Post1_early (:,:,3) = l.Data.(SORT)(7).Early_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,2);
L_Post1_early  = nanmean(L_Post1_early ,3);

L_Post1_late (:,:,1) = l.Data.(SORT)(5).Late_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,3);%Full
L_Post1_late (:,:,2) = l.Data.(SORT)(6).Late_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,3);
L_Post1_late (:,:,3) = l.Data.(SORT)(7).Late_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,3);
L_Post1_late  = nanmean(L_Post1_late ,3);

Y_Post1_full(:,:,1) = y.Data.(SORT)(5).Full_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,1);%Full
Y_Post1_full(:,:,2) = y.Data.(SORT)(6).Full_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,1);
Y_Post1_full(:,:,3) = y.Data.(SORT)(7).Full_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,1);
Y_Post1_full = nanmean(Y_Post1_full,3);

Y_Post1_early (:,:,1) = y.Data.(SORT)(5).Early_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,2);%Full
Y_Post1_early (:,:,2) = y.Data.(SORT)(6).Early_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,2);
Y_Post1_early (:,:,3) = y.Data.(SORT)(7).Early_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,2);
Y_Post1_early  = nanmean(Y_Post1_early ,3);

Y_Post1_late (:,:,1) = y.Data.(SORT)(5).Late_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,3);%Full
Y_Post1_late (:,:,2) = y.Data.(SORT)(6).Late_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,3);
Y_Post1_late (:,:,3) = y.Data.(SORT)(7).Late_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,3);
Y_Post1_late  = nanmean(Y_Post1_late ,3);

%%% Post2
C_Post2_full(:,:,1) = c.Data.(SORT)(8).Full_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,1);%Full
C_Post2_full(:,:,2) = c.Data.(SORT)(9).Full_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,1);
C_Post2_full(:,:,3) = c.Data.(SORT)(10).Full_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,1);
C_Post2_full(:,:,4) = c.Data.(SORT)(11).Full_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,1);
C_Post2_full = nanmean(C_Post2_full,3);

C_Post2_early (:,:,1) = c.Data.(SORT)(8).Early_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,2);%Full
C_Post2_early (:,:,2) = c.Data.(SORT)(9).Early_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,2);
C_Post2_early (:,:,3) = c.Data.(SORT)(10).Early_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,2);
C_Post2_early (:,:,4) = c.Data.(SORT)(11).Early_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,2);
C_Post2_early  = nanmean(C_Post2_early ,3);

C_Post2_late (:,:,1) = c.Data.(SORT)(8).Late_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,3);%Full
C_Post2_late (:,:,2) = c.Data.(SORT)(9).Late_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,3);
C_Post2_late (:,:,3) = c.Data.(SORT)(10).Late_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,3);
C_Post2_late (:,:,4) = c.Data.(SORT)(11).Late_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,3);
C_Post2_late  = nanmean(C_Post2_late ,3);

L_Post2_full(:,:,1) = l.Data.(SORT)(8).Full_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,1);%Full
L_Post2_full(:,:,2) = l.Data.(SORT)(9).Full_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,1);
L_Post2_full(:,:,3) = l.Data.(SORT)(10).Full_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,1);
L_Post2_full(:,:,4) = l.Data.(SORT)(11).Full_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,1);
L_Post2_full = nanmean(L_Post2_full,3);

L_Post2_early (:,:,1) = l.Data.(SORT)(8).Early_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,2);%Full
L_Post2_early (:,:,2) = l.Data.(SORT)(9).Early_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,2);
L_Post2_early (:,:,3) = l.Data.(SORT)(10).Early_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,2);
L_Post2_early (:,:,4) = l.Data.(SORT)(11).Early_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,2);
L_Post2_early  = nanmean(L_Post2_early ,3);

L_Post2_late (:,:,1) = l.Data.(SORT)(8).Late_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,3);%Full
L_Post2_late (:,:,2) = l.Data.(SORT)(9).Late_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,3);
L_Post2_late (:,:,3) = l.Data.(SORT)(10).Late_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,3);
L_Post2_late (:,:,4) = l.Data.(SORT)(11).Late_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,3);
L_Post2_late  = nanmean(L_Post2_late ,3);

Y_Post2_full(:,:,1) = y.Data.(SORT)(8).Full_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,1);%Full
Y_Post2_full(:,:,2) = y.Data.(SORT)(9).Full_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,1);
Y_Post2_full(:,:,3) = y.Data.(SORT)(10).Full_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,1);
Y_Post2_full(:,:,4) = y.Data.(SORT)(11).Full_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,1);
Y_Post2_full = nanmean(Y_Post2_full,3);

Y_Post2_early (:,:,1) = y.Data.(SORT)(8).Early_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,2);%Full
Y_Post2_early (:,:,2) = y.Data.(SORT)(9).Early_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,2);
Y_Post2_early (:,:,3) = y.Data.(SORT)(10).Early_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,2);
Y_Post2_early (:,:,4) = y.Data.(SORT)(11).Early_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,2);
Y_Post2_early  = nanmean(Y_Post2_early ,3);

Y_Post2_late (:,:,1) = y.Data.(SORT)(8).Late_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,3);%Full
Y_Post2_late (:,:,2) = y.Data.(SORT)(9).Late_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,3);
Y_Post2_late (:,:,3) = y.Data.(SORT)(10).Late_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,3);
Y_Post2_late (:,:,4) = y.Data.(SORT)(11).Late_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,3);
Y_Post2_late  = nanmean(Y_Post2_late ,3);


%%%Post 9-11
% C_Post3_full(:,:,1) = c.Data.(SORT)(9).Full_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,1);%Full
% C_Post3_full(:,:,2) = c.Data.(SORT)(10).Full_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,1);
% C_Post3_full(:,:,3) = c.Data.(SORT)(11).Full_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,1);
% C_Post3_full = nanmean(C_Post3_full,3);
% 
% C_Post3_early (:,:,1) = c.Data.(SORT)(9).Early_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,2);%Full
% C_Post3_early (:,:,2) = c.Data.(SORT)(10).Early_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,2);
% C_Post3_early (:,:,3) = c.Data.(SORT)(11).Early_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,2);
% C_Post3_early  = nanmean(C_Post3_early ,3);
% 
% C_Post3_late (:,:,1) = c.Data.(SORT)(9).Late_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,3);%Full
% C_Post3_late (:,:,2) = c.Data.(SORT)(10).Late_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,3);
% C_Post3_late (:,:,3) = c.Data.(SORT)(11).Late_RMS_AVREC(:,BF_C+BFlim)./CNorm(:,:,3);
% C_Post3_late  = nanmean(C_Post3_late ,3);
% 
% L_Post3_full(:,:,1) = l.Data.(SORT)(9).Full_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,1);%Full
% L_Post3_full(:,:,2) = l.Data.(SORT)(10).Full_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,1);
% L_Post3_full(:,:,3) = l.Data.(SORT)(11).Full_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,1);
% L_Post3_full = nanmean(L_Post3_full,3);
% 
% L_Post3_early (:,:,1) = l.Data.(SORT)(9).Early_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,2);%Full
% L_Post3_early (:,:,2) = l.Data.(SORT)(10).Early_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,2);
% L_Post3_early (:,:,3) = l.Data.(SORT)(11).Early_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,2);
% L_Post3_early  = nanmean(L_Post3_early ,3);
% 
% L_Post3_late (:,:,1) = l.Data.(SORT)(9).Late_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,3);%Full
% L_Post3_late (:,:,2) = l.Data.(SORT)(10).Late_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,3);
% L_Post3_late (:,:,3) = l.Data.(SORT)(11).Late_RMS_AVREC(:,BF_L+BFlim)./LNorm(:,:,3);
% L_Post3_late  = nanmean(L_Post3_late ,3);
% 
% Y_Post3_full(:,:,1) = y.Data.(SORT)(9).Full_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,1);%Full
% Y_Post3_full(:,:,2) = y.Data.(SORT)(10).Full_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,1);
% Y_Post3_full(:,:,3) = y.Data.(SORT)(11).Full_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,1);
% Y_Post3_full = nanmean(Y_Post3_full,3);
% 
% Y_Post3_early (:,:,1) = y.Data.(SORT)(9).Early_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,2);%Full
% Y_Post3_early (:,:,2) = y.Data.(SORT)(10).Early_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,2);
% Y_Post3_early (:,:,3) = y.Data.(SORT)(11).Early_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,2);
% Y_Post3_early  = nanmean(Y_Post3_early ,3);
% 
% Y_Post3_late (:,:,1) = y.Data.(SORT)(9).Late_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,3);%Full
% Y_Post3_late (:,:,2) = y.Data.(SORT)(10).Late_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,3);
% Y_Post3_late (:,:,3) = y.Data.(SORT)(11).Late_RMS_AVREC(:,BF_Y+BFlim)./YNorm(:,:,3);
% Y_Post3_late  = nanmean(Y_Post3_late ,3);
% close all


% rmanova
clear Pre_O L_O Post1_O Post2_O C_O Opto_O
% L = [3 1];
% Pre_O = teg_repeated_measures_ANOVA([vertcat(C_Pre_full, nan(5,1)) L_Pre_full...
%     vertcat(Y_Pre_full, nan(8,1))],L,{'Group','Measurement'});
% L = [3 1];
% L_O = teg_repeated_measures_ANOVA([vertcat(C_L_full, nan(5,1)) L_L_full...
%     vertcat(Y_L_full, nan(8,1))],L,{'Group','Measurement'});
% L = [3 1];
% Post1_O = teg_repeated_measures_ANOVA([vertcat(C_Post1_full, nan(5,1)) L_Post1_full...
%     vertcat(Y_Post1_full, nan(8,1))],L,{'Group','Measurement'});
% L = [3 1];
% Post2_O = teg_repeated_measures_ANOVA([vertcat(C_Post2_full, nan(5,1)) L_Post2_full...
%     vertcat(Y_Post2_full, nan(8,1))],L,{'Group','Measurement'});
% 
% L = [1 4];
% C_O = teg_repeated_measures_ANOVA([C_Pre_full C_L_full C_Post1_full C_Post2_full],L,{'Group','Measurement'});
% L = [1 4];
% Opto_O = teg_repeated_measures_ANOVA([L_Pre_full L_L_full L_Post1_full L_Post2_full],L,{'Group','Measurement'});
% L = [1 4];
% Y_O = teg_repeated_measures_ANOVA([Y_Pre_full Y_L_full Y_Post1_full Y_Post2_full],L,{'Group','Measurement'});

%  keyboard
figure
subplot(2,3,1)
bar([1 4 7 10 ],[nanmean(C_Pre_full) nanmean(C_L_full) nanmean(C_Post1_full) nanmean(C_Post2_full) ],0.3)
hold on
bar([2 5 8 11 ],[nanmean(L_Pre_full) nanmean(L_L_full) nanmean(L_Post1_full) nanmean(L_Post2_full) ],0.3)
bar([3 6 9 12 ],[nanmean(Y_Pre_full) nanmean(Y_L_full) nanmean(Y_Post1_full) nanmean(Y_Post2_full) ],0.3)

errorbar([1 4 7 10 ],[nanmean(C_Pre_full) nanmean(C_L_full) nanmean(C_Post1_full) nanmean(C_Post2_full) ],...
    [nanstd(C_Pre_full)/sqrt(sum(~isnan(C_Pre_full)))...
    nanstd(C_L_full)/sqrt(sum(~isnan(C_L_full)))...
    nanstd(C_Post1_full)/sqrt(sum(~isnan(C_Post1_full)))...
    nanstd(C_Post2_full)/sqrt(sum(~isnan(C_Post2_full)))...
%     nanstd(C_Post3_full)/sqrt(sum(~isnan(C_Post3_full)))...
    ],'.', 'LineWidth', 2)

errorbar([2 5 8 11 ],[nanmean(L_Pre_full) nanmean(L_L_full) nanmean(L_Post1_full) nanmean(L_Post2_full) ],...
    [nanstd(L_Pre_full)/sqrt(sum(~isnan(L_Pre_full)))...
    nanstd(L_L_full)/sqrt(sum(~isnan(L_L_full)))...
    nanstd(L_Post1_full)/sqrt(sum(~isnan(L_Post1_full)))...
    nanstd(L_Post2_full)/sqrt(sum(~isnan(L_Post2_full)))...
%     nanstd(L_Post3_full)/sqrt(sum(~isnan(L_Post3_full)))...
    ],'.', 'LineWidth', 2)

errorbar([3 6 9 12 ],[nanmean(Y_Pre_full) nanmean(Y_L_full) nanmean(Y_Post1_full) nanmean(Y_Post2_full)],...
    [nanstd(Y_Pre_full)/sqrt(sum(~isnan(Y_Pre_full)))...
    nanstd(Y_L_full)/sqrt(sum(~isnan(Y_L_full)))...
    nanstd(Y_Post1_full)/sqrt(sum(~isnan(Y_Post1_full)))...
    nanstd(Y_Post2_full)/sqrt(sum(~isnan(Y_Post2_full)))...
%     nanstd(Y_Post3_full)/sqrt(sum(~isnan(Y_Post3_full)))...
    ],'.', 'LineWidth', 2)

xticks([2 5 8 11 14])
xticklabels({'Pre 1-3','Opto','Post 1-3','Post 4-7'})
ylim([0.8 1.4])
set(gca,'fontsize',12,'FontWeight','bold')
title ('Full RMS AVREC')
legend('Control','Opto','YFP','Location','best','AutoUpdate','off')


% rmanova
clear Pre_O L_O Post1_O Post2_O C_O Opto_O
% L = [3 1];
% Pre_O = teg_repeated_measures_ANOVA([vertcat(C_Pre_early, nan(5,1)) L_Pre_early...
%     vertcat(Y_Pre_early, nan(8,1))],L,{'Group','Measurement'});
% L = [3 1];
% L_O = teg_repeated_measures_ANOVA([vertcat(C_L_early, nan(5,1)) L_L_early...
%     vertcat(Y_L_early, nan(8,1))],L,{'Group','Measurement'});
% L = [3 1];
% Post1_O = teg_repeated_measures_ANOVA([vertcat(C_Post1_early, nan(5,1)) L_Post1_early...
%     vertcat(Y_Post1_early, nan(8,1))],L,{'Group','Measurement'});
% L = [3 1];
% Post2_O = teg_repeated_measures_ANOVA([vertcat(C_Post2_early, nan(5,1)) L_Post2_early...
%     vertcat(Y_Post2_early, nan(8,1))],L,{'Group','Measurement'});
% 
% L = [1 4];
% C_O = teg_repeated_measures_ANOVA([C_Pre_early C_L_early C_Post1_early C_Post2_early],L,{'Group','Measurement'});
% L = [1 4];
% Opto_O = teg_repeated_measures_ANOVA([L_Pre_early L_L_early L_Post1_early L_Post2_early],L,{'Group','Measurement'});
% L = [1 4];
% Y_O = teg_repeated_measures_ANOVA([Y_Pre_early Y_L_early Y_Post1_early Y_Post2_early],L,{'Group','Measurement'});

%  keyboard

subplot(2,3,2)
bar([1 4 7 10 ],[nanmean(C_Pre_early) nanmean(C_L_early) nanmean(C_Post1_early) nanmean(C_Post2_early) ],0.3)
hold on
bar([2 5 8 11 ],[nanmean(L_Pre_early) nanmean(L_L_early) nanmean(L_Post1_early) nanmean(L_Post2_early) ],0.3)
bar([3 6 9 12 ],[nanmean(Y_Pre_early) nanmean(Y_L_early) nanmean(Y_Post1_early) nanmean(Y_Post2_early) ],0.3)

errorbar([1 4 7 10 ],[nanmean(C_Pre_early) nanmean(C_L_early) nanmean(C_Post1_early) nanmean(C_Post2_early) ],...
    [nanstd(C_Pre_early)/sqrt(sum(~isnan(C_Pre_early)))...
    nanstd(C_L_early)/sqrt(sum(~isnan(C_L_early)))...
    nanstd(C_Post1_early)/sqrt(sum(~isnan(C_Post1_early)))...
    nanstd(C_Post2_early)/sqrt(sum(~isnan(C_Post2_early)))...
%     nanstd(C_Post3_early)/sqrt(sum(~isnan(C_Post3_early)))...
    ],'.', 'LineWidth', 2)

errorbar([2 5 8 11 ],[nanmean(L_Pre_early) nanmean(L_L_early) nanmean(L_Post1_early) nanmean(L_Post2_early) ],...
    [nanstd(L_Pre_early)/sqrt(sum(~isnan(L_Pre_early)))...
    nanstd(L_L_early)/sqrt(sum(~isnan(L_L_early)))...
    nanstd(L_Post1_early)/sqrt(sum(~isnan(L_Post1_early)))...
    nanstd(L_Post2_early)/sqrt(sum(~isnan(L_Post2_early)))...
%     nanstd(L_Post3_early)/sqrt(sum(~isnan(L_Post3_early)))...
    ],'.', 'LineWidth', 2)

errorbar([3 6 9 12 ],[nanmean(Y_Pre_early) nanmean(Y_L_early) nanmean(Y_Post1_early) nanmean(Y_Post2_early) ],...
    [nanstd(Y_Pre_early)/sqrt(sum(~isnan(Y_Pre_early)))...
    nanstd(Y_L_early)/sqrt(sum(~isnan(Y_L_early)))...
    nanstd(Y_Post1_early)/sqrt(sum(~isnan(Y_Post1_early)))...
    nanstd(Y_Post2_early)/sqrt(sum(~isnan(Y_Post2_early)))...
%     nanstd(Y_Post3_early)/sqrt(sum(~isnan(Y_Post3_early)))...
    ],'.', 'LineWidth', 2)
% legend('Control','Opto','YFP','Location','best','AutoUpdate','off')

% if isfield(Post1_O,'groups')
%         try
%             if Post1_O.groups(1:2) == [2 1]
%                 plot([7.05 7.95],[1.2 1.2],'LineWidth', 2,'Color','black')
%                 plot(7.5,1.23,'k*')
%             end
%         catch
%         end
% 
%         try 
%             if Post1_O.groups(3:4) == [3 2]
%                 plot([8.05 8.95],[1.2 1.2],'LineWidth', 2,'Color','black')
%                 plot(8.5,1.23,'k*')
%             end
%         catch
%         end
%     
% end

% if isfield(Post2_O,'groups')
%         try
%             if Post1_O.groups(1:2) == [2 1]
%                 plot([10.05 10.95],[1.2 1.2],'LineWidth', 2,'Color','black')
%                 plot(10.5,1.23,'k*')
%             end
%         catch
%         end
% 
%         try 
%             if Post2_O.groups(3:4) == [3 2]
%                 plot([11.05 11.95],[1.2 1.2],'LineWidth', 2,'Color','black')
%                 plot(11.5,1.23,'k*')
%             end
%         catch
%         end
%     
% end
% if isfield(Opto_O,'groups')
%         try
%             if Opto_O.groups(1:2) == [3 1]
%                 plot([2 2 8 8],[1.25 1.3 1.3 1.25],'LineWidth', 2,'Color','red')
%                 plot(5,1.33,'r*')
%             end
%         catch
%         end
% 
%    
% end
% 
% 
xticks([2 5 8 11 ])
xticklabels({'Pre 1-3','Opto','Post 1-3','Post 4-7'})
ylim([0.8 1.4])
set(gca,'fontsize',12,'FontWeight','bold')
title ('early RMS AVREC')

% clear Pre_O L_O Post1_O Post2_O C_O Opto_O
% L = [3 1];
% Pre_O = teg_repeated_measures_ANOVA([vertcat(C_Pre_late, nan(5,1)) L_Pre_late...
%     vertcat(Y_Pre_late, nan(8,1))],L,{'Group','Measurement'});
% L = [3 1];
% L_O = teg_repeated_measures_ANOVA([vertcat(C_L_late, nan(5,1)) L_L_late...
%     vertcat(Y_L_late, nan(8,1))],L,{'Group','Measurement'});
% L = [3 1];
% Post1_O = teg_repeated_measures_ANOVA([vertcat(C_Post1_late, nan(5,1)) L_Post1_late...
%     vertcat(Y_Post1_late, nan(8,1))],L,{'Group','Measurement'});
% L = [3 1];
% Post2_O = teg_repeated_measures_ANOVA([vertcat(C_Post2_late, nan(5,1)) L_Post2_late...
%     vertcat(Y_Post2_late, nan(8,1))],L,{'Group','Measurement'});
% 
% L = [1 4];
% C_O = teg_repeated_measures_ANOVA([C_Pre_late C_L_late C_Post1_late C_Post2_late],L,{'Group','Measurement'});
% L = [1 4];
% Opto_O = teg_repeated_measures_ANOVA([L_Pre_late L_L_late L_Post1_late L_Post2_late],L,{'Group','Measurement'});
% L = [1 4];
% Y_O = teg_repeated_measures_ANOVA([Y_Pre_late Y_L_late Y_Post1_late Y_Post2_late],L,{'Group','Measurement'});

%  keyboard

subplot(2,3,3)
bar([1 4 7 10 ],[nanmean(C_Pre_late) nanmean(C_L_late) nanmean(C_Post1_late) nanmean(C_Post2_late) ],0.3)
hold on
bar([2 5 8 11 ],[nanmean(L_Pre_late) nanmean(L_L_late) nanmean(L_Post1_late) nanmean(L_Post2_late) ],0.3)
bar([3 6 9 12 ],[nanmean(Y_Pre_late) nanmean(Y_L_late) nanmean(Y_Post1_late) nanmean(Y_Post2_late) ],0.3)

errorbar([1 4 7 10 ],[nanmean(C_Pre_late) nanmean(C_L_late) nanmean(C_Post1_late) nanmean(C_Post2_late) ],...
    [nanstd(C_Pre_late)/sqrt(sum(~isnan(C_Pre_late)))...
    nanstd(C_L_late)/sqrt(sum(~isnan(C_L_late)))...
    nanstd(C_Post1_late)/sqrt(sum(~isnan(C_Post1_late)))...
    nanstd(C_Post2_late)/sqrt(sum(~isnan(C_Post2_late)))...
%     nanstd(C_Post3_late)/sqrt(sum(~isnan(C_Post3_late)))...
    ],'.', 'LineWidth', 2)

errorbar([2 5 8 11 ],[nanmean(L_Pre_late) nanmean(L_L_late) nanmean(L_Post1_late) nanmean(L_Post2_late) ],...
    [nanstd(L_Pre_late)/sqrt(sum(~isnan(L_Pre_late)))...
    nanstd(L_L_late)/sqrt(sum(~isnan(L_L_late)))...
    nanstd(L_Post1_late)/sqrt(sum(~isnan(L_Post1_late)))...
    nanstd(L_Post2_late)/sqrt(sum(~isnan(L_Post2_late)))...
%     nanstd(L_Post3_late)/sqrt(sum(~isnan(L_Post3_late)))...
    ],'.', 'LineWidth', 2)

errorbar([3 6 9 12 ],[nanmean(Y_Pre_late) nanmean(Y_L_late) nanmean(Y_Post1_late) nanmean(Y_Post2_late) ],...
    [nanstd(Y_Pre_late)/sqrt(sum(~isnan(Y_Pre_late)))...
    nanstd(Y_L_late)/sqrt(sum(~isnan(Y_L_late)))...
    nanstd(Y_Post1_late)/sqrt(sum(~isnan(Y_Post1_late)))...
    nanstd(Y_Post2_late)/sqrt(sum(~isnan(Y_Post2_late)))...
%     nanstd(Y_Post3_late)/sqrt(sum(~isnan(Y_Post3_late)))...
    ],'.', 'LineWidth', 2)

xticks([2 5 8 11 ])
xticklabels({'Pre 1-3','Opto','Post 1-3','Post 4-7'})
ylim([0.8 1.4])
set(gca,'fontsize',12,'FontWeight','bold')
title ('late RMS AVREC')
% legend('Control','Opto','YFP','Location','best','AutoUpdate','off')


%%keyboard

%%%%%%%%%%%%RELRES
CNorm = nan(length(c.Data.names),3);
for i1 = 1:3
    CNorm(:,i1,1) =  c.Data.(SORT)(i1).Full_RMS_RELRES(:,BF_C+BFlim);%full
    CNorm(:,i1,2) =  c.Data.(SORT)(i1).Early_RMS_RELRES(:,BF_C+BFlim);%early
    CNorm(:,i1,3) =  c.Data.(SORT)(i1).Late_RMS_RELRES(:,BF_C+BFlim);%late
end
CNorm = nanmean(CNorm,2);

LNorm = nan(length(l.Data.names),3);
for i1 = 1:3
    LNorm(:,i1,1) =  l.Data.(SORT)(i1).Full_RMS_RELRES(:,BF_L+BFlim);%full
    LNorm(:,i1,2) =  l.Data.(SORT)(i1).Early_RMS_RELRES(:,BF_L+BFlim);%early
    LNorm(:,i1,3) =  l.Data.(SORT)(i1).Late_RMS_RELRES(:,BF_L+BFlim);%late
end
LNorm = nanmean(LNorm,2);

YNorm = nan(length(y.Data.names),3);
for i1 = 1:3
    YNorm(:,i1,1) =  y.Data.(SORT)(i1).Full_RMS_RELRES(:,BF_Y+BFlim);%full
    YNorm(:,i1,2) =  y.Data.(SORT)(i1).Early_RMS_RELRES(:,BF_Y+BFlim);%early
    YNorm(:,i1,3) =  y.Data.(SORT)(i1).Late_RMS_RELRES(:,BF_Y+BFlim);%late
end

YNorm = nanmean(YNorm,2);
%Pre 1-3
C_Pre_full(:,:,1) = c.Data.(SORT)(1).Full_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,1);%Full
C_Pre_full(:,:,2) = c.Data.(SORT)(2).Full_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,1);
C_Pre_full(:,:,3) = c.Data.(SORT)(3).Full_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,1);
C_Pre_full = nanmean(nanmean(C_Pre_full,3),2);

C_Pre_early (:,:,1) = c.Data.(SORT)(1).Early_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,2);%Full
C_Pre_early (:,:,2) = c.Data.(SORT)(2).Early_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,2);
C_Pre_early (:,:,3) = c.Data.(SORT)(3).Early_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,2);
C_Pre_early  = nanmean(C_Pre_early ,3);

C_Pre_late (:,:,1) = c.Data.(SORT)(1).Late_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,3);%Full
C_Pre_late (:,:,2) = c.Data.(SORT)(2).Late_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,3);
C_Pre_late (:,:,3) = c.Data.(SORT)(3).Late_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,3);
C_Pre_late  = nanmean(C_Pre_late ,3);

L_Pre_full(:,:,1) = l.Data.(SORT)(1).Full_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,1);%Full
L_Pre_full(:,:,2) = l.Data.(SORT)(2).Full_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,1);
L_Pre_full(:,:,3) = l.Data.(SORT)(3).Full_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,1);
L_Pre_full = nanmean(L_Pre_full,3);

L_Pre_early (:,:,1) = l.Data.(SORT)(1).Early_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,2);%Full
L_Pre_early (:,:,2) = l.Data.(SORT)(2).Early_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,2);
L_Pre_early (:,:,3) = l.Data.(SORT)(3).Early_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,2);
L_Pre_early  = nanmean(L_Pre_early ,3);

L_Pre_late (:,:,1) = l.Data.(SORT)(1).Late_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,3);%Full
L_Pre_late (:,:,2) = l.Data.(SORT)(2).Late_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,3);
L_Pre_late (:,:,3) = l.Data.(SORT)(3).Late_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,3);
L_Pre_late  = nanmean(L_Pre_late ,3);

Y_Pre_full(:,:,1) = y.Data.(SORT)(1).Full_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,1);%Full
Y_Pre_full(:,:,2) = y.Data.(SORT)(2).Full_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,1);
Y_Pre_full(:,:,3) = y.Data.(SORT)(3).Full_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,1);
Y_Pre_full = nanmean(Y_Pre_full,3);

Y_Pre_early (:,:,1) = y.Data.(SORT)(1).Early_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,2);%Full
Y_Pre_early (:,:,2) = y.Data.(SORT)(2).Early_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,2);
Y_Pre_early (:,:,3) = y.Data.(SORT)(3).Early_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,2);
Y_Pre_early  = nanmean(Y_Pre_early ,3);

Y_Pre_late (:,:,1) = y.Data.(SORT)(1).Late_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,3);%Full
Y_Pre_late (:,:,2) = y.Data.(SORT)(2).Late_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,3);
Y_Pre_late (:,:,3) = y.Data.(SORT)(3).Late_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,3);
Y_Pre_late  = nanmean(Y_Pre_late ,3);

%Laser
C_L_full = c.Data.(SORT)(4).Full_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,1);%Full
C_L_early = c.Data.(SORT)(4).Early_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,2);%early
C_L_late = c.Data.(SORT)(4).Late_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,3);%late

L_L_full = l.Data.(SORT)(4).Full_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,1);%Full
L_L_early = l.Data.(SORT)(4).Early_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,2);%early
L_L_late = l.Data.(SORT)(4).Late_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,3);%late

Y_L_full = y.Data.(SORT)(4).Full_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,1);%Full
Y_L_early = y.Data.(SORT)(4).Early_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,2);%early
Y_L_late = y.Data.(SORT)(4).Late_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,3);%late


%Post 1-3
C_Post1_full(:,:,1) = c.Data.(SORT)(5).Full_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,1);%Full
C_Post1_full(:,:,2) = c.Data.(SORT)(6).Full_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,1);
C_Post1_full(:,:,3) = c.Data.(SORT)(7).Full_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,1);
C_Post1_full = nanmean(C_Post1_full,3);

C_Post1_early (:,:,1) = c.Data.(SORT)(5).Early_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,2);%Full
C_Post1_early (:,:,2) = c.Data.(SORT)(6).Early_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,2);
C_Post1_early (:,:,3) = c.Data.(SORT)(7).Early_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,2);
C_Post1_early  = nanmean(C_Post1_early ,3);

C_Post1_late (:,:,1) = c.Data.(SORT)(5).Late_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,3);%Full
C_Post1_late (:,:,2) = c.Data.(SORT)(6).Late_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,3);
C_Post1_late (:,:,3) = c.Data.(SORT)(7).Late_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,3);
C_Post1_late  = nanmean(C_Post1_late ,3);

L_Post1_full(:,:,1) = l.Data.(SORT)(5).Full_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,1);%Full
L_Post1_full(:,:,2) = l.Data.(SORT)(6).Full_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,1);
L_Post1_full(:,:,3) = l.Data.(SORT)(7).Full_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,1);
L_Post1_full = nanmean(L_Post1_full,3);

L_Post1_early (:,:,1) = l.Data.(SORT)(5).Early_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,2);%Full
L_Post1_early (:,:,2) = l.Data.(SORT)(6).Early_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,2);
L_Post1_early (:,:,3) = l.Data.(SORT)(7).Early_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,2);
L_Post1_early  = nanmean(L_Post1_early ,3);

L_Post1_late (:,:,1) = l.Data.(SORT)(5).Late_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,3);%Full
L_Post1_late (:,:,2) = l.Data.(SORT)(6).Late_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,3);
L_Post1_late (:,:,3) = l.Data.(SORT)(7).Late_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,3);
L_Post1_late  = nanmean(L_Post1_late ,3);

Y_Post1_full(:,:,1) = y.Data.(SORT)(5).Full_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,1);%Full
Y_Post1_full(:,:,2) = y.Data.(SORT)(6).Full_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,1);
Y_Post1_full(:,:,3) = y.Data.(SORT)(7).Full_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,1);
Y_Post1_full = nanmean(Y_Post1_full,3);

Y_Post1_early (:,:,1) = y.Data.(SORT)(5).Early_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,2);%Full
Y_Post1_early (:,:,2) = y.Data.(SORT)(6).Early_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,2);
Y_Post1_early (:,:,3) = y.Data.(SORT)(7).Early_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,2);
Y_Post1_early  = nanmean(Y_Post1_early ,3);

Y_Post1_late (:,:,1) = y.Data.(SORT)(5).Late_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,3);%Full
Y_Post1_late (:,:,2) = y.Data.(SORT)(6).Late_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,3);
Y_Post1_late (:,:,3) = y.Data.(SORT)(7).Late_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,3);
Y_Post1_late  = nanmean(Y_Post1_late ,3);

%%% Post2
C_Post2_full(:,:,1) = c.Data.(SORT)(8).Full_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,1);%Full
C_Post2_full(:,:,2) = c.Data.(SORT)(9).Full_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,1);
C_Post2_full(:,:,3) = c.Data.(SORT)(10).Full_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,1);
C_Post2_full(:,:,4) = c.Data.(SORT)(11).Full_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,1);
C_Post2_full = nanmean(C_Post2_full,3);

C_Post2_early (:,:,1) = c.Data.(SORT)(8).Early_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,2);%Full
C_Post2_early (:,:,2) = c.Data.(SORT)(9).Early_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,2);
C_Post2_early (:,:,3) = c.Data.(SORT)(10).Early_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,2);
C_Post2_early (:,:,4) = c.Data.(SORT)(11).Early_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,2);
C_Post2_early  = nanmean(C_Post2_early ,3);

C_Post2_late (:,:,1) = c.Data.(SORT)(8).Late_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,3);%Full
C_Post2_late (:,:,2) = c.Data.(SORT)(9).Late_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,3);
C_Post2_late (:,:,3) = c.Data.(SORT)(10).Late_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,3);
C_Post2_late (:,:,4) = c.Data.(SORT)(11).Late_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,3);
C_Post2_late  = nanmean(C_Post2_late ,3);

L_Post2_full(:,:,1) = l.Data.(SORT)(8).Full_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,1);%Full
L_Post2_full(:,:,2) = l.Data.(SORT)(9).Full_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,1);
L_Post2_full(:,:,3) = l.Data.(SORT)(10).Full_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,1);
L_Post2_full(:,:,4) = l.Data.(SORT)(11).Full_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,1);
L_Post2_full = nanmean(L_Post2_full,3);

L_Post2_early (:,:,1) = l.Data.(SORT)(8).Early_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,2);%Full
L_Post2_early (:,:,2) = l.Data.(SORT)(9).Early_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,2);
L_Post2_early (:,:,3) = l.Data.(SORT)(10).Early_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,2);
L_Post2_early (:,:,4) = l.Data.(SORT)(11).Early_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,2);
L_Post2_early  = nanmean(L_Post2_early ,3);

L_Post2_late (:,:,1) = l.Data.(SORT)(8).Late_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,3);%Full
L_Post2_late (:,:,2) = l.Data.(SORT)(9).Late_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,3);
L_Post2_late (:,:,3) = l.Data.(SORT)(10).Late_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,3);
L_Post2_late (:,:,4) = l.Data.(SORT)(11).Late_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,3);
L_Post2_late  = nanmean(L_Post2_late ,3);

Y_Post2_full(:,:,1) = y.Data.(SORT)(8).Full_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,1);%Full
Y_Post2_full(:,:,2) = y.Data.(SORT)(9).Full_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,1);
Y_Post2_full(:,:,3) = y.Data.(SORT)(10).Full_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,1);
Y_Post2_full(:,:,4) = y.Data.(SORT)(11).Full_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,1);
Y_Post2_full = nanmean(Y_Post2_full,3);

Y_Post2_early (:,:,1) = y.Data.(SORT)(8).Early_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,2);%Full
Y_Post2_early (:,:,2) = y.Data.(SORT)(9).Early_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,2);
Y_Post2_early (:,:,3) = y.Data.(SORT)(10).Early_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,2);
Y_Post2_early (:,:,4) = y.Data.(SORT)(11).Early_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,2);
Y_Post2_early  = nanmean(Y_Post2_early ,3);

Y_Post2_late (:,:,1) = y.Data.(SORT)(8).Late_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,3);%Full
Y_Post2_late (:,:,2) = y.Data.(SORT)(9).Late_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,3);
Y_Post2_late (:,:,3) = y.Data.(SORT)(10).Late_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,3);
Y_Post2_late (:,:,4) = y.Data.(SORT)(11).Late_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,3);
Y_Post2_late  = nanmean(Y_Post2_late ,3);


% %%%Post 9-11
% C_Post3_full(:,:,1) = c.Data.(SORT)(9).Full_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,1);%Full
% C_Post3_full(:,:,2) = c.Data.(SORT)(10).Full_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,1);
% C_Post3_full(:,:,3) = c.Data.(SORT)(11).Full_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,1);
% C_Post3_full = nanmean(C_Post3_full,3);
% 
% C_Post3_early (:,:,1) = c.Data.(SORT)(9).Early_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,2);%Full
% C_Post3_early (:,:,2) = c.Data.(SORT)(10).Early_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,2);
% C_Post3_early (:,:,3) = c.Data.(SORT)(11).Early_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,2);
% C_Post3_early  = nanmean(C_Post3_early ,3);
% 
% C_Post3_late (:,:,1) = c.Data.(SORT)(9).Late_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,3);%Full
% C_Post3_late (:,:,2) = c.Data.(SORT)(10).Late_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,3);
% C_Post3_late (:,:,3) = c.Data.(SORT)(11).Late_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,3);
% C_Post3_late  = nanmean(C_Post3_late ,3);
% 
% L_Post3_full(:,:,1) = l.Data.(SORT)(9).Full_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,1);%Full
% L_Post3_full(:,:,2) = l.Data.(SORT)(10).Full_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,1);
% L_Post3_full(:,:,3) = l.Data.(SORT)(11).Full_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,1);
% L_Post3_full = nanmean(L_Post3_full,3);
% 
% L_Post3_early (:,:,1) = l.Data.(SORT)(9).Early_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,2);%Full
% L_Post3_early (:,:,2) = l.Data.(SORT)(10).Early_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,2);
% L_Post3_early (:,:,3) = l.Data.(SORT)(11).Early_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,2);
% L_Post3_early  = nanmean(L_Post3_early ,3);
% 
% L_Post3_late (:,:,1) = l.Data.(SORT)(9).Late_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,3);%Full
% L_Post3_late (:,:,2) = l.Data.(SORT)(10).Late_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,3);
% L_Post3_late (:,:,3) = l.Data.(SORT)(11).Late_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,3);
% L_Post3_late  = nanmean(L_Post3_late ,3);
% 
% Y_Post3_full(:,:,1) = y.Data.(SORT)(9).Full_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,1);%Full
% Y_Post3_full(:,:,2) = y.Data.(SORT)(10).Full_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,1);
% Y_Post3_full(:,:,3) = y.Data.(SORT)(11).Full_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,1);
% Y_Post3_full = nanmean(Y_Post3_full,3);
% 
% Y_Post3_early (:,:,1) = y.Data.(SORT)(9).Early_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,2);%Full
% Y_Post3_early (:,:,2) = y.Data.(SORT)(10).Early_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,2);
% Y_Post3_early (:,:,3) = y.Data.(SORT)(11).Early_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,2);
% Y_Post3_early  = nanmean(Y_Post3_early ,3);
% 
% Y_Post3_late (:,:,1) = y.Data.(SORT)(9).Late_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,3);%Full
% Y_Post3_late (:,:,2) = y.Data.(SORT)(10).Late_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,3);
% Y_Post3_late (:,:,3) = y.Data.(SORT)(11).Late_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,3);
% Y_Post3_late  = nanmean(Y_Post3_late ,3);
clear Pre_O L_O Post1_O Post2_O C_O Opto_O
% L = [3 1];
% Pre_O = teg_repeated_measures_ANOVA([vertcat(C_Pre_full, nan(5,1)) L_Pre_full...
%     vertcat(Y_Pre_full, nan(8,1))],L,{'Group','Measurement'});
% L = [3 1];
% L_O = teg_repeated_measures_ANOVA([vertcat(C_L_full, nan(5,1)) L_L_full...
%     vertcat(Y_L_full, nan(8,1))],L,{'Group','Measurement'});
% L = [3 1];
% Post1_O = teg_repeated_measures_ANOVA([vertcat(C_Post1_full, nan(5,1)) L_Post1_full...
%     vertcat(Y_Post1_full, nan(8,1))],L,{'Group','Measurement'});
% L = [3 1];
% Post2_O = teg_repeated_measures_ANOVA([vertcat(C_Post2_full, nan(5,1)) L_Post2_full...
%     vertcat(Y_Post2_full, nan(8,1))],L,{'Group','Measurement'});
% 
% L = [1 4];
% C_O = teg_repeated_measures_ANOVA([C_Pre_full C_L_full C_Post1_full C_Post2_full],L,{'Group','Measurement'});
% L = [1 4];
% Opto_O = teg_repeated_measures_ANOVA([L_Pre_full L_L_full L_Post1_full L_Post2_full],L,{'Group','Measurement'});
% L = [1 4];
% Y_O = teg_repeated_measures_ANOVA([Y_Pre_full Y_L_full Y_Post1_full Y_Post2_full],L,{'Group','Measurement'});

%  keyboard

subplot(2,3,4)
bar([1 4 7 10 ],[nanmean(C_Pre_full) nanmean(C_L_full) nanmean(C_Post1_full) nanmean(C_Post2_full) ],0.3)
hold on
bar([2 5 8 11 ],[nanmean(L_Pre_full) nanmean(L_L_full) nanmean(L_Post1_full) nanmean(L_Post2_full) ],0.3)
bar([3 6 9 12 ],[nanmean(Y_Pre_full) nanmean(Y_L_full) nanmean(Y_Post1_full) nanmean(Y_Post2_full) ],0.3)

errorbar([1 4 7 10 ],[nanmean(C_Pre_full) nanmean(C_L_full) nanmean(C_Post1_full) nanmean(C_Post2_full) ],...
    [nanstd(C_Pre_full)/sqrt(sum(~isnan(C_Pre_full)))...
    nanstd(C_L_full)/sqrt(sum(~isnan(C_L_full)))...
    nanstd(C_Post1_full)/sqrt(sum(~isnan(C_Post1_full)))...
    nanstd(C_Post2_full)/sqrt(sum(~isnan(C_Post2_full)))...
%     nanstd(C_Post3_full)/sqrt(sum(~isnan(C_Post3_full)))...
    ],'.', 'LineWidth', 2)

errorbar([2 5 8 11 ],[nanmean(L_Pre_full) nanmean(L_L_full) nanmean(L_Post1_full) nanmean(L_Post2_full) ],...
    [nanstd(L_Pre_full)/sqrt(sum(~isnan(L_Pre_full)))...
    nanstd(L_L_full)/sqrt(sum(~isnan(L_L_full)))...
    nanstd(L_Post1_full)/sqrt(sum(~isnan(L_Post1_full)))...
    nanstd(L_Post2_full)/sqrt(sum(~isnan(L_Post2_full)))...
%     nanstd(L_Post3_full)/sqrt(sum(~isnan(L_Post3_full)))...
    ],'.', 'LineWidth', 2)

errorbar([3 6 9 12 ],[nanmean(Y_Pre_full) nanmean(Y_L_full) nanmean(Y_Post1_full) nanmean(Y_Post2_full) ],...
    [nanstd(Y_Pre_full)/sqrt(sum(~isnan(Y_Pre_full)))...
    nanstd(Y_L_full)/sqrt(sum(~isnan(Y_L_full)))...
    nanstd(Y_Post1_full)/sqrt(sum(~isnan(Y_Post1_full)))...
    nanstd(Y_Post2_full)/sqrt(sum(~isnan(Y_Post2_full)))...
%     nanstd(Y_Post3_full)/sqrt(sum(~isnan(Y_Post3_full)))...
    ],'.', 'LineWidth', 2)
% 
% if isfield(Opto_O,'groups')
%         try
%             if Opto_O.groups(1:2) == [2 1]
%                 plot([2 2 5 5],[1.5 1.55 1.55 1.5],'LineWidth', 2,'Color','red')
%                 plot(3.5,1.58,'r*')
%             end
%         catch
%         end
% 
%    
% end

xticks([2 5 8 11 ])
xticklabels({'Pre 1-3','Opto','Post 1-3','Post 4-7'})
ylim([0.75 1.65])
set(gca,'fontsize',12,'FontWeight','bold')
title ('Full RMS RELRES')
% legend('Control','Opto','YFP','Location','best')

clear Pre_O L_O Post1_O Post2_O C_O Opto_O
% L = [3 1];
% Pre_O = teg_repeated_measures_ANOVA([vertcat(C_Pre_early, nan(5,1)) L_Pre_early...
%     vertcat(Y_Pre_early, nan(8,1))],L,{'Group','Measurement'});
% L = [3 1];
% L_O = teg_repeated_measures_ANOVA([vertcat(C_L_early, nan(5,1)) L_L_early...
%     vertcat(Y_L_early, nan(8,1))],L,{'Group','Measurement'});
% L = [3 1];
% Post1_O = teg_repeated_measures_ANOVA([vertcat(C_Post1_early, nan(5,1)) L_Post1_early...
%     vertcat(Y_Post1_early, nan(8,1))],L,{'Group','Measurement'});
% L = [3 1];
% Post2_O = teg_repeated_measures_ANOVA([vertcat(C_Post2_early, nan(5,1)) L_Post2_early...
%     vertcat(Y_Post2_early, nan(8,1))],L,{'Group','Measurement'});
% 
% L = [1 4];
% C_O = teg_repeated_measures_ANOVA([C_Pre_early C_L_early C_Post1_early C_Post2_early],L,{'Group','Measurement'});
% L = [1 4];
% Opto_O = teg_repeated_measures_ANOVA([L_Pre_early L_L_early L_Post1_early L_Post2_early],L,{'Group','Measurement'});
% L = [1 4];
% Y_O = teg_repeated_measures_ANOVA([Y_Pre_early Y_L_early Y_Post1_early Y_Post2_early],L,{'Group','Measurement'});
% 
%  keyboard


subplot(2,3,5)
bar([1 4 7 10 ],[nanmean(C_Pre_early) nanmean(C_L_early) nanmean(C_Post1_early) nanmean(C_Post2_early) ],0.3)
hold on
bar([2 5 8 11 ],[nanmean(L_Pre_early) nanmean(L_L_early) nanmean(L_Post1_early) nanmean(L_Post2_early) ],0.3)
bar([3 6 9 12 ],[nanmean(Y_Pre_early) nanmean(Y_L_early) nanmean(Y_Post1_early) nanmean(Y_Post2_early) ],0.3)

errorbar([1 4 7 10 ],[nanmean(C_Pre_early) nanmean(C_L_early) nanmean(C_Post1_early) nanmean(C_Post2_early) ],...
    [nanstd(C_Pre_early)/sqrt(sum(~isnan(C_Pre_early)))...
    nanstd(C_L_early)/sqrt(sum(~isnan(C_L_early)))...
    nanstd(C_Post1_early)/sqrt(sum(~isnan(C_Post1_early)))...
    nanstd(C_Post2_early)/sqrt(sum(~isnan(C_Post2_early)))...
%     nanstd(C_Post3_early)/sqrt(sum(~isnan(C_Post3_early)))...
    ],'.', 'LineWidth', 2)

errorbar([2 5 8 11 ],[nanmean(L_Pre_early) nanmean(L_L_early) nanmean(L_Post1_early) nanmean(L_Post2_early) ],...
    [nanstd(L_Pre_early)/sqrt(sum(~isnan(L_Pre_early)))...
    nanstd(L_L_early)/sqrt(sum(~isnan(L_L_early)))...
    nanstd(L_Post1_early)/sqrt(sum(~isnan(L_Post1_early)))...
    nanstd(L_Post2_early)/sqrt(sum(~isnan(L_Post2_early)))...
%     nanstd(L_Post3_early)/sqrt(sum(~isnan(L_Post3_early)))...
    ],'.', 'LineWidth', 2)

errorbar([3 6 9 12 ],[nanmean(Y_Pre_early) nanmean(Y_L_early) nanmean(Y_Post1_early) nanmean(Y_Post2_early) ],...
    [nanstd(Y_Pre_early)/sqrt(sum(~isnan(Y_Pre_early)))...
    nanstd(Y_L_early)/sqrt(sum(~isnan(Y_L_early)))...
    nanstd(Y_Post1_early)/sqrt(sum(~isnan(Y_Post1_early)))...
    nanstd(Y_Post2_early)/sqrt(sum(~isnan(Y_Post2_early)))...
%     nanstd(Y_Post3_early)/sqrt(sum(~isnan(Y_Post3_early)))...
    ],'.', 'LineWidth', 2)

xticks([2 5 8 11 ])
xticklabels({'Pre 1-3','Opto','Post 1-3','Post 4-7'})
ylim([0.75 1.65])
set(gca,'fontsize',12,'FontWeight','bold')
title ('early RMS RELRES')
% legend('Control','Opto','YFP','Location','best')

clear Pre_O L_O Post1_O Post2_O C_O Opto_O
% L = [3 1];
% Pre_O = teg_repeated_measures_ANOVA([vertcat(C_Pre_late, nan(5,1)) L_Pre_late...
%     vertcat(Y_Pre_late, nan(8,1))],L,{'Group','Measurement'});
% L = [3 1];
% L_O = teg_repeated_measures_ANOVA([vertcat(C_L_late, nan(5,1)) L_L_late...
%     vertcat(Y_L_late, nan(8,1))],L,{'Group','Measurement'});
% L = [3 1];
% Post1_O = teg_repeated_measures_ANOVA([vertcat(C_Post1_late, nan(5,1)) L_Post1_late...
%     vertcat(Y_Post1_late, nan(8,1))],L,{'Group','Measurement'});
% L = [3 1];
% Post2_O = teg_repeated_measures_ANOVA([vertcat(C_Post2_late, nan(5,1)) L_Post2_late...
%     vertcat(Y_Post2_late, nan(8,1))],L,{'Group','Measurement'});
% 
% L = [1 4];
% C_O = teg_repeated_measures_ANOVA([C_Pre_late C_L_late C_Post1_late C_Post2_late],L,{'Group','Measurement'});
% L = [1 4];
% Opto_O = teg_repeated_measures_ANOVA([L_Pre_late L_L_late L_Post1_late L_Post2_late],L,{'Group','Measurement'});
% L = [1 4];
% Y_O = teg_repeated_measures_ANOVA([Y_Pre_late Y_L_late Y_Post1_late Y_Post2_late],L,{'Group','Measurement'});
% 
%  keyboard


subplot(2,3,6)
bar([1 4 7 10 ],[nanmean(C_Pre_late) nanmean(C_L_late) nanmean(C_Post1_late) nanmean(C_Post2_late) ],0.3)
hold on
bar([2 5 8 11 ],[nanmean(L_Pre_late) nanmean(L_L_late) nanmean(L_Post1_late) nanmean(L_Post2_late) ],0.3)
bar([3 6 9 12 ],[nanmean(Y_Pre_late) nanmean(Y_L_late) nanmean(Y_Post1_late) nanmean(Y_Post2_late) ],0.3)

errorbar([1 4 7 10 ],[nanmean(C_Pre_late) nanmean(C_L_late) nanmean(C_Post1_late) nanmean(C_Post2_late) ],...
    [nanstd(C_Pre_late)/sqrt(sum(~isnan(C_Pre_late)))...
    nanstd(C_L_late)/sqrt(sum(~isnan(C_L_late)))...
    nanstd(C_Post1_late)/sqrt(sum(~isnan(C_Post1_late)))...
    nanstd(C_Post2_late)/sqrt(sum(~isnan(C_Post2_late)))...
%     nanstd(C_Post3_late)/sqrt(sum(~isnan(C_Post3_late)))...
    ],'.', 'LineWidth', 2)

errorbar([2 5 8 11 ],[nanmean(L_Pre_late) nanmean(L_L_late) nanmean(L_Post1_late) nanmean(L_Post2_late) ],...
    [nanstd(L_Pre_late)/sqrt(sum(~isnan(L_Pre_late)))...
    nanstd(L_L_late)/sqrt(sum(~isnan(L_L_late)))...
    nanstd(L_Post1_late)/sqrt(sum(~isnan(L_Post1_late)))...
    nanstd(L_Post2_late)/sqrt(sum(~isnan(L_Post2_late)))...
%     nanstd(L_Post3_late)/sqrt(sum(~isnan(L_Post3_late)))...
    ],'.', 'LineWidth', 2)

errorbar([3 6 9 12 ],[nanmean(Y_Pre_late) nanmean(Y_L_late) nanmean(Y_Post1_late) nanmean(Y_Post2_late) ],...
    [nanstd(Y_Pre_late)/sqrt(sum(~isnan(Y_Pre_late)))...
    nanstd(Y_L_late)/sqrt(sum(~isnan(Y_L_late)))...
    nanstd(Y_Post1_late)/sqrt(sum(~isnan(Y_Post1_late)))...
    nanstd(Y_Post2_late)/sqrt(sum(~isnan(Y_Post2_late)))...
%     nanstd(Y_Post3_late)/sqrt(sum(~isnan(Y_Post3_late)))...
    ],'.', 'LineWidth', 2)

% 
% if isfield(Y_O,'groups')
%         try
%             if Y_O.groups(1:2) == [3 1]
%                 plot([3 3 9 9],[1.5 1.55 1.55 1.5],'LineWidth', 2,'Color','yellow')
%                 plot(6,1.58,'y*')
%             end
%         catch
%         end
% 
%    
% end




xticks([2 5 8 11 14])
xticklabels({'Pre 1-3','Opto','Post 1-3','Post 4-7'})
ylim([0.75 1.65])
set(gca,'fontsize',12,'FontWeight','bold')
title ('late RMS RELRES')
legend('Control','Opto','YFP','Location','best')



%%keyboard
%%%%%%%%%%%%RELRES
CNorm = nan(length(c.Data.names),3);
for i1 = 1:3
    CNorm(:,i1,1) =  c.Data.(SORT)(i1).Full_RMS_RELRES(:,BF_C+BFlim);%full
    CNorm(:,i1,2) =  c.Data.(SORT)(i1).Early_RMS_RELRES(:,BF_C+BFlim);%early
    CNorm(:,i1,3) =  c.Data.(SORT)(i1).Late_RMS_RELRES(:,BF_C+BFlim);%late
end
CNorm = nanmean(CNorm,2);

LNorm = nan(length(l.Data.names),3);
for i1 = 1:3
    LNorm(:,i1,1) =  l.Data.(SORT)(i1).Full_RMS_RELRES(:,BF_L+BFlim);%full
    LNorm(:,i1,2) =  l.Data.(SORT)(i1).Early_RMS_RELRES(:,BF_L+BFlim);%early
    LNorm(:,i1,3) =  l.Data.(SORT)(i1).Late_RMS_RELRES(:,BF_L+BFlim);%late
end
LNorm = nanmean(LNorm,2);

YNorm = nan(length(y.Data.names),3);
for i1 = 1:3
    YNorm(:,i1,1) =  y.Data.(SORT)(i1).Full_RMS_RELRES(:,BF_Y+BFlim);%full
    YNorm(:,i1,2) =  y.Data.(SORT)(i1).Early_RMS_RELRES(:,BF_Y+BFlim);%early
    YNorm(:,i1,3) =  y.Data.(SORT)(i1).Late_RMS_RELRES(:,BF_Y+BFlim);%late
end

YNorm = nanmean(YNorm,2);
%Pre 1-3
C_Pre_full(:,:,1) = c.Data.(SORT)(1).Full_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,1);%Full
C_Pre_full(:,:,2) = c.Data.(SORT)(2).Full_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,1);
C_Pre_full(:,:,3) = c.Data.(SORT)(3).Full_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,1);
C_Pre_full = nanmean(nanmean(C_Pre_full,3),2);

C_Pre_early (:,:,1) = c.Data.(SORT)(1).Early_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,2);%Full
C_Pre_early (:,:,2) = c.Data.(SORT)(2).Early_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,2);
C_Pre_early (:,:,3) = c.Data.(SORT)(3).Early_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,2);
C_Pre_early  = nanmean(C_Pre_early ,3);

C_Pre_late (:,:,1) = c.Data.(SORT)(1).Late_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,3);%Full
C_Pre_late (:,:,2) = c.Data.(SORT)(2).Late_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,3);
C_Pre_late (:,:,3) = c.Data.(SORT)(3).Late_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,3);
C_Pre_late  = nanmean(C_Pre_late ,3);

L_Pre_full(:,:,1) = l.Data.(SORT)(1).Full_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,1);%Full
L_Pre_full(:,:,2) = l.Data.(SORT)(2).Full_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,1);
L_Pre_full(:,:,3) = l.Data.(SORT)(3).Full_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,1);
L_Pre_full = nanmean(L_Pre_full,3);

L_Pre_early (:,:,1) = l.Data.(SORT)(1).Early_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,2);%Full
L_Pre_early (:,:,2) = l.Data.(SORT)(2).Early_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,2);
L_Pre_early (:,:,3) = l.Data.(SORT)(3).Early_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,2);
L_Pre_early  = nanmean(L_Pre_early ,3);

L_Pre_late (:,:,1) = l.Data.(SORT)(1).Late_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,3);%Full
L_Pre_late (:,:,2) = l.Data.(SORT)(2).Late_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,3);
L_Pre_late (:,:,3) = l.Data.(SORT)(3).Late_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,3);
L_Pre_late  = nanmean(L_Pre_late ,3);

Y_Pre_full(:,:,1) = y.Data.(SORT)(1).Full_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,1);%Full
Y_Pre_full(:,:,2) = y.Data.(SORT)(2).Full_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,1);
Y_Pre_full(:,:,3) = y.Data.(SORT)(3).Full_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,1);
Y_Pre_full = nanmean(Y_Pre_full,3);

Y_Pre_early (:,:,1) = y.Data.(SORT)(1).Early_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,2);%Full
Y_Pre_early (:,:,2) = y.Data.(SORT)(2).Early_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,2);
Y_Pre_early (:,:,3) = y.Data.(SORT)(3).Early_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,2);
Y_Pre_early  = nanmean(Y_Pre_early ,3);

Y_Pre_late (:,:,1) = y.Data.(SORT)(1).Late_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,3);%Full
Y_Pre_late (:,:,2) = y.Data.(SORT)(2).Late_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,3);
Y_Pre_late (:,:,3) = y.Data.(SORT)(3).Late_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,3);
Y_Pre_late  = nanmean(Y_Pre_late ,3);

%Laser
C_L_full = c.Data.(SORT)(4).Full_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,1);%Full
C_L_early = c.Data.(SORT)(4).Early_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,2);%early
C_L_late = c.Data.(SORT)(4).Late_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,3);%late

L_L_full = l.Data.(SORT)(4).Full_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,1);%Full
L_L_early = l.Data.(SORT)(4).Early_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,2);%early
L_L_late = l.Data.(SORT)(4).Late_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,3);%late

Y_L_full = y.Data.(SORT)(4).Full_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,1);%Full
Y_L_early = y.Data.(SORT)(4).Early_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,2);%early
Y_L_late = y.Data.(SORT)(4).Late_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,3);%late


%Post 1-3
C_Post1_full(:,:,1) = c.Data.(SORT)(5).Full_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,1);%Full
C_Post1_full(:,:,2) = c.Data.(SORT)(6).Full_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,1);
C_Post1_full(:,:,3) = c.Data.(SORT)(7).Full_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,1);
C_Post1_full = nanmean(C_Post1_full,3);

C_Post1_early (:,:,1) = c.Data.(SORT)(5).Early_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,2);%Full
C_Post1_early (:,:,2) = c.Data.(SORT)(6).Early_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,2);
C_Post1_early (:,:,3) = c.Data.(SORT)(7).Early_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,2);
C_Post1_early  = nanmean(C_Post1_early ,3);

C_Post1_late (:,:,1) = c.Data.(SORT)(5).Late_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,3);%Full
C_Post1_late (:,:,2) = c.Data.(SORT)(6).Late_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,3);
C_Post1_late (:,:,3) = c.Data.(SORT)(7).Late_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,3);
C_Post1_late  = nanmean(C_Post1_late ,3);

L_Post1_full(:,:,1) = l.Data.(SORT)(5).Full_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,1);%Full
L_Post1_full(:,:,2) = l.Data.(SORT)(6).Full_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,1);
L_Post1_full(:,:,3) = l.Data.(SORT)(7).Full_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,1);
L_Post1_full = nanmean(L_Post1_full,3);

L_Post1_early (:,:,1) = l.Data.(SORT)(5).Early_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,2);%Full
L_Post1_early (:,:,2) = l.Data.(SORT)(6).Early_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,2);
L_Post1_early (:,:,3) = l.Data.(SORT)(7).Early_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,2);
L_Post1_early  = nanmean(L_Post1_early ,3);

L_Post1_late (:,:,1) = l.Data.(SORT)(5).Late_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,3);%Full
L_Post1_late (:,:,2) = l.Data.(SORT)(6).Late_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,3);
L_Post1_late (:,:,3) = l.Data.(SORT)(7).Late_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,3);
L_Post1_late  = nanmean(L_Post1_late ,3);

Y_Post1_full(:,:,1) = y.Data.(SORT)(5).Full_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,1);%Full
Y_Post1_full(:,:,2) = y.Data.(SORT)(6).Full_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,1);
Y_Post1_full(:,:,3) = y.Data.(SORT)(7).Full_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,1);
Y_Post1_full = nanmean(Y_Post1_full,3);

Y_Post1_early (:,:,1) = y.Data.(SORT)(5).Early_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,2);%Full
Y_Post1_early (:,:,2) = y.Data.(SORT)(6).Early_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,2);
Y_Post1_early (:,:,3) = y.Data.(SORT)(7).Early_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,2);
Y_Post1_early  = nanmean(Y_Post1_early ,3);

Y_Post1_late (:,:,1) = y.Data.(SORT)(5).Late_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,3);%Full
Y_Post1_late (:,:,2) = y.Data.(SORT)(6).Late_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,3);
Y_Post1_late (:,:,3) = y.Data.(SORT)(7).Late_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,3);
Y_Post1_late  = nanmean(Y_Post1_late ,3);

%%% Post2
C_Post2_full(:,:,1) = c.Data.(SORT)(8).Full_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,1);%Full
C_Post2_full(:,:,2) = c.Data.(SORT)(9).Full_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,1);
C_Post2_full(:,:,3) = c.Data.(SORT)(10).Full_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,1);
C_Post2_full(:,:,4) = c.Data.(SORT)(11).Full_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,1);
C_Post2_full = nanmean(C_Post2_full,3);

C_Post2_early (:,:,1) = c.Data.(SORT)(8).Early_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,2);%Full
C_Post2_early (:,:,2) = c.Data.(SORT)(9).Early_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,2);
C_Post2_early (:,:,3) = c.Data.(SORT)(10).Early_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,2);
C_Post2_early (:,:,4) = c.Data.(SORT)(11).Early_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,2);
C_Post2_early  = nanmean(C_Post2_early ,3);

C_Post2_late (:,:,1) = c.Data.(SORT)(8).Late_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,3);%Full
C_Post2_late (:,:,2) = c.Data.(SORT)(9).Late_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,3);
C_Post2_late (:,:,3) = c.Data.(SORT)(10).Late_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,3);
C_Post2_late (:,:,4) = c.Data.(SORT)(11).Late_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,3);
C_Post2_late  = nanmean(C_Post2_late ,3);

L_Post2_full(:,:,1) = l.Data.(SORT)(8).Full_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,1);%Full
L_Post2_full(:,:,2) = l.Data.(SORT)(9).Full_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,1);
L_Post2_full(:,:,3) = l.Data.(SORT)(10).Full_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,1);
L_Post2_full(:,:,4) = l.Data.(SORT)(11).Full_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,1);
L_Post2_full = nanmean(L_Post2_full,3);

L_Post2_early (:,:,1) = l.Data.(SORT)(8).Early_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,2);%Full
L_Post2_early (:,:,2) = l.Data.(SORT)(9).Early_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,2);
L_Post2_early (:,:,3) = l.Data.(SORT)(10).Early_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,2);
L_Post2_early (:,:,4) = l.Data.(SORT)(11).Early_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,2);
L_Post2_early  = nanmean(L_Post2_early ,3);

L_Post2_late (:,:,1) = l.Data.(SORT)(8).Late_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,3);%Full
L_Post2_late (:,:,2) = l.Data.(SORT)(9).Late_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,3);
L_Post2_late (:,:,3) = l.Data.(SORT)(10).Late_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,3);
L_Post2_late (:,:,4) = l.Data.(SORT)(11).Late_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,3);
L_Post2_late  = nanmean(L_Post2_late ,3);

Y_Post2_full(:,:,1) = y.Data.(SORT)(8).Full_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,1);%Full
Y_Post2_full(:,:,2) = y.Data.(SORT)(9).Full_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,1);
Y_Post2_full(:,:,3) = y.Data.(SORT)(10).Full_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,1);
Y_Post2_full(:,:,4) = y.Data.(SORT)(11).Full_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,1);
Y_Post2_full = nanmean(Y_Post2_full,3);

Y_Post2_early (:,:,1) = y.Data.(SORT)(8).Early_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,2);%Full
Y_Post2_early (:,:,2) = y.Data.(SORT)(9).Early_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,2);
Y_Post2_early (:,:,3) = y.Data.(SORT)(10).Early_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,2);
Y_Post2_early (:,:,4) = y.Data.(SORT)(11).Early_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,2);
Y_Post2_early  = nanmean(Y_Post2_early ,3);

Y_Post2_late (:,:,1) = y.Data.(SORT)(8).Late_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,3);%Full
Y_Post2_late (:,:,2) = y.Data.(SORT)(9).Late_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,3);
Y_Post2_late (:,:,3) = y.Data.(SORT)(10).Late_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,3);
Y_Post2_late (:,:,4) = y.Data.(SORT)(11).Late_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,3);
Y_Post2_late  = nanmean(Y_Post2_late ,3);


% %%%Post 9-11
% C_Post3_full(:,:,1) = c.Data.(SORT)(9).Full_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,1);%Full
% C_Post3_full(:,:,2) = c.Data.(SORT)(10).Full_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,1);
% C_Post3_full(:,:,3) = c.Data.(SORT)(11).Full_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,1);
% C_Post3_full = nanmean(C_Post3_full,3);
% 
% C_Post3_early (:,:,1) = c.Data.(SORT)(9).Early_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,2);%Full
% C_Post3_early (:,:,2) = c.Data.(SORT)(10).Early_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,2);
% C_Post3_early (:,:,3) = c.Data.(SORT)(11).Early_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,2);
% C_Post3_early  = nanmean(C_Post3_early ,3);
% 
% C_Post3_late (:,:,1) = c.Data.(SORT)(9).Late_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,3);%Full
% C_Post3_late (:,:,2) = c.Data.(SORT)(10).Late_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,3);
% C_Post3_late (:,:,3) = c.Data.(SORT)(11).Late_RMS_RELRES(:,BF_C+BFlim)./CNorm(:,:,3);
% C_Post3_late  = nanmean(C_Post3_late ,3);
% 
% L_Post3_full(:,:,1) = l.Data.(SORT)(9).Full_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,1);%Full
% L_Post3_full(:,:,2) = l.Data.(SORT)(10).Full_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,1);
% L_Post3_full(:,:,3) = l.Data.(SORT)(11).Full_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,1);
% L_Post3_full = nanmean(L_Post3_full,3);
% 
% L_Post3_early (:,:,1) = l.Data.(SORT)(9).Early_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,2);%Full
% L_Post3_early (:,:,2) = l.Data.(SORT)(10).Early_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,2);
% L_Post3_early (:,:,3) = l.Data.(SORT)(11).Early_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,2);
% L_Post3_early  = nanmean(L_Post3_early ,3);
% 
% L_Post3_late (:,:,1) = l.Data.(SORT)(9).Late_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,3);%Full
% L_Post3_late (:,:,2) = l.Data.(SORT)(10).Late_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,3);
% L_Post3_late (:,:,3) = l.Data.(SORT)(11).Late_RMS_RELRES(:,BF_L+BFlim)./LNorm(:,:,3);
% L_Post3_late  = nanmean(L_Post3_late ,3);
% 
% Y_Post3_full(:,:,1) = y.Data.(SORT)(9).Full_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,1);%Full
% Y_Post3_full(:,:,2) = y.Data.(SORT)(10).Full_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,1);
% Y_Post3_full(:,:,3) = y.Data.(SORT)(11).Full_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,1);
% Y_Post3_full = nanmean(Y_Post3_full,3);
% 
% Y_Post3_early (:,:,1) = y.Data.(SORT)(9).Early_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,2);%Full
% Y_Post3_early (:,:,2) = y.Data.(SORT)(10).Early_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,2);
% Y_Post3_early (:,:,3) = y.Data.(SORT)(11).Early_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,2);
% Y_Post3_early  = nanmean(Y_Post3_early ,3);
% 
% Y_Post3_late (:,:,1) = y.Data.(SORT)(9).Late_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,3);%Full
% Y_Post3_late (:,:,2) = y.Data.(SORT)(10).Late_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,3);
% Y_Post3_late (:,:,3) = y.Data.(SORT)(11).Late_RMS_RELRES(:,BF_Y+BFlim)./YNorm(:,:,3);
% Y_Post3_late  = nanmean(Y_Post3_late ,3);
clear Pre_O L_O Post1_O Post2_O C_O Opto_O
% L = [3 1];
% Pre_O = teg_repeated_measures_ANOVA([vertcat(C_Pre_full, nan(5,1)) L_Pre_full...
%     vertcat(Y_Pre_full, nan(8,1))],L,{'Group','Measurement'});
% L = [3 1];
% L_O = teg_repeated_measures_ANOVA([vertcat(C_L_full, nan(5,1)) L_L_full...
%     vertcat(Y_L_full, nan(8,1))],L,{'Group','Measurement'});
% L = [3 1];
% Post1_O = teg_repeated_measures_ANOVA([vertcat(C_Post1_full, nan(5,1)) L_Post1_full...
%     vertcat(Y_Post1_full, nan(8,1))],L,{'Group','Measurement'});
% L = [3 1];
% Post2_O = teg_repeated_measures_ANOVA([vertcat(C_Post2_full, nan(5,1)) L_Post2_full...
%     vertcat(Y_Post2_full, nan(8,1))],L,{'Group','Measurement'});
% 
% L = [1 4];
% C_O = teg_repeated_measures_ANOVA([C_Pre_full C_L_full C_Post1_full C_Post2_full],L,{'Group','Measurement'});
% L = [1 4];
% Opto_O = teg_repeated_measures_ANOVA([L_Pre_full L_L_full L_Post1_full L_Post2_full],L,{'Group','Measurement'});
% L = [1 4];
% Y_O = teg_repeated_measures_ANOVA([Y_Pre_full Y_L_full Y_Post1_full Y_Post2_full],L,{'Group','Measurement'});
% keyboard
figure
subplot(1,3,1)
bar([1 4 7 10 ],[nanmean(C_Pre_full) nanmean(C_L_full) nanmean(C_Post1_full) nanmean(C_Post2_full) ],0.3)
hold on
bar([2 5 8 11 ],[nanmean(L_Pre_full) nanmean(L_L_full) nanmean(L_Post1_full) nanmean(L_Post2_full) ],0.3)
bar([3 6 9 12 ],[nanmean(Y_Pre_full) nanmean(Y_L_full) nanmean(Y_Post1_full) nanmean(Y_Post2_full) ],0.3)

errorbar([1 4 7 10 ],[nanmean(C_Pre_full) nanmean(C_L_full) nanmean(C_Post1_full) nanmean(C_Post2_full) ],...
    [nanstd(C_Pre_full)/sqrt(sum(~isnan(C_Pre_full)))...
    nanstd(C_L_full)/sqrt(sum(~isnan(C_L_full)))...
    nanstd(C_Post1_full)/sqrt(sum(~isnan(C_Post1_full)))...
    nanstd(C_Post2_full)/sqrt(sum(~isnan(C_Post2_full)))...
%     nanstd(C_Post3_full)/sqrt(sum(~isnan(C_Post3_full)))...
    ],'.', 'LineWidth', 2)

errorbar([2 5 8 11 ],[nanmean(L_Pre_full) nanmean(L_L_full) nanmean(L_Post1_full) nanmean(L_Post2_full) ],...
    [nanstd(L_Pre_full)/sqrt(sum(~isnan(L_Pre_full)))...
    nanstd(L_L_full)/sqrt(sum(~isnan(L_L_full)))...
    nanstd(L_Post1_full)/sqrt(sum(~isnan(L_Post1_full)))...
    nanstd(L_Post2_full)/sqrt(sum(~isnan(L_Post2_full)))...
%     nanstd(L_Post3_full)/sqrt(sum(~isnan(L_Post3_full)))...
    ],'.', 'LineWidth', 2)

errorbar([3 6 9 12 ],[nanmean(Y_Pre_full) nanmean(Y_L_full) nanmean(Y_Post1_full) nanmean(Y_Post2_full) ],...
    [nanstd(Y_Pre_full)/sqrt(sum(~isnan(Y_Pre_full)))...
    nanstd(Y_L_full)/sqrt(sum(~isnan(Y_L_full)))...
    nanstd(Y_Post1_full)/sqrt(sum(~isnan(Y_Post1_full)))...
    nanstd(Y_Post2_full)/sqrt(sum(~isnan(Y_Post2_full)))...
%     nanstd(Y_Post3_full)/sqrt(sum(~isnan(Y_Post3_full)))...
    ],'.', 'LineWidth', 2)

xticks([2 5 8 11 ])
xticklabels({'Pre 1-3','Opto','Post 1-3','Post 4-7'})
ylim([0.75 2])
set(gca,'fontsize',12,'FontWeight','bold')
title ('Full RMS RELRES')
legend('Control','Opto','YFP','Location','best')



subplot(1,3,2)
bar([1 4 7 10 ],[nanmean(C_Pre_early) nanmean(C_L_early) nanmean(C_Post1_early) nanmean(C_Post2_early) ],0.3)
hold on
bar([2 5 8 11 ],[nanmean(L_Pre_early) nanmean(L_L_early) nanmean(L_Post1_early) nanmean(L_Post2_early) ],0.3)
bar([3 6 9 12 ],[nanmean(Y_Pre_early) nanmean(Y_L_early) nanmean(Y_Post1_early) nanmean(Y_Post2_early) ],0.3)

errorbar([1 4 7 10 ],[nanmean(C_Pre_early) nanmean(C_L_early) nanmean(C_Post1_early) nanmean(C_Post2_early) ],...
    [nanstd(C_Pre_early)/sqrt(sum(~isnan(C_Pre_early)))...
    nanstd(C_L_early)/sqrt(sum(~isnan(C_L_early)))...
    nanstd(C_Post1_early)/sqrt(sum(~isnan(C_Post1_early)))...
    nanstd(C_Post2_early)/sqrt(sum(~isnan(C_Post2_early)))...
%     nanstd(C_Post3_early)/sqrt(sum(~isnan(C_Post3_early)))...
    ],'.', 'LineWidth', 2)

errorbar([2 5 8 11 ],[nanmean(L_Pre_early) nanmean(L_L_early) nanmean(L_Post1_early) nanmean(L_Post2_early) ],...
    [nanstd(L_Pre_early)/sqrt(sum(~isnan(L_Pre_early)))...
    nanstd(L_L_early)/sqrt(sum(~isnan(L_L_early)))...
    nanstd(L_Post1_early)/sqrt(sum(~isnan(L_Post1_early)))...
    nanstd(L_Post2_early)/sqrt(sum(~isnan(L_Post2_early)))...
%     nanstd(L_Post3_early)/sqrt(sum(~isnan(L_Post3_early)))...
    ],'.', 'LineWidth', 2)

errorbar([3 6 9 12 ],[nanmean(Y_Pre_early) nanmean(Y_L_early) nanmean(Y_Post1_early) nanmean(Y_Post2_early) ],...
    [nanstd(Y_Pre_early)/sqrt(sum(~isnan(Y_Pre_early)))...
    nanstd(Y_L_early)/sqrt(sum(~isnan(Y_L_early)))...
    nanstd(Y_Post1_early)/sqrt(sum(~isnan(Y_Post1_early)))...
    nanstd(Y_Post2_early)/sqrt(sum(~isnan(Y_Post2_early)))...
%     nanstd(Y_Post3_early)/sqrt(sum(~isnan(Y_Post3_early)))...
    ],'.', 'LineWidth', 2)

xticks([2 5 8 11 ])
xticklabels({'Pre 1-3','Opto','Post 1-3','Post 4-7'})
ylim([0.75 2])
set(gca,'fontsize',12,'FontWeight','bold')
title ('early RMS RELRES')
legend('Control','Opto','YFP','Location','best')


subplot(1,3,3)
bar([1 4 7 10 ],[nanmean(C_Pre_late) nanmean(C_L_late) nanmean(C_Post1_late) nanmean(C_Post2_late) ],0.3)
hold on
bar([2 5 8 11 ],[nanmean(L_Pre_late) nanmean(L_L_late) nanmean(L_Post1_late) nanmean(L_Post2_late) ],0.3)
bar([3 6 9 12 ],[nanmean(Y_Pre_late) nanmean(Y_L_late) nanmean(Y_Post1_late) nanmean(Y_Post2_late) ],0.3)

errorbar([1 4 7 10 ],[nanmean(C_Pre_late) nanmean(C_L_late) nanmean(C_Post1_late) nanmean(C_Post2_late) ],...
    [nanstd(C_Pre_late)/sqrt(sum(~isnan(C_Pre_late)))...
    nanstd(C_L_late)/sqrt(sum(~isnan(C_L_late)))...
    nanstd(C_Post1_late)/sqrt(sum(~isnan(C_Post1_late)))...
    nanstd(C_Post2_late)/sqrt(sum(~isnan(C_Post2_late)))...
%     nanstd(C_Post3_late)/sqrt(sum(~isnan(C_Post3_late)))...
    ],'.', 'LineWidth', 2)

errorbar([2 5 8 11 ],[nanmean(L_Pre_late) nanmean(L_L_late) nanmean(L_Post1_late) nanmean(L_Post2_late) ],...
    [nanstd(L_Pre_late)/sqrt(sum(~isnan(L_Pre_late)))...
    nanstd(L_L_late)/sqrt(sum(~isnan(L_L_late)))...
    nanstd(L_Post1_late)/sqrt(sum(~isnan(L_Post1_late)))...
    nanstd(L_Post2_late)/sqrt(sum(~isnan(L_Post2_late)))...
%     nanstd(L_Post3_late)/sqrt(sum(~isnan(L_Post3_late)))...
    ],'.', 'LineWidth', 2)

errorbar([3 6 9 12 ],[nanmean(Y_Pre_late) nanmean(Y_L_late) nanmean(Y_Post1_late) nanmean(Y_Post2_late) ],...
    [nanstd(Y_Pre_late)/sqrt(sum(~isnan(Y_Pre_late)))...
    nanstd(Y_L_late)/sqrt(sum(~isnan(Y_L_late)))...
    nanstd(Y_Post1_late)/sqrt(sum(~isnan(Y_Post1_late)))...
    nanstd(Y_Post2_late)/sqrt(sum(~isnan(Y_Post2_late)))...
%     nanstd(Y_Post3_late)/sqrt(sum(~isnan(Y_Post3_late)))...
    ],'.', 'LineWidth', 2)

xticks([2 5 8 11 14])
xticklabels({'Pre 1-3','Opto','Post 1-3','Post 4-7'})
ylim([0.75 2])
set(gca,'fontsize',12,'FontWeight','bold')
title ('late RMS RELRES')
legend('Control','Opto','YFP','Location','best')

% keyboard
clear 
M = [1, 4, 9]
load('K:\CSD_dynamic_analysis\DATA\Output\Output_Input_HighP_7post_Data_Threshold_25_Zscore_0_binned_1_mirror_0.mat');
BF = Data.BF_Pos;
BF_lim = 2;
SORT = 'ST_based';
% Sinks = {'I_IIL'  'InfE'  'VIL'};
Sinks = {'I_IIL' 'IVE'  'VaE' 'VbE' 'VIE' 'InfE' 'VIL'};%    
c = {'b','g','r','y','c'};

LNorm = nan(length(Data.names),3, length(Sinks));
for i2 = 1:length(Sinks)
    
    for i1 = 1:3
        LNorm(:,i1,i2) =  Data.(SORT)(i1).SinkPeakAmp.(Sinks{i2})(:,BF);
    end
end
LNorm = nanmean(LNorm,2);

 figure
for i1 = 1:length(Sinks)   
    
    subplot(length(Sinks),1,i1)
    
    hold on
    for i2 = 1:length(M)
    Dummy = Data.(SORT)(M(i2)).SinkPeakAmp.(Sinks{i1})(:,BF-BF_lim:BF+BF_lim)./LNorm(:,:,i1); 
    SEM = nanstd(Dummy)./sqrt(sum(~isnan(Dummy)));
    shadedErrorBar([],nanmean(Dummy),SEM,'lineprops',{ c{i2},'markerfacecolor',c{i2}});
    end
    title(Sinks{i1})
    if i1 == 1
    legend('Pre','','','','Opto','','','','Post','','','','Location','best')
    end    
    xticks([1 2 3 4 5])
    xticklabels({'- 3/4','- 1/2','BF','+ 1/2','+3/4'})
    
end

% %%%keyboard
figure

for i1 = 1:length(Sinks)   
    
    subplot(length(Sinks),1,i1)
    
    hold on
    
   %%%Pre 1-3
    Dummy(:,:,1) = Data.(SORT)(1).SinkPeakAmp.(Sinks{i1})(:,BF-BF_lim:BF+BF_lim)./LNorm(:,:,i1); 
    Dummy(:,:,2) = Data.(SORT)(2).SinkPeakAmp.(Sinks{i1})(:,BF-BF_lim:BF+BF_lim)./LNorm(:,:,i1);
    Dummy(:,:,3) = Data.(SORT)(3).SinkPeakAmp.(Sinks{i1})(:,BF-BF_lim:BF+BF_lim)./LNorm(:,:,i1);
    Dummymean = nanmean(Dummy,3);
    STD = vertcat(Dummy(:,:,1), Dummy(:,:,2), Dummy(:,:,3));    
    SEM = nanstd(STD)./sqrt(sum(~isnan(STD)));  
%     shadedErrorBar([],nanmean(Dummymean),SEM,'lineprops',{ c{1},'markerfacecolor',c{1}});
    plot(nanmean(Dummymean),'LineWidth',2)
    clear Dummy
    
    Dummy(:,:,1) = Data.(SORT)(4).SinkPeakAmp.(Sinks{i1})(:,BF-BF_lim:BF+BF_lim)./LNorm(:,:,i1); 
    Dummymean = nanmean(Dummy,3);
    STD = Dummy(:,:,1);    
    SEM = nanstd(STD)./sqrt(sum(~isnan(STD)));  
%     shadedErrorBar([],nanmean(Dummymean),SEM,'lineprops',{ c{2},'markerfacecolor',c{2}});
    plot(nanmean(Dummymean),'LineWidth',2)
    clear Dummy
    
    %%%Post 1-3
    Dummy(:,:,1) = Data.(SORT)(5).SinkPeakAmp.(Sinks{i1})(:,BF-BF_lim:BF+BF_lim)./LNorm(:,:,i1); 
    Dummy(:,:,2) = Data.(SORT)(6).SinkPeakAmp.(Sinks{i1})(:,BF-BF_lim:BF+BF_lim)./LNorm(:,:,i1);
    Dummy(:,:,3) = Data.(SORT)(7).SinkPeakAmp.(Sinks{i1})(:,BF-BF_lim:BF+BF_lim)./LNorm(:,:,i1);
    Dummymean = nanmean(Dummy,3);
    STD = vertcat(Dummy(:,:,1), Dummy(:,:,2), Dummy(:,:,3));    
    SEM = nanstd(STD)./sqrt(sum(~isnan(STD))); 
     plot(nanmean(Dummymean),'LineWidth',2)
%     shadedErrorBar([],nanmean(Dummymean),SEM,'lineprops',{ c{3},'markerfacecolor',c{3}});
    clear Dummy
    
    %%%Post 3-5
    Dummy(:,:,1) = Data.(SORT)(8).SinkPeakAmp.(Sinks{i1})(:,BF-BF_lim:BF+BF_lim)./LNorm(:,:,i1); 
    Dummy(:,:,2) = Data.(SORT)(9).SinkPeakAmp.(Sinks{i1})(:,BF-BF_lim:BF+BF_lim)./LNorm(:,:,i1);
%     Dummy(:,:,3) = Data.(SORT)(9).SinkPeakAmp.(Sinks{i1})(:,BF-BF_lim:BF+BF_lim)./LNorm(:,:,i1);
    Dummymean = nanmean(Dummy,3);
    STD = vertcat(Dummy(:,:,1), Dummy(:,:,2));    
    SEM = nanstd(STD)./sqrt(sum(~isnan(STD)));  
%     shadedErrorBar([],nanmean(Dummymean),SEM,'lineprops',{ c{4},'markerfacecolor',c{4}});
    plot(nanmean(Dummymean),'LineWidth',2)
    clear Dummy
    
    %%%Post 3-5
    Dummy(:,:,1) = Data.(SORT)(10).SinkPeakAmp.(Sinks{i1})(:,BF-BF_lim:BF+BF_lim)./LNorm(:,:,i1); 
    Dummy(:,:,2) = Data.(SORT)(11).SinkPeakAmp.(Sinks{i1})(:,BF-BF_lim:BF+BF_lim)./LNorm(:,:,i1);
%     Dummy(:,:,3) = Data.(SORT)(11).SinkPeakAmp.(Sinks{i1})(:,BF-BF_lim:BF+BF_lim)./LNorm(:,:,i1);
    Dummymean = nanmean(Dummy,3);
    STD = vertcat(Dummy(:,:,1), Dummy(:,:,2));    
    SEM = nanstd(STD)./sqrt(sum(~isnan(STD)));  
%     shadedErrorBar([],nanmean(Dummymean),SEM,'lineprops',{ c{5},'markerfacecolor',c{5}});
    plot(nanmean(Dummymean),'LineWidth',2)
    clear Dummy
   
    title(Sinks{i1})
    if i1 == 1
    legend('Pre 1-3','Opto','Post 1-3','Post 4-5','Post 6-7','Location','best')
    end    
    xticks([1 2 3 4 5])
    xticklabels({'- 3/4','- 1/2','BF','+ 1/2','+3/4'})
    
end


clear


% %%%keyboard
M = [1, 4, 9]
SORT = 'ST_based';
Sinks = {'I_IIL' 'IVE'  'VaE' 'VbE' 'VIE' 'InfE' 'VIL'};
% Sinks = {'IVE' 'VaE'  'InfE'};%'I_IIL'  'VIL'    
BFlim= 1;
para = 'SinkPeakAmp'; 

c = load('K:\CSD_dynamic_analysis\DATA\Output\Output_Input_Control_7post_n7_Data_Threshold_25_Zscore_0_binned_1_mirror_1.mat');
BF_C = c.Data.BF_Pos;
l = load('K:\CSD_dynamic_analysis\DATA\Output\Output_Input_HighP_7post_Data_Threshold_25_Zscore_0_binned_1_mirror_1.mat');
BF_L = l.Data.BF_Pos;
y = load('K:\CSD_dynamic_analysis\DATA\Output\Output_Input_YFP_7post_Data_Threshold_25_Zscore_0_binned_1_mirror_1.mat');
BF_Y = y.Data.BF_Pos;

CNorm = nan(length(c.Data.names),3,length(Sinks));
LNorm = nan(length(l.Data.names),3,length(Sinks));
YNorm = nan(length(y.Data.names),3,length(Sinks));

for i2 = 1:length(Sinks)   

    for i1 = 1:3
        CNorm(:,i1,i2) =  c.Data.(SORT)(i1).(para).(Sinks{i2})(:,BF_C+BFlim);%full 
        LNorm(:,i1,i2) =  l.Data.(SORT)(i1).(para).(Sinks{i2})(:,BF_L+BFlim);%full
        YNorm(:,i1,i2) =  y.Data.(SORT)(i1).(para).(Sinks{i2})(:,BF_Y+BFlim);%full 
    end

end
CNorm = nanmean(CNorm,2);
LNorm = nanmean(LNorm,2);
YNorm = nanmean(YNorm,2);

figure
for i1 = 1:length(Sinks)
%Pre
C_Pre_full = c.Data.(SORT)(M(1)).(para).(Sinks{i1})(:,BF_C+BFlim)./CNorm(:,:,i1);%Full

L_Pre_full = l.Data.(SORT)(M(1)).(para).(Sinks{i1})(:,BF_L+BFlim)./LNorm(:,:,i1);%Full

Y_Pre_full = y.Data.(SORT)(M(1)).(para).(Sinks{i1})(:,BF_Y+BFlim)./YNorm(:,:,i1);%Full

%Laser
C_L_full = c.Data.(SORT)(M(2)).(para).(Sinks{i1})(:,BF_C+BFlim)./CNorm(:,:,i1);%Full

L_L_full = l.Data.(SORT)(M(2)).(para).(Sinks{i1})(:,BF_L+BFlim)./LNorm(:,:,i1);%Full

Y_L_full = y.Data.(SORT)(M(2)).(para).(Sinks{i1})(:,BF_Y+BFlim)./YNorm(:,:,i1);%Full

%Post
C_Post_full = c.Data.(SORT)(M(3)).(para).(Sinks{i1})(:,BF_C+BFlim)./CNorm(:,:,i1);%Full

L_Post_full = l.Data.(SORT)(M(3)).(para).(Sinks{i1})(:,BF_L+BFlim)./LNorm(:,:,i1);%Full

Y_Post_full = y.Data.(SORT)(M(3)).(para).(Sinks{i1})(:,BF_Y+BFlim)./YNorm(:,:,i1);%Full


subplot(length(Sinks),1,i1)
bar([1 4 7],[nanmean(C_Pre_full) nanmean(C_L_full) nanmean(C_Post_full)],0.3)
hold on
bar([2 5 8],[nanmean(L_Pre_full) nanmean(L_L_full) nanmean(L_Post_full)],0.3)
bar([3 6 9],[nanmean(Y_Pre_full) nanmean(Y_L_full) nanmean(Y_Post_full)],0.3)

errorbar([1 4 7],[nanmean(C_Pre_full) nanmean(C_L_full) nanmean(C_Post_full)],...
    [nanstd(C_Pre_full)/sqrt(sum(~isnan(C_Pre_full)))...
    nanstd(C_L_full)/sqrt(sum(~isnan(C_L_full)))...
    nanstd(C_Post_full)/sqrt(sum(~isnan(C_Post_full)))],'.','LineWidth',2)

errorbar([2 5 8],[nanmean(L_Pre_full) nanmean(L_L_full) nanmean(L_Post_full)],...
    [nanstd(L_Pre_full)/sqrt(sum(~isnan(L_Pre_full)))...
    nanstd(L_L_full)/sqrt(sum(~isnan(L_L_full)))...
    nanstd(L_Post_full)/sqrt(sum(~isnan(L_Post_full)))],'.','LineWidth',2)

errorbar([3 6 9],[nanmean(Y_Pre_full) nanmean(Y_L_full) nanmean(Y_Post_full)],...
    [nanstd(Y_Pre_full)/sqrt(sum(~isnan(Y_Pre_full)))...
    nanstd(Y_L_full)/sqrt(sum(~isnan(Y_L_full)))...
    nanstd(Y_Post_full)/sqrt(sum(~isnan(Y_Post_full)))],'.','LineWidth',2)

xticks([2 5 8])
xticklabels({'Pre','Opto','Post'})
set(gca,'fontsize',12,'FontWeight','bold')
ylim([0.5 1.5])

title (['normalized Peak Amps Layer' Sinks{i1}])
legend('Control','Opto','YFP','Location','best')
end

% keyboard

figure
for i1 = 1:length(Sinks)
    
%Pre 1-3

    Dummy(:,:,1) = c.Data.(SORT)(1).(para).(Sinks{i1})(:,BF_C+BFlim)./CNorm(:,:,i1); 
    Dummy(:,:,2) = c.Data.(SORT)(2).(para).(Sinks{i1})(:,BF_C+BFlim)./CNorm(:,:,i1);
    Dummy(:,:,3) = c.Data.(SORT)(3).(para).(Sinks{i1})(:,BF_C+BFlim)./CNorm(:,:,i1);
    C_Pre_full = nanmean(Dummy,3);    
    clear Dummy
    
    Dummy(:,:,1) = l.Data.(SORT)(1).(para).(Sinks{i1})(:,BF_L+BFlim)./LNorm(:,:,i1); 
    Dummy(:,:,2) = l.Data.(SORT)(2).(para).(Sinks{i1})(:,BF_L+BFlim)./LNorm(:,:,i1);
    Dummy(:,:,3) = l.Data.(SORT)(3).(para).(Sinks{i1})(:,BF_L+BFlim)./LNorm(:,:,i1);
    L_Pre_full = nanmean(Dummy,3);    
    clear Dummy

    Dummy(:,:,1) = y.Data.(SORT)(1).(para).(Sinks{i1})(:,BF_Y+BFlim)./YNorm(:,:,i1); 
    Dummy(:,:,2) = y.Data.(SORT)(2).(para).(Sinks{i1})(:,BF_Y+BFlim)./YNorm(:,:,i1);
    Dummy(:,:,3) = y.Data.(SORT)(3).(para).(Sinks{i1})(:,BF_Y+BFlim)./YNorm(:,:,i1);
    Y_Pre_full = nanmean(Dummy,3);    
    clear Dummy


%Laser
C_L_full = c.Data.(SORT)(4).(para).(Sinks{i1})(:,BF_C+BFlim)./CNorm(:,:,i1);%Full

L_L_full = l.Data.(SORT)(4).(para).(Sinks{i1})(:,BF_L+BFlim)./LNorm(:,:,i1);%Full

Y_L_full = y.Data.(SORT)(4).(para).(Sinks{i1})(:,BF_Y+BFlim)./YNorm(:,:,i1);%Full

%Post 1-3

    Dummy(:,:,1) = c.Data.(SORT)(5).(para).(Sinks{i1})(:,BF_C+BFlim)./CNorm(:,:,i1); 
    Dummy(:,:,2) = c.Data.(SORT)(6).(para).(Sinks{i1})(:,BF_C+BFlim)./CNorm(:,:,i1);
    Dummy(:,:,3) = c.Data.(SORT)(7).(para).(Sinks{i1})(:,BF_C+BFlim)./CNorm(:,:,i1);
    C_Post1_full = nanmean(Dummy,3);    
    clear Dummy
    
    Dummy(:,:,1) = l.Data.(SORT)(5).(para).(Sinks{i1})(:,BF_L+BFlim)./LNorm(:,:,i1); 
    Dummy(:,:,2) = l.Data.(SORT)(6).(para).(Sinks{i1})(:,BF_L+BFlim)./LNorm(:,:,i1);
    Dummy(:,:,3) = l.Data.(SORT)(7).(para).(Sinks{i1})(:,BF_L+BFlim)./LNorm(:,:,i1);
    L_Post1_full = nanmean(Dummy,3);    
    clear Dummy

    Dummy(:,:,1) = y.Data.(SORT)(5).(para).(Sinks{i1})(:,BF_Y+BFlim)./YNorm(:,:,i1); 
    Dummy(:,:,2) = y.Data.(SORT)(6).(para).(Sinks{i1})(:,BF_Y+BFlim)./YNorm(:,:,i1);
    Dummy(:,:,3) = y.Data.(SORT)(7).(para).(Sinks{i1})(:,BF_Y+BFlim)./YNorm(:,:,i1);
    Y_Post1_full = nanmean(Dummy,3);    
    clear Dummy

%Post 3-5

    Dummy(:,:,1) = c.Data.(SORT)(8).(para).(Sinks{i1})(:,BF_C+BFlim)./CNorm(:,:,i1); 
    Dummy(:,:,2) = c.Data.(SORT)(9).(para).(Sinks{i1})(:,BF_C+BFlim)./CNorm(:,:,i1);
    Dummy(:,:,3) = c.Data.(SORT)(10).(para).(Sinks{i1})(:,BF_C+BFlim)./CNorm(:,:,i1);
    Dummy(:,:,4) = c.Data.(SORT)(11).(para).(Sinks{i1})(:,BF_C+BFlim)./CNorm(:,:,i1);
    C_Post2_full = nanmean(Dummy,3);    
    clear Dummy
    
    Dummy(:,:,1) = l.Data.(SORT)(8).(para).(Sinks{i1})(:,BF_L+BFlim)./LNorm(:,:,i1); 
    Dummy(:,:,2) = l.Data.(SORT)(9).(para).(Sinks{i1})(:,BF_L+BFlim)./LNorm(:,:,i1);
    Dummy(:,:,3) = l.Data.(SORT)(10).(para).(Sinks{i1})(:,BF_L+BFlim)./LNorm(:,:,i1);
    Dummy(:,:,4) = l.Data.(SORT)(11).(para).(Sinks{i1})(:,BF_L+BFlim)./LNorm(:,:,i1);
    L_Post2_full = nanmean(Dummy,3);    
    clear Dummy

    Dummy(:,:,1) = y.Data.(SORT)(8).(para).(Sinks{i1})(:,BF_Y+BFlim)./YNorm(:,:,i1); 
    Dummy(:,:,2) = y.Data.(SORT)(9).(para).(Sinks{i1})(:,BF_Y+BFlim)./YNorm(:,:,i1);
    Dummy(:,:,3) = y.Data.(SORT)(10).(para).(Sinks{i1})(:,BF_Y+BFlim)./YNorm(:,:,i1);
    Dummy(:,:,4) = y.Data.(SORT)(11).(para).(Sinks{i1})(:,BF_Y+BFlim)./YNorm(:,:,i1);
    Y_Post2_full = nanmean(Dummy,3);    
    clear Dummy

% %Post 5-7
% 
%     Dummy(:,:,1) = c.Data.(SORT)(9).SinkPeakAmp.(Sinks{i1})(:,BF_C+BFlim)./CNorm(:,:,i1); 
%     Dummy(:,:,2) = c.Data.(SORT)(10).SinkPeakAmp.(Sinks{i1})(:,BF_C+BFlim)./CNorm(:,:,i1);
%     Dummy(:,:,3) = c.Data.(SORT)(11).SinkPeakAmp.(Sinks{i1})(:,BF_C+BFlim)./CNorm(:,:,i1);
%     C_Post3_full = nanmean(Dummy,3);    
%     clear Dummy
%     
%     Dummy(:,:,1) = l.Data.(SORT)(9).SinkPeakAmp.(Sinks{i1})(:,BF_L+BFlim)./LNorm(:,:,i1); 
%     Dummy(:,:,2) = l.Data.(SORT)(10).SinkPeakAmp.(Sinks{i1})(:,BF_L+BFlim)./LNorm(:,:,i1);
%     Dummy(:,:,3) = l.Data.(SORT)(11).SinkPeakAmp.(Sinks{i1})(:,BF_L+BFlim)./LNorm(:,:,i1);
%     L_Post3_full = nanmean(Dummy,3);    
%     clear Dummy
% 
%     Dummy(:,:,1) = y.Data.(SORT)(9).SinkPeakAmp.(Sinks{i1})(:,BF_Y+BFlim)./YNorm(:,:,i1); 
%     Dummy(:,:,2) = y.Data.(SORT)(10).SinkPeakAmp.(Sinks{i1})(:,BF_Y+BFlim)./YNorm(:,:,i1);
%     Dummy(:,:,3) = y.Data.(SORT)(11).SinkPeakAmp.(Sinks{i1})(:,BF_Y+BFlim)./YNorm(:,:,i1);
%     Y_Post3_full = nanmean(Dummy,3);    
%     clear Dummy




subplot(length(Sinks),1,i1)
bar([1 4 7 10 ],[nanmean(C_Pre_full) nanmean(C_L_full) nanmean(C_Post1_full) nanmean(C_Post2_full) ],0.3)
hold on
bar([2 5 8 11 ],[nanmean(L_Pre_full) nanmean(L_L_full) nanmean(L_Post1_full) nanmean(L_Post2_full) ],0.3)
bar([3 6 9 12 ],[nanmean(Y_Pre_full) nanmean(Y_L_full) nanmean(Y_Post1_full) nanmean(Y_Post2_full) ],0.3)

errorbar([1 4 7 10 ],[nanmean(C_Pre_full) nanmean(C_L_full) nanmean(C_Post1_full) nanmean(C_Post2_full)],...
    [nanstd(C_Pre_full)/sqrt(sum(~isnan(C_Pre_full)))...
    nanstd(C_L_full)/sqrt(sum(~isnan(C_L_full)))...
    nanstd(C_Post1_full)/sqrt(sum(~isnan(C_Post1_full)))...
    nanstd(C_Post2_full)/sqrt(sum(~isnan(C_Post2_full)))...
%     nanstd(C_Post3_full)/sqrt(sum(~isnan(C_Post3_full)))...
    ],'.','LineWidth',2)

errorbar([2 5 8 11 ],[nanmean(L_Pre_full) nanmean(L_L_full) nanmean(L_Post1_full) nanmean(L_Post2_full) ],...
    [nanstd(L_Pre_full)/sqrt(sum(~isnan(L_Pre_full)))...
    nanstd(L_L_full)/sqrt(sum(~isnan(L_L_full)))...
    nanstd(L_Post1_full)/sqrt(sum(~isnan(L_Post1_full)))...
    nanstd(C_Post2_full)/sqrt(sum(~isnan(C_Post2_full)))...
%     nanstd(C_Post3_full)/sqrt(sum(~isnan(C_Post3_full)))...
    ],'.','LineWidth',2)

errorbar([3 6 9 12 ],[nanmean(Y_Pre_full) nanmean(Y_L_full) nanmean(Y_Post1_full) nanmean(Y_Post2_full) ],...
    [nanstd(Y_Pre_full)/sqrt(sum(~isnan(Y_Pre_full)))...
    nanstd(Y_L_full)/sqrt(sum(~isnan(Y_L_full)))...
    nanstd(Y_Post1_full)/sqrt(sum(~isnan(Y_Post1_full)))...
    nanstd(Y_Post2_full)/sqrt(sum(~isnan(Y_Post2_full)))...
%     nanstd(Y_Post3_full)/sqrt(sum(~isnan(Y_Post3_full)))...
    ],'.','LineWidth',2)

xticks([2 5 8 11 ])
xticklabels({'Pre 1-3','Opto','Post 1-3','Post 4-7'})
set(gca,'fontsize',12,'FontWeight','bold')
ylim([0.5 1.5])

title (['normalized ' para ' Layer' Sinks{i1}])
legend('Control','Opto','YFP','Location','best')
end

% %%%keyboard



%%% Fig 6?
clear
M = [1, 4, 9]
SORT = 'ST_based';
BFlim= 2;

c = load('K:\CSD_dynamic_analysis\DATA\Output\Output_Input_Control_7post_n7_Data_Threshold_25_Zscore_0_binned_1_mirror_0.mat');
BF_C = c.Data.BF_Pos;
l = load('K:\CSD_dynamic_analysis\DATA\Output\Output_Input_HighP_7post_Data_Threshold_25_Zscore_0_binned_1_mirror_0.mat');
BF_L = l.Data.BF_Pos;
y = load('K:\CSD_dynamic_analysis\DATA\Output\Output_Input_YFP_7post_Data_Threshold_25_Zscore_0_binned_1_mirror_0.mat');
BF_Y = y.Data.BF_Pos;

% CNorm = nan(length(c.Data.names),3);
% for i1 = 1:3
%     CNorm(:,i1,1) =  c.Data.(SORT)(i1).Full_RMS_AVREC(:,BF_C+BFlim);%full
%     CNorm(:,i1,2) =  c.Data.(SORT)(i1).Early_RMS_AVREC(:,BF_C+BFlim);%early
%     CNorm(:,i1,3) =  c.Data.(SORT)(i1).Late_RMS_AVREC(:,BF_C+BFlim);%late
% end
% CNorm = nanmean(CNorm,2);
% 
% LNorm = nan(length(l.Data.names),3);
% for i1 = 1:3
%     LNorm(:,i1,1) =  l.Data.(SORT)(i1).Full_RMS_AVREC(:,BF_L+BFlim);%full
%     LNorm(:,i1,2) =  l.Data.(SORT)(i1).Early_RMS_AVREC(:,BF_L+BFlim);%early
%     LNorm(:,i1,3) =  l.Data.(SORT)(i1).Late_RMS_AVREC(:,BF_L+BFlim);%late
% end
% LNorm = nanmean(LNorm,2);
% 
% YNorm = nan(length(y.Data.names),3);
% for i1 = 1:3
%     YNorm(:,i1,1) =  y.Data.(SORT)(i1).Full_RMS_AVREC(:,BF_Y+BFlim);%full
%     YNorm(:,i1,2) =  y.Data.(SORT)(i1).Early_RMS_AVREC(:,BF_Y+BFlim);%early
%     YNorm(:,i1,3) =  y.Data.(SORT)(i1).Late_RMS_AVREC(:,BF_Y+BFlim);%late
% end
% YNorm = nanmean(YNorm,2);

%Pre
C_Pre_full_AVREC = c.Data.(SORT)(M(1)).Full_RMS_AVREC(:,BF_C-BFlim:BF_C+BFlim);%./CNorm(:,:,1);%Full
C_Pre_early_AVREC = c.Data.(SORT)(M(1)).Early_RMS_AVREC(:,BF_C-BFlim:BF_C+BFlim);%./CNorm(:,:,2);%early
C_Pre_late_AVREC = c.Data.(SORT)(M(1)).Late_RMS_AVREC(:,BF_C-BFlim:BF_C+BFlim);%./CNorm(:,:,3);%late

L_Pre_full_AVREC = l.Data.(SORT)(M(1)).Full_RMS_AVREC(:,BF_C-BFlim:BF_L+BFlim);%./LNorm(:,:,1);%Full
L_Pre_early_AVREC = l.Data.(SORT)(M(1)).Early_RMS_AVREC(:,BF_C-BFlim:BF_L+BFlim);%./LNorm(:,:,2);%early
L_Pre_late_AVREC = l.Data.(SORT)(M(1)).Late_RMS_AVREC(:,BF_C-BFlim:BF_L+BFlim);%./LNorm(:,:,3);%late

Y_Pre_full_AVREC = y.Data.(SORT)(M(1)).Full_RMS_AVREC(:,BF_C-BFlim:BF_Y+BFlim);%./YNorm(:,:,1);%Full
Y_Pre_early_AVREC = y.Data.(SORT)(M(1)).Early_RMS_AVREC(:,BF_C-BFlim:BF_Y+BFlim);%./YNorm(:,:,2);%early
Y_Pre_late_AVREC = y.Data.(SORT)(M(1)).Late_RMS_AVREC(:,BF_C-BFlim:BF_Y+BFlim);%./YNorm(:,:,3);%late

%Laser
C_L_full_AVREC = c.Data.(SORT)(M(2)).Full_RMS_AVREC(:,BF_C-BFlim:BF_C+BFlim);%./CNorm(:,:,1);%Full
C_L_early_AVREC = c.Data.(SORT)(M(2)).Early_RMS_AVREC(:,BF_C-BFlim:BF_C+BFlim);%./CNorm(:,:,2);%early
C_L_late_AVREC = c.Data.(SORT)(M(2)).Late_RMS_AVREC(:,BF_C-BFlim:BF_C+BFlim);%./CNorm(:,:,3);%late

L_L_full_AVREC = l.Data.(SORT)(M(2)).Full_RMS_AVREC(:,BF_C-BFlim:BF_L+BFlim);%./LNorm(:,:,1);%Full
L_L_early_AVREC = l.Data.(SORT)(M(2)).Early_RMS_AVREC(:,BF_C-BFlim:BF_L+BFlim);%./LNorm(:,:,2);%early
L_L_late_AVREC = l.Data.(SORT)(M(2)).Late_RMS_AVREC(:,BF_C-BFlim:BF_L+BFlim);%./LNorm(:,:,3);%late

Y_L_full_AVREC = y.Data.(SORT)(M(2)).Full_RMS_AVREC(:,BF_C-BFlim:BF_Y+BFlim);%./YNorm(:,:,1);%Full
Y_L_early_AVREC = y.Data.(SORT)(M(2)).Early_RMS_AVREC(:,BF_C-BFlim:BF_Y+BFlim);%./YNorm(:,:,2);%early
Y_L_late_AVREC = y.Data.(SORT)(M(2)).Late_RMS_AVREC(:,BF_C-BFlim:BF_Y+BFlim);%./YNorm(:,:,3);%late


%Post
C_Post_full_AVREC = c.Data.(SORT)(M(3)).Full_RMS_AVREC(:,BF_C+BFlim);%./CNorm(:,:,1);%Full
C_Post_early_AVREC = c.Data.(SORT)(M(3)).Early_RMS_AVREC(:,BF_C+BFlim);%./CNorm(:,:,2);%early
C_Post_late_AVREC = c.Data.(SORT)(M(3)).Late_RMS_AVREC(:,BF_C+BFlim);%./CNorm(:,:,3);%late

L_Post_full_AVREC = l.Data.(SORT)(M(3)).Full_RMS_AVREC(:,BF_L+BFlim);%./LNorm(:,:,1);%Full
L_Post_early_AVREC = l.Data.(SORT)(M(3)).Early_RMS_AVREC(:,BF_L+BFlim);%./LNorm(:,:,2);%early
L_Post_late_AVREC = l.Data.(SORT)(M(3)).Late_RMS_AVREC(:,BF_L+BFlim);%./LNorm(:,:,3);%late

Y_Post_full_AVREC = y.Data.(SORT)(M(3)).Full_RMS_AVREC(:,BF_Y+BFlim);%./YNorm(:,:,1);%Full
Y_Post_early_AVREC = y.Data.(SORT)(M(3)).Early_RMS_AVREC(:,BF_Y+BFlim);%./YNorm(:,:,2);%early
Y_Post_late_AVREC = y.Data.(SORT)(M(3)).Late_RMS_AVREC(:,BF_Y+BFlim);%./YNorm(:,:,3);%late


%Pre
C_Pre_full_RELRES = c.Data.(SORT)(M(1)).Full_RMS_RELRES(:,BF_C-BFlim:BF_C+BFlim);%./CNorm(:,:,1);%Full
C_Pre_early_RELRES = c.Data.(SORT)(M(1)).Early_RMS_RELRES(:,BF_C-BFlim:BF_C+BFlim);%./CNorm(:,:,2);%early
C_Pre_late_RELRES = c.Data.(SORT)(M(1)).Late_RMS_RELRES(:,BF_C-BFlim:BF_C+BFlim);%./CNorm(:,:,3);%late

L_Pre_full_RELRES = l.Data.(SORT)(M(1)).Full_RMS_RELRES(:,BF_C-BFlim:BF_L+BFlim);%./LNorm(:,:,1);%Full
L_Pre_early_RELRES = l.Data.(SORT)(M(1)).Early_RMS_RELRES(:,BF_C-BFlim:BF_L+BFlim);%./LNorm(:,:,2);%early
L_Pre_late_RELRES = l.Data.(SORT)(M(1)).Late_RMS_RELRES(:,BF_C-BFlim:BF_L+BFlim);%./LNorm(:,:,3);%late

Y_Pre_full_RELRES = y.Data.(SORT)(M(1)).Full_RMS_RELRES(:,BF_C-BFlim:BF_Y+BFlim);%./YNorm(:,:,1);%Full
Y_Pre_early_RELRES = y.Data.(SORT)(M(1)).Early_RMS_RELRES(:,BF_C-BFlim:BF_Y+BFlim);%./YNorm(:,:,2);%early
Y_Pre_late_RELRES = y.Data.(SORT)(M(1)).Late_RMS_RELRES(:,BF_C-BFlim:BF_Y+BFlim);%./YNorm(:,:,3);%late

%Laser
C_L_full_RELRES = c.Data.(SORT)(M(2)).Full_RMS_RELRES(:,BF_C-BFlim:BF_C+BFlim);%./CNorm(:,:,1);%Full
C_L_early_RELRES = c.Data.(SORT)(M(2)).Early_RMS_RELRES(:,BF_C-BFlim:BF_C+BFlim);%./CNorm(:,:,2);%early
C_L_late_RELRES = c.Data.(SORT)(M(2)).Late_RMS_RELRES(:,BF_C-BFlim:BF_C+BFlim);%./CNorm(:,:,3);%late

L_L_full_RELRES = l.Data.(SORT)(M(2)).Full_RMS_RELRES(:,BF_C-BFlim:BF_L+BFlim);%./LNorm(:,:,1);%Full
L_L_early_RELRES = l.Data.(SORT)(M(2)).Early_RMS_RELRES(:,BF_C-BFlim:BF_L+BFlim);%./LNorm(:,:,2);%early
L_L_late_RELRES = l.Data.(SORT)(M(2)).Late_RMS_RELRES(:,BF_C-BFlim:BF_L+BFlim);%./LNorm(:,:,3);%late

Y_L_full_RELRES = y.Data.(SORT)(M(2)).Full_RMS_RELRES(:,BF_C-BFlim:BF_Y+BFlim);%./YNorm(:,:,1);%Full
Y_L_early_RELRES = y.Data.(SORT)(M(2)).Early_RMS_RELRES(:,BF_C-BFlim:BF_Y+BFlim);%./YNorm(:,:,2);%early
Y_L_late_RELRES = y.Data.(SORT)(M(2)).Late_RMS_RELRES(:,BF_C-BFlim:BF_Y+BFlim);%./YNorm(:,:,3);%late


%Post
C_Post_full_RELRES = c.Data.(SORT)(M(3)).Full_RMS_RELRES(:,BF_C+BFlim);%./CNorm(:,:,1);%Full
C_Post_early_RELRES = c.Data.(SORT)(M(3)).Early_RMS_RELRES(:,BF_C+BFlim);%./CNorm(:,:,2);%early
C_Post_late_RELRES = c.Data.(SORT)(M(3)).Late_RMS_RELRES(:,BF_C+BFlim);%./CNorm(:,:,3);%late

L_Post_full_RELRES = l.Data.(SORT)(M(3)).Full_RMS_RELRES(:,BF_L+BFlim);%./LNorm(:,:,1);%Full
L_Post_early_RELRES = l.Data.(SORT)(M(3)).Early_RMS_RELRES(:,BF_L+BFlim);%./LNorm(:,:,2);%early
L_Post_late_RELRES = l.Data.(SORT)(M(3)).Late_RMS_RELRES(:,BF_L+BFlim);%./LNorm(:,:,3);%late

Y_Post_full_RELRES = y.Data.(SORT)(M(3)).Full_RMS_RELRES(:,BF_Y+BFlim);%./YNorm(:,:,1);%Full
Y_Post_early_RELRES = y.Data.(SORT)(M(3)).Early_RMS_RELRES(:,BF_Y+BFlim);%./YNorm(:,:,2);%early
Y_Post_late_RELRES = y.Data.(SORT)(M(3)).Late_RMS_RELRES(:,BF_Y+BFlim);%./YNorm(:,:,3);%late

figure
subplot(1,3,1)
I = ~isnan(C_Pre_full_AVREC) & ~isnan(C_Pre_full_RELRES);
scatter(C_Pre_full_RELRES(I),C_Pre_full_AVREC(I),'MarkerEdgeColor',[0 0 1])

hold on
[rho,pval] = corr(C_Pre_full_RELRES(I),C_Pre_full_AVREC(I),'Type','Pearson');                    
                    fitvars = polyfit(C_Pre_full_RELRES(I),C_Pre_full_AVREC(I) , 1);
                    m = fitvars(1);
                    c = fitvars(2);                   
                    
                    Xnan2 =vertcat(0, C_Pre_full_RELRES(I));
                    Y_Fit = polyval(fitvars, Xnan2);
                    plot(Xnan2,Y_Fit,'Color','b');  
%                     message = sprintf(['y = ' num2str(m) '*x +' num2str(c) ...
%                         '\nlinear Correlation (Pearson) R = ' num2str(rho) ' p = ' num2str(pval)]);
%                     text(min(C_Pre_full_RELRES(I)),max(Y_Fit),message);


I = ~isnan(L_Pre_full_AVREC) & ~isnan(L_Pre_full_RELRES);
scatter(L_Pre_full_RELRES(I),L_Pre_full_AVREC(I),'MarkerEdgeColor',[0 1 0])

[rho,pval] = corr(L_Pre_full_RELRES(I),L_Pre_full_AVREC(I),'Type','Pearson');                    
                    fitvars = polyfit(L_Pre_full_RELRES(I),L_Pre_full_AVREC(I) , 1);
                    m = fitvars(1);
                    c = fitvars(2);                   
                    
                    Xnan2 =vertcat(0, L_Pre_full_RELRES(I));
                    Y_Fit = polyval(fitvars, Xnan2);
                    plot(Xnan2,Y_Fit,'Color','g');  


I = ~isnan(Y_Pre_full_AVREC) & ~isnan(Y_Pre_full_RELRES);
scatter(Y_Pre_full_RELRES(I),Y_Pre_full_AVREC(I),'MarkerEdgeColor',[1 0 0])


[rho,pval] = corr(Y_Pre_full_RELRES(I),Y_Pre_full_AVREC(I),'Type','Pearson');                    
                    fitvars = polyfit(Y_Pre_full_RELRES(I),Y_Pre_full_AVREC(I) , 1);
                    m = fitvars(1);
                    c = fitvars(2);                   
                    
                    Xnan2 =vertcat(0, Y_Pre_full_RELRES(I));
                    Y_Fit = polyval(fitvars, Xnan2);
                    plot(Xnan2,Y_Fit,'Color','r');  
                    
xlabel('Full RMS RELRES')
ylabel('Full RMS AVREC')
title('Full')
legend('Control','','Opto','','YFP','')




subplot(1,3,2)
I = ~isnan(C_Pre_early_AVREC) & ~isnan(C_Pre_early_RELRES);
scatter(C_Pre_early_RELRES(I),C_Pre_early_AVREC(I),'MarkerEdgeColor',[0 0 1])

hold on
[rho,pval] = corr(C_Pre_early_RELRES(I),C_Pre_early_AVREC(I),'Type','Pearson');                    
                    fitvars = polyfit(C_Pre_early_RELRES(I),C_Pre_early_AVREC(I) , 1);
                    m = fitvars(1);
                    c = fitvars(2);                   
                    
                    Xnan2 =vertcat(0, C_Pre_early_RELRES(I));
                    Y_Fit = polyval(fitvars, Xnan2);
                    plot(Xnan2,Y_Fit,'Color','b');  
%                     message = sprintf(['y = ' num2str(m) '*x +' num2str(c) ...
%                         '\nlinear Correlation (Pearson) R = ' num2str(rho) ' p = ' num2str(pval)]);
%                     text(min(C_Pre_early_RELRES(I)),max(Y_Fit),message);


I = ~isnan(L_Pre_early_AVREC) & ~isnan(L_Pre_early_RELRES);
scatter(L_Pre_early_RELRES(I),L_Pre_early_AVREC(I),'MarkerEdgeColor',[0 1 0])

[rho,pval] = corr(L_Pre_early_RELRES(I),L_Pre_early_AVREC(I),'Type','Pearson');                    
                    fitvars = polyfit(L_Pre_early_RELRES(I),L_Pre_early_AVREC(I) , 1);
                    m = fitvars(1);
                    c = fitvars(2);                   
                    
                    Xnan2 =vertcat(0, L_Pre_early_RELRES(I));
                    Y_Fit = polyval(fitvars, Xnan2);
                    plot(Xnan2,Y_Fit,'Color','g');  


I = ~isnan(Y_Pre_early_AVREC) & ~isnan(Y_Pre_early_RELRES);
scatter(Y_Pre_early_RELRES(I),Y_Pre_early_AVREC(I),'MarkerEdgeColor',[1 0 0])


[rho,pval] = corr(Y_Pre_early_RELRES(I),Y_Pre_early_AVREC(I),'Type','Pearson');                    
                    fitvars = polyfit(Y_Pre_early_RELRES(I),Y_Pre_early_AVREC(I) , 1);
                    m = fitvars(1);
                    c = fitvars(2);                   
                    
                    Xnan2 =vertcat(0, Y_Pre_early_RELRES(I));
                    Y_Fit = polyval(fitvars, Xnan2);
                    plot(Xnan2,Y_Fit,'Color','r');  
                    
xlabel('early RMS RELRES')
ylabel('early RMS AVREC')
title('early')
legend('Control','','Opto','','YFP','')



subplot(1,3,3)
I = ~isnan(C_Pre_late_AVREC) & ~isnan(C_Pre_late_RELRES);
scatter(C_Pre_late_RELRES(I),C_Pre_late_AVREC(I),'MarkerEdgeColor',[0 0 1])

hold on
[rho,pval] = corr(C_Pre_late_RELRES(I),C_Pre_late_AVREC(I),'Type','Pearson');                    
                    fitvars = polyfit(C_Pre_late_RELRES(I),C_Pre_late_AVREC(I) , 1);
                    m = fitvars(1);
                    c = fitvars(2);                   
%                     %%%keyboard
                    Xnan2 =vertcat(0, C_Pre_late_RELRES(I));
                    Y_Fit = polyval(fitvars, Xnan2);
                    plot(Xnan2,Y_Fit,'Color','b');  
%                     message = sprintf(['y = ' num2str(m) '*x +' num2str(c) ...
%                         '\nlinear Correlation (Pearson) R = ' num2str(rho) ' p = ' num2str(pval)]);
%                     text(min(C_Pre_late_RELRES(I)),max(Y_Fit),message);


I = ~isnan(L_Pre_late_AVREC) & ~isnan(L_Pre_late_RELRES);
scatter(L_Pre_late_RELRES(I),L_Pre_late_AVREC(I),'MarkerEdgeColor',[0 1 0])

[rho,pval] = corr(L_Pre_late_RELRES(I),L_Pre_late_AVREC(I),'Type','Pearson');                    
                    fitvars = polyfit(L_Pre_late_RELRES(I),L_Pre_late_AVREC(I) , 1);
                    m = fitvars(1);
                    c = fitvars(2);                   
                    
                    Xnan2 =vertcat(0, L_Pre_late_RELRES(I));
                    Y_Fit = polyval(fitvars, Xnan2);
                    plot(Xnan2,Y_Fit,'Color','g');  


I = ~isnan(Y_Pre_late_AVREC) & ~isnan(Y_Pre_late_RELRES);
scatter(Y_Pre_late_RELRES(I),Y_Pre_late_AVREC(I),'MarkerEdgeColor',[1 0 0])


[rho,pval] = corr(Y_Pre_late_RELRES(I),Y_Pre_late_AVREC(I),'Type','Pearson');                    
                    fitvars = polyfit(Y_Pre_late_RELRES(I),Y_Pre_late_AVREC(I) , 1);
                    m = fitvars(1);
                    c = fitvars(2);                   
                    
                    Xnan2 =vertcat(0, Y_Pre_late_RELRES(I));
                    Y_Fit = polyval(fitvars, Xnan2);
                    plot(Xnan2,Y_Fit,'Color','r');  
                    
xlabel('late RMS RELRES')
ylabel('late RMS AVREC')
title('late')
legend('Control','','Opto','','YFP','')

 clear
%  close all


SORT = 'ST_based';
Para =  'SinkPeakAmp';% ,'SinkRMS', 'tempSinkRMS''SinkDur', 'SinkPeakLate','Sinkonset'  
% Sinks = {'I_IIL'    'VIL'};%'IVE'  'VaE' 'InfE'Sinks = {'I_IIL'  'IVE'  'VaE' 'InfE'  'VIL'};
Sinks = {'I_IIL' 'IVE'  'VaE' 'VbE' 'VIE' 'InfE' 'VIL'};
Order = {'l', 'c','y' };%,'y'


c = load('K:\CSD_dynamic_analysis\DATA\Output\Output_Input_Control_7post_n7_Data_Threshold_25_Zscore_0_binned_1_mirror_1.mat');
BF_c = c.Data.BF_Pos;
l = load('K:\CSD_dynamic_analysis\DATA\Output\Output_Input_HighP_7post_Data_Threshold_25_Zscore_0_binned_1_mirror_1.mat');
BF_l = l.Data.BF_Pos;
y = load('K:\CSD_dynamic_analysis\DATA\Output\Output_Input_YFP_7post_Data_Threshold_25_Zscore_0_binned_1_mirror_1.mat');
BF_y = y.Data.BF_Pos;

nl = length(l.Data.names);

try
  ny = length(y.Data.names);  
catch
  ny = 0;
end
 try
  nc = length(c.Data.names);  
catch
  nc = 0;
 end   

 
nAnimal = sum([nl,ny, nc]);
PreMatrix = nan (nAnimal,11,length(Sinks),3,length(Order));% #animals, 11 entries, #Sinks, BF/near/non, l/C/y

n = 0;

for i1 = 1:length(Order)% l, c, x
    for i2 = 1:11
        for i3 = 1:eval(['n' Order{i1}])% animals
            for i4 = 1:length(Sinks)% sinks
                X =eval([Order{i1} '.Data.' (SORT) '(i2).' (Para) '.' Sinks{i4}]);
                BF = eval(['BF_'  Order{i1} ]);
                PreMatrix(1:eval(['n' Order{i1}]),i2,i4,1,i1) = X(:,BF);
                PreMatrix(1:eval(['n' Order{i1}]),i2,i4,2,i1) = X(:,BF+1);
                PreMatrix(1:eval(['n' Order{i1}]),i2,i4,3,i1) = X(:,BF+2);
            end
        end
    end
end

Pre = PreMatrix(:,1:3,:,:,:);
Pre=nanmean(Pre,2);
figure('name',Para)
Farbe = {'r','b','g'};

for i1 = 1:length(Sinks)
    for i2 = 1:3
%  %%keyboard
        if i1 == 1
        subplot(length(Sinks),3,i2)
        elseif i1 == 2
        subplot(length(Sinks),3,i2+(3))
        else
        subplot(length(Sinks),3,i2+(3)*(i1-1))    
        end
        hold on
        common = [];
        for i3=1:length(Order)
           X= PreMatrix(:,:,i1,i2,i3)./Pre(:,:,i1,i2,i3);
           if i3 == 1
               S1 = X;
           else
               S2 = X;
           end
           M =nanmean(X);
           SEM = nanstd(X)./sqrt(nansum(~isnan(X)));
           shadedErrorBar([],M, SEM,'lineprops',{Farbe{i3},'markerfacecolor',Farbe{i3}})
          
           common = vertcat(common, X(:,1:3));
        end
            L =[2 11];  
            S3 = [S1, S2];
            O = teg_repeated_measures_ANOVA(S3,L,{'Group','Measurement'});
            [p,n]=ttest(S1,S2);
%         [stat,pperm,hperm,statperm]= elbe_permtest_005(S1,S2,1,'raw',1000,0,0.05,0);
%         SigPos=find(pperm <= 0.025 | pperm >= 1.025);
        %keyboard
        
T1 = [O.labels{1} ' F(' num2str(O.R(1,2)) ','  num2str(O.R(1,3)) ') = ' num2str(O.R(1,1)) ' p = ' num2str(O.R(1,4))];
T2 = [O.labels{2} ' F(' num2str(O.R(2,2)) ','  num2str(O.R(2,3)) ') = ' num2str(O.R(2,1)) ' p = ' num2str(O.R(2,4))];
T3 = [O.labels{3} ' F(' num2str(O.R(3,2)) ','  num2str(O.R(3,3)) ') = ' num2str(O.R(3,1)) ' p = ' num2str(O.R(3,4))];
T =strvcat(T1, T2, T3);
t = text(1.1,1.1, T,'Fontsize',10);
% s = t.Fontsize;
% t.Fontsize = 6;

        ALLSTD = nanmean(nanstd(common'));
        %ALLSTD = nanmean(nanstd(common'))/sqrt(sum(~isnan(nanmean(common'))));
        ALLMean = nanmean(nanmean(common'));
        plot([1 11],[1 1],':k');
        plot([1:11],p*1.7,'k*')
%         plot([1 11],[ALLMean-ALLSTD ALLMean-ALLSTD],':k');
%         if isfield(O,'groups')
%              keyboard
%         end
%         plot([SigPos],[ones(size(SigPos))*1.5],'k*')
        ylim([0.5 2])
        xlim([1 11])
        if i2 ==1
        ylabel(Sinks{i1})
        end        

        xticks([1 2 3 4 5 6 7 8 9 10 11])
        xticklabels({'Pre(n-2)','Pre(n-1)', 'Pre(n)','Laser', 'Post 1','Post 2', 'Post 3','Post 4','Post 5','Post 6','Post 7'})
        xtickangle(270)
        if i1 == 1 && i2 == 1
            legend({'L','','','','C','','','','r','','',''},'FontSize',6,'Location','best')
        end
        
        if  i2 == 1
            title('BF')
        elseif  i2 == 2
            title('± near BF')       
        else
            title('± non BF')  
        end

    end
end

clear


load('K:\CSD_dynamic_analysis\DATA\Output\Output_Input_ALL_OPTO_7post_Data_Threshold_25_Zscore_0_binned_1_mirror_1.mat')
BF =Data.BF_Pos;
SORT = 'ST_based'; 
Para = 'Full_RMS_RELRES';
M = [1, 4, 9]
IDs = Data.names;
LP = [];
 for i1=1:length(IDs)
    X = GOTs_LPs(IDs{i1});
    LP = [LP X];
 end
LP = LP';

Norm = [];

for i1 = 1:3
    Norm(:,i1,1) = Data.(SORT)(i1).(Para)(:,BF); 
    Norm(:,i1,2) = Data.(SORT)(i1).(Para)(:,BF+1);
    Norm(:,i1,3) = Data.(SORT)(i1).(Para)(:,BF+2);
end

Norm = nanmean(Norm,2);

Pre_BFNorm = Data.(SORT)(M (1)).(Para)(:,BF)./Norm(:,:,1);
Pre_nearBFNorm = Data.(SORT)(M (1)).(Para)(:,BF)./Norm(:,:,2);
Pre_nonBFNorm = Data.(SORT)(M (1)).(Para)(:,BF)./Norm(:,:,3);

L_BFNorm = Data.(SORT)(M (2)).(Para)(:,BF)./Norm(:,:,1);
L_nearBFNorm = Data.(SORT)(M (2)).(Para)(:,BF)./Norm(:,:,2);
L_nonBFNorm = Data.(SORT)(M (2)).(Para)(:,BF)./Norm(:,:,3);

Post_BFNorm = Data.(SORT)(M (3)).(Para)(:,BF)./Norm(:,:,1);
Post_nearBFNorm = Data.(SORT)(M (3)).(Para)(:,BF)./Norm(:,:,2);
Post_nonBFNorm = Data.(SORT)(M (3)).(Para)(:,BF)./Norm(:,:,3);

figure
subplot(2,3,1)
hold on
I = ~isnan(LP) &  ~isnan(Pre_BFNorm);
scatter(LP(I), Pre_BFNorm(I),'MarkerEdgeColor',[1 0 0])

[rho,pval] = corr(LP(I), Pre_BFNorm(I),'Type','Pearson');                    
                    fitvars = polyfit(LP(I), Pre_BFNorm(I) , 1);
                    m = fitvars(1);
                    c = fitvars(2);                   
                    
                    Xnan2 =vertcat(0, LP);
                    Y_Fit = polyval(fitvars, Xnan2);
                    plot(Xnan2,Y_Fit,'Color','r'); 
                    
I = ~isnan(LP) &  ~isnan(L_BFNorm);                   
scatter(LP(I), L_BFNorm(I),'MarkerEdgeColor',[0 1 0])

[rho,pval] = corr(LP(I), L_BFNorm(I),'Type','Pearson');                    
                    fitvars = polyfit(LP(I), L_BFNorm(I) , 1);
                    m = fitvars(1);
                    c = fitvars(2);                   
                    
                    Xnan2 =vertcat(0, LP);
                    Y_Fit = polyval(fitvars, Xnan2);
                    plot(Xnan2,Y_Fit,'Color','g');   
                    
I = ~isnan(LP) &  ~isnan(Post_BFNorm);                    
scatter(LP(I), Post_BFNorm(I),'MarkerEdgeColor',[0 0 1])

[rho,pval] = corr(LP(I), Post_BFNorm(I),'Type','Pearson');                    
                    fitvars = polyfit(LP(I), Post_BFNorm(I) , 1);
                    m = fitvars(1);
                    c = fitvars(2);                   
                    
                    Xnan2 =vertcat(0, LP);
                    Y_Fit = polyval(fitvars, Xnan2);
                    plot(Xnan2,Y_Fit,'Color','b');  
                    
legend('Pre','','Opto','','Post','','Location','best')
 title('BF')                                   
                    
subplot(2,3,2)
hold on
I = ~isnan(LP) &  ~isnan(Pre_nearBFNorm);
scatter(LP(I), Pre_nearBFNorm(I),'MarkerEdgeColor',[1 0 0])

[rho,pval] = corr(LP(I), Pre_nearBFNorm(I),'Type','Pearson');                    
                    fitvars = polyfit(LP(I), Pre_nearBFNorm(I) , 1);
                    m = fitvars(1);
                    c = fitvars(2);                   
                    
                    Xnan2 =vertcat(0, LP);
                    Y_Fit = polyval(fitvars, Xnan2);
                    plot(Xnan2,Y_Fit,'Color','r'); 
                    
I = ~isnan(LP) &  ~isnan(L_nearBFNorm);                   
scatter(LP(I), L_nearBFNorm(I),'MarkerEdgeColor',[0 1 0])

[rho,pval] = corr(LP(I), L_nearBFNorm(I),'Type','Pearson');                    
                    fitvars = polyfit(LP(I), L_nearBFNorm(I) , 1);
                    m = fitvars(1);
                    c = fitvars(2);                   
                    
                    Xnan2 =vertcat(0, LP);
                    Y_Fit = polyval(fitvars, Xnan2);
                    plot(Xnan2,Y_Fit,'Color','g');   
                    
I = ~isnan(LP) &  ~isnan(Post_nearBFNorm);                    
scatter(LP(I), Post_nearBFNorm(I),'MarkerEdgeColor',[0 0 1])

[rho,pval] = corr(LP(I), Post_nearBFNorm(I),'Type','Pearson');                    
                    fitvars = polyfit(LP(I), Post_nearBFNorm(I) , 1);
                    m = fitvars(1);
                    c = fitvars(2);                   
                    
                    Xnan2 =vertcat(0, LP);
                    Y_Fit = polyval(fitvars, Xnan2);
                    plot(Xnan2,Y_Fit,'Color','b');  
                    
legend('Pre','','Opto','','Post','','Location','best')
 title('nearBF')                       

subplot(2,3,3)
hold on
I = ~isnan(LP) &  ~isnan(Pre_nonBFNorm);
scatter(LP(I), Pre_nonBFNorm(I),'MarkerEdgeColor',[1 0 0])

[rho,pval] = corr(LP(I), Pre_nonBFNorm(I),'Type','Pearson');                    
                    fitvars = polyfit(LP(I), Pre_nonBFNorm(I) , 1);
                    m = fitvars(1);
                    c = fitvars(2);                   
                    
                    Xnan2 =vertcat(0, LP);
                    Y_Fit = polyval(fitvars, Xnan2);
                    plot(Xnan2,Y_Fit,'Color','r'); 
                    
I = ~isnan(LP) &  ~isnan(L_nonBFNorm);                   
scatter(LP(I), L_nonBFNorm(I),'MarkerEdgeColor',[0 1 0])

[rho,pval] = corr(LP(I), L_nonBFNorm(I),'Type','Pearson');                    
                    fitvars = polyfit(LP(I), L_nonBFNorm(I) , 1);
                    m = fitvars(1);
                    c = fitvars(2);                   
                    
                    Xnan2 =vertcat(0, LP);
                    Y_Fit = polyval(fitvars, Xnan2);
                    plot(Xnan2,Y_Fit,'Color','g');   
                    
I = ~isnan(LP) &  ~isnan(Post_nonBFNorm);                    
scatter(LP(I), Post_nonBFNorm(I),'MarkerEdgeColor',[0 0 1])

[rho,pval] = corr(LP(I), Post_nonBFNorm(I),'Type','Pearson');                    
                    fitvars = polyfit(LP(I), Post_nonBFNorm(I) , 1);
                    m = fitvars(1);
                    c = fitvars(2);                   
                    
                    Xnan2 =vertcat(0, LP);
                    Y_Fit = polyval(fitvars, Xnan2);
                    plot(Xnan2,Y_Fit,'Color','b');  
                    
legend('Pre','','Opto','','Post','','Location','best')
 title('nonBF')    
 
 %%%% Z-Score on Pres
Pre_BFNorm = (Pre_BFNorm-nanmean(Pre_BFNorm))/nanstd(Pre_BFNorm);
Pre_nearBFNorm = (Pre_nearBFNorm-nanmean(Pre_nearBFNorm))/nanstd(Pre_nearBFNorm);
Pre_nonBFNorm = (Pre_nonBFNorm-nanmean(Pre_nonBFNorm))/nanstd(Pre_nonBFNorm);

L_BFNorm = (L_BFNorm-nanmean(Pre_BFNorm))/nanstd(Pre_BFNorm);
L_nearBFNorm = (L_nearBFNorm-nanmean(Pre_nearBFNorm))/nanstd(Pre_nearBFNorm);
L_nonBFNorm = (L_nonBFNorm-nanmean(Pre_nonBFNorm))/nanstd(Pre_nonBFNorm);

Post_BFNorm = (Post_BFNorm-nanmean(Pre_BFNorm))/nanstd(Pre_BFNorm);
Post_nearBFNorm = (Post_nearBFNorm-nanmean(Pre_nearBFNorm))/nanstd(Pre_nearBFNorm);
Post_nonBFNorm = (Post_nonBFNorm-nanmean(Pre_nonBFNorm))/nanstd(Pre_nonBFNorm);
 
subplot(2,3,4)
hold on
I = ~isnan(LP) &  ~isnan(Pre_BFNorm);
scatter(LP(I), Pre_BFNorm(I),'MarkerEdgeColor',[1 0 0])

[rho,pval] = corr(LP(I), Pre_BFNorm(I),'Type','Pearson');                    
                    fitvars = polyfit(LP(I), Pre_BFNorm(I) , 1);
                    m = fitvars(1);
                    c = fitvars(2);                   
                    
                    Xnan2 =vertcat(0, LP);
                    Y_Fit = polyval(fitvars, Xnan2);
                    plot(Xnan2,Y_Fit,'Color','r'); 
                    
I = ~isnan(LP) &  ~isnan(L_BFNorm);                   
scatter(LP(I), L_BFNorm(I),'MarkerEdgeColor',[0 1 0])

[rho,pval] = corr(LP(I), L_BFNorm(I),'Type','Pearson');                    
                    fitvars = polyfit(LP(I), L_BFNorm(I) , 1);
                    m = fitvars(1);
                    c = fitvars(2);                   
                    
                    Xnan2 =vertcat(0, LP);
                    Y_Fit = polyval(fitvars, Xnan2);
                    plot(Xnan2,Y_Fit,'Color','g');   
                    
I = ~isnan(LP) &  ~isnan(Post_BFNorm);                    
scatter(LP(I), Post_BFNorm(I),'MarkerEdgeColor',[0 0 1])

[rho,pval] = corr(LP(I), Post_BFNorm(I),'Type','Pearson');                    
                    fitvars = polyfit(LP(I), Post_BFNorm(I) , 1);
                    m = fitvars(1);
                    c = fitvars(2);                   
                    
                    Xnan2 =vertcat(0, LP);
                    Y_Fit = polyval(fitvars, Xnan2);
                    plot(Xnan2,Y_Fit,'Color','b');  
                    
legend('Pre','','Opto','','Post','','Location','best')
 title('Z-scored BF')                                   
                    
subplot(2,3,5)
hold on
I = ~isnan(LP) &  ~isnan(Pre_nearBFNorm);
scatter(LP(I), Pre_nearBFNorm(I),'MarkerEdgeColor',[1 0 0])

[rho,pval] = corr(LP(I), Pre_nearBFNorm(I),'Type','Pearson');                    
                    fitvars = polyfit(LP(I), Pre_nearBFNorm(I) , 1);
                    m = fitvars(1);
                    c = fitvars(2);                   
                    
                    Xnan2 =vertcat(0, LP);
                    Y_Fit = polyval(fitvars, Xnan2);
                    plot(Xnan2,Y_Fit,'Color','r'); 
                    
I = ~isnan(LP) &  ~isnan(L_nearBFNorm);                   
scatter(LP(I), L_nearBFNorm(I),'MarkerEdgeColor',[0 1 0])

[rho,pval] = corr(LP(I), L_nearBFNorm(I),'Type','Pearson');                    
                    fitvars = polyfit(LP(I), L_nearBFNorm(I) , 1);
                    m = fitvars(1);
                    c = fitvars(2);                   
                    
                    Xnan2 =vertcat(0, LP);
                    Y_Fit = polyval(fitvars, Xnan2);
                    plot(Xnan2,Y_Fit,'Color','g');   
                    
I = ~isnan(LP) &  ~isnan(Post_nearBFNorm);                    
scatter(LP(I), Post_nearBFNorm(I),'MarkerEdgeColor',[0 0 1])

[rho,pval] = corr(LP(I), Post_nearBFNorm(I),'Type','Pearson');                    
                    fitvars = polyfit(LP(I), Post_nearBFNorm(I) , 1);
                    m = fitvars(1);
                    c = fitvars(2);                   
                    
                    Xnan2 =vertcat(0, LP);
                    Y_Fit = polyval(fitvars, Xnan2);
                    plot(Xnan2,Y_Fit,'Color','b');  
                    
legend('Pre','','Opto','','Post','','Location','best')
 title('zscored nearBF')                       

subplot(2,3,6)
hold on
I = ~isnan(LP) &  ~isnan(Pre_nonBFNorm);
scatter(LP(I), Pre_nonBFNorm(I),'MarkerEdgeColor',[1 0 0])

[rho,pval] = corr(LP(I), Pre_nonBFNorm(I),'Type','Pearson');                    
                    fitvars = polyfit(LP(I), Pre_nonBFNorm(I) , 1);
                    m = fitvars(1);
                    c = fitvars(2);                   
                    
                    Xnan2 =vertcat(0, LP);
                    Y_Fit = polyval(fitvars, Xnan2);
                    plot(Xnan2,Y_Fit,'Color','r'); 
                    
I = ~isnan(LP) &  ~isnan(L_nonBFNorm);                   
scatter(LP(I), L_nonBFNorm(I),'MarkerEdgeColor',[0 1 0])

[rho,pval] = corr(LP(I), L_nonBFNorm(I),'Type','Pearson');                    
                    fitvars = polyfit(LP(I), L_nonBFNorm(I) , 1);
                    m = fitvars(1);
                    c = fitvars(2);                   
                    
                    Xnan2 =vertcat(0, LP);
                    Y_Fit = polyval(fitvars, Xnan2);
                    plot(Xnan2,Y_Fit,'Color','g');   
                    
I = ~isnan(LP) &  ~isnan(Post_nonBFNorm);                    
scatter(LP(I), Post_nonBFNorm(I),'MarkerEdgeColor',[0 0 1])

[rho,pval] = corr(LP(I), Post_nonBFNorm(I),'Type','Pearson');                    
                    fitvars = polyfit(LP(I), Post_nonBFNorm(I) , 1);
                    m = fitvars(1);
                    c = fitvars(2);                   
                    
                    Xnan2 =vertcat(0, LP);
                    Y_Fit = polyval(fitvars, Xnan2);
                    plot(Xnan2,Y_Fit,'Color','b');  
                    
legend('Pre','','Opto','','Post','','Location','best')
 title('zscored nonBF')     