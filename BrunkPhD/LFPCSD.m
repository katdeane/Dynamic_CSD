clear
dbstop if error
home = cd;
addpath(genpath(home));

channels =...
    [32; 11; 31; 12; 30; 13; 29; 14; 17;  1; 19; 10; 18;  2; 28;  9; 20;  3; 27;  8; 21;  4; 23;  7; 22; 15; 25;  6; 24;  5; 26; 16];




load('K:\CSD_dynamic_analysis\raw\GOT049_0014.mat')
          
BL =(Header.t_pre)*1000;          
[FP, ~, CSD, ~, ~, ~, ~, ~,~,~, ~,~] =...
    SingleTrialCSD_00(SWEEP,channels,BL); 



figure,...
imagesc(CSD{7}),colormap(jet),caxis([-.0005 .0005])
hold on


%%% LFP
LP =avgFP(:,:,7);
LP2 = ones(size(LP));
Chan = 1:32;

LP2 = LP2.*Chan'; % equidistant Chans
LP2 = (LP2 + LP*4);% 4 fold amplitude of LFP
% figure, 
plot(flipud(LP2)','Color','black')
set(gca,'Ydir','reverse')
ylim([0 33])
