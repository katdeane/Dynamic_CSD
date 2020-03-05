function prettySpectralPlots(rel2BFin)
%% Pretty Spectral Plots
% Input:    evokeAll.mat from Andrew Curran's wavelet analysis
% Output:   Scatter plots of target sinks for peak latency/amplitude and 
%           best waveform.

%% INIT
% 
% clear
cd('D:\MyCode\Dynamic_CSD_Analysis');
home = pwd;
addpath(genpath(home));

folder = 'D:\MyCode\Dynamic_CSD_Analysis\AndrewSpectralData\Data\';
savFold = 'D:\MyCode\Dynamic_CSD_Analysis\AndrewSpectralData\Pictures';
dataFile = 'evokAll.mat';

cd(home); cd Data;
load([folder dataFile], 'evokAll') % titled 'evokAll'

% rel2BFin = 0;
evokAll = evokAll(evokAll.rel2Bf==rel2BFin,:);
pos = [421 44.2 904 738.4];

cd(home); cd(savFold);
para = {'preStim' 'offSink' 'onSink'};
sink = {'I_IIE' 'IVE' 'VbE' 'VIaE'};


%% PLOT EVOKE DAT

for isink = 1:length(sink)
    
    %pull out layer sink of all other vectors %I_IIE IVE VbE VIaE
%     animalchunk = evokAll.animal(strcmp(sink{isink},evokAll.layer)); % if
%     %we want to visualize animals rather than groups, uncomment this line
    conditionchunk = evokAll.condition(strcmp(sink{isink},evokAll.layer));
    rel2SinkCenterchunk = evokAll.rel2SinkCenter(strcmp(sink{isink},evokAll.layer));
    
    
    
    for ipara = 1:length(para)
        
        %pull out parameter
        PT_all = vertcat(evokAll.peakTime(:,:).(para{ipara})); %preStim offSink onSink
        PS_all = vertcat(evokAll.peakStrength(:,:).(para{ipara}));
        BW_all = vertcat(evokAll.bestWavelet(:,:).(para{ipara}));
        %pull out sink
        PT_chunk = PT_all(strcmp(sink{isink},evokAll.layer));
        PS_chunk = PS_all(strcmp(sink{isink},evokAll.layer));
        BW_chunk = BW_all(strcmp(sink{isink},evokAll.layer));
        
        % IVE PreStim BW vs PT
        clear g
        
        g = gramm('x',PT_chunk,'y',BW_chunk,'color',conditionchunk, ...
            'subset',(rel2SinkCenterchunk>-2 & rel2SinkCenterchunk<2)); % this means the center 5 channels are in the plot!
        g.set_names('x','Peak Time','y','','row','')
        g.set_title([sink{isink} ' ' para{ipara} ' Best Wavelet vs Peak Time Evoked Data ' num2str(rel2BFin)])
        g.geom_point();
        g.set_color_options('map','matlab');
        g.axe_property('YScale','log')
        f=figure; set(gcf,'Position',pos);
        g.draw();
        g.export('file_name',[sink{isink} ' ' para{ipara} ' Best Wavelet vs Peak Time Evoked Data ' num2str(rel2BFin)], 'file_type','pdf');
        g.export('file_name',[sink{isink} ' ' para{ipara} ' Best Wavelet vs Peak Time Evoked Data ' num2str(rel2BFin)], 'file_type','png');
%         saveas(f,[savFold [sink{isink} '_BWvPTevok_' para{ipara} '.emf']])
        
        %BW vs PS
        clear g
        
        g = gramm('x',BW_chunk,'y',PS_chunk,'color',conditionchunk,...
            'subset',(rel2SinkCenterchunk>-2 & rel2SinkCenterchunk<2));
        g.set_names('x','Best Freq','y','','row','')
        g.geom_point();
        g.set_color_options('map','matlab');
        g.axe_property('XScale','log')
        g.set_title([sink{isink} ' ' para{ipara} ' Best Wavelet vs Peak Strength Evoked Data ' num2str(rel2BFin)])
        f=figure; set(gcf,'Position',pos);
        g.draw();
        g.export('file_name',[sink{isink} ' ' para{ipara} ' Best Wavelet vs Peak Strength Evoked Data ' num2str(rel2BFin)], 'file_type','pdf');
        g.export('file_name',[sink{isink} ' ' para{ipara} ' Best Wavelet vs Peak Strength Evoked Data ' num2str(rel2BFin)], 'file_type','png');
%         saveas(f,[savFold [sink{isink} '_BFvPSevok_' para{ipara} '.emf']])
        
        %PT vs PS
        g = gramm('x',PT_chunk,'y',PS_chunk,'color',conditionchunk,...
            'subset',(rel2SinkCenterchunk>-2 & rel2SinkCenterchunk<2));
        g.set_names('x','Peak Time','y','','row','')
        g.set_title([sink{isink} ' ' para{ipara} ' Peak Time vs Peak Strength Evoked Data ' num2str(rel2BFin)])
        g.geom_point();
        g.set_color_options('map','matlab');
        f=figure; set(gcf,'Position',pos);
        g.draw();
        g.export('file_name',[sink{isink} ' ' para{ipara} ' Peak Time vs Peak Strength Evoked Data ' num2str(rel2BFin)], 'file_type','pdf');
        g.export('file_name',[sink{isink} ' ' para{ipara} ' Peak Time vs Peak Strength Evoked Data ' num2str(rel2BFin)], 'file_type','png');
%         saveas(f,[savFold [sink{isink} '_PTvPSevok_' para{ipara} '.emf']])
        
        
        % Best Wavelet stats and boxplots
        BWket = BW_chunk(strcmp('AnesthetizedPre',conditionchunk));
        BWawa = BW_chunk(strcmp('Awake10dB',conditionchunk));
        BWmus = BW_chunk(strcmp('Muscimol',conditionchunk));
        
        Group = [ ones(size(BWket)); 2 * ones(size(BWawa)); 3 * ones(size(BWmus))];
        BWfig = figure;
        boxplot([BWket; BWawa; BWmus],Group,1)
        set(gca,'XTickLabel',{'Ket','Awake','Mus'})
        title([sink{isink} ' ' para{ipara} ' Best Wavelet Frequency ' num2str(rel2BFin)])
        
        [abwH,abwP,abwCI,abwStat] = ttest2(BWket, BWawa);
        AvK_BWCohensd = iMakeCohensD(BWket,BWawa);
        if rel2BFin < 2 && rel2BFin > -2
            [mbwH,mbwP,mbwCI,mbwStat] = ttest(BWket, BWmus);
            MvK_BWCohensd = iMakeCohensD(BWket,BWmus);
        end
        
        h = gcf;
        set(h, 'PaperType', 'A4');
        set(h, 'PaperOrientation', 'landscape');
        set(h, 'PaperUnits', 'centimeters');
        savefig([sink{isink} ' ' para{ipara} ' Best Wavelet Frequency ' num2str(rel2BFin)])
        saveas(gcf, [sink{isink} ' ' para{ipara} ' Best Wavelet Frequency ' num2str(rel2BFin) '.pdf'])
        
        % Peak Strength stats and boxplots
        PSket = PS_chunk(strcmp('AnesthetizedPre',conditionchunk));
        PSawa = PS_chunk(strcmp('Awake10dB',conditionchunk));
        PSmus = PS_chunk(strcmp('Muscimol',conditionchunk));
        
        PSfig = figure;
        boxplot([PSket; PSawa; PSmus],Group,1)
        set(gca,'XTickLabel',{'Ket','Awake','Mus'})
        title([sink{isink} ' ' para{ipara} ' Peak Strength at BW ' num2str(rel2BFin)])
        
        h = gcf;
        set(h, 'PaperType', 'A4');
        set(h, 'PaperOrientation', 'landscape');
        set(h, 'PaperUnits', 'centimeters');
        savefig([sink{isink} ' ' para{ipara} ' Peak Strength at BW ' num2str(rel2BFin)])
        saveas(gcf, [sink{isink} ' ' para{ipara} ' Peak Strength at BW ' num2str(rel2BFin) '.pdf'])
        
        [apsH,apsP,apsCI,apsStat] = ttest2(PSket, PSawa);
        AvK_PSCohensd = iMakeCohensD(PSket,PSawa);
        if rel2BFin < 2 && rel2BFin > -2
            [mpsH,mpsP,mpsCI,mpsStat] = ttest(PSket, PSmus);
            MvK_PSCohensd = iMakeCohensD(PSket,PSmus);
        end
      
        % Peak Time stats and boxplots
        PTket = PT_chunk(strcmp('AnesthetizedPre',conditionchunk));
        PTawa = PT_chunk(strcmp('Awake10dB',conditionchunk));
        PTmus = PT_chunk(strcmp('Muscimol',conditionchunk));
        
        PTfig = figure;
        boxplot([PTket; PTawa; PTmus],Group,1)
        set(gca,'XTickLabel',{'Ket','Awake','Mus'})
        title([sink{isink} ' ' para{ipara} ' Peak Time at BW ' num2str(rel2BFin) '.pdf'])
        
        h = gcf;
        set(h, 'PaperType', 'A4');
        set(h, 'PaperOrientation', 'landscape');
        set(h, 'PaperUnits', 'centimeters');
        savefig([sink{isink} ' ' para{ipara} ' Peak Time at BW ' num2str(rel2BFin)])
        saveas(gcf, [sink{isink} ' ' para{ipara} ' Peak Time at BW ' num2str(rel2BFin) '.pdf'])
        
        [aptH,aptP,aptCI,aptStat] = ttest2(PTket, PTawa);
        AvK_PTCohensd = iMakeCohensD(PTket,PTawa);
        if rel2BFin < 2 && rel2BFin > -2
            [mptH,mptP,mptCI,mptStat] = ttest(PTket, PTmus);
            MvK_PTCohensd = iMakeCohensD(PTket,PTmus);
        end
        
        if rel2BFin < 2 && rel2BFin > -2
            save([sink{isink} '_' para{ipara} '_BWstats_' num2str(rel2BFin) '.mat'],'abwH', 'abwP', 'abwCI',...
                'abwStat','AvK_BWCohensd','mbwH', 'mbwP', 'mbwCI', 'mbwStat','MvK_BWCohensd')
            save([sink{isink} '_' para{ipara} '_PSstats_' num2str(rel2BFin) '.mat'],'apsH', 'apsP', 'apsCI',...
                'apsStat','AvK_PSCohensd','mpsH', 'mpsP', 'mpsCI', 'mpsStat','MvK_PSCohensd')
            save([sink{isink} '_' para{ipara} '_PTstats_' num2str(rel2BFin) '.mat'],'aptH', 'aptP', 'aptCI',...
                'aptStat','AvK_PTCohensd','mptH', 'mptP', 'mptCI', 'mptStat','MvK_PTCohensd')
        else
            save([sink{isink} '_' para{ipara} '_BWstats_' num2str(rel2BFin) '.mat'],'abwH', 'abwP', 'abwCI',...
                'abwStat','AvK_BWCohensd')
            save([sink{isink} '_' para{ipara} '_PSstats_' num2str(rel2BFin) '.mat'],'apsH', 'apsP', 'apsCI',...
                'apsStat','AvK_PSCohensd')
            save([sink{isink} '_' para{ipara} '_PTstats_' num2str(rel2BFin) '.mat'],'aptH', 'aptP', 'aptCI',...
                'aptStat','AvK_PTCohensd')
        end
        
    end
    
    close all
    
end
%% PLOT INDUCED DAT
% 
% g = gramm('x',bfDat.peakTime_indDat,'y',bfDat.bestWavelet_indDat,'color',bfDat.animal,...
%      'column',bfDat.condition,'row',bfDat.rel2SinkCenter,'subset',(bfDat.rel2SinkCenter~=-4 & bfDat.rel2SinkCenter~=4));
% g.set_names('x','Peak Time','y','','row','')
% g.set_title('Best Freq vs Peak Time Induced Data')
% g.geom_point();
% f=figure; set(gcf,'Position',pos);
% g.draw();
% saveas(f,[savFold 'BFvPTevok.emf'])
% 
% g = gramm('x',bfDat.bestWavelet_indDat,'y',bfDat.peakStrength_indDat,'color',bfDat.animal,...
%      'column',bfDat.condition,'row',bfDat.rel2SinkCenter,'subset',(bfDat.rel2SinkCenter~=-4 & bfDat.rel2SinkCenter~=4));
% g.set_names('x','Best Freq','y','','row','')
% g.geom_point();
% g.set_title('Best Freq vs Peak Strength Induced Data')
% f=figure; set(gcf,'Position',pos);
% g.draw();
% 
% g = gramm('x',bfDat.peakTime_indDat,'y',bfDat.peakStrength_indDat,'color',bfDat.animal,...
%      'column',bfDat.condition,'row',bfDat.rel2SinkCenter,'subset',(bfDat.rel2SinkCenter~=-4 & bfDat.rel2SinkCenter~=4));
% g.set_names('x','Peak Time','y','','row','')
% g.set_title('Peak Time vs Peak Strength Induced Data')
% g.geom_point();
% f=figure; set(gcf,'Position',pos);
% g.draw();

%% TEST PLOT DIFFERENCE
% 
% g = gramm('x',bfDat.peakTime_indDat-bfDat.peakTime_evokDat,'y',bfDat.bestWavelet_indDat-bfDat.bestWavelet_evokDat,'color',bfDat.animal,...
%      'column',bfDat.condition,'row',bfDat.rel2SinkCenter,'subset',(bfDat.rel2SinkCenter~=-4 & bfDat.rel2SinkCenter~=4));
% g.set_names('x','Peak Time','y','','row','')
% g.set_title('Best Freq vs Peak Time Induced minus Evoked Data')
% g.geom_point();
% figure; set(gcf,'Position',pos);
% g.draw();
% 
% g = gramm('x',bfDat.bestWavelet_indDat-bfDat.bestWavelet_evokDat,'y',bfDat.peakStrength_indDat-bfDat.peakStrength_evokDat,'color',bfDat.animal,...
%      'column',bfDat.condition,'row',bfDat.rel2SinkCenter,'subset',(bfDat.rel2SinkCenter~=-4 & bfDat.rel2SinkCenter~=4));
% g.set_names('x','Best Freq','y','','row','')
% g.geom_point();
% g.set_title('Best Freq vs Peak Strength Induced minus Evoked Data')
% figure; set(gcf,'Position',pos);
% g.draw();