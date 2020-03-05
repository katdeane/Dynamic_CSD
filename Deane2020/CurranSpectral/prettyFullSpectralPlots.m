%% INIT

folder = 'D:\MyCode\Dynamic_CSD_Analysis\AndrewSpectralData\Data\';

dataFile = 'evok_ind_fullWavSpect.mat'; 

load([folder dataFile]) % titled 'allDat'

bfDat = allDat(allDat.rel2Bf==0,:);
pos = [421 44.2 904 738.4];

%% PLOTS Ind Wavs
pos = [421 44.2 904 738.4];
bfDat = allDat((allDat.rel2Bf == 0),:);
varNames = bfDat.Properties.VariableNames';
lims = [5 10 3 3 2]*1e-04;
for ii = 1:5
    figure('Name',varNames{ii+5})
    set(gcf,'Position',pos);
    for jj = 1:7
        subDat = bfDat(bfDat.rel2SinkCenter == jj-4,:);
        [G,TID] = findgroups(subDat.condition);
        subplot(3,3,jj)
        hold on;
        splitapply(@(x) plot(bfDat.freqs_indDatWavs(1,:),mean(x,1,'omitnan')),subDat.(varNames{ii+5}),G)
        %xlim([0 200])
        ylim([0,lims(ii)])
        set(gca,'XScale','log')
        title(sprintf('rel2SC: %d',jj-5));
        hold off
    end
    legend(TID)
end

%% PLOTS Evok waves
pos = [421 44.2 904 738.4];
bfDat = allDat((allDat.rel2Bf == 0),:);
varNames = bfDat.Properties.VariableNames';
lims = [5 10 3 3 2]*1e-04;
for ii = 1:5
    figure('Name',varNames{ii+5})
    set(gcf,'Position',pos);
    for jj = 1:9
        subDat = bfDat(bfDat.rel2SinkCenter == jj-5,:);
        [G,TID] = findgroups(subDat.condition);
        subplot(3,3,jj)
        hold on;
        splitapply(@(x) plot(bfDat.freqs_indDatWavs(1,:),mean(x,1,'omitnan')),subDat.(varNames{ii+5}),G)
        %xlim([0 200])
        ylim([0,lims(ii)])
        set(gca,'XScale','log')
        title(sprintf('rel2SC: %d',jj-5));
        hold off
    end
    legend(TID)
end

%% PLOT rel2BF = 0, rel2SC = 0, ONSINK for EV. IND. and DIFF

bfDat = allDat((allDat.rel2Bf == 0),:);
varNames = bfDat.Properties.VariableNames';
figure('Name',varNames{7})
set(gcf,'Position',pos);
for jj = 1:5
    subDat = bfDat(bfDat.rel2SinkCenter == jj-3,:);
    [G,TID] = findgroups(subDat.condition);
    subplot(3,2,jj)
    hold on;
    splitapply(@(x) plot(bfDat.freqs_evokDatWavs(1,:),mean(x,1,'omitnan')),subDat.(varNames{7}),G)
    set(gca,'XScale','log')
    title(sprintf('rel2SC: %d',jj-3));
    hold off
end
legend(TID)

figure('Name',varNames{7})
set(gcf,'Position',pos);
for jj = 1:5
    subDat = bfDat(bfDat.rel2SinkCenter == jj-3,:);
    [G,TID] = findgroups(subDat.condition);
    subplot(3,2,jj)
    hold on;
    x = subDat.(varNames{7})(G==1,:);
    plotshaded(bfDat.freqs_evokDatWavs(1,:),[mean(x,1,'omitnan')-std(x,1,'omitnan'); mean(x,1,'omitnan'); mean(x,1,'omitnan')+std(x,1,'omitnan')],'b')
    x = subDat.(varNames{7})(G==2,:);
    plotshaded(bfDat.freqs_evokDatWavs(1,:),[mean(x,1,'omitnan')-std(x,1,'omitnan'); mean(x,1,'omitnan'); mean(x,1,'omitnan')+std(x,1,'omitnan')],'r')
    x = subDat.(varNames{7})(G==3,:);
    plotshaded(bfDat.freqs_evokDatWavs(1,:),[mean(x,1,'omitnan')-std(x,1,'omitnan'); mean(x,1,'omitnan'); mean(x,1,'omitnan')+std(x,1,'omitnan')],'g')
    set(gca,'XScale','log')
    title(sprintf('rel2SC: %d',jj-3));
    hold off
end
leg = repelem(TID,2);
legend(leg)

%% GRAMM PLOT
figure('Name',varNames{7})
set(gcf,'Position',pos);
g = gramm(  'x',bfDat.freqs_evokDatWavs(1,:),...
            'y',bfDat.(varNames{7}),...
            'color',bfDat.condition,...
            'column',bfDat.rel2SinkCenter,...
            'subset',any((bfDat.rel2SinkCenter == [-2:2])'));
g.facet_wrap(bfDat.rel2SinkCenter,'ncols',2,'scale','independent')
g.stat_summary( 'type','quartile',...
                'geom','area')
g.set_names('x','Freq Hz','y','PSD','column','rel2SinkCenter')
g.axe_property('xlim',[5 160],'XScale','log')
g.draw();


% figure('Name',varNames{12})
% set(gcf,'Position',pos);
% for jj = 1:9
%     subDat = bfDat(bfDat.rel2SinkCenter == jj-5,:);
%     [G,TID] = findgroups(subDat.condition);
%     subplot(3,3,jj)
%     hold on;
%     splitapply(@(x) plot(bfDat.freqs_evokDatWavs(1,:),mean(x,1,'omitnan')),subDat.(varNames{12}),G)
%     set(gca,'XScale','log')
%     title(sprintf('rel2SC: %d',jj-5));
%     hold off
% end
% legend(TID)
% 
% figure('Name',[varNames{12} ' minus ' varNames{7}])
% set(gcf,'Position',pos);
% for jj = 1:9
%     subDat = bfDat(bfDat.rel2SinkCenter == jj-5,:);
%     [G,TID] = findgroups(subDat.condition);
%     subplot(3,3,jj)
%     hold on;
%     splitapply(@(x) plot(bfDat.freqs_evokDatWavs(1,:),mean(x,1,'omitnan')),subDat.(varNames{12})-subDat.(varNames{7}),G)
%     set(gca,'XScale','log')
%     title(sprintf('rel2SC: %d',jj-5));
%     hold off
% end
% legend(TID)





