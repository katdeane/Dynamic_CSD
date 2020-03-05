function GMOScal_LayerSelect(layer, rel2BFin)
% GMOScal gets the averaged scalogram for each condition at layer 4
% Uses data produced by CSD_allLayers_scalogram from Andrew Curran, named
% here scalogramsfull.mat
% Outputs 3 figures with subplots for 3 slightly different visualizations:
%   - scaled to full scale (minimum found value, maximum found value)
% 03.06.2019 added function to specify layer under consideration

%% INIT 
cd('D:\MyCode\Dynamic_CSD_Analysis');
home = pwd;
addpath(genpath(home));

cd AndrewSpectralData; cd Data;

load('scalogramsfull.mat')% var = wtTable

varNames = unique(wtTable.layer);
if ~any(contains(varNames,layer)) && ~strcmp(layer, 'ALL')
    error('Please enter a valid layer, KATRINA')
end
params.startTime = -0.2; % seconds
params.limit = 600;

%check if layer or full mat needed
if ~strcmp(layer, 'ALL')
    %Pull out layer 
    wt2 = wtTable(contains(wtTable.layer,layer),:);
    wtTable = wt2;
else
    % if condition is all, it needs to still pull out only early sinks
    % because there is no time distinction here (i.e. = VIE == VIL)
    wt2 = wtTable(contains(wtTable.layer,'E'),:);
    wtTable = wt2;
end

%Pull out conditions and limit them to 600ms
awake = table2cell(wtTable(contains(wtTable.condition,'Awake')&wtTable.rel2Bf==rel2BFin,1));
awake =  cellfun(@(x) x(:,1:params.limit),awake,'UniformOutput',false);

anest = table2cell(wtTable(contains(wtTable.condition,'Anesth')&wtTable.rel2Bf==rel2BFin,1));
anest =  cellfun(@(x) x(:,1:params.limit),anest,'UniformOutput',false);

musc = table2cell(wtTable(contains(wtTable.condition,'Muscimol')&wtTable.rel2Bf==rel2BFin,1));
musc =  cellfun(@(x) x(:,1:params.limit),musc,'UniformOutput',false);

%Stack the individual animals (animal#x54x600)
for ii = 1:length(awake)
    awake2(ii,:,:)=awake{ii};
end

for ii = 1:length(anest)
    anest2(ii,:,:)=anest{ii};
end

for ii = 1:length(musc)
    musc2(ii,:,:)=musc{ii};
end

cd(home); cd AndrewSpectralData; cd Pictures;
%% Scaled for total (Minimum to Maximum caxis values between all three conditions)
%close all
[X,Y]=meshgrid(wtTable.freq{1},params.startTime*1000:(params.limit-201));
awake2 = abs(awake2);
%awakeFig = figure;
completeScaleFig = figure('Name',['Complete Scale ' layer ' ' num2str(rel2BFin)],'Position',[-1070 500 1065 400]); 
awakeFig = subplot(131);
surf(Y,X,squeeze(mean(awake2,1))','EdgeColor','None'); view(2);
set(gca,'YScale','log'); title('Awake')
yticks([0 10 20 30 40 50 60 80 100 200 300 500])
colorbar
clim = get(gca,'clim');

anest2 = abs(anest2);
%ketFig = figure; 
ketFig = subplot(132);
surf(Y,X,squeeze(mean(anest2,1))','EdgeColor','None'); view(2);
set(gca,'YScale','log'); title('Ketamine')
yticks([0 10 20 30 40 50 60 80 100 200 300 500])
colorbar
clim = [clim; get(gca,'clim')];

musc2 = abs(musc2);
%muscFig = figure; 
muscFig = subplot(133);
surf(Y,X,squeeze(mean(musc2,1))','EdgeColor','None'); view(2);
set(gca,'YScale','log'); title('Muscimol')
yticks([0 10 20 30 40 50 60 80 100 200 300 500])
colorbar
clim = [clim; get(gca,'clim')];

newC = [min(clim(:)) max(clim(:))];

%figure(awakeFig); 
set(awakeFig,'Clim',newC);colorbar;
%figure(ketFig); 
set(ketFig,'Clim',newC);colorbar;
%figure(muscFig); 
set(muscFig,'Clim',newC);colorbar;

h = gcf;
h.Renderer = 'Painters';
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(['Complete Scale ' layer ' ' num2str(rel2BFin)])
saveas(gcf, ['Complete Scale ' layer ' ' num2str(rel2BFin) '.pdf'])
%% Scaled for awake condition (considers also Musc condition)

%close all
[X,Y]=meshgrid(wtTable.freq{1},params.startTime*1000:(params.limit-201));
awake2 = abs(awake2);
anestScaleFig = figure('Name', ['Awake Scale ' layer ' ' num2str(rel2BFin)],'Position',[-1070 500 1065 400]); 
awakeFig = subplot(131);
surf(Y,X,squeeze(mean(awake2,1))','EdgeColor','None'); view(2);
set(gca,'YScale','log'); title('Awake')
colorbar
clim = get(gca,'clim');

anest2 = abs(anest2);
%ketFig = figure; 
ketFig = subplot(132);
surf(Y,X,squeeze(mean(anest2,1))','EdgeColor','None'); view(2);
set(gca,'YScale','log'); title('Ketamine')
colorbar
%clim = [clim; get(gca,'clim')];

musc2 = abs(musc2);
%muscFig = figure; 
muscFig = subplot(133);
surf(Y,X,squeeze(mean(musc2,1))','EdgeColor','None'); view(2);
set(gca,'YScale','log'); title('Muscimol')
colorbar
clim = [clim; get(gca,'clim')];

newC = [min(clim(:)) max(clim(:))];

%figure(awakeFig); 
set(awakeFig,'Clim',newC);colorbar;
%figure(ketFig); 
set(ketFig,'Clim',newC);colorbar;
%figure(muscFig); 
set(muscFig,'Clim',newC);colorbar;

h = gcf;
h.Renderer = 'Painters';
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(['Awake Scale ' layer ' ' num2str(rel2BFin)])
saveas(gcf, ['Complete Scale ' layer ' ' num2str(rel2BFin) '.pdf'])
%% Log scale
%close all
[X,Y]=meshgrid(wtTable.freq{1},params.startTime*1000:(params.limit-201));
logScaleFig = figure('Name', ['Log Scale ' layer ' ' num2str(rel2BFin)],'Position',[-1070 500 1065 400]); 
awakeFig = subplot(131);
surf(Y,X,log10(squeeze(mean(awake2,1))'),'EdgeColor','None'); view(2); %LOG10 BEING USED
set(gca,'YScale','log'); title('Awake')
colorbar
clim = get(gca,'clim');

anest2 = abs(anest2);
ketFig = subplot(132); 
surf(Y,X,log10(squeeze(mean(anest2,1))'),'EdgeColor','None'); view(2);
set(gca,'YScale','log'); title('Ketamine')
colorbar
clim = [clim; get(gca,'clim')];

musc2 = abs(musc2);
muscFig = subplot(133); 
surf(Y,X,log10(squeeze(mean(musc2,1))'),'EdgeColor','None'); view(2);
set(gca,'YScale','log'); title('Muscimol')
colorbar
clim = [clim; get(gca,'clim')];

newC = [min(clim(:)) max(clim(:))];

%figure(awakeFig); 
set(awakeFig,'Clim',newC);colorbar;
%figure(ketFig); 
set(ketFig,'Clim',newC);colorbar;
%figure(muscFig); 
set(muscFig,'Clim',newC);colorbar;

h = gcf;
h.Renderer = 'Painters';
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(['Log Scale ' layer ' ' num2str(rel2BFin)])
saveas(gcf, ['Complete Scale ' layer ' ' num2str(rel2BFin) '.pdf'])
%% Save
% saveAllFigs('D:\MyCode\Dynamic_CSD_Analysis\AndrewSpectralData\Pictures');
