function [sig_mass, pVal, permMean, permSTD] = Perm_makegramm(nperms,obsInputMASS, permInputMASS, obsInputMAX, permInputMAX,...
    tickorder, ticklen, Ticknames, plotname)

% This function exists to make gramm plots with the given permutation test
% data

% INPUT:
% obsInput and permInput of clustermass and max t test depend on what part of the CSD is being run
% tickorder, ticklen, Ticknames are set by main script header and must be
% changed there if other frequency bins are required
% plotname needs to be according to which layer is being run to save out
% the grammplots

% OUTPUT: 
% allmassmat go to main script to calculate final stats
% .png and .pdf are exported of all gramm plots in the PermuteCSDs figure
% folder

pthresh = nperms*0.0055; %for sig checks -- Bonferroni correction for 9 tests

allmass = horzcat(obsInputMASS, permInputMASS);
allmassmat = cell2mat(allmass);

allmax = horzcat(obsInputMAX, permInputMAX);
allmaxmat = cell2mat(allmax);
% tick names,
freq_bin = repmat(tickorder, [1,length(allmass)/ticklen]);
% and group names
permgroup = {'Permutation'}; obsgroup = {'Observed'};
permrepmat = repmat(permgroup, [1,(length(allmass)-ticklen)]);
obsrepmat = repmat(obsgroup, [1,ticklen]);
obsVSperm = horzcat(obsrepmat, permrepmat);

% Cluster Mass plot
clear g
g=gramm('x',freq_bin,'y',allmassmat,'color',obsVSperm);
g.stat_summary('type','std','geom','area'); %mean and std shown
g.set_layout_options('Position',[0 0 0.7 0.7],...
    'legend_pos',[0.71 0.66 0.2 0.2],... %We detach the legend from the plot and move it to the top right
    'margin_height',[0.1 0.1],...
    'margin_width',[0.1 0.1],...
    'redraw',false);
g.set_names('x','Frequency Bin','y','Cluster Mass','color','Group');
% g.axe_property('ylim',[0 5000])
g.set_color_options('map','brewer_dark');
g.axe_property('XTickLabel',Ticknames)
g.draw();
g.export('file_name',['Cluster ' plotname], 'file_type','png');
g.export('file_name',['Cluster ' plotname], 'file_type','pdf');

close all;

% Check Significance of full clustermass
acreshape = reshape(allmassmat(ticklen+1:end),ticklen,1000); %pull out all values per tick
comp = cell2mat(obsInputMASS)'; %turn observed into matrix
sig_mass = sum(acreshape>comp,2); %compare how many instances the sum of the permutation is greater than observed
pVal = sig_mass/nperms;
permMean = mean(acreshape,2);
permSTD = std(acreshape,0,2);

sig_mass(sig_mass <= pthresh) = true; %if it's not more than 50, it's siginificant for that frequency
sig_mass(sig_mass > pthresh) = false;


% BF histogram of permutation distribution
f = figure('Name',['BFhist ' plotname]);
hist(gca,acreshape(4,:),20)
g=gca;
hold on;
plot([comp(4),comp(4)],[0,g.YLim(2)],'r');
saveas(f,['BFhist ' plotname])
close all