function [g] = grammStarSignAdd(g,group1,group2,sig,offset,bracketFlag)
%grammStarSignAdd easily add significance stars to gramm plot
%   This is experimental code for adding stars into gramm plots (author: Andrew Curran)

% g being a gramm boxplot object
% group1 and group2 indicating the two box plots you want the star between. If you have only x vs y, then one number for each group is enough, if there are groups you gotta give a 2N vector of [x group]
% sig is just your sig value. 0.05 or whatever, itll do NS if its too high
% offset is axes ordinate based offset from the tip of the highest whisker between the two groups
% offset now dodges from absolute 
% bracketFlag allows the choice of having small dangling vert lines from
% the edges of the main sig line
% For now, assume boxplot

% Ok, adding in multipass functionality. Needs three rounds of dev:
% 1. specify heights of each sigstar separately
% 2. heights automatically bump upwards based on first offset. 
% 3. heights only bump if no interference (would be difficult)...

%Known issues: boxplot only. Interference with intermediate boxplots very
%easy to accidentally cause. 

nSubGroups = size(g.results.stat_boxplot,1);
%Convoluted solution to prevent trickery from a couple empty groups
nMainGroups = max(cellfun(@numel,{g.results.stat_boxplot.upper_whisker_handle}));

for i = 1:length(sig) % Repeat for several sig levels
    thisGroup1 = group1(i,:);
    thisGroup2 = group2(i,:);
    thisSig = sig(i);
    thisOffset = offset(i);
    % check input against number of x points
    if length(thisGroup1) == 2 %Indicative of double positioning (x and group)
        group1x = thisGroup1(1);
        group1g = thisGroup1(2);
        group2x = thisGroup2(1);
        group2g = thisGroup2(2);
    elseif length(thisGroup1) == 1 %one position, assume no grouping
        group1x = thisGroup1(1);
        group1g = 1;
        group2x = thisGroup2(1);
        group2g = 1;
    else
        error('Improper number of elements in grouping variable!')
    end
    
    plottedX1 = g.results.stat_boxplot(group1g).upper_whisker_handle(group1x).XData(end);
    plottedX2 = g.results.stat_boxplot(group2g).upper_whisker_handle(group2x).XData(end);
    
    plottedY1 = max(g.results.stat_boxplot(group1g).upper_whisker_handle(group1x).YData(end));
    plottedY2 = max(g.results.stat_boxplot(group2g).upper_whisker_handle(group2x).YData(end));
    
    plottedY = max([plottedY1 plottedY2]);
    
    % Figure out star string
    if thisSig > 0.05
        t = 'NS';
        tOffset = 100;
    elseif thisSig <= 0.05 && thisSig > 0.01
        t = '*';
        tOffset = 50;
    elseif thisSig <= 0.01 && thisSig > 0.001
        t = '**';
        tOffset = 50;
    elseif thisSig <= 0.001
        t = '***';
        tOffset = 50;
    else
        error('Did you try using a number for your significance value?')
    end
    
    % Draw the line
    L = line(g.facet_axes_handles,...
        [plottedX1 plottedX2],...
        [plottedY plottedY]+thisOffset,...
        'Color','k');
    
    if nargin > 5 && bracketFlag == 1
        line(g.facet_axes_handles,...
            [plottedX1 plottedX1],...
            [plottedY-50 plottedY]+thisOffset,...
            'Color','k');
        line(g.facet_axes_handles,...
            [plottedX2 plottedX2],...
            [plottedY-50 plottedY]+thisOffset,...
            'Color','k');
    end
    
    textHeight = plottedY+thisOffset+tOffset;
    T = text(g.facet_axes_handles,...
        mean([plottedX1 plottedX2]),...
        textHeight,...
        t,...
        'HorizontalAlignment','center');
    
    if textHeight > g.facet_axes_handles.YLim(2)
        g.facet_axes_handles.YLim(2) = textHeight+10;
    end
    
end

end

