function PlotAnimalAvrec(avrec, names)
%plots individual animal AVRECs and average group AVREC, called by
%snglAVREC_BF and avgAVREC
   %BF
    h = figure('Name','sAvgRec - BF');
    average = [];
    for i3 = 1:length(names)
        line = avrec.(names{i3}).BF;
        line = line(1:600); %GXL have 600ms and GKD have 800!
        plot(line);
        average = [average line'];
        hold on
    end
    
    average = nanmean(average,2)';
    plot(average,'k','LineWidth',2);
    allnames = vertcat(names,{'Average'});
    legend(allnames)
    
    set(h, 'PaperType', 'A4');
    set(h, 'PaperOrientation', 'landscape');
    set(h, 'PaperUnits', 'centimeters');
    savefig(h, 'sAvgRec - BF','compact')
    close (h)
    
    %BF -1
    h = figure('Name','sAvgRec - BF -1');
    average = [];
    for i3 = 1:length(names)
        line = avrec.(names{i3}).min_one;
        line = line(1:600); %GXL have 600ms and GKD have 800!
        plot(line);
        average = [average line'];
        hold on
    end
    
    average = nanmean(average,2)';
    plot(average,'k','LineWidth',2);
    legend(allnames)
    
    set(h, 'PaperType', 'A4');
    set(h, 'PaperOrientation', 'landscape');
    set(h, 'PaperUnits', 'centimeters');
    savefig(h, 'sAvgRec - BF -1','compact')
    close (h)
    
    %BF -2
    h = figure('Name','sAvgRec - BF -2');
    average = [];
    for i3 = 1:length(names)
        line = avrec.(names{i3}).min_two;
        line = line(1:600); %GXL have 600ms and GKD have 800!
        plot(line);
        average = [average line'];
        hold on
    end
    
    average = nanmean(average,2)';
    plot(average,'k','LineWidth',2);
    legend(allnames)
    
    set(h, 'PaperType', 'A4');
    set(h, 'PaperOrientation', 'landscape');
    set(h, 'PaperUnits', 'centimeters');
    savefig(h, 'sAvgRec - BF -2','compact')
    close (h)
    
     %BF -3
    h = figure('Name','sAvgRec - BF -3');
    average = [];
    for i3 = 1:length(names)
        line = avrec.(names{i3}).min_three;
        line = line(1:600); %GXL have 600ms and GKD have 800!
        plot(line);
        average = [average line'];
        hold on
    end
    
    average = nanmean(average,2)';
    plot(average,'k','LineWidth',2);
    legend(allnames)
    
    set(h, 'PaperType', 'A4');
    set(h, 'PaperOrientation', 'landscape');
    set(h, 'PaperUnits', 'centimeters');
    savefig(h, 'sAvgRec - BF -3','compact')
    close (h)
   
    %BF +1
    h = figure('Name','sAvgRec - BF +1');
    average = [];
    for i3 = 1:length(names)
        line = avrec.(names{i3}).plus_one;
        line = line(1:600); %GXL have 600ms and GKD have 800!
        plot(line);
        average = [average line'];
        hold on
    end
    
    average = nanmean(average,2)';
    plot(average,'k','LineWidth',2);
    legend(allnames)
    
    set(h, 'PaperType', 'A4');
    set(h, 'PaperOrientation', 'landscape');
    set(h, 'PaperUnits', 'centimeters');
    savefig(h, 'sAvgRec - BF +1','compact')
    close (h)
    
    %BF +2
    h = figure('Name','sAvgRec - BF +2');
    average = [];
    for i3 = 1:length(names)
        line = avrec.(names{i3}).plus_two;
        line = line(1:600); %GXL have 600ms and GKD have 800!
        plot(line);
        average = [average line'];
        hold on
    end
    
    average = nanmean(average,2)';
    plot(average,'k','LineWidth',2);
    legend(allnames)
    
    set(h, 'PaperType', 'A4');
    set(h, 'PaperOrientation', 'landscape');
    set(h, 'PaperUnits', 'centimeters');
    savefig(h, 'sAvgRec - BF +2','compact')
    close (h)
    
     %BF +3
    h = figure('Name','sAvgRec - BF +3');
    average = [];
    for i3 = 1:length(names)
        line = avrec.(names{i3}).plus_three;
        line = line(1:600); %GXL have 600ms and GKD have 800!
        plot(line);
        average = [average line'];
        hold on
    end
    
    average = nanmean(average,2)';
    plot(average,'k','LineWidth',2);
    legend(allnames)
    
    set(h, 'PaperType', 'A4');
    set(h, 'PaperOrientation', 'landscape');
    set(h, 'PaperUnits', 'centimeters');
    savefig(h, 'sAvgRec - BF +3','compact')
    close (h)
    close all;