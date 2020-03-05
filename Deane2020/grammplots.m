%% Start
clear all
cd('D:\MyCode\Dynamic_CSD_Analysis');
warning('OFF');
dbstop if error

home = pwd;
addpath(genpath(home));

%for teg_repeated measures_ANOVA
levels = [2,7];
varinames = {'Groups','Frequencies'};
rowofnans = NaN(1,7);
srowofnans = [{NaN(1,50)} {NaN(1,50)} {NaN(1,50)} {NaN(1,50)} {NaN(1,50)} {NaN(1,50)} {NaN(1,50)}];
ticks = {'a-3' 'b-2' 'c-1' 'dBF' 'e+1' 'f+2' 'g+3'};

% Order = {'IVE','IVL','I_IIE','I_IIL', 'VaE','VaL','VbE','VbL','VIE','VIL'};
Order = {'IVE','IVL','I_IIE','I_IIL', 'VaE','VaL','VbE','VbL','VIaE','VIaL','VIbE','VIbL'};
Parameter = {'SinkRMS','SinkPeakAmp','SinkPeakLate','Sinkonset'};
SParameter = {'SingleSinkRMS','SingleSinkPeakAmp','SingleSinkPeakLat'};


%% Load in the appropriate files
cd DATA;cd output;
load('AnesthetizedPre_Data.m_Threshold_0.25_Zscore_0_binned_1_mirror_0.mat')
Anesthetized = Data; clear Data; 
load('Awake10dB_Data.m_Threshold_0.25_Zscore_0_binned_1_mirror_0.mat')
Awake = Data; clear Data;
load('Muscimol_Data.m_Threshold_0.25_Zscore_0_binned_1_mirror_0.mat')
Muscimol = Data; clear Data;

cd(home);cd figs;
% mkdir('Group Plots Gramm'); cd('Group Plots Gramm')
mkdir('Group Plots Gramm'); cd('Group Plots Gramm'); mkdir('AVG'); mkdir('SINGLE'); addpath(genpath(home));

%% Sink Loop
for isink = 1:length(Order)
    %% Avg Trials Loop
    cd(home); cd figs; cd('Group Plots Gramm'); cd('AVG');
    for ipara = 1:length(Parameter)
        %% 3 Groups
        An_data = vertcat(Anesthetized.GS_based.(Parameter{ipara}).(Order{isink})(:,5:11));
        An_datamean = nanmean(An_data,1);
        An_datasem = nanstd(An_data,1)/sqrt(sum(~isnan(An_datamean)));
        
        Aw_data = vertcat(Awake.GS_based.(Parameter{ipara}).(Order{isink})(:,4:10),rowofnans, rowofnans, rowofnans);
        Aw_datamean = nanmean(Aw_data,1);
        Aw_datasem = nanstd(Aw_data,1)/sqrt(sum(~isnan(Aw_datamean)));
        
        M_data = vertcat(Muscimol.GS_based.(Parameter{ipara}).(Order{isink})(:,5:11));
        M_datamean = nanmean(M_data,1);
        M_datasem = nanstd(M_data,1)/sqrt(sum(~isnan(M_datamean)));
        
        % Create Appropriate Structure
        plotdata = struct;
        
        tone = ticks';
        Angrammdata = An_data(1,:)';
        Mugrammdata = M_data(1,:)';
        Awgrammdata = Aw_data(1,:)';
        
        Angroup = {'Anesthetized' 'Anesthetized' 'Anesthetized' 'Anesthetized' 'Anesthetized' 'Anesthetized' 'Anesthetized'}; %7 to match the amount of ticks
        An = Angroup';
        Mugroup = {'Muscimol' 'Muscimol' 'Muscimol' 'Muscimol' 'Muscimol' 'Muscimol' 'Muscimol'};
        Mu = Mugroup';
        Awgroup = {'Awake' 'Awake' 'Awake' 'Awake' 'Awake' 'Awake' 'Awake'};
        Aw = Awgroup';
        
        for istack = 1:size(An_data,1)-1
            
            tone = vertcat(tone, ticks');
            Angrammdata = vertcat(Angrammdata, An_data(istack+1,:)');
            Mugrammdata = vertcat(Mugrammdata, M_data(istack+1,:)');
            Awgrammdata = vertcat(Awgrammdata, Aw_data(istack+1,:)');
            An = vertcat(An, Angroup');
            Mu = vertcat(Mu, Mugroup');
            Aw = vertcat(Aw, Awgroup');
            
        end
        
        plotdata.tone = vertcat(tone, tone, tone);
        plotdata.data = vertcat(Angrammdata, Mugrammdata, Awgrammdata);
        plotdata.group = vertcat(An, Mu, Aw);

        clear g
        figure('Position',[100 100 650 550]);
        
        %Tuning curve
        g(1,1)=gramm('x',plotdata.tone,'y',plotdata.data, 'color', plotdata.group); %,'y',cars.Acceleration,'color',cars.Cylinders,'subset',cars.Cylinders~=3 & cars.Cylinders~=5
        % g(1,1).geom_point(); %to see each point, removed the YLim (highest points are closer to 2.5
        g(1,1).stat_summary('type','sem','geom','errorbar'); %mean and sem shown
        g(1,1).stat_summary('type','sem','geom','point'); %mean and sem shown
        g(1,1).stat_summary('type','sem','geom','area'); %mean and sem shown
        g(1,1).set_layout_options('Position',[0 0 .8 1],...
            'legend_pos',[0.6 0.68 0.2 0.2],... %We detach the legend from the plot and move it to the top right
            'margin_height',[0.1 0.02],...
            'margin_width',[0.1 0.02],...
            'redraw',false);
        g(1,1).set_text_options('base_size',12)
        g(1,1).axe_property('Ygrid','on'); %,'YLim',[0 0.0015]
        g(1,1).set_names('x','Tone','y','mV','color','Group');
        g(1,1).set_color_options('map','matlab');
        g(1,1).axe_property('XTickLabel',{'-3', '-2', '-1', 'BF', '+1', '+2', '+3'})
          
        %side histogram
        g(2,1)=gramm('x',plotdata.data,'color',plotdata.group);
        g(2,1).set_layout_options('Position',[0.8 0 0.2 1],...
            'legend',false,...
            'margin_height',[0.1 0.02],...
            'margin_width',[0.02 0.05],...
            'redraw',false);
        g(2,1).set_text_options('base_size',12)
        g(2,1).set_names('x','');
        g(2,1).stat_bin('geom','overlaid_bar','fill','transparent'); %histogram
        g(2,1).coord_flip();
        g(2,1).axe_property('XTickLabel','');
        g(2,1).set_color_options('map','matlab');
        
        g.set_title([(Order{isink}) ' ' (Parameter{ipara})]);
        g.draw();
        g.export('file_name',[(Order{isink}) (Parameter{ipara})], 'file_type','pdf');
        g.export('file_name',[(Order{isink}) '_' (Parameter{ipara})], 'file_type','png');
        close all;
        
        %% Awake vs Anesthetized
        fullmean = nanmean(vertcat(An_data,Aw_data),1)';
        
        % Create Appropriate Structure
        plotdata = struct;
        
        tone = ticks';
        Angrammdata = An_data(1,:)';
        Awgrammdata = Aw_data(1,:)';
        
        Angroup = {'Anesthetized' 'Anesthetized' 'Anesthetized' 'Anesthetized' 'Anesthetized' 'Anesthetized' 'Anesthetized'}; %7 to match the amount of ticks
        An = Angroup';
        Awgroup = {'Awake' 'Awake' 'Awake' 'Awake' 'Awake' 'Awake' 'Awake'};
        Aw = Awgroup';
        
        for istack = 1:size(An_data,1)-1
            
            tone = vertcat(tone, ticks');
            Angrammdata = vertcat(Angrammdata, An_data(istack+1,:)');
            Awgrammdata = vertcat(Awgrammdata, Aw_data(istack+1,:)');
            An = vertcat(An, Angroup');
            Aw = vertcat(Aw, Awgroup');
            
        end
        
        plotdata.tone = vertcat(tone, tone);
        plotdata.data = vertcat(Angrammdata, Awgrammdata);
        plotdata.group = vertcat(An, Aw);
        
        
        clear g
        figure('Position',[100 100 650 650]);
        
        g(1,1)=gramm('x',ticks','y',fullmean);
        g(1,1).geom_line();
        g(1,1).set_layout_options('Position',[0 0.8 0.8 0.2],... %Set the position in the figure (as in standard 'Position' axe property)
            'legend',false,... % No need to display legend for side histograms
            'margin_height',[0.02 0.05],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
            'margin_width',[0.1 0.02],...
            'redraw',false); %We deactivate automatic redrawing/resizing so that the axes stay aligned according to the margin options
        g(1,1).set_text_options('base_size',12)
        g(1,1).set_names('x','','y','Mean');
        g(1,1).axe_property('XTickLabel',''); % We deactivate the ticks
        g(1,1).set_line_options('base_size',3)
        g(1,1).set_color_options('map',[0 0 0])
        
        
        g(2,1)=gramm('x',plotdata.tone,'y',plotdata.data, 'color', plotdata.group); %,'y',cars.Acceleration,'color',cars.Cylinders,'subset',cars.Cylinders~=3 & cars.Cylinders~=5
        % g(2,1).geom_point(); %to see each point, removed the YLim (highest points are closer to 2.5
        g(2,1).stat_summary('type','sem','geom','errorbar'); %mean and sem shown
        g(2,1).stat_summary('type','sem','geom','point'); %mean and sem shown
        g(2,1).stat_summary('type','sem','geom','area'); %mean and sem shown
        g(2,1).set_layout_options('Position',[0 0 0.8 0.8],...
            'legend_pos',[0.75 0.75 0.2 0.2],... %We detach the legend from the plot and move it to the top right
            'margin_height',[0.1 0.02],...
            'margin_width',[0.1 0.02],...
            'redraw',false);
        g(2,1).set_text_options('base_size',12)
        g(2,1).axe_property('Ygrid','on'); %,'YLim',[0 0.0015]
        g(2,1).set_names('x','Tone','y','mV','color','Group');
        g(2,1).set_color_options('map','matlab');
        g(2,1).axe_property('XTickLabel',{'-3', '-2', '-1', 'BF', '+1', '+2', '+3'})
        
        
        g(3,1)=gramm('x',plotdata.data,'color',plotdata.group);
        g(3,1).set_layout_options('Position',[0.8 0 0.2 0.8],...
            'legend',false,...
            'margin_height',[0.1 0.02],...
            'margin_width',[0.02 0.05],...
            'redraw',false);
        g(3,1).set_text_options('base_size',12)
        g(3,1).set_names('x','');
        g(3,1).stat_bin('geom','overlaid_bar','fill','transparent'); %histogram
        g(3,1).coord_flip();
        g(3,1).axe_property('XTickLabel','');
        g(3,1).set_color_options('map','matlab');
        
        
        g.set_title([(Order{isink}) ' ' (Parameter{ipara})]);
        g.draw();
        g.export('file_name',['Awake_' (Order{isink}) '_' (Parameter{ipara})], 'file_type','png');
        g.export('file_name',['Awake_' (Order{isink}) (Parameter{ipara})], 'file_type','pdf');
        close all;
        
        %% Muscimol vs Anesthetized
        fullmean = nanmean(vertcat(An_data,M_data),1)';
        
        % Create Appropriate Structure
        plotdata = struct;
        
        tone = ticks';
        Angrammdata = An_data(1,:)';
        Mugrammdata = M_data(1,:)';
        
        Angroup = {'Anesthetized' 'Anesthetized' 'Anesthetized' 'Anesthetized' 'Anesthetized' 'Anesthetized' 'Anesthetized'}; %7 to match the amount of ticks
        An = Angroup';
        Mugroup = {'Muscimol' 'Muscimol' 'Muscimol' 'Muscimol' 'Muscimol' 'Muscimol' 'Muscimol'};
        Mu = Mugroup';
        
        for istack = 1:size(An_data,1)-1;
            
            tone = vertcat(tone, ticks');
            Angrammdata = vertcat(Angrammdata, An_data(istack+1,:)');
            Mugrammdata = vertcat(Mugrammdata, M_data(istack+1,:)');
            An = vertcat(An, Angroup');
            Mu = vertcat(Mu, Mugroup');
            
        end
        
        plotdata.tone = vertcat(tone, tone);
        plotdata.data = vertcat(Angrammdata, Mugrammdata);
        plotdata.group = vertcat(An, Mu);
        
        
        clear g
        figure('Position',[100 100 650 650]);
        
        g(1,1)=gramm('x',ticks','y',fullmean);
        g(1,1).geom_line();
        g(1,1).set_layout_options('Position',[0 0.8 0.8 0.2],... %Set the position in the figure (as in standard 'Position' axe property)
            'legend',false,... % No need to display legend for side histograms
            'margin_height',[0.02 0.05],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
            'margin_width',[0.1 0.02],...
            'redraw',false); %We deactivate automatic redrawing/resizing so that the axes stay aligned according to the margin options
        g(1,1).set_text_options('base_size',12)
        g(1,1).set_names('x','','y','Mean');
        g(1,1).axe_property('XTickLabel',''); % We deactivate the ticks
        g(1,1).set_line_options('base_size',3)
        g(1,1).set_color_options('map',[0 0 0])
        
        
        g(2,1)=gramm('x',plotdata.tone,'y',plotdata.data, 'color', plotdata.group); %,'y',cars.Acceleration,'color',cars.Cylinders,'subset',cars.Cylinders~=3 & cars.Cylinders~=5
        % g(2,1).geom_point(); %to see each point, removed the YLim (highest points are closer to 2.5
        g(2,1).stat_summary('type','sem','geom','errorbar'); %mean and sem shown
        g(2,1).stat_summary('type','sem','geom','point'); %mean and sem shown
        g(2,1).stat_summary('type','sem','geom','area'); %mean and sem shown
        g(2,1).set_layout_options('Position',[0 0 0.8 0.8],...
            'legend_pos',[0.75 0.75 0.2 0.2],... %We detach the legend from the plot and move it to the top right
            'margin_height',[0.1 0.02],...
            'margin_width',[0.1 0.02],...
            'redraw',false);
        g(2,1).set_text_options('base_size',12)
        g(2,1).axe_property('Ygrid','on'); %,'YLim',[0 0.0015]
        g(2,1).set_names('x','Tone','y','mV','color','Group');
        g(2,1).set_color_options('map','matlab');
        g(2,1).axe_property('XTickLabel',{'-3', '-2', '-1', 'BF', '+1', '+2', '+3'})
        
        
        g(3,1)=gramm('x',plotdata.data,'color',plotdata.group);
        g(3,1).set_layout_options('Position',[0.8 0 0.2 0.8],...
            'legend',false,...
            'margin_height',[0.1 0.02],...
            'margin_width',[0.02 0.05],...
            'redraw',false);
        g(3,1).set_text_options('base_size',12)
        g(3,1).set_names('x','');
        g(3,1).stat_bin('geom','overlaid_bar','fill','transparent'); %histogram
        g(3,1).coord_flip();
        g(3,1).axe_property('XTickLabel','');
        g(3,1).set_color_options('map','matlab');
        
        
        g.set_title([(Order{isink}) ' ' (Parameter{ipara})]);
        g.draw();
        g.export('file_name',['Musc_' (Order{isink}) '_' (Parameter{ipara})], 'file_type','png');
        g.export('file_name',['Musc_' (Order{isink}) (Parameter{ipara})], 'file_type','pdf');
        close all;
    end
    
    cd(home); cd figs; cd('Group Plots Gramm'); cd('SINGLE');
    %% Single Trials Loop
    for ispara = 1:length(SParameter)
        %% 3 Groups
        An_data = vertcat(Anesthetized.singleGS_based.(SParameter{ispara}).(Order{isink})(:,5:11)); %, rowofnans adds extra row if needed to make dimensions match
        shapeAn = cellfun(@(x) x', An_data,'UniformOutput', false); %so 50 rows per stimulus
        expandAn = cell2mat(shapeAn); %so all 50 rows are laid out in matrix
        
        Aw_data = vertcat(Awake.singleGS_based.(SParameter{ispara}).(Order{isink})(:,4:10),srowofnans,srowofnans,srowofnans);
        shapeAw = cellfun(@(x) x', Aw_data,'UniformOutput', false);
        expandAw = cell2mat(shapeAw);
        
        M_data = vertcat(Muscimol.singleGS_based.(SParameter{ispara}).(Order{isink})(:,5:11));
        shapeM = cellfun(@(x) x', M_data,'UniformOutput', false);
        expandM = cell2mat(shapeM);
        
        % Create Appropriate Structure
        plotdata = struct;
        
        tone = ticks';
        Angrammdata = expandAn(1,:)';
        Mugrammdata = expandM(1,:)';
        Awgrammdata = expandAw(1,:)';
        
        Angroup = {'Anesthetized' 'Anesthetized' 'Anesthetized' 'Anesthetized' 'Anesthetized' 'Anesthetized' 'Anesthetized'}; %7 to match the amount of ticks
        An = Angroup';
        Mugroup = {'Muscimol' 'Muscimol' 'Muscimol' 'Muscimol' 'Muscimol' 'Muscimol' 'Muscimol'};
        Mu = Mugroup';
        Awgroup = {'Awake' 'Awake' 'Awake' 'Awake' 'Awake' 'Awake' 'Awake'};
        Aw = Awgroup';
        
        for istack = 1:size(expandAn,1)-1
            
            tone = vertcat(tone, ticks');
            Angrammdata = vertcat(Angrammdata, expandAn(istack+1,:)');
            Mugrammdata = vertcat(Mugrammdata, expandM(istack+1,:)');
            Awgrammdata = vertcat(Awgrammdata, expandAw(istack+1,:)');
            An = vertcat(An, Angroup');
            Mu = vertcat(Mu, Mugroup');
            Aw = vertcat(Aw, Awgroup');
            
        end
        
        plotdata.tone = vertcat(tone, tone, tone);
        plotdata.data = vertcat(Angrammdata, Mugrammdata, Awgrammdata);
        plotdata.group = vertcat(An, Mu, Aw);

        clear g
        figure('Position',[100 100 650 550]);
        
        %Tuning curve
        g(1,1)=gramm('x',plotdata.tone,'y',plotdata.data, 'color', plotdata.group); %,'y',cars.Acceleration,'color',cars.Cylinders,'subset',cars.Cylinders~=3 & cars.Cylinders~=5
        % g(1,1).geom_point(); %to see each point, removed the YLim (highest points are closer to 2.5
        g(1,1).stat_summary('type','sem','geom','errorbar'); %mean and sem shown
        g(1,1).stat_summary('type','sem','geom','point'); %mean and sem shown
        g(1,1).stat_summary('type','sem','geom','area'); %mean and sem shown
        g(1,1).set_layout_options('Position',[0 0 .8 1],...
            'legend_pos',[0.6 0.68 0.2 0.2],... %We detach the legend from the plot and move it to the top right
            'margin_height',[0.1 0.02],...
            'margin_width',[0.1 0.02],...
            'redraw',false);
        g(1,1).set_text_options('base_size',12)
        g(1,1).axe_property('Ygrid','on'); %'YLim',[0 0.0025]
        g(1,1).set_names('x','Tone','y','mV','color','Group');
        g(1,1).set_color_options('map','matlab');
        g(1,1).axe_property('XTickLabel',{'-3', '-2', '-1', 'BF', '+1', '+2', '+3'})
          
        %side histogram
        g(2,1)=gramm('x',plotdata.data,'color',plotdata.group);
        g(2,1).set_layout_options('Position',[0.8 0 0.2 1],...
            'legend',false,...
            'margin_height',[0.1 0.02],...
            'margin_width',[0.02 0.05],...
            'redraw',false);
        g(2,1).set_text_options('base_size',12)
        g(2,1).set_names('x','');
        g(2,1).stat_bin('geom','overlaid_bar','fill','transparent'); %histogram
        g(2,1).coord_flip();
        g(2,1).axe_property('XTickLabel','');
        g(2,1).set_color_options('map','matlab');
        
        g.set_title([(Order{isink}) ' ' (SParameter{ispara})]);
        g.draw();
        g.export('file_name',[(Order{isink}) (SParameter{ispara})], 'file_type','pdf');
        g.export('file_name',[(Order{isink}) '_' (SParameter{ispara})], 'file_type','png');
        close all;
        
        %% Awake vs Anesthetized
        fullmean = nanmean(vertcat(expandAn,expandAw),1)';
        
        % Create Appropriate Structure
        plotdata = struct;
        
        tone = ticks';
        Angrammdata = expandAn(1,:)';
        Awgrammdata = expandAw(1,:)';
        
        Angroup = {'Anesthetized' 'Anesthetized' 'Anesthetized' 'Anesthetized' 'Anesthetized' 'Anesthetized' 'Anesthetized'}; %7 to match the amount of ticks
        An = Angroup';
        Awgroup = {'Awake' 'Awake' 'Awake' 'Awake' 'Awake' 'Awake' 'Awake'};
        Aw = Awgroup';
        
        for istack = 1:size(An_data,1)-1
            
            tone = vertcat(tone, ticks');
            Angrammdata = vertcat(Angrammdata, expandAn(istack+1,:)');
            Awgrammdata = vertcat(Awgrammdata, expandAw(istack+1,:)');
            An = vertcat(An, Angroup');
            Aw = vertcat(Aw, Awgroup');
            
        end
        
        plotdata.tone = vertcat(tone, tone);
        plotdata.data = vertcat(Angrammdata, Awgrammdata);
        plotdata.group = vertcat(An, Aw);
        
        
        clear g
        figure('Position',[100 100 650 650]);
        
        g(1,1)=gramm('x',ticks','y',fullmean);
        g(1,1).geom_line();
        g(1,1).set_layout_options('Position',[0 0.8 0.8 0.2],... %Set the position in the figure (as in standard 'Position' axe property)
            'legend',false,... % No need to display legend for side histograms
            'margin_height',[0.02 0.05],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
            'margin_width',[0.1 0.02],...
            'redraw',false); %We deactivate automatic redrawing/resizing so that the axes stay aligned according to the margin options
        g(1,1).set_text_options('base_size',12)
        g(1,1).set_names('x','','y','Mean');
        g(1,1).axe_property('XTickLabel',''); % We deactivate the ticks
        g(1,1).set_line_options('base_size',3)
        g(1,1).set_color_options('map',[0 0 0])
        
        
        g(2,1)=gramm('x',plotdata.tone,'y',plotdata.data, 'color', plotdata.group); %,'y',cars.Acceleration,'color',cars.Cylinders,'subset',cars.Cylinders~=3 & cars.Cylinders~=5
        % g(2,1).geom_point(); %to see each point, removed the YLim (highest points are closer to 2.5
        g(2,1).stat_summary('type','sem','geom','errorbar'); %mean and sem shown
        g(2,1).stat_summary('type','sem','geom','point'); %mean and sem shown
        g(2,1).stat_summary('type','sem','geom','area'); %mean and sem shown
        g(2,1).set_layout_options('Position',[0 0 0.8 0.8],...
            'legend_pos',[0.75 0.75 0.2 0.2],... %We detach the legend from the plot and move it to the top right
            'margin_height',[0.1 0.02],...
            'margin_width',[0.1 0.02],...
            'redraw',false);
        g(2,1).set_text_options('base_size',12)
        g(2,1).axe_property('Ygrid','on'); %,'YLim',[0 0.0015]
        g(2,1).set_names('x','Tone','y','mV','color','Group');
        g(2,1).set_color_options('map','matlab');
        g(2,1).axe_property('XTickLabel',{'-3', '-2', '-1', 'BF', '+1', '+2', '+3'})
        
        
        g(3,1)=gramm('x',plotdata.data,'color',plotdata.group);
        g(3,1).set_layout_options('Position',[0.8 0 0.2 0.8],...
            'legend',false,...
            'margin_height',[0.1 0.02],...
            'margin_width',[0.02 0.05],...
            'redraw',false);
        g(3,1).set_text_options('base_size',12)
        g(3,1).set_names('x','');
        g(3,1).stat_bin('geom','overlaid_bar','fill','transparent'); %histogram
        g(3,1).coord_flip();
        g(3,1).axe_property('XTickLabel','');
        g(3,1).set_color_options('map','matlab');
        
        
        g.set_title([(Order{isink}) ' ' (SParameter{ispara})]);
        g.draw();
        g.export('file_name',['Awake_' (Order{isink}) '_' (SParameter{ispara})], 'file_type','png');
        g.export('file_name',['Awake_' (Order{isink}) (SParameter{ispara})], 'file_type','pdf');
        close all;
        
        %% Muscimol vs Anesthetized
        fullmean = nanmean(vertcat(expandAn,expandM),1)';
        
        % Create Appropriate Structure
        plotdata = struct;
        
        tone = ticks';
        Angrammdata = expandAn(1,:)';
        Mugrammdata = expandM(1,:)';
        
        Angroup = {'Anesthetized' 'Anesthetized' 'Anesthetized' 'Anesthetized' 'Anesthetized' 'Anesthetized' 'Anesthetized'}; %7 to match the amount of ticks
        An = Angroup';
        Mugroup = {'Muscimol' 'Muscimol' 'Muscimol' 'Muscimol' 'Muscimol' 'Muscimol' 'Muscimol'};
        Mu = Mugroup';
        
        for istack = 1:size(expandAn,1)-1
            
            tone = vertcat(tone, ticks');
            Angrammdata = vertcat(Angrammdata, expandAn(istack+1,:)');
            Mugrammdata = vertcat(Mugrammdata, expandM(istack+1,:)');
            An = vertcat(An, Angroup');
            Mu = vertcat(Mu, Mugroup');
            
        end
        
        plotdata.tone = vertcat(tone, tone);
        plotdata.data = vertcat(Angrammdata, Mugrammdata);
        plotdata.group = vertcat(An, Mu);
        
        
        clear g
        figure('Position',[100 100 650 650]);
        
        g(1,1)=gramm('x',ticks','y',fullmean);
        g(1,1).geom_line();
        g(1,1).set_layout_options('Position',[0 0.8 0.8 0.2],... %Set the position in the figure (as in standard 'Position' axe property)
            'legend',false,... % No need to display legend for side histograms
            'margin_height',[0.02 0.05],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
            'margin_width',[0.1 0.02],...
            'redraw',false); %We deactivate automatic redrawing/resizing so that the axes stay aligned according to the margin options
        g(1,1).set_text_options('base_size',12)
        g(1,1).set_names('x','','y','Mean');
        g(1,1).axe_property('XTickLabel',''); % We deactivate the ticks
        g(1,1).set_line_options('base_size',3)
        g(1,1).set_color_options('map',[0 0 0])
        
        
        g(2,1)=gramm('x',plotdata.tone,'y',plotdata.data, 'color', plotdata.group); %,'y',cars.Acceleration,'color',cars.Cylinders,'subset',cars.Cylinders~=3 & cars.Cylinders~=5
        % g(2,1).geom_point(); %to see each point, removed the YLim (highest points are closer to 2.5
        g(2,1).stat_summary('type','sem','geom','errorbar'); %mean and sem shown
        g(2,1).stat_summary('type','sem','geom','point'); %mean and sem shown
        g(2,1).stat_summary('type','sem','geom','area'); %mean and sem shown
        g(2,1).set_layout_options('Position',[0 0 0.8 0.8],...
            'legend_pos',[0.75 0.75 0.2 0.2],... %We detach the legend from the plot and move it to the top right
            'margin_height',[0.1 0.02],...
            'margin_width',[0.1 0.02],...
            'redraw',false);
        g(2,1).set_text_options('base_size',12)
        g(2,1).axe_property('Ygrid','on'); %,'YLim',[0 0.0015]
        g(2,1).set_names('x','Tone','y','mV','color','Group');
        g(2,1).set_color_options('map','matlab');
        g(2,1).axe_property('XTickLabel',{'-3', '-2', '-1', 'BF', '+1', '+2', '+3'})
        
        
        g(3,1)=gramm('x',plotdata.data,'color',plotdata.group);
        g(3,1).set_layout_options('Position',[0.8 0 0.2 0.8],...
            'legend',false,...
            'margin_height',[0.1 0.02],...
            'margin_width',[0.02 0.05],...
            'redraw',false);
        g(3,1).set_text_options('base_size',12)
        g(3,1).set_names('x','');
        g(3,1).stat_bin('geom','overlaid_bar','fill','transparent'); %histogram
        g(3,1).coord_flip();
        g(3,1).axe_property('XTickLabel','');
        g(3,1).set_color_options('map','matlab');
        
        
        g.set_title([(Order{isink}) ' ' (SParameter{ispara})]);
        g.draw();
        g.export('file_name',['Musc_' (Order{isink}) '_' (SParameter{ispara})], 'file_type','png');
        g.export('file_name',['Musc_' (Order{isink}) (SParameter{ispara})], 'file_type','pdf');
        close all;
    end
end