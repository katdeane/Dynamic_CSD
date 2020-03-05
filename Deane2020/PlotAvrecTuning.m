function PlotAvrecTuning(feature, Name, plottune, home)

for i1 = 1:length(feature)
% tuning curves
clear g
figure();

g(1,1)=gramm('x',plottune.frqz,'y',plottune.(feature{i1}), 'color', plottune.tgroup); %,'y',cars.Acceleration,'color',cars.Cylinders,'subset',cars.Cylinders~=3 & cars.Cylinders~=5
g(1,1).stat_summary('type','sem','geom','errorbar'); %mean and sem shown
g(1,1).stat_summary('type','sem','geom','point'); %mean and sem shown
g(1,1).stat_summary('type','sem','geom','area'); %mean and sem shown
g(1,1).set_text_options('base_size',12) 
if contains(feature{i1},'_norm')
    g(1,1).set_names('x','Tone','y','%','color','Group');
    g(1,1).axe_property('YLim',[.25 1.5]);
elseif strcmp(feature{i1},'peaklate')
    g(1,1).set_names('x','Tone','y','ms','color','Group');
    g(1,1).axe_property('YLim',[220 270]);
elseif strcmp(feature{i1},'peakamp')
    g(1,1).set_names('x','Tone','y','mV/mm²','color','Group');
    g(1,1).axe_property('YLim',[0 0.003]);
else
    g(1,1).set_names('x','Tone','y','mV/mm²','color','Group');
    g(1,1).axe_property('YLim',[0.0002 0.0015]);
end
g(1,1).set_color_options('map','matlab');
g(1,1).axe_property('XTickLabel',{'-2', '-1', 'BF', '+1', '+2'})

g.set_title(Name{i1});
g.draw();
g.export('file_name',Name{i1}, 'file_type','pdf');
g.export('file_name',Name{i1}, 'file_type','png');
close all;

% box plots, single feature
clear g
figure();

g=gramm('x',plottune.tgroup,'y',plottune.(feature{i1}), 'color', plottune.tgroup);
g.facet_grid(plottune.frqz,[]);
g.stat_boxplot('notch',true);
if strcmp(feature{i1},'peaklate')
    g(1,1).set_names('x','Group','y','Latency (ms)','color','Group');
elseif strcmp(feature{i1},'peakamp')
    g(1,1).set_names('x','Group','y','Amplituded (mV/mm²)','color','Group');
else
    g(1,1).set_names('x','Tone','y','RMS (mV/mm²)','color','Group');

end
g(1,1).set_color_options('map','matlab');

g.set_title([Name{i1} 'boxplot']);
g.draw();
g.export('file_name',[Name{i1} 'boxplot'], 'file_type','pdf');
g.export('file_name',[Name{i1} 'boxplot'], 'file_type','png');
close all;


end

%all frequencies amp vs latency
clear g
figure();

g=gramm('x',plottune.peaklate,'y',plottune.peakamp, 'color', plottune.animal, 'column', plottune.tgroup);
g.facet_grid(plottune.frqz,[]);
g.set_point_options('base_size',2);
g.geom_point();
% g.set_color_options('map','matlab');

g.set_title('Peak Latency against Amp');
g.draw();
g.export('file_name','Peak Latency against Amp', 'file_type','pdf');
g.export('file_name','Peak Latency against Amp', 'file_type','png');
close all;



BFlogical = strcmp('d BF',plottune.frqz);
min2logical = strcmp('b -2',plottune.frqz);

BF_avrec_late = plottune.peaklate(BFlogical);
BF_avrec_amp = plottune.peakamp(BFlogical);
BF_avrec_group = plottune.tgroup(BFlogical);
BF_avrec_animal = plottune.animal(BFlogical);

m2_avrec_late = plottune.peaklate(min2logical);
m2_avrec_amp = plottune.peakamp(min2logical);
m2_avrec_group = plottune.tgroup(min2logical);
m2_avrec_animal = plottune.animal(min2logical);


%% peak latency against peak amplitute
clear g
figure();

g=gramm('x',BF_avrec_late,'y',BF_avrec_amp, 'color', BF_avrec_group);
g.geom_point();
% g.axe_property('YLim',[-1 0.5],'XLim',[200 300]);
g.set_color_options('map','matlab');

g.set_names('x','Peak Latency','y','Peak Amplitude','color','Group');
g.set_title('Avrec Peak Latency against Amp BF');
g.draw();
g.export('file_name','Peak Latency against Amp BF', 'file_type','pdf');
g.export('file_name','Peak Latency against Amp BF', 'file_type','png');
close all;

clear g
figure();

g=gramm('x',m2_avrec_late,'y',m2_avrec_amp, 'color', m2_avrec_group);
g.geom_point();
g.set_color_options('map','matlab');

g.set_names('x','Peak Latency','y','Peak Amplitude','color','Group');
g.set_title('Avrec Peak Latency against Amp BF-2');
g.draw();
g.export('file_name','Peak Latency against Amp BF-2', 'file_type','pdf');
g.export('file_name','Peak Latency against Amp BF-2', 'file_type','png');
close all;

cd(home);cd DATA; cd avrec_compare;
save('AvrecPlotData_single.mat','plottune')


clear g
figure();

g=gramm('x',BF_avrec_late,'y',BF_avrec_amp, 'color', BF_avrec_animal);
g.geom_point();
% g.axe_property('YLim',[-1 0.5],'XLim',[200 300]);
% g.set_color_options('map','matlab');

g.set_names('x','Peak Latency','y','Peak Amplitude','color','Group');
g.set_title('Avrec Peak Latency against Amp BF');
g.draw();
g.export('file_name','Peak Latency against Amp BF Animal', 'file_type','pdf');
g.export('file_name','Peak Latency against Amp BF Animal', 'file_type','png');
close all;

clear g
figure();

g=gramm('x',m2_avrec_late,'y',m2_avrec_amp, 'color', m2_avrec_animal);
g.geom_point();
% g.set_color_options('map','matlab');

g.set_names('x','Peak Latency','y','Peak Amplitude','color','Group');
g.set_title('Avrec Peak Latency against Amp BF-2');
g.draw();
g.export('file_name','Peak Latency against Amp BF-2 Animal', 'file_type','pdf');
g.export('file_name','Peak Latency against Amp BF-2 Animal', 'file_type','png');
close all;

cd(home);cd DATA; cd avrec_compare;
save('AvrecPlotData_single.mat','plottune')