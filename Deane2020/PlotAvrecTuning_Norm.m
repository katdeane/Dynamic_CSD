function PlotAvrecTuning_Norm(feature, Name, plottune)

for i1 = 1:length(feature)
% tuning curves
clear g
figure();

g(1,1)=gramm('x',plottune.frqzMir,'y',plottune.(feature{i1}), 'color', plottune.tgroup); %,'y',cars.Acceleration,'color',cars.Cylinders,'subset',cars.Cylinders~=3 & cars.Cylinders~=5
g(1,1).stat_summary('type','std','geom','errorbar'); %mean and std shown
g(1,1).stat_summary('type','std','geom','point'); %mean and std shown
g(1,1).stat_summary('type','std','geom','area'); %mean and std shown
g(1,1).set_text_options('base_size',12) 

g(1,1).set_names('x','Tone','y','%','color','Group');
g(1,1).axe_property('YLim',[.25 1.5]);

g(1,1).set_color_options('map','matlab');
g(1,1).axe_property('XTickLabel',{'BF', '+/-1', '+/-2'})

g.set_title(Name{i1});
g.draw();
g.export('file_name',Name{i1}, 'file_type','pdf');
g.export('file_name',Name{i1}, 'file_type','png');
close all;

end