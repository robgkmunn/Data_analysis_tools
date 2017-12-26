

dir = table;
slope = vertcat(dat.slope,dat.slope_sq);
%slope_comp = vertcat(dat.slope_x_sq,dat.slope_y_sq)
dir.slope = slope;
direction = cell((length(dat.slope)*2),1);
direction(1:length(dat.slope),1) = {'Baseline'};
direction((length(dat.slope)+1):(length(dat.slope)*2),1) = {'Expansion'};
dir.dir = direction;
neuron(1:length(dat.slope),1) = 1:length(dat.slope);
neuron((length(dat.slope)+1):(length(dat.slope)*2),1) = 1:length(dat.slope);
dir.neuron = neuron;

statarray = grpstats(dir.slope,dir.dir);
means = statarray';
ax = {'Baseline','Expansion'};
ax = ax';
mean = cell2table(ax);
mean.mean = means';
mean.group = ones(2,1);


g1 = gramm('x',dir.dir,'y',dir.slope,'group',dir.neuron);
g1.set_title('Slope (Both Axes) - Cells that pass in Both')
g1.geom_line('alpha',0.5)
g1.set_line_options('base_size',1.5)
g1.set_color_options('map',[0 0 0])
g1.set_order_options('x',0)
g1.set_names('x','Axis','y','Slope of Speed/Firing Rate Function')
g1.no_legend()
g1.draw()
g1.update('color',dir.dir)
g1.geom_point('alpha',0.5)
g1.set_point_options('base_size',8)
g1.set_color_options('map',[1 0.8398 0
                            1 0.5469 0])
g1.no_legend()
g1.draw()
g1.update('x',mean.ax,'y',mean.mean,'group',mean.group)
g1.set_order_options('x',0)
g1.geom_point('alpha',0.5)
g1.set_color_options('map','brewer')
g1.draw()
g1.update('x',mean.ax,'y',mean.mean,'group',mean.group)
g1.set_order_options('x',0)
g1.geom_line()
g1.set_line_options('base_size',3,'style',{':' '-' '--' '-.'})
g1.set_color_options('map',[1 0.569 0])
set([g1.facet_axes_handles],'LineWidth',2)
g1.set_text_options('base_size',14)
%g1.axe_property('YLim',[-0.02 0.23]);
g1.draw()

dir = table;
slope = vertcat(dat.slope,dat.slope_sq);
%slope_comp = vertcat(selectdat.slope_x_sq,selectdat.slope_y_sq)
dir.slope = slope;
direction = cell((length(dat.slope)*2),1);
direction(1:length(dat.slope),1) = {'Baseline'};
direction((length(dat.slope)+1):(length(dat.slope)*2),1) = {'Expansion'};
dir.dir = direction;
neuron(1:length(dat.slope),1) = 1:length(dat.slope);
neuron((length(dat.slope)+1):(length(dat.slope)*2),1) = 1:length(dat.slope);
dir.neuron = neuron;

%%%%%% MAKE ARRAY OF CELLS THAT PASS IN BOTH OF AND SQ AND SCATTERPLOT %%%

rows = dat.speedScore > 0.044 & dat.speedScore_sq > 0.044;
selectdat = dat(rows,:);

g2 = gramm('x',selectdat.mean_speed,'y',selectdat.mean_speed_sq);
g2.geom_point('alpha',0.5)
g2.stat_glm()
g2.set_title('Mean Running Speed','FontSize',12)
g2.stat_cornerhist('nbins',8,'geom','line')
g2.set_color_options('map',[1 0.5469 0])
g2.set_names('x','Mean Running Speed in Baseline (cm/s)','y','Mean Running Speed in Expansion (Hz)')
g2.set_point_options('base_size',8)
g2.axe_property('LineWidth', 1.5)
g2.geom_abline('intercept',0,'slope',1,'style','r:')
g2.set_text_options('base_size',14)
%g2.axe_property('XLim',[-0.05  0.2])
%g2.axe_property('YLim',[-0.05  0.2])
g2.draw()


%%%%%%%%%%%%%%%% SLOPE PLOT WITH ONLY CELLS THAT PASS IN BOTH %%%%%%%%%%%
dir = table;
slope = vertcat(selectdat.slope,selectdat.slope_sq);
%slope_comp = vertcat(selectdat.slope_x_sq,selectdat.slope_y_sq)
dir.slope = slope;
direction = cell((length(selectdat.slope)*2),1);
direction(1:length(selectdat.slope),1) = {'Baseline'};
direction((length(selectdat.slope)+1):(length(selectdat.slope)*2),1) = {'Expansion'};
dir.dir = direction;
neuron(1:length(selectdat.slope),1) = 1:length(selectdat.slope);
neuron((length(selectdat.slope)+1):(length(selectdat.slope)*2),1) = 1:length(selectdat.slope);
dir.neuron = neuron;


statarray = grpstats(dir.slope,dir.dir);
means = statarray';
ax = {'Baseline','Expansion'};
ax = ax';
mean = cell2table(ax);
mean.mean = means';
mean.group = ones(2,1);

figure()
g1 = gramm('x',dir.dir,'y',dir.slope,'group',dir.neuron);
g1.set_title('Slope (Both Axes)')
g1.geom_line('alpha',0.5)
g1.set_line_options('base_size',1.5)
g1.set_color_options('map',[0 0 0])
g1.set_order_options('x',0)
g1.set_names('x','Condition','y','Slope of Speed/Firing Rate Function')
g1.no_legend()
g1.draw()
g1.update('color',dir.dir)
g1.geom_point('alpha',0.5)
g1.set_point_options('base_size',8)
g1.set_color_options('map',[1 0.5469 0
                            1 0.5469 0])
g1.no_legend()
g1.draw()
g1.update('x',mean.ax,'y',mean.mean,'group',mean.group)
g1.set_order_options('x',0)
g1.geom_point('alpha',0.5)
g1.set_color_options('map','brewer')
g1.draw()
g1.update('x',mean.ax,'y',mean.mean,'group',mean.group)
g1.set_order_options('x',0)
g1.geom_line()
g1.set_line_options('base_size',3,'style',{':' '-' '--' '-.'})
g1.set_color_options('map',[1 0.569 0])
set([g1.facet_axes_handles],'LineWidth',2)
g1.set_text_options('base_size',14)
%g1.axe_property('YLim',[-0.02 0.23]);
g1.draw()
