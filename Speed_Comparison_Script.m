%%%%%% Custom script for plotting wildtype against trip8b speed cell data
%%%%%% from Munn....Giocomo et.al 2018
%%%%%% Robert Munn, PhD, Stanford University 26 March 2018

%%%% Expects two .xlsx sheets with columns of values for the speed score,
%%%% slope, speed etc data. Uses the excellent GRAMM library by Pier Morel to do the
%%%% plotting: https://github.com/piermorel/gramm

data = uigetfile('*.xlsx','Choose WT Data');
wt = readtable(data);
clear data
data = uigetfile('*.xlsx','Choose T8B Data');
t8b = readtable(data);

% make new tables containing only cells that pass the speed thresholdin
% both the open field baseline and tyhe compression in both the wildtype
% and trip8b groups
rows_t8b = t8b.speedScore >= 0.044 & t8b.speedScore_sq >= 0.044;
t8b_both = t8b(rows_t8b,:);

rows_wt = wt.speedScore >= 0.044 & wt.speedScore_sq >= 0.044;
wt_both = wt(rows_wt,:);

% make direction array 2x the length of t8b array in order to plot the axes
% against each other; also concatenate the axes data into one array
direction(1:length(t8b_both.slope),1) = {'Static Axis'}; direction(length(t8b_both.slope):length(t8b_both.slope)*2,1) = {'Compressed Axis'};
condition(1:length(t8b_both.slope),1) = {'Baseline'}; condition(length(t8b_both.slope):length(t8b_both.slope)*2,1) = {'Compression'};
t8b_slopes = vertcat(t8b_both.slope,t8b_both.slope_sq);
t8b_slopes_bas = vertcat(t8b_both.slope_x,t8b_both.slope_y);
t8b_slopes_sq = vertcat(t8b_both.slope_x_sq,t8b_both.slope_y_sq);
t8b_score = vertcat(t8b_both.speedScore,t8b_both.speedScore_sq);

% basic plots of trip8b speed cells

figure() % Overall speed score differences
g = gramm('x',condition,'y',t8b_score);
g.geom_jitter('dodge',0.7)
g.stat_violin('half','true','fill','edge')
g.stat_boxplot('width',0.1)
g.set_names('x',{},'y','Speed Score')
g.coord_flip()
g.set_color_options('map',[0.19 0.64 0.33])
g.set_text_options('font','Helvetica','base_size',12,'legend_title_scaling',0.75,'legend_scaling',0.75,'title_scaling',1)
g.axe_property('LineWidth', 1.5)
g.draw()
set([g.results.stat_boxplot.outliers_handle],'visible','off')
set([g.results.stat_boxplot.box_handle],'FaceColor','none')
[p_score,h_score,stat_score] = signrank(t8b_both.speedScore,t8b_both.speedScore_sq);

figure() %Baseline slopes in each axis
g1 = gramm('x',direction,'y',t8b_slopes_bas);
g1.geom_jitter('dodge',0.7)
g1.stat_violin('half','true','fill','edge')
g1.stat_boxplot('width',0.1)
g1.set_names('x',{},'y','Slope of the Running Speed/Firing Rate Fit')
g1.coord_flip()
g1.set_color_options('map',[0.19 0.64 0.33])
g1.set_text_options('font','Helvetica','base_size',12,'legend_title_scaling',0.75,'legend_scaling',0.75,'title_scaling',1)
g1.axe_property('LineWidth', 1.5)
g1.draw()
set([g1.results.stat_boxplot.outliers_handle],'visible','off')
set([g1.results.stat_boxplot.box_handle],'FaceColor','none')
[p_t8b_baseline_slopes,h_t8b_baseline_slopes,stats_t8b_baseline_slopes] = signrank(t8b_both.slope_x,t8b_both.slope_y);

figure() %Compressed slopes in each axis
g2 = gramm('x',direction,'y',t8b_slopes_sq);
g2.geom_jitter('dodge',0.7)
g2.stat_violin('half','true','fill','edge')
g2.stat_boxplot('width',0.1)
g2.set_names('x',{},'y','Slope of the Running Speed/Firing Rate Fit')
g2.coord_flip()
g2.set_color_options('map',[0.19 0.64 0.33])
g2.set_text_options('font','Helvetica','base_size',12,'legend_title_scaling',0.75,'legend_scaling',0.75,'title_scaling',1)
g2.axe_property('LineWidth', 1.5)
g2.draw()
set([g2.results.stat_boxplot.outliers_handle],'visible','off')
set([g2.results.stat_boxplot.box_handle],'FaceColor','none')
[p_t8b_compression_slopes,h_t8b_compression_slopes,stats_t8b_compression_slopes] = signrank(t8b_both.slope_x_sq,t8b_both.slope_y_sq);

figure() % Overall slopes in the baseline vs. compression
g3 = gramm('x',condition,'y',t8b_slopes);
g3.geom_jitter('dodge',0.7)
g3.stat_violin('half','true','fill','edge')
g3.stat_boxplot('width',0.1)
g3.set_names('x',{},'y','Slope of the Firing Rate/Running Speed Fit')
g3.coord_flip()
g3.set_color_options('map',[0.19 0.64 0.33])
g3.set_text_options('font','Helvetica','base_size',12,'legend_title_scaling',0.75,'legend_scaling',0.75,'title_scaling',1)
g3.axe_property('LineWidth', 1.5)
g3.draw()
set([g3.results.stat_boxplot.outliers_handle],'visible','off')
set([g3.results.stat_boxplot.box_handle],'FaceColor','none')

% Make arrays combining wt and t8b data and a genotype var the length of
% each vertical array
allslope = vertcat(wt_both.slope,t8b_both.slope);   allslope_sq = vertcat(wt_both.slope_sq,t8b_both.slope_sq);
allslope_x = vertcat(wt_both.slope_x,t8b_both.slope_x);   allslope_x_sq = vertcat(wt_both.slope_x_sq,t8b_both.slope_x_sq);
allslope_y = vertcat(wt_both.slope_y,t8b_both.slope_y);   allslope_y_sq = vertcat(wt_both.slope_y_sq,t8b_both.slope_y_sq);
allspeed = vertcat(wt_both.speedScore,t8b_both.speedScore); allspeed_sq = vertcat(wt_both.speedScore_sq,t8b_both.speedScore_sq);
wt_both.slope_diff = wt_both.slope_sq - wt_both.slope;  wt_both.score_diff = wt_both.speedScore_sq - wt_both.speedScore;
t8b_both.slope_diff = t8b_both.slope_sq - t8b_both.slope;   t8b_both.score_diff = t8b_both.speedScore_sq - t8b_both.speedScore;
allslope_diff = vertcat(wt_both.slope_diff,t8b_both.slope_diff);    allspeed_diff = vertcat(wt_both.score_diff,t8b_both.score_diff);
genotype(1:length(wt_both.slope),1) = {'Wildtype'}; genotype(length(wt_both.slope):length(wt_both.slope)+length(t8b_both.slope),1) = {'Trip8b'};

figure() % Differences in overall slope between wt and t8b
g4 = gramm('x',genotype,'y',allslope_diff,'color',genotype);
g4.geom_jitter()
g4.stat_violin('half','true','fill','edge')
g4.stat_boxplot('width',0.1)
g4.set_color_options('map',[0.19 0.51 0.74
    0.19 0.64 0.33])
g4.set_names('x',{},'y','Difference in Slope Between Baseline and Compression')
g4.set_order_options('color',{'Wildtype','Trip8b'},'x',-1)
g4.set_limit_extra([0,0.2],[0.1 0])
g4.set_text_options('font','Helvetica','base_size',12,'legend_title_scaling',0.75,'legend_scaling',0.75,'title_scaling',1)
g4.axe_property('LineWidth', 1.5)
g4.no_legend()
g4.coord_flip()
g4.draw()
set([g4.results.stat_boxplot.outliers_handle],'visible','off')
set([g4.results.stat_boxplot.box_handle],'FaceColor','none')
[h_overallslope_diff,p_overallslope_diff,stat_overallslope_diff] = kstest2(t8b_both.slope_diff,wt_both.slope_diff);

figure() % Differences in overall speed score between wt and t8b
g4 = gramm('x',genotype,'y',allspeed_diff,'color',genotype);
g4.geom_jitter()
g4.stat_violin('half','true','fill','edge')
g4.stat_boxplot('width',0.1)
g4.set_color_options('map',[0.19 0.51 0.74
    0.19 0.64 0.33])
g4.set_names('x',{},'y','Difference in Speed Score (Baseline - Compression)')
g4.set_order_options('color',{'Wildtype','Trip8b'},'x',-1)
g4.set_limit_extra([0,0.2],[0.1 0])
g4.set_text_options('font','Helvetica','base_size',12,'legend_title_scaling',0.75,'legend_scaling',0.75,'title_scaling',1)
g4.axe_property('LineWidth', 1.5)
g4.no_legend()
g4.coord_flip()
g4.draw()
set([g4.results.stat_boxplot.outliers_handle],'visible','off')
set([g4.results.stat_boxplot.box_handle],'FaceColor','none')
[h_overallscore_diff,p_overallscore_diff,stat_overallscore_diff] = kstest2(t8b_both.score_diff,wt_both.score_diff);

figure() % Compare axes in the baseline (t8b)
g5 = gramm('x',t8b_both.slope_x,'y',t8b_both.slope_y);
g5.geom_point()
g5.geom_abline('style','r--')
g5.set_color_options('map',[0.19 0.64 0.33])
g5.set_names('x','Slope (Static Axis)','y','Slope (Compressed Axis)')
g5.stat_cornerhist('nbins',10)
g5.set_text_options('font','Helvetica','base_size',12,'legend_title_scaling',0.75,'legend_scaling',0.75,'title_scaling',1)
g5.axe_property('LineWidth', 1.5,'XLim',[-0.05 0.2],'Ylim',[-0.05 0.2])
g5.draw()
[p_axes_baseline,h_axes_baseline,stat_axes_baseline] = signrank(t8b_both.slope_x,t8b_both.slope_y);

figure() % Compare axes in the compression (t8b)
g5 = gramm('x',t8b_both.slope_x_sq,'y',t8b_both.slope_y_sq);
g5.geom_point()
g5.geom_abline('style','r--')
g5.set_color_options('map',[0.19 0.64 0.33])
g5.set_names('x','Slope (Static Axis)','y','Slope (Compressed Axis)')
g5.stat_cornerhist('nbins',10)
g5.set_text_options('font','Helvetica','base_size',12,'legend_title_scaling',0.75,'legend_scaling',0.75,'title_scaling',1)
g5.axe_property('LineWidth', 1.5,'XLim',[-0.05 0.2],'YLim',[-0.05 0.2])
g5.draw()
[p_axes_compression,h_axes_compression,stat_axes_compression] = signrank(t8b_both.slope_x_sq,t8b_both.slope_y_sq);


both_comp_ax = vertcat(wt_both.slope_y_sq,t8b_both.slope_y_sq);
both_stat_ax = vertcat(wt_both.slope_x_sq,t8b_both.slope_x_sq);

figure()
g6 = gramm('x',both_stat_ax,'y',both_comp_ax,'color',genotype);
g6.geom_point()
g6.geom_abline('style','r--')
g6.stat_cornerhist('nbins',10)
g6.set_order_options('color',{'Wildtype','Trip8b'})
g6.set_color_options('map',[0.19 0.51 0.74
    0.19 0.64 0.33])
g6.axe_property('LineWidth',1.5)
g6.set_names('x','Speed Slope (Static Axis)','y','Speed Slope (Compressed Axis)')
g6.draw()

allrunspd = vertcat(wt_both.mean_speed,t8b_both.mean_speed);


figure() % Differences in overall mean running speed
g7 = gramm('x',genotype,'y',allrunspd,'color',genotype);
g7.geom_jitter()
g7.stat_boxplot('width',1.2)
g7.set_color_options('map',[0.19 0.51 0.74
    0.19 0.64 0.33])
g7.set_names('x',{},'y','Running Speed (cm/s)')
g7.set_order_options('color',{'Wildtype','Trip8b'},'x',-1)
%g7.set_limit_extra([0,0.2],[0.1 0])
g7.set_text_options('font','Helvetica','base_size',12,'legend_title_scaling',0.75,'legend_scaling',0.75,'title_scaling',1)
g7.axe_property('LineWidth', 1.5)
g7.no_legend()
g7.draw()
set([g7.results.stat_boxplot.outliers_handle],'visible','off')
[h_overallspeed_diff,p_overallspeed_diff,stat_overallspeed_diff] = kstest2(t8b_both.mean_speed,wt_both.mean_speed);

probt8b = fitdist(t8b_both.slope_diff,'normal'); % make cumulative probability distributions
x = [0 0.025 0.05 0.075 0.1 0.125 0.15 0.175 0.2 0.25 0.25];
tripcdf = cdf(probt8b,x);
probwt = fitdist(wt_both.slope_diff,'normal');
wtcdf = cdf(probwt,x);

figure()
g8 = gramm('x',x,'y',wtcdf);
g8.geom_line()
g8.set_color_options('map',[0.19 0.51 0.74])
g8.set_text_options('font','Helvetica','base_size',12,'legend_title_scaling',0.75,'legend_scaling',0.75,'title_scaling',1)
g8.axe_property('LineWidth', 1.5)
g8.set_names('x','Difference in Slope (Baseline - Compression)','y','Cumulative Probability');
g8.draw()
g8.update('x',x,'y',tripcdf);
g8.set_color_options('map',[0.19 0.64 0.33]);
g8.geom_line()
g8.draw()

figure()
g9 = gramm('x',allslope_diff,'color',genotype);
g9.stat_bin('normalization','cdf','nbins',100,'geom','stairs')
g9.set_color_options('map',[0.19 0.51 0.74
    0.19 0.64 0.33])
g9.set_names('x','Difference in Slope (Baseline - Compression)','y','Cumulative Probability')
g9.axe_property('LineWidth',1.5)
g9.set_text_options('font','Helvetica','base_size',12,'legend_title_scaling',0.75,'legend_scaling',0.75,'title_scaling',1)
g9.draw()

figure()

