%%%%%% Custom script for plotting wildtype against trip8b speed cell data
%%%%%% from Munn....Giocomo et.al 2018
%%%%%% Robert Munn, PhD, Stanford University 26 March 2018

%%%% Expects two .xlsx sheets with columns of values for the speed score,
%%%% slope, speed etc data

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
t8b_slopes = vertcat(t8b_both.slope_x,t8b_both.slope_y);
t8b_slopes_sq = vertcat(t8b_both.slope_x_sq,t8b_both.slope_y_sq);

figure() % basic plots of trip8b speed cells
g1 = gramm('x',direction,'y',t8b_slopes);
g1.geom_jitter()
g1.stat_violin('half','true')
g1.set_names('x',{},'y','Slope of the Running Speed/Firing Rate Fit')
g1.coord_flip()
g1.set_color_options('map',[0.19 0.64 0.33])
g1.draw()
[p_t8b_baseline_slopes,h_t8b_baseline_slopes,stats_t8b_baseline_slopes] = signrank(t8b_both.slope_x,t8b_both.slope_y);

figure()
g2 = gramm('x',direction,'y',t8b_slopes_sq);
g2.geom_jitter()
g2.stat_violin('half','true')
g2.set_names('x',{},'y','Slope of the Running Speed/Firing Rate Fit')
g2.coord_flip()
g2.set_color_options('map',[0.19 0.64 0.33])
g2.draw()
[p_t8b_compression_slopes,h_t8b_compression_slopes,stats_t8b_compression_slopes] = signrank(t8b_both.slope_x_sq,t8b_both.slope_y_sq);

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

figure() 
g1 = gramm('x',genotype,'y',allslope_diff,'color',genotype);
g1.geom_jitter()
g1.stat_violin('half','left')
g1.coord_flip()
g1.draw()