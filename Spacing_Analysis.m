
clear all
file = uigetfile('*.mat','Select the analysis workspace "Spacing Analysis"');
load(file);
%%%%%%%%%%%%%% Set up arrays for plotting %%%%%%%%%%%%%
xShift_wt = xShift_wt';yShift_wt = yShift_wt';
xShift_t8b = xShift_t8b';yShift_t8b = yShift_t8b';
groups(1:length(lambda_wt),1) = {'Wildtype'};
groups(length(groups):(length(groups)+length(lambda_t8b)),1) = {'Trip8b'};
lambda = vertcat(lambda_wt,lambda_t8b);
rho = vertcat(rho_max_wt,rho_max_t8b);
xShift = vertcat(xShift_wt,xShift_t8b);
yShift = vertcat(yShift_wt,yShift_t8b);
spacing_wt = vertcat(spacing_ratio_wt(:,1),spacing_ratio_wt(:,2));
spacing_t8b = vertcat(spacing_ratio_t8b(:,1),spacing_ratio_t8b(:,2));
spacing_wt_groups(1:length(lambda_wt),1) = {'Baseline'};
spacing_wt_groups(length(spacing_wt_groups):length(spacing_wt_groups)+length(lambda_wt),1) = {'Compression'};
spacing_t8b_groups(1:length(lambda_t8b),1) = {'Baseline'};
spacing_t8b_groups(length(spacing_t8b_groups):length(spacing_t8b_groups)+length(spacing_ratio_t8b(:,2)),1) = {'Compression'};
spacing_wt_genotype(1:2*length(lambda_wt),1) = {'Wildtype'};
spacing_t8b_genotype(1:2*length(lambda_t8b),1) = {'Trip8b'}
spacing = vertcat(spacing_wt,spacing_t8b);
sp_groups = vertcat(spacing_wt_groups,spacing_t8b_groups);
sp_genotype = vertcat(spacing_wt_genotype,spacing_t8b_genotype);
s = table(spacing,sp_groups,sp_genotype);
s.Properties.VariableNames = {'Spacing_ratio','group','genotype'}
%%%%%%%%%%%% plot some figures %%%%%%%%%%%%%%%%%%%
figure()
g = gramm('x',lambda,'color',groups);
g.stat_bin('nbins',400,'normalization','cdf','geom','line')
g.set_color_options('map',[0.1922 0.6392 0.3294
    0.4196 0.6824 0.8392])
g.set_names('x','Stretch Factor (lambda)','y','Cumulative Probability','color','')
g.set_line_options('base_size',3)
g.set_text_options('base_size',14)
g.draw

figure()
g1 = gramm('x',lambda,'color',groups);
g1.stat_bin('nbins',8,'geom','line')
g1.set_color_options('map',[0.1922 0.6392 0.3294
    0.4196 0.6824 0.8392])
g1.set_names('x','Stretch Factor (lambda)','y','Number of Neurons','color','')
g1.set_line_options('base_size',3)
g1.set_text_options('base_size',14)
g1.draw

figure()
g2 = gramm('x',rho,'color',groups);
g2.stat_bin('nbins',100,'normalization','cdf','geom','line')
g2.set_color_options('map',[0.1922 0.6392 0.3294
    0.4196 0.6824 0.8392])
g2.set_names('x','Map to map correlation (Rho)','y','Cumulative Probability','color','')
g2.set_line_options('base_size',3)
g2.set_text_options('base_size',14)
g2.draw

figure()
g3 = gramm('x',xShift,'color',groups);
g3.stat_bin('nbins',100,'normalization','cdf','geom','line')
g3.set_color_options('map',[0.1922 0.6392 0.3294
    0.4196 0.6824 0.8392])
g3.set_names('x','Optimum X Axis Shift','y','Cumulative Probability','color','')
g3.set_line_options('base_size',3)
g3.set_text_options('base_size',14)
g3.draw

figure()
g4 = gramm('x',yShift,'color',groups);
g4.stat_bin('nbins',100,'normalization','cdf','geom','line')
g4.set_color_options('map',[0.1922 0.6392 0.3294
    0.4196 0.6824 0.8392])
g4.set_names('x','Optimum Y Axis Shift','y','Cumulative Probability','color','')
g4.set_line_options('base_size',3)
g4.set_text_options('base_size',14)
g4.draw

figure()
g5 = gramm('x',s.group,'y',s.Spacing_ratio,'color',s.genotype);
g5.stat_boxplot()
g5.set_color_options('map',[0.1922 0.6392 0.3294
    0.4196 0.6824 0.8392])
g5.set_names('x','','y','Min/Max Spacing Ratio','color','')
g5.set_line_options('base_size',3)
g5.set_text_options('base_size',14)
g5.draw()


%%%% Do a stat %%%%%%%%%%%

[h_lambda,p_lambda,ks_lambda] = kstest2(lambda_wt,lambda_t8b);
[h_rho,p_rho,ks_rho] = kstest2(rho_max_wt,rho_max_t8b);
[h_xShift,p_xShift,ks_xShift] = kstest2(xShift_wt,xShift_t8b);
[h_yShift,p_yShift,ks_yShift] = kstest2(yShift_wt,yShift_t8b);
[p_wt_spacing,h_wt_spacing,stats_wt_spacing] = ranksum(spacing_ratio_wt(:,1),spacing_ratio_wt(:,2));
[p_t8b_spacing,h_t8b_spacing,stats_t8b_spacing] = ranksum(spacing_ratio_t8b(:,1),spacing_ratio_t8b(:,2));