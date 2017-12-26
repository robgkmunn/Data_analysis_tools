function [outfile] = Conjunctive_Analysis

table = uigetfile('*xlsx');
A = readtable(table);
numCell = numel(A.Tetrode);

for n = 1:numCell % do all the open field sessions
    posfile{n} = strcat(A.Sessions{n},'_pos.mat');
    spikefile{n} = strcat(A.Sessions{n},'_T',num2str(A.Tetrode(n)),'C',num2str(A.Unit(n)),'.mat');
    direction{n} = A.Squish_Direction{n};
    Box(n) = A.Box(n);
    OF{1,1}='Animal';OF{1,2}='date';OF{1,3}='mvl';OF{1,4}='mv_arg';OF{1,5}='peak_rate';OF{1,6}='pref_angle';OF{1,7}='stability';OF{1,8}='speedScore';OF{1,9}='speedScore_x';OF{1,10}='speedScore_y';OF{1,11}='Intercept';OF{1,12}='slope';OF{1,13}='intercept_x';OF{1,14}='slope_x';OF{1,15} ='intercept_y';OF{1,16}='slope_y';OF{1,17}='speedStability';
    OF{1,18} = 'Speed_Stability_x';OF{1,19} = 'Speed_Stability_y';OF{1,20} = 'mean_speed';OF{1,21} = 'mean_speed_x';OF{1,22} = 'mean_speed_y';OF{1,23} = 'mean_fr';OF{1,24} = 'mean_fr_x';OF{1,25} = 'mean_fr_y';OF{1,26} = 'dynamic_range';OF{1,27} = 'dynamic_range_x';OF{1,28} = 'dynamic_range_y';OF{1,29} = 'lin_intercept';OF{1,30} = 'exp_intercept';OF{1,31} = 'lin_slope';OF{1,32} = 'exp_slope';OF{1,33} = 'lin_intercept_x';OF{1,34} = 'exp_intercept_x';OF{1,35} = 'lin_slope_x';OF{1,36} = 'exp_slope_x';OF{1,37} = 'lin_intercept_y';OF{1,38} = 'exp_intercept_y';OF{1,39} = 'lin_slope_y';OF{1,40} = 'exp_slope_y'; 
    [OF{n+1,1},OF{n+1,2},OF{n+1,3},OF{n+1,4},OF{n+1,5},OF{n+1,6},OF{n+1,7}] = head_direction(posfile{n},spikefile{n},Box(n),direction{n});
    [OF{n+1,8},OF{n+1,9},OF{n+1,10},OF{n+1,11},OF{n+1,12},OF{n+1,13},OF{n+1,14},OF{n+1,15},OF{n+1,16},OF{n+1,17},OF{n+1,18},OF{n+1,19},OF{n+1,20},OF{n+1,21},OF{n+1,22},OF{n+1,23},OF{n+1,24},OF{n+1,25},OF{n+1,26},OF{n+1,27},OF{n+1,28},OF{n+1,29},OF{n+1,30},OF{n+1,31},OF{n+1,32},OF{n+1,33},OF{n+1,34},OF{n+1,35},OF{n+1,36},OF{n+1,37},OF{n+1,38},OF{n+1,39},OF{n+1,40},speed_bins(:,n)] = speedScore_stability_rm_v2(spikefile{n},posfile{n},direction{n},Box(n));
end

mean_wt_speed = mean(speed_bins,2);sd_wt_speed = std(speed_bins,[],2);
mean_wt_speed = mean_wt_speed'; sd_wt_speed = sd_wt_speed';
%OF{1,41} = 'binned speed'; OF{1,42} = 'stdev speed';
%OF{2:end,41} = mean_wt_speed; OF{2:end,42} = sd_wt_speed;
OF_dat = OF(2:end,:);
Open_Field = cell2table(OF_dat);
Open_Field.Properties.VariableNames = OF(1,:);

for ns = 1:numCell % do all the squish sessions
    posfile_SQ{ns} = strcat(A.Sessions_1{ns},'_pos.mat');
    spikefile_SQ{ns} = strcat(A.Sessions_1{ns},'_T',num2str(A.Tetrode_1(ns)),'C',num2str(A.Unit_1(ns)),'.mat');
    direction_SQ{ns} = A.Squish_Direction_1{ns};
    Box_SQ(ns) = A.Box_1(ns);
    SQ{1,1}='Animal_sq';SQ{1,2}='date_sq';SQ{1,3}='mvl_sq';SQ{1,4}='mv_arg_sq';SQ{1,5}='peak_rate_sq';SQ{1,6}='pref_angle_sq';SQ{1,7}='stability_sq';SQ{1,8}='speedScore_sq';SQ{1,9}='speedScore_x_sq';SQ{1,10}='speedScore_y_sq';SQ{1,11}='Intercept_sq';SQ{1,12}='slope_sq';SQ{1,13}='intercept_x_sq';SQ{1,14}='slope_x_sq';SQ{1,15} ='intercept_y_sq';SQ{1,16}='slope_y_sq';SQ{1,17}='speedStability_sq';
    SQ{1,18} = 'Speed_Stability_x_sq';SQ{1,19} = 'Speed_Stability_y_sq';SQ{1,20} = 'mean_speed_sq';SQ{1,21} = 'mean_speed_x_sq';SQ{1,22} = 'mean_speed_y_sq';SQ{1,23} = 'mean_fr_sq';SQ{1,24} = 'mean_fr_x_sq';SQ{1,25} = 'mean_fr_y_sq';SQ{1,26} = 'dynamic_range_sq';SQ{1,27} = 'dynamic_range_x_sq';SQ{1,28} = 'dynamic_range_y_sq';SQ{1,29} = 'lin_intercept_sq';SQ{1,30} = 'exp_intercept_sq';SQ{1,31} = 'lin_slope_sq';SQ{1,32} = 'exp_slope_sq';SQ{1,33} = 'lin_intercept_x_sq';SQ{1,34} = 'exp_intercept_x_sq';SQ{1,35} = 'lin_slope_x_sq';SQ{1,36} = 'exp_slope_x_sq';SQ{1,37} = 'lin_intercept_y_sq';SQ{1,38} = 'exp_intercept_y_sq';SQ{1,39} = 'lin_slope_y_sq';SQ{1,40} = 'exp_slope_y_sq'; 
    [SQ{ns+1,1},SQ{ns+1,2},SQ{ns+1,3},SQ{ns+1,4},SQ{ns+1,5},SQ{ns+1,6},SQ{ns+1,7}] = head_direction(posfile_SQ{ns},spikefile_SQ{ns},Box_SQ(ns),direction_SQ{ns});
    [SQ{ns+1,8},SQ{ns+1,9},SQ{ns+1,10},SQ{ns+1,11},SQ{ns+1,12},SQ{ns+1,13},SQ{ns+1,14},SQ{ns+1,15},SQ{ns+1,16},SQ{ns+1,17},SQ{ns+1,18},SQ{ns+1,19},SQ{ns+1,20},SQ{ns+1,21},SQ{ns+1,22},SQ{ns+1,23},SQ{ns+1,24},SQ{ns+1,25},SQ{ns+1,26},SQ{ns+1,27},SQ{ns+1,28},SQ{ns+1,29},SQ{ns+1,30},SQ{ns+1,31},SQ{ns+1,32},SQ{ns+1,33},SQ{ns+1,34},SQ{ns+1,35},SQ{ns+1,36},SQ{ns+1,37},SQ{ns+1,38},SQ{ns+1,39},SQ{ns+1,40}] = speedScore_stability_rm_v2(spikefile_SQ{ns},posfile_SQ{ns},direction_SQ{ns},Box_SQ(ns));
end

SQ_dat = SQ(2:end,:);
Squish = cell2table(SQ_dat);
Squish.Properties.VariableNames = SQ(1,:);
dat = [Open_Field,Squish];
writetable(dat,'Output.xlsx','WriteVariableNames',1,'Filetype','spreadsheet')
filename = 'Output.mat';
save(filename,'dat');
keyboard
change_MVL = A.MVL - A.MVL_1;
change_Score = A.Score - A.Score_1;
change_Slope = A.Slope - A.Slope_1;

mvl = vertcat(dat.mvl,dat.mvl_SQ);
group(1:47,1) = 1; group(48:94,1) = 2;

g = gramm('x',group,'y',mvl,'color',group);
g.geom_point();
g.geom_line();
g.draw

[r,p] = corr(change_MVL,change_Score);
[r1,p1] = corr(change_MVL,change_Slope);

plot(change_MVL,change_Score,'ko')
plot(change_MVL,change_Slope,'rO')

% line-group comparisons
figure()
compareMVL(1,:) = dat.mvl;  compareMVL(2,:) = dat.mvl_SQ;
p = plot(compareMVL);
set(p,'LineWidth',2);
set(p,'Color',[0.6602 0.6602 0.6602]);
set(p,'Marker','O');
set(p,'MarkerEdge','k');
axis([ 0.75 2.25 0 1]);
ylabel('Mean Vector Length')
xticks([1 2])
xticklabels({'Open Field','Compression'});
ax = gca;
ax.FontName = 'Myriad Pro';
ax.FontSize = 12;
meanmvl(1,1) = mean(dat.mvl);
meanmvl(2,1) = mean(dat.mvl_SQ);
hold on 
p1 = plot(meanmvl);
set(p1,'Linewidth',2,'Color','r','Marker','O');
hold off

figure()
compareangle(1,:) = dat.pref_angle;  compareangle(2,:) = dat.pref_angle_SQ;
p = plot(compareangle);
set(p,'LineWidth',2);
set(p,'Color',[0.6602 0.6602 0.6602]);
set(p,'Marker','O');
set(p,'MarkerEdge','k');
axis([ 0.75 2.25 0 400]);
ylabel('Preferred Phase Angle (Degrees)')
xticks([1 2])
xticklabels({'Open Field','Compression'});
ax = gca;
ax.FontName = 'Myriad Pro';
ax.FontSize = 12;

figure()
compareslope(1,:) = dat.slope;  compareslope(2,:) = dat.slope_SQ;
p = plot(compareslope);
set(p,'LineWidth',2);
set(p,'Color',[0.6602 0.6602 0.6602]);
set(p,'Marker','O');
set(p,'MarkerEdge','k');
axis([ 0.75 2.25 -0.05  0.15]);
ylabel('Slope')
xticks([1 2])
xticklabels({'Open Field','Compression'});
ax = gca;
ax.FontName = 'Myriad Pro';
ax.FontSize = 12;
meanslope(1,1) = mean(dat.slope);
meanslope(2,1) = mean(dat.slope_SQ);
hold on 
p1 = plot(meanslope);
set(p1,'Linewidth',2,'Color','r','Marker','O');
hold off

figure()
compareslopex(1,:) = dat.slope_x;  compareslopex(2,:) = dat.slope_x_SQ;
p = plot(compareslopex);
set(p,'LineWidth',2);
set(p,'Color',[0.6602 0.6602 0.6602]);
set(p,'Marker','O');
set(p,'MarkerEdge','k');
axis([ 0.75 2.25 -0.05  0.15]);
ylabel('Slope (X axis)')
xticks([1 2])
xticklabels({'Open Field','Compression'});
ax = gca;
ax.FontName = 'Myriad Pro';
ax.FontSize = 12;
meanslopex(1,1) = mean(dat.slope_x);
meanslopex(2,1) = mean(dat.slope_x_SQ);
hold on 
p1 = plot(meanslopex);
set(p1,'Linewidth',2,'Color','r','Marker','O');
hold off

figure()
compareslopey(1,:) = dat.slope_y;  compareslopey(2,:) = dat.slope_y_SQ;
p = plot(compareslopey);
set(p,'LineWidth',2);
set(p,'Color',[0.6602 0.6602 0.6602]);
set(p,'Marker','O');
set(p,'MarkerEdge','k');
axis([ 0.75 2.25 -0.05  0.15]);
ylabel('Slope (Y axis)')
xticks([1 2])
xticklabels({'Open Field','Compression'});
ax = gca;
ax.FontName = 'Myriad Pro';
ax.FontSize = 12;
meanslopey(1,1) = mean(dat.slope_y);
meanslopey(2,1) = mean(dat.slope_y_SQ);
hold on 
p1 = plot(meanslopey);
set(p1,'Linewidth',2,'Color','r','Marker','O');
hold off

figure()
compareintercept(1,:) = dat.Intercept;  compareintercept(2,:) = dat.Intercept_SQ;
p = plot(compareintercept);
set(p,'LineWidth',2);
set(p,'Color',[0.6602 0.6602 0.6602]);
set(p,'Marker','O');
set(p,'MarkerEdge','k');
axis([ 0.75 2.25 0 7]);
ylabel('Intercept')
xticks([1 2])
xticklabels({'Open Field','Compression'});
ax = gca;
ax.FontName = 'Myriad Pro';
ax.FontSize = 12;
meanintercept(1,1) = mean(dat.Intercept);
meanintercept(2,1) = mean(dat.Intercept_SQ);
hold on 
p1 = plot(meanintercept);
set(p1,'Linewidth',2,'Color','r','Marker','O');
hold off

end


