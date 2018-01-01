function [] = directional_controls_v1_2

%This script performs the following checks on head direction/speed cell
%data
%
%1. Amount of time sampling each of 45-315 (N),  315-225 (W), 225-135 (S),
%135-45 (E) directional bins
%
%2. Speed distribution
%
%3. Comparison of speed distributions and angular sampling between a
%baseline and spatially modified condition. Required fields in the .csv are
%pointers to each individual cell's database files (from database maker) as
%"Sessions" - pointers to the tetrode and unit number (under headings
%"Tetrode" and "Unit". 
% Also checks which wall has moved (if applicable) under
% "Squish_Direction". Requires a heading "Box" for scale information.

%
% V1.0 Robert Munn, Stanford University. 8 August 2017
% V1.1 Change from traditional spreadsheet data input to table
%      Also added check to ensure correct orientation 29 Dec 2017
% V1.2 Clean up code and add comments. Various analyses are added
%      (commented) that can be individually uncommented. 1 Jan 2017

table = uigetfile('*.csv');
A = readtable(table);
numCell = numel(A.Tetrode);
spreadsheet = input('What shall we call the output file?','s');
posfile = cell(length(A.Tetrode),1);
spikefile = cell(length(A.Tetrode),1);
direction = cell(length(A.Tetrode),1);
numbins = 36;
for n = 1:length(A.Tetrode)
    posfile{n} = strcat(A.Sessions{n},'_pos.mat');
    spikefile{n} = strcat(A.Sessions{n},'_T',num2str(A.Tetrode(n)),'C',num2str(A.Unit(n)),'.mat');
    direction{n} = A.Squish_Direction{n};
    Box(n) = A.Box(n);
end

directional_array = cell(length(A.Tetrode),numbins); % preallocate arrays
direction_bins = ones(length(A.Tetrode),numbins);

for n = 1:length(A.Tetrode)
edges = [2:2:80];
hdMatrix = [];
load(posfile{n});
load(spikefile{n});
% load the position and spike times

posy = -posy;
posy2 = -posy2; % flip db output to be not a mirror image

if strcmp(direction{n},"West") % if the stable wall (cue wall) is west % rotate 90 degrees ccw
    posx_new = -posy; posx2_new = -posy2; posy = posx; posy2 = posx2;
    posx = posx_new; posx2 = posx2_new;
    disp("did a rotation, yo")
    
else
    disp("Stayed the same")
end

    minX = nanmin(posx);
    maxX = nanmax(posx);
    minY = nanmin(posy);
    maxY = nanmax(posy);
    xLength = maxX - minX;
    yLength = maxY - minY;
    sLength = max([xLength, yLength]);
    scale = Box(n) / sLength;
    posx = posx * scale;
    posy = posy * scale;
    posx2 = posx2 * scale;
    posy2 = posy2 * scale;
%%%%%% Speed filter %%%%%%

% Low speed threshold. Segments of the path where the mouse moves slower
% than this threshold will be removed.
p.lowSpeedThreshold = 2.5; % [cm/s]

% High speed threshold. Segments of the path where the mouse moves faster
% than this threshold will be removed.
p.highSpeedThreshold = 50; % [cm/s]

posy_y = posy;
posy_y2 = posy2;
posx_x = posx;
posx_x2 = posx2;

if p.lowSpeedThreshold > 0 || p.highSpeedThreshold > 0
    % Calculate the speed of the mouse, sample by sample
    speed = speed2D(posx,posy,post);
    speedy = speed1D_rob(posy,post);
    speedx = speed1D_rob(posx,post);
    
    if p.lowSpeedThreshold > 0 && p.highSpeedThreshold > 0
        bad_ind = find(speed < p.lowSpeedThreshold | speed > p.highSpeedThreshold);
    elseif p.lowSpeedThreshold > 0 && p.highSpeedThreshold == 0
        bad_ind = find(speed < p.lowSpeedThreshold );
    else
        bad_ind = find(speed > p.highSpeedThreshold );
    end
    
        if p.lowSpeedThreshold > 0 && p.highSpeedThreshold > 0
        bad_indy = find(speedy < p.lowSpeedThreshold | speedy > p.highSpeedThreshold);
    elseif p.lowSpeedThreshold > 0 && p.highSpeedThreshold == 0
        bad_indy = find(speedy < p.lowSpeedThreshold );
    else
        bad_indy = find(speedy > p.highSpeedThreshold );
        end
    
    if p.lowSpeedThreshold > 0 && p.highSpeedThreshold > 0
        bad_indx = find(speedx < p.lowSpeedThreshold | speedx > p.highSpeedThreshold);
    elseif p.lowSpeedThreshold > 0 && p.highSpeedThreshold == 0
        bad_indx = find(speedx < p.lowSpeedThreshold );
    else
        bad_indx = find(speedx > p.highSpeedThreshold );
    end
    
    % Remove the segments that have too high or too low speed
    posx(bad_ind) = NaN;posx_x(bad_indx) = NaN;
    posy(bad_ind) = NaN;posy_y(bad_indy) = NaN;
    posx2(bad_ind) = NaN;posx_x2(bad_indx) = NaN;
    posy2(bad_ind) = NaN;posy_y2(bad_indy) = NaN;
    
    
    filtered_speed = speed2D(posx,posy,post);
    filtered_speedy = speed1D_rob(posy_y,post);
    filtered_speedx = speed1D_rob(posx_x,post);
end
%%%%%%%%%%%%% SPEED CONTROLS %%%%%%%%%%%%%%%%%%%%%
[speed_count] = histcounts(filtered_speed,edges); % Omnibus
speed_array(n,:) = speed_count;
speed_weight = speed_array(n,:);
speed_dist = speed_weight.*edges(1,1:39);
speed_dist = sum(speed_dist(:,:));
speed_out(1,1:39) = edges(1,1:39);
speed_out(n+1,1:39) = speed_array(n,:);

[speed_countx] = histcounts(filtered_speedx,edges); % X Direction
speed_arrayx(n,:) = speed_countx;
speed_weightx = speed_arrayx(n,:);
speed_distx = speed_weightx.*edges(1,1:39);
speed_distx = sum(speed_distx(:,:));
speed_outx(1,1:39) = edges(1,1:39);
speed_outx(n+1,1:39) = speed_arrayx(n,:);

[speed_county] = histcounts(filtered_speedy,edges); %Y direction
speed_arrayy(n,:) = speed_county;
speed_weighty = speed_arrayy(n,:);
speed_disty = speed_weighty.*edges(1,1:39);
speed_disty = sum(speed_disty(:,:));
speed_outy(1,1:39) = edges(1,1:39);
speed_outy(n+1,1:39) = speed_arrayy(n,:);

%%%%%%%%%%%%% HEAD DIRECTION CONTROLS %%%%%%%%%%%%%%
hdDir =  atan2(posy2-posy,posx2-posx);
hdDir(isnan(hdDir)) = []; %take out the nans
%hdDir = smooth(hdDir,2,'rlowess'); %smooth the data
pos = find(hdDir > 0);
neg = find(hdDir < 0);
positive(:,n) = (length(pos)/length(hdDir))*100/1;
negative(:,n) = (length(neg)/length(hdDir))*100/1;

hdDir(hdDir < 0) = hdDir(hdDir<0)+2*pi;
[dir(1:length(hdDir),n)] = discretize(hdDir,4);  
if dir(:,n) == 0
    dir = [];
end
q = length(dir(:,n));
q1 = sum(dir(:,n) == 1); q2 = sum(dir(:,n) == 2); q3 = sum(dir(:,n) == 3); q4 = sum(dir(:,n) == 4);

divider = q1+q2+q3+q4;
percq1 = (q1/divider)*100/1; quadrants(n,1) = percq1(1,1);
percq2 = (q2/divider)*100/1; quadrants(n,2) = percq2(1,1);
percq3 = (q3/divider)*100/1; quadrants(n,3) = percq3(1,1);
percq4 = (q4/divider)*100/1; quadrants(n,4) = percq4(1,1);

failsafe = sum(quadrants(n,:)); % check to make sure this calculation is done properly
if failsafe >= 99.99 & failsafe<= 100.001
    disp 'yay'
else
    disp 'ruh roh raggy'
end

Npart1 = find(hdDir > 5.49779);
Npart2 = find(hdDir <= 0.785398);
if size(Npart1,1) > 1 % check whether this comes out as rows or columns and orthogonalize
    N = cat(1,Npart1,Npart2);
elseif size(Npart1,2) > 1
    N = cat(2,Npart1,Npart2);
end
E = find(hdDir > 0.785398 & hdDir <= 2.35619);
S = find(hdDir > 2.35619 & hdDir <= 3.92699);
W = find(hdDir > 3.92699 & hdDir <= 5.49779);
totallen = length(N)+length(E)+length(S)+length(W);
North(n,:) = (length(N)/totallen)*100/1; East(n,:) = (length(E)/totallen)*100/1;South(n,:) = (length(S)/totallen)*100/1;West(n,:) = (length(W)/totallen)*100/1;
quad(n,1) = North(n,:); quad(n,2) = East(n,:); quad(n,3) = South(n,:);quad(n,4) = West(n,:);
if North(n,:)+East(n,:)+South(n,:)+West(n,:) >= 99.999 & North(n,:)+East(n,:)+South(n,:)+West(n,:) <= 100.0001
    disp 'it worked (% quadrants add up to 100)'
else
    disp 'you messed this up somewhere (% quadrants do not add up to 100)'
end

maxlength(1,n) = length(hdDir);
dir_array{:,n} = hdDir;
x = [1:2*pi];
d = linspace(0,2*pi,length(hdDir)); % make an array between 0 and 2pi with n equally spaced steps
pd = makedist('Uniform','lower',0,'upper',(2*pi));
cd = cdf(pd,d);
[f,x] = ecdf(hdDir); %make contin
[direction_count(:,n), binning] = histcounts(hdDir,20); % count the number of directional epochs in each of 20 bins
mean_circular_direction = circ_mean(hdDir); %compute mean direction
if mean_circular_direction<0
    mean_circular_direction = mean_circular_direction+(2*pi);
else
end

vector = circ_r(hdDir);% compute directional vector
rayleigh = circ_rtest(hdDir);
[pval, m] = circ_otest(hdDir);

h = hist(hdDir,numbins);
if size(h,2) < numbins
    h(1,end:numbins) = NaN;
else
end
direction_bins(n,1:numbins) = h;
for t = 1:numbins
directional_array{n,t} = h(t);
end

end


%%%%%%%%%%%%%%%%ANALYSIS AND PLOTTING%%%%%%%%%%%%%%%
meanpos = mean(positive);
meanneg = mean(negative);
errorpos = std(positive);
errorneg = std(negative);

posneg = [meanpos meanneg];
posnegerror = [errorpos errorneg];
posnegnames = {'0-180','180-360'};
figure()
hold on
b1 = bar(posneg);
e1 = errorbar(posneg,posnegerror, '.');
set(gca,'XTick',1:4,'XTickLabel',posnegnames)
e.Color = 'black';
b.FaceAlpha = 0.6;
ylabel('Mean percent time');
hold off

figure(); hold on
for k = 1:numCell
 colorVec = jet(k);
 plot(edges(1,1:39),speed_array(k,:), 'Color', colorVec(k,:));
 xlabel('cm/s')
 ylabel('time spent in speed epoch')
 title('running speed distributions')
end
 hold off

%for s = 1:numCell
 %   if length(dir_array{1,s})< max(maxlength)
  %       pad = max(maxlength);
   %      dir_array{1,s}(end:pad,1) = NaN;
   % end
%end

%new_array = dir_array{1}
%ca = cell(max(maxlength),numCell);
%for row = 1:size(new_array,1)
   % ca{row,1} = new_array(row)
%end
    
%hd_array = NaN(max(maxlength),numCell);  %pulls directional values out of cell array and puts them in a matrix by preallocating a n x n matrix of NaNs. EXTREMELY TIME CONSUMING
%for g = 1:numCell
%hd_array(:,g) = dir_array{1,g}
%end
%mean_hd = total_hd/numCell
binning(~0) = [];
thing = direction_bins(:,1:4);temp = direction_bins(:,33:36);
nord = horzcat(thing,temp);
west = direction_bins(:,24:31);
sud = direction_bins(:,15:22);
ost = direction_bins(:,6:13);

NS = horzcat(nord,sud);
EW = horzcat(west,ost);
NS = mean(NS,2);
EW = mean(EW,2);
[pvalue,hypothesis,statssignrank] = signrank(NS,EW)
dir_plotter = horzcat(NS,EW);
figure()
nbp = notBoxPlot(dir_plotter);
ax = gca
ax.XTickLabels = {"North/South","East/West"};
ylabel('Number of Time Bins Facing Direction');
%%%%%%%%%%%%% speed plotting %%%%%%%%%%%%%%%%%%%%%
mean_speed = mean(speed_array); % plot overall mean speed distribution
stdev_speed = std(speed_array,0,1);
sterror_speed = stdev_speed/(sqrt(length(speed_array)));
two_sterror_speed = sterror_speed*2;
figure()
hold on
color = [0 0.4275 0.1725]
s = shadedErrorBar(edges(1:1:39),mean_speed,two_sterror_speed,{'-','markerfacecolor',[0.1922 0.6392 0.3294]},1)
 s.patch.FaceColor = [color];
 s.mainLine.Color = [color];
 s.edge(1,1).Color = [color];
 s.edge(1,2).Color = [color];
xlabel('speed (cm/s)'),ylabel('number of time bins');
title('Mean running speed');
axis([0 50 0 (max(mean_speed)+2000)])
hold off

mean_speedx = mean(speed_arrayx); % plot mean speed distribution in the X direction
stdev_speedx = std(speed_arrayx,0,1);
sterror_speedx = stdev_speedx/(sqrt(length(speed_arrayx)));
two_sterror_speedx = sterror_speedx*2;
figure()
hold on
s1 = shadedErrorBar(edges(1:1:39),mean_speedx,two_sterror_speedx,{'-','markerfacecolor',[0.1922 0.6392 0.3294]},1)
 s1.patch.FaceColor = [color];
 s1.mainLine.Color = [color];
 s1.edge(1,1).Color = [color];
 s1.edge(1,2).Color = [color];
xlabel('speed (cm/s)'),ylabel('number of time bins');
title('Mean running speed in the X axis');
axis([0 50 0 (max(mean_speedx)+2000)])
hold off

mean_speedy = mean(speed_arrayy); % plot mean speed distribution in the Y direction
stdev_speedy = std(speed_arrayy,0,1);
sterror_speedy = stdev_speedy/(sqrt(length(speed_arrayy)));
two_sterror_speedy = sterror_speedy*2;
figure()
hold on
s2 = shadedErrorBar(edges(1:1:39),mean_speedy,two_sterror_speedy,{'-','markerfacecolor',[0.1922 0.6392 0.3294]},1)
 s2.patch.FaceColor = [color];
 s2.mainLine.Color = [color];
 s2.edge(1,1).Color = [color];
 s2.edge(1,2).Color = [color];
xlabel('speed (cm/s)'),ylabel('number of time bins');
title('Mean running speed in the Y axis');
axis([0 50 0 (max(mean_speedy)+2000)])
hold off

meanN = mean(North);meanE = mean(East);meanS = mean(South);meanW = mean(West);
centered_quads = [meanN meanE meanS meanW];
errorN = std(North)/sqrt(numCell); errorE = std(East)/sqrt(numCell); errorS = std(South)/sqrt(numCell); errorW = std(West)/sqrt(numCell);
centered_error = 2*[errorN errorE errorS errorW];

axisnames = {'N','S','W','E'};
%%%%%%%%%%%%%%% direction plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure()
%hold on
%b = bar(eighthnames,mean_eighth);
%e = errorbar(eighthnames,mean_eighth,erroreighth, '.');
%e.Color = 'black';
%b.FaceAlpha = 0.6;
%ylabel('Mean percent time facing');
%hold off
figure()
hold on
b = bar(centered_quads);
e = errorbar(centered_quads,centered_error, '.');
set(gca,'XTick',1:4,'XTickLabel',axisnames)
e.Color = 'black';
b.FaceAlpha = 0.6;
ylabel('Mean percent time');
hold off
%figure()
%hold on
%b = bar(mean_NSEW);
%e = errorbar(mean_NSEW,error_NSEW, '.');
%e.Color = 'black';
%b.FaceAlpha = 0.6;
%ylabel('Mean percent time');
%hold off
%plusstdev = mean_direction + stdev_direction;
%minusstdev = mean_direction - stdev_direction;
%sterr_direction = stdev_direction/(sqrt(length(numCell)));
%meandir = mean(dir);
%pol = histogram(meandir);

%speed_test = [speed_arrayx speed_arrayy];
%speed_groups(1:length(speed_arrayx(1,:))) = 1;
%speed_groups(length(speed_arrayx(1,:)):((length(speed_arrayx(1,:))+1)+(length(speed_arrayx(1,:))-1))) = 2;
%[p_quad,tbl_quad,stats_quad] = anova1(quad);
%[cquad,mquad,hquad,nmsquad] = multcompare(stats_quad);
%[p_speed,tbl_speed,stats_speed] = anova1(speed_test,speed_groups);
%[cspeed,mspeed,hspeed,nmsspeed] = multcompare(stats_speed);
%[p_kw_quadrant,tbl_kw_quad,stats_kw_quadrant] = kruskalwallis(quadrant);
%[p_kw_percent,tbl_kw_percent,stats_kw_percent] = kruskalwallis(percent);
%[p_kw_speed,tbl_kw_speed,stats_kw_speed] = kruskalwallis(speed_test,speed_groups);
speedName = [spreadsheet,'binned_speed_array','.xlsx'];
fileName = [spreadsheet,'_binned_direction_array','.xlsx'];
xlswrite(fileName,directional_array);
xlswrite(speedName,speed_out);