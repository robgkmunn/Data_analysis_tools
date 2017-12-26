function [speedScore, speedScore_x, speedScore_y, intercept, slope, intercept_x, slope_x, intercept_y, slope_y, speedStability, speedStability_x, speedStability_y, mean_speed, mean_speed_x, mean_speed_y, mean_fr, mean_fr_x, mean_fr_y, dynamic_range, dynamic_range_x, dynamic_range_y, lin_intercept, exp_intercept, lin_slope, exp_slope, lin_intercept_x, exp_intercept_x, lin_slope_x, exp_slope_x, lin_intercept_y, exp_intercept_y, lin_slope_y, exp_slope_y,speed_bins] = speedScore_stability_rm_v2(spikefile,posfile,direction,Box)
load(spikefile);
load(posfile);
makeplot = 1;
makeplotx = 1;
makeploty = 1;
makeheatplot = 0;
posy = -posy;
posy2 = -posy2;

%%%%%% Scale the coordinates using the shape information %%%%%%

% rotate box so that direction is correct

if strcmp(direction,"West") % if the stable wall (cue wall) is west % rotate 90 degrees ccw
    posx_new = -posy; posx2_new = -posy2; posy = posx; posy2 = posx2;
    posx = posx_new; posx2 = posx2_new;
    disp("did a rotation, yo")
    rotation = 'No Rotation'
else
    rotation = 'Rotated'
end

    minX = nanmin(posx); maxX = nanmax(posx);
    minY = nanmin(posy); maxY = nanmax(posy);
    xLength = maxX - minX; yLength = maxY - minY; sLength = max([xLength, yLength]);
    scale = Box / sLength;
  
if xLength > yLength+40 
    posx = posx * scale; posy = posy * scale;
    posx2 = posx2 * scale; posy2 = posy2 * scale;
   disp('This was a Squish or Expansion Baseline')
   Squish = 'True'
else
    posx = posx * scale; posy = posy * scale;
    posx2 = posx2 * scale; posy2 = posy2 * scale;
   disp('This was a normal 1x1m box')
   Squish = 'False'
end
     
 minX_fix = nanmin(posx); maxX_fix = nanmax(posx); % sanity check
 minY_fix = nanmin(posy); maxY_fix = nanmax(posy); 
 xLength_fix = maxX_fix - minX_fix; yLength_fix = maxY_fix - minY_fix;

% take out NaN's and replace them with neighboring values
positions = {posx, posy};

for k = 1:2
    pos_temp = positions{k};
    nan_ind = find(isnan(pos_temp));
    for m = 1:numel(nan_ind)
        if nan_ind(m) - 1 == 0
            temp = find(~isnan(pos_temp),1,'first');
            pos_temp(nan_ind(m)) = pos_temp(temp);
        else
            pos_temp(nan_ind(m)) = pos_temp(nan_ind(m)-1);
        end
    end
    
    positions{k} = pos_temp;
    
end

boxSize = Box;
posx = positions{1}; posy = positions{2};

% compute speed at every time point
velx = diff([posx(1); posx]); vely = diff([posy(1); posy]); dt = 0.02;
speed = sqrt(velx.^2+vely.^2)/dt;
speed(speed > 100) = NaN;
speedx = abs(velx/dt);
speedx(speedx > 100) = NaN;
speedy = abs(vely/dt);
speedy(speedy > 100) = NaN;

binWidth = 2;
speedVec = 5:binWidth:50;
sessionInd = round(linspace(1,numel(post),5));

h = hist(cellTS,post); %the tracking sampling rate is 50Hz and thus 20ms frames. The speed has the same                                                                            %sampling rate
fr = gauss_smoothing(h,20)*50;  % see attachment.
select = speed > 2;
speed(isnan(speed)) = interp1(find(~isnan(speed)), speed(~isnan(speed)), find(isnan(speed)), 'pchip'); % interpolate NaNs
speed = gauss_smoothing(speed,10); % gaussian kernel smoothing

selectx = speedx > 2;
speedx(isnan(speedx)) = interp1(find(~isnan(speedx)), speedx(~isnan(speedx)), find(isnan(speedx)), 'pchip'); % interpolate NaNs
speedx = gauss_smoothing(speedx,10); % gaussian kernel smoothing

selecty = speedy > 2;
speedy(isnan(speedy)) = interp1(find(~isnan(speedy)), speedy(~isnan(speedy)), find(isnan(speedy)), 'pchip'); % interpolate NaNs
speedy = gauss_smoothing(speedy,10); % gaussian kernel smoothing

meanFR = numel(cellTS)/max(post); % mean firing rate over whole session

% Bin speed for tuning curve
numBins = 10; %define the number of bins to use for speed
maxSpeed = 100; % cm/s
speedAxis = linspace(0,maxSpeed,numBins+1); %create values for binning
speed_tuning = nan(numBins,1); %make array for tuning curve

for z = 1:numBins
    indices = find(speed >= speedAxis(z) & speed < speedAxis(z+1)); % times when the animal runs at a certain speed
    speed_tuning(z) = nanmean(h(indices))*50;   
    number_spikes(z) = nansum(h(indices))*50;
end

speedAxisx = linspace(0,maxSpeed,numBins+1); %create values for binning %%%%%%%%%%%%%%%%%%%%%%% X AXIS
speed_tuningx = nan(numBins,1); %make array for tuning curve

for zx = 1:numBins
    indicesx = find(speedx >= speedAxisx(zx) & speedx < speedAxisx(zx+1)); % times when the animal runs at a certain speed
    speed_tuningx(zx) = nanmean(h(indicesx))*50;   
    number_spikes(zx) = nansum(h(indicesx))*50;
end

speedAxisy = linspace(0,maxSpeed,numBins+1); %create values for binning %%%%%%%%%%%%%%%%%%%%%%% Y AXIS
speed_tuningy = nan(numBins,1); %make array for tuning curve

for zy = 1:numBins
    indicesy = find(speedy >= speedAxisy(zy) & speed < speedAxisy(zy+1)); % times when the animal runs at a certain speed
    speed_tuningy(zy) = nanmean(h(indicesy))*50;   
    number_spikes(zy) = nansum(h(indicesy))*50;
end


% calculate speed score by correlating instantaneous firing rate and speed
speed = reshape(speed,numel(post),1);
speedx = reshape(speedx,numel(post),1);
speedy = reshape(speedy,numel(post),1);
fr = reshape(fr,numel(post),1);
speedFilt = speed(select);
speedxFilt = speedx(selectx);
speedyFilt = speedy(selecty);
frFilt = fr(select);
frFiltx = fr(selectx);
frFilty = fr(selecty);
speedScore =corr(speedFilt,frFilt);
speedScore_x =corr(speedxFilt,frFiltx);
speedScore_y =corr(speedyFilt,frFilty);
slope = [ones(sum(select),1) speedFilt]\frFilt;
intercept = slope(1);
slope = slope(2);
slope_x = [ones(sum(selectx),1) speedxFilt]\frFiltx;
intercept_x = slope_x(1);
slope_x = slope_x(2);
slope_y = [ones(sum(selecty),1) speedyFilt]\frFilty;
intercept_y = slope_y(1);
slope_y = slope_y(2);
mean_speed = mean(speedFilt);
max_speed = max(speedFilt)
mean_speed_x = mean(speedxFilt);
mean_speed_y = mean(speedyFilt);
mean_fr = nanmean(frFilt);
mean_fr_x = nanmean(frFiltx);
mean_fr_y = nanmean(frFilty);
dynamic_range = (max(frFilt) - min(frFilt))/mean_fr;
dynamic_range_x = (max(frFiltx) - min(frFiltx))/mean_fr_x;
dynamic_range_y = (max(frFilty) - min(frFilty))/mean_fr_y;
speed_bins = hist(speedFilt);
% compute stability
% divide session into 4 and correlate
split_session = strsplit(posfile,'\');
animal = split_session{9};
sesh = split_session{10};
sesh = strsplit(sesh,'.');
sesh = sesh{1};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PLOTTING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if makeheatplot
array = [speedFilt,frFilt]; % plot occupancy heatmap for each cell
scatter(speedFilt,frFilt)
map = hist3(array,[10 10]);
map1 = map';
map1(size(map1,1) + 1, size(map1,2) + 1) = 0;
xb = linspace(min(array(:,1)),max(array(:,1)),size(map,1)+1);
yb = linspace(min(array(:,2)),max(array(:,2)),size(map,1)+1);
map1(map1 == Inf) = 0;
mapmax = max(map1);
    map1 = map1./mapmax
    map1 = log(map1)
h = pcolor(xb,yb,map1);
ylabel('Firing Rate (Hz)');
xlabel('Speed (cm/s)')
h.LineWidth = 1.5;
axisscale = min(map1(map1 > 0))./max(map1);
r = refline(slope,intercept);
    r.Color = 'b'; r.LineStyle = '--'; r.LineWidth = 5;
formatspec = ("%s %s Score: %0.2f  Slope: %0.2f");
title(sprintf(formatspec,animal,sesh,speedScore,slope));
t = sprintf('%s','%s',animal,sesh);
saveas(h,t,'svg');
end

if makeplot
    figure()
    edges = 0:5:max(speedFilt);
    edges = edges(1:(length(edges)-1));
    [~,bin] = histc(speedFilt,0:5:max(speedFilt));
    frMeans = nan(max(bin),1);
    frErrs = nan(max(bin),1);
    for i = 1:max(bin)
        frMeans(i) = mean(frFilt(bin==i));
        frErrs(i) = std(frFilt(bin==i))./sqrt(max(bin));
    end
    eb = shadedErrorBar(edges,frMeans,frErrs,'lineProps','-b','patchSaturation',0.15)
    axb = gca; axb.FontName = 'Helvetica'; axb.FontSize = 14;
    xlabel('Running speed (cm/s)');
    ylabel('Mean Firing Rate (Hz)');
    formatspec = ("%s %s Score: %0.2f  Slope: %0.2f");
    title(sprintf(formatspec,animal,sesh,speedScore,slope));
    t = sprintf('%s',animal,sesh);
    saveas(gcf,t,'svg');
end

if makeplotx
    figure()
    edgesx = 0:5:max(speedxFilt);
    edgesx = edgesx(1:(length(edgesx)-1));
    [~,binx] = histc(speedxFilt,0:5:max(speedxFilt));
    frMeansx = nan(max(binx),1);
    frErrsx = nan(max(binx),1);
    for i = 1:max(binx)
        frMeansx(i) = mean(frFiltx(binx==i));
        frErrsx(i) = std(frFiltx(binx==i))./sqrt(max(binx));
    end
    ebx = shadedErrorBar(edgesx,frMeansx,frErrsx,'lineProps','-b','patchSaturation',0.15)
    xlabel('Running speed (cm/s)');
    ylabel('Mean Firing Rate (Hz)');
    formatspec = ("%s %s Score: %0.2f  Slope: %0.2f");
    title(sprintf(formatspec,animal,sesh,speedScore,slope));
    t = sprintf('%s %s X direction',animal,sesh);
    saveas(gcf,t,'svg');
end

if makeploty
    figure()
    edgesy = 0:5:max(speedyFilt);
    edgesy = edgesy(1:(length(edgesy)-1));
    [~,bin] = histc(speedyFilt,0:5:max(speedyFilt));
    frMeansy = nan(max(bin),1);
    frErrsy = nan(max(bin),1);
    for i = 1:max(bin)
        frMeansy(i) = mean(frFilty(bin==i));
        frErrsy(i) = std(frFilty(bin==i))./sqrt(max(bin));
    end
    eby = shadedErrorBar(edgesy,frMeansy,frErrsy,'lineProps','-b','patchSaturation',0.15)
    xlabel('Running speed (cm/s)');
    ylabel('Mean Firing Rate (Hz)');
    formatspec = ("%s %s Score: %0.2f  Slope: %0.2f");
    title(sprintf(formatspec,animal,sesh,speedScore,slope));
    t = sprintf('%s %s Y Direction',animal,sesh);
    saveas(gcf,t,'svg');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fr_speed_maps = nan(4,numel(speedVec)-1);
for n = 1:4
    speed_n = speed(sessionInd(n):sessionInd(n+1));
    fr_n = fr(sessionInd(n):sessionInd(n+1));
    for j = 1:numel(speedVec)-1
        start = speedVec(j); stop = speedVec(j+1);
        fr_speed_maps(n,j) = nanmean(fr_n((speed_n > start & speed_n < stop)));
    end
    % take away any nan's
    [~,badCol] = find(isnan(fr_speed_maps));
    fr_speed_maps = fr_speed_maps(:,setdiff(1:numel(speedVec)-1,unique(badCol)));  
end
% compute rmse for each comparison
session1 = [1 1 1 2 2 3]; session2 = [2 3 4 3 4 4];
correlations = nan(6,1);
for m = 1:6
    correlations(m) = corr(fr_speed_maps(session1(m),:)',fr_speed_maps(session2(m),:)');
end
speedStability  = nanmean(correlations);

%%%%%%%%%%%%%%%% X AXIS STABILITY %%%%%%%%%%%%%%%%%%
fr_speed_maps_x = nan(4,numel(speedVec)-1);
for n = 1:4
    speed_x = speedx(sessionInd(n):sessionInd(n+1));
    fr_x = fr(sessionInd(n):sessionInd(n+1));
    for j = 1:numel(speedVec)-1
        start = speedVec(j); stop = speedVec(j+1);
        fr_speed_maps_x(n,j) = nanmean(fr_x((speed_x > start & speed_x < stop)));
    end
    % take away any nan's
    [~,badCol_x] = find(isnan(fr_speed_maps_x));
    fr_speed_maps_x = fr_speed_maps_x(:,setdiff(1:numel(speedVec)-1,unique(badCol_x)));  
end
% compute rmse for each comparison
session1 = [1 1 1 2 2 3]; session2 = [2 3 4 3 4 4];
correlations_x = nan(6,1);
for m = 1:6
    correlations_x(m) = corr(fr_speed_maps_x(session1(m),:)',fr_speed_maps_x(session2(m),:)');
end
speedStability_x  = nanmean(correlations_x);


%%%%%%%%%%%%% Y AXIS STABILITY %%%%%%%%%%%%%%%%%%%%%
fr_speed_maps_y = nan(4,numel(speedVec)-1);
for n = 1:4
    speed_y = speedy(sessionInd(n):sessionInd(n+1));
    fr_y = fr(sessionInd(n):sessionInd(n+1));
    for j = 1:numel(speedVec)-1
        start = speedVec(j); stop = speedVec(j+1);
        fr_speed_maps_y(n,j) = nanmean(fr_y((speed_y > start & speed_y < stop)));
    end
    % take away any nan's
    [~,badColy] = find(isnan(fr_speed_maps_y));
    fr_speed_maps_y = fr_speed_maps_y(:,setdiff(1:numel(speedVec)-1,unique(badColy)));  
end
% compute rmse for each comparison
session1 = [1 1 1 2 2 3]; session2 = [2 3 4 3 4 4];
correlations_y = nan(6,1);
for m = 1:6
    correlations_y(m) = corr(fr_speed_maps_y(session1(m),:)',fr_speed_maps_y(session2(m),:)');
end
speedStability_y  = nanmean(correlations_y);
nonNegativeFiringRate = 1;
spikeCounts = hist(cellTS,post);
[lin_slope,lin_intercept,exp_slope,exp_intercept] = fitting(speed,fr);
[lin_slope_x,lin_intercept_x,exp_slope_x,exp_intercept_x] = fitting(speed_x,fr_x);
[lin_slope_y,lin_intercept_y,exp_slope_y,exp_intercept_y] = fitting(speed_y,fr_y);
end

function [lin_slope,lin_intercept,exp_slope,exp_intercept] = fitting(speed,fr)
quadfit = polyfit(speed,fr,2);
linfit = polyfit(speed,fr,1);
exp_slope = quadfit(2);
exp_intercept = quadfit(3);
lin_slope = linfit(1);
lin_intercept = linfit(2);
end
