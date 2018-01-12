function [map,rotate] = map_maker_squish(file1, file2, squish, boxSize, rotate)

global sampleTime minBinTime minX minY numBinsX numBinsY xBinWidth yBinWidth 

global smoothingMode xLength yLength

%Loading test file
load(file1);
load(file2);

posy = -posy; %dbmaker is mirror image of tint
posy2 = -posy2;


if squish == 1 %this is a squish session
    
    % check which side is squished
    minX = nanmin(posx);
    maxX = nanmax(posx);
    minY = nanmin(posy);
    maxY = nanmax(posy);
    xLength = maxX - minX;
    yLength = maxY - minY;
    
    % if squish occurred along x-axis
    if yLength > xLength + 20 % add 20cm to the squish axis to ensure it's correct
       rotate_squish = 0;   
    % if squish occurred along y-axis    
    elseif xLength > yLength + 20 % 
        rotate_squish = 1;
    else
        disp('this is not a squish session')
        return
    end
if rotate_squish == 1
    posx_new = -posy; posx2_new = -posy2; posy = posx; posy2 = posx2;
    posx = posx_new; posx2 = posx2_new;
    disp("rotated the squish map")
else
end
end

if rotate == 1
    posx_new = -posy; posx2_new = -posy2; posy = posx; posy2 = posx2;
    posx = posx_new; posx2 = posx2_new;
    disp("rotated the baseline map, because the squish was rotated")
else
end

%scales box size and bin size - ASSUMES SHAPE IS A BOX
%max and min x and y coordinates from data
minX = nanmin(posx);
maxX = nanmax(posx);
minY = nanmin(posy);
maxY = nanmax(posy);

%scaled box size
xLength = maxX - minX;
yLength = maxY - minY;
scaleFactor = boxSize/(max(xLength,yLength));
posx = posx * scaleFactor;
posy = posy * scaleFactor;

minX = nanmin(posx);
maxX = nanmax(posx);
minY = nanmin(posy);
maxY = nanmax(posy);

xLength = maxX - minX;
yLength = maxY - minY;

%why are we doing this?
% xDiff = abs(minX)-x_size/2;
% posx = posx+xDiff*ones(size(posx));
% 
% yDiff = abs(minY)-y_size/2;
% posy = posy+yDiff*ones(size(posy));


% Number of bins in each direction of the map
binSize = 2; %cm

%bin width ( will vary if the scale is changed) - increase bin width if there
%are a lot of white squares in the figures
xBinWidth = binSize;
yBinWidth = binSize;

numBinsX = ceil(xLength/xBinWidth);
numBinsY = ceil(yLength/yBinWidth);

%Aligning spike data with position data
post22 = post*10;
TS2 = cellTS*10;
post2 = round(post22);
TS = round(TS2);
inds = find(ismember(post2, TS));
spkx = posx(inds);
spky = posy(inds);

sampleTime = .02;
minBinTime = 0.020;
smoothingMode= 0;

[map,~] = RateMapAdaptSmooth(spkx,spky,posx,posy);

%flip the map along the y-axis so that the images from matlab match the
%images from tint!
%map = flipud(map);

%Not currently being taken into account: if the wall was moved from the West inward, the images need to
%both be rotated by 180 degrees in order for the "anchoring" to still be on
% the left.
if squish ==1
if rotate_squish == 1
    rotate = 1;
else
    rotate = 0
end
end
%if rotate_map == 1
%    map = rot90(map,1);
%elseif rotate_map == 0
%    map = rot180(map,1);
%end

return