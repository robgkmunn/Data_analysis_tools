function [rho_max, yShift, xShift, lambda] = squish_factor_translation_RM_v3(file1, file2, file3, file4,boxSize)

% file1 = Spike file, open field
% file2 = Position file, open field
% file3 = Spike file, squish
% file4 = Position file, squish

rotate = 0; %start out assuming that no rotation is needed

[map_squish, rotate] = map_maker_squish(file3,file4,1,boxSize,rotate);
%run map maker on the squish session first to determine if a rotation
%occured. if so, will need to also rotate the baseline map!
[map_no_squish, rotate] = map_maker_squish(file1,file2,0,boxSize,rotate);

% find maps of squish and non-squish sessions
%map_no_squish = map_maker(file1,file2,0);
%map_squish = map_maker(file3,file4,1);

% find dimensions of maps
[x1,y1] = size(map_squish);
[x2,y2] = size(map_no_squish);

% check to make sure number of rows is the same (sanity check)
if x1 ~= x2
 map_squish = map_squish(1:x2,1:y1);
end

x_maxShift = 10;
y_maxShift = 10;

%if y is negative, the squish map has shifted "down" relative to the
%original.
%if y is positive, the squish map has shifted "up" relative to the
%original.
% if x is negative, the squish map has shifted "left" relative to the
% original
%if x is positive, the squish map has shifted "right" relative to the
%original

% rho = zeros(2*x_maxShift+1,2*y_maxShift+1,y2+10-y1+1); %initialize a 2D matrix (number of y-shifts, # of stretches (difference in columns between maps))
%this will need to be ammended in order to take into account x-shifts as
%well.

rho = zeros(2*x_maxShift+1,2*y_maxShift+1,(y2-y1+1));


for x = -x_maxShift:1:x_maxShift
    %start by shifting the squish map left or right relative to the
    %original map.
    
    for y = -y_maxShift:1:y_maxShift
        % also shift the map in the up or down direction relative to the
        % original map
        
        %for s = y1:y2
        for s = y1:y2 % let the stretch go an extra 20 units
            
            map_resize = imresize(map_squish, [size(map_squish,1), s]); %resize squished map to be same number of rows as original, and new number of columns
            if x >= 0 && y >= 0
                map_squish_to_compare = map_resize(1:size(map_resize,1)-y, 1:min(size(map_no_squish,2)-x,size(map_resize,2)));
                small_map_no_squish = map_no_squish(y+1:end,x+1:min(x + size(map_resize,2),size(map_no_squish,2)));
                
            elseif x <= 0 && y >= 0
                map_squish_to_compare = map_resize(1:size(map_resize,1)-y,abs(x)+1:size(map_resize,2));
                small_map_no_squish = map_no_squish(y+1:end,1:min(size(map_resize,2)- abs(x),size(map_no_squish,2)));
                
            elseif x >= 0 && y <= 0
                
                map_squish_to_compare = map_resize(abs(y)+1:end,1:min(size(map_no_squish,2) - x,size(map_resize,2)));
                small_map_no_squish = map_no_squish(1:end-abs(y),x+1:min(size(map_no_squish,2),x+size(map_resize,2)));
                
            elseif x <= 0 && y <= 0
                
                map_squish_to_compare = map_resize(abs(y)+1:end,abs(x)+1:end);
                small_map_no_squish = map_no_squish(1:size(map_no_squish,1)-abs(y),1:size(map_resize,2)-abs(x));
                
            end
           
            %calculate correlation coefficient
            X = corr2(map_squish_to_compare,small_map_no_squish);
            index1 = y + y_maxShift + 1;
            index2 = x + x_maxShift + 1;
            index3 = s-y1+1;
            rho(index1,index2,index3) = X;
        end
    end
end
xBinWidth = 2;
yBinWidth  = 2;
[rho_max, linearIndexesOfMaxes] = max(rho(:));
[yShiftIndex,xShiftIndex,stretchIndex] = ind2sub(size(rho),linearIndexesOfMaxes);
lambda = (stretchIndex-1)/(y2-y1); %normalized squish factor
stretch_units = stretchIndex - 1;

yShift = -y_maxShift + yShiftIndex - 1; %in bins
xShift = -x_maxShift + xShiftIndex - 1; %in bins
xShift = xBinWidth*xShift;
yShift = yBinWidth*yShift; %in cm
shift = sqrt(xShift^2 + yShift^2); %in cm

disp(lambda);
disp(yShift);
disp(xShift);
disp(rho_max);


map_resize = imresize(map_squish, [size(map_squish,1), size(map_squish,2) + stretch_units]);
x = xShift;
y = yShift;
    if x >= 0 && y >= 0
                map_squish_to_compare = map_resize(1:size(map_resize,1)-y, 1:min(size(map_no_squish,2)-x,size(map_resize,2)));
                small_map_no_squish = map_no_squish(y+1:end,x+1:min(x + size(map_resize,2),size(map_no_squish,2)));
                small_map_no_squish_rowStart = y+1;
                small_map_no_squish_rowStop = size(map_no_squish,1);
                small_map_no_squish_colStart = x + 1;
                small_map_no_squish_colStop = min(x + size(map_resize,2),size(map_no_squish,2));
                
            elseif x <= 0 && y >= 0
                map_squish_to_compare = map_resize(1:size(map_resize,1)-y,abs(x)+1:size(map_resize,2));  
                small_map_no_squish = map_no_squish(y+1:end,1:min(size(map_resize,2)- abs(x),size(map_no_squish,2)));
                small_map_no_squish_rowStart = y+1;
                small_map_no_squish_rowStop = size(map_no_squish,1);
                small_map_no_squish_colStart = 1;
                small_map_no_squish_colStop = min(size(map_resize,2)- abs(x),size(map_no_squish,2));
                
            elseif x >= 0 && y <= 0
                map_squish_to_compare = map_resize(abs(y)+1:end,1:min(size(map_no_squish,2) - x,size(map_resize,2)));
                small_map_no_squish = map_no_squish(1:end-abs(y),x+1:min(size(map_no_squish,2),x+size(map_resize,2)));
                small_map_no_squish_rowStart = 1;
                small_map_no_squish_rowStop = size(map_no_squish,1) - abs(y);
                small_map_no_squish_colStart = x+1;
                small_map_no_squish_colStop = min(size(map_no_squish,2),x+size(map_resize,2));
            elseif x <= 0 && y <= 0
                map_squish_to_compare = map_resize(abs(y)+1:end,abs(x)+1:end);
                small_map_no_squish = map_no_squish(1:size(map_no_squish,1)-abs(y),1:size(map_resize,2)-abs(x));
                small_map_no_squish_rowStart = 1;
                small_map_no_squish_rowStop = size(map_no_squish,1)-abs(y);
                small_map_no_squish_colStart = 1;
                small_map_no_squish_colStop = size(map_resize,2)-abs(x);
    end
    
fig1 = figure(1);
set(fig1,'Position',[0,0,800,800]);
ax1 = subplot(2,2,1);
imagesc(map_no_squish)
%recall which portion of the original map is being compared and draw a box around it.
hold on
col = [small_map_no_squish_colStart; small_map_no_squish_colStop; small_map_no_squish_colStop; small_map_no_squish_colStart; small_map_no_squish_colStart];
row = [small_map_no_squish_rowStart; small_map_no_squish_rowStart; small_map_no_squish_rowStop; small_map_no_squish_rowStop; small_map_no_squish_rowStart];
plot(col,row,'-k','lineWidth',3)
title('Original')
set(gca,'fontsize',20)
hold on
ax2 = subplot(2,2,2);
imagesc(map_squish)
title('Squish')
set(gca,'fontsize',20)    
ax3 = subplot(2,2,3);
imagesc(map_squish_to_compare)
title('Squish + Shift Map')
set(gca,'fontsize',20)
axis([0 50 0 50]);
linkaxes([ax1,ax2,ax3],'xy')
title(['l = ',num2str(lambda, '%0.2f'), ' r = ', num2str(rho_max, '%0.2f'), '  xShift = ', num2str(x), '  yShift = ', num2str(y)])
name = strsplit(file1,'\');
try
name = [name{1,end-2},name{1,end-1},name{1,end}];
catch
    name = file1;
end
fileName = [name,'_gridsquish_and_10xy_shift.jpg'];
saveas(gcf,fileName,'jpg');
end