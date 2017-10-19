% This program grabs the relevant EEG files from Axona recordings and then computes speed filtered theta frequency and amplitude overall
% and separately for N, S, E & W directions
% by Robert Munn, Stanford University 29 August 2017

function directional_filtering_1_2(spreadsheet,sheet)

[~,~,data] = xlsread(spreadsheet,sheet);
Session = data(2:end,find(strcmp(data(1,:),'Session(s)') == 1));

%Tetrode = data(2:end,find(strcmp(data(1,:),'Tetrode') == 1));
%Unit = data(2:end,find(strcmp(data(1,:),'Unit') == 1));
Box = data(2:end,find(strcmp(data(1,:),'Box') == 1));
%stable_wall = data(2:end,find(strcmp(data(1,:),'stable_wall') == 1));
squish = data(2:end,find(strcmp(data(1,:), 'Squish') == 1));
numCell = length(Box);
Theta_x = zeros(120050,numCell);
Theta_y = zeros(120050,numCell);
boxSize = 100;
% find spike and pos files
for n = 1:numCell
    
    directory = Session{n};
    posfile = strcat(directory,'_pos.mat');
   % spikefile = strcat(directory,'_T',num2str(Tetrode{n}),'C',num2str(Unit{n}),'.mat');
    eegfile = strcat(directory,'_eeg.mat');

  if ~exist(posfile)
      disp([posfile,' does not exist']);
  end
  
 % if ~exist(spikefile)
 %     disp([spikefile,' does not exist']);
%  end
  
  if ~exist(eegfile)
      disp(['yo eeg file','_',eegfile,'_','be missing, yo']);
  end

pos = load(posfile); % load the position 
%spike = load(spikefile); % load the spike times
eeg = load(eegfile); % load the eeg file

posx = pos.posx;
posy = pos.posy;
posy2 = pos.posy2;
posx2 = pos.posx2;
post = pos.post;
%cellTS = spike.cellTS;
recSystem = pos.recSystem;
eegdata = eeg.eegData;
eeg = eegdata{1,1};
sampling = eegdata{1,2};
eeg_filetype = eegdata{1,3};

% Bin width for the speed in theta anlysis. Theta frequency for each speed 
% bin is calculated.
p.binWidthSpeed = 2; % [cm/s]

% Set the maximum speed bin value for the speed analysis. Threshold is used
% both for the theta frequency and the intrinsic frequency analysis.
p.maximumSpeedForPlot = 50; % [cm/s]

% Low speed threshold. Segments of the path where the rat moves slower
% than this threshold will be removed. Set it to zero (0) to keep
% everything. Value in centimeters per second.
p.lowSpeedThreshold = 2.5; % [cm/s]

% High speed threshold. Segments of the path where the rat moves faster
% than this threshold will be removed. Set it to zero (0) to keep
% everything. Value in centimeters per second.
p.highSpeedThreshold = 50; % [cm/s]

% Bandpass limits for the theta band filter.
% The passbands (Fp1 Fp2) and stop bands (Fs1 Fs2) are defined as
%                       ------                   
%                     /           \
%                    /             \
%                   /               \
%                  /                 \
%   ----------                  ----------------- 
%           Fs1  Fp1       Fp2  Fs2
%
% p.thetaBand = [Fs1, Fp1, Fp2, Fs2];
p.thetaBand = [2, 4, 12, 14];

% Bin width for the autocorrelogram (ISI)
p.binWidth = 0.005; % [Sec]

% Length of the autocorrelogram. 
p.diagramLength = 0.500; % [Sec]

% Offset used when calculating the linear acceleration. Number sets the
% number of samples to look back and ahead when calculating the
% acceleration.
p.accelrationSampleOffset = 3;

posy = -posy;
posy2 = -posy2; % flip db output to be not a mirror image

for f = 1:numCell
if squish{f,1} == 1 % this means that it was a squish (assume squish from right - west wall is stable)
    
    minX = nanmin(posx); maxX = nanmax(posx);
    minY = nanmin(posy); maxY = nanmax(posy);
    xLength = maxX - minX; yLength = maxY - minY; sLength = max([xLength, yLength]);
    scale = boxSize / sLength;
    % x-direction is half of y-direction
    posx = posx * scale; posy = posy * scale;
    posx2 = posx2 * scale; posy2 = posy2 * scale;
    
else
    minX = nanmin(posx); maxX = nanmax(posx);
    minY = nanmin(posy); maxY = nanmax(posy);
    xLength = maxX - minX; yLength = maxY - minY; sLength = max([xLength, yLength]);
    scale = boxSize / sLength;
    posx = posx * scale; posy = posy * scale;
    posx2 = posx2 * scale; posy2 = posy2 * scale;
end
end

% change x to y if the squish seems to be coming from the y direction - squish is always x
if yLength < xLength 
    posx_new = posy; posx2_new = posy2; posy = posx; posy2 = posx2;
    posx = posx_new; posx2 = posx2_new;
end
    
% compute head direction sampling
hdDir = atan2(posy2-posy,posx2-posx)+pi/2;
hdDir(hdDir < 0) = hdDir(hdDir<0)+2*pi;

%%%%%% Speed filter %%%%%%

% Low speed threshold. Segments of the path where the mouse moves slower
% than this threshold will be removed.
p.lowSpeedThreshold = 2.5; % [cm/s]

% High speed threshold. Segments of the path where the mouse moves faster
% than this threshold will be removed.
p.highSpeedThreshold = 50; % [cm/s]

if p.lowSpeedThreshold > 0 || p.highSpeedThreshold > 0
% Calculate the speed of the mouse, sample by sample
    speed = speed2D(posx,posy,post);
    
    if p.lowSpeedThreshold > 0 && p.highSpeedThreshold > 0
        bad_ind = find(speed < p.lowSpeedThreshold | speed > p.highSpeedThreshold);
    elseif p.lowSpeedThreshold > 0 && p.highSpeedThreshold == 0
        bad_ind = find(speed < p.lowSpeedThreshold );
    else
        bad_ind = find(speed > p.highSpeedThreshold );
    end
 
% Remove the segments that have too high or too low speed
    posx(bad_ind) = NaN;
    posy(bad_ind) = NaN;
    
    posx2(bad_ind) = NaN;
    posy2(bad_ind) = NaN;    
  speed(speed < 2) = NaN; %actually take out those speed variables
  speed(speed > 50) = NaN;
% find times when trajectories are consistently X
%or Y by considering the change in head direction. Reject
%all those positions that are not X or Y. the windows are  X (N pi/6 - 11pi/6,
%S 7pi/6 - 5pi/6) Y (E, pi/3 - 2pi/3; W, 5pi/3 - 4pi/3)
post_south = post;
post_north = post;
post_east = post;
post_west = post;

bad_dir_south = find(hdDir > 7*pi/6 | hdDir < 5*pi/6);
post_south(bad_dir_south) = NaN;
bad_dir_north = find(hdDir > 11*pi/6 | hdDir < pi/6);
post_north(bad_dir_north) = NaN;
bad_dir_east = find(hdDir < pi/3 | hdDir > 2*pi/3);
post_east(bad_dir_east) = NaN;
bad_dir_west = find(hdDir < 4*pi/3 | hdDir > 5*pi/3);
post_west(bad_dir_west) = NaN;

%one dimensional speed
speed_x = speed1D(posx,post);
speed_x(speed_x < 2.5) = NaN;
speed_x(speed_x > 50) = NaN;
speed_y = speed1D(posy,post);
speed_x(speed_x < 2.5) = NaN;
speed_x(speed_x > 50) = NaN;
meanspeed_x = nanmean(speed_x);
meanspeed_y = nanmean(speed_y);
meanspeed_both = nanmean(speed);
decimated_eeg = decimate(eeg,5); % decmimate EEG to make it the same length as pos vars and speed
%bandpass fast fourier transform the eeg signal

theta = fftbandpass(eeg,sampling,2,4,12,14);

%indexes where the peaks of theta are
  [~, peakInd] = thetaPhase(theta);
    
% Calculate instantaneous frequency and amplitude
 thetaProperties = thetaFrequency(theta, peakInd, sampling);
% thetaProperties = thetaFiltering(thetaProperties);

   % Calculate the theta frequency and amplitude for each position sample
  [posTheta,posAmplitude]  = thetaFrequencyPosSample(post,thetaProperties);
 [posTheta_s,posAmplitude_s]  = thetaFrequencyPosSample(post_south,thetaProperties);
 [posTheta_n,posAmplitude_n]  = thetaFrequencyPosSample(post_north,thetaProperties);
 [posTheta_e,posAmplitude_e]  = thetaFrequencyPosSample(post_east,thetaProperties);
 [posTheta_w,posAmplitude_w]  = thetaFrequencyPosSample(post_west,thetaProperties);

%find the NaN positions in one array (i.e W) and replace with corresponding values from the opposite direction (i.e E) making two directional vectors (X and Y)
  x = isnan(posTheta_n);                 
  posTheta_n(x) = posTheta_s(x);
  y = isnan(posTheta_w);
  posTheta_w(y) = posTheta_e(y);
  
  xamp = isnan(posAmplitude_n);
  posAmplitude_n(xamp) = posAmplitude_s(xamp);
  yamp = isnan(posAmplitude_w);
  posAmplitude_w(yamp) = posAmplitude_e(yamp);
  
  Theta_x = posTheta_n;
  Theta_y = posTheta_w;
  
  Amplitude_x = posAmplitude_n;
  Amplitude_y = posAmplitude_w;
  
  % make running speed/theta frequency bins for each of the X and Y axes
  % and perform linear regression on them
 [speedAxis_x, thetaFreqBin_x, thetaAmpBin_x] =  speedThetaBinning(speed_x, Theta_x, Amplitude_x, p);
 [speedAxis_y, thetaFreqBin_y, thetaAmpBin_y] =  speedThetaBinning(speed_y, Theta_y, Amplitude_y, p);
 [R2_freq_x, slope_freq_x, intercept_freq_x] = regressionAnalysis(speedAxis_x, thetaFreqBin_x);
 R_freq_x = sqrt(R2_freq_x);
 [R2_freq_y, slope_freq_y, intercept_freq_y] = regressionAnalysis(speedAxis_y, thetaFreqBin_y);
 R_freq_y = sqrt(R2_freq_y);
 figure();
 posTheta(posTheta > 12) = NaN; % take out instantaneous theta frequencies greater than 12Hz
scatter(speed,posTheta);
keyboard
 % store the values in an array
 X_slope_frequency(1,n) = slope_freq_x;
 Y_slope_frequency(1,n) = slope_freq_y;
 X_intercept_frequency(1,n) = intercept_freq_x;
 Y_intercept_frequency(1,n) = intercept_freq_y;
 R_val_x_freq(1,n) = R_freq_x;
 R_val_y_freq(1,n) = R_freq_y;
%array_test = isnan(Theta_x);

%for o = 1:length(Theta_x) %
   % index_first = find(array_test,1,'first');
    %find(array_test,1)
%end
 [theta_x_consec] = findConseqSeq(Theta_x,50);
 [theta_y_consec] = findConseqSeq(Theta_y,50);



 Theta_x(isnan(theta_x_consec))=[]; %take out all the NaNs
 Theta_y(isnan(theta_y_consec))=[];
 

 Theta_x_dir{1,n} = Theta_x;
 Theta_y_dir{1,n} = Theta_y;

  X_theta_freq(1,n) = mean(Theta_x);
  Y_theta_freq(1,n) = mean(Theta_y);
 
 X_theta_freq_stdev(1,n) = std(Theta_x);
 Y_theta_freq_stdev(1,n) = std(Theta_y);
 
% store all that junk in a cell array
 Splittar = strsplit(directory,'\');
 Animalname = Splittar{1,6};
 sesh = Splittar{1,end};
 session_st = strcat(Animalname,'_',sesh);
 theta_array(1,1:14) = {'Animal','theta frequency in the x axis','theta frequency in the y axis','standard deviation theta x','standard deviation theta y','Frequency slope X','Frequency slope Y','Frequency intercept X','Frequency intercept Y','R value of the fit (x)','R value of the fit (y)','Speed (X)','Speed (Y)', 'Overall Speed'};
 theta_array{n+1,1} = (session_st);
 theta_array{n+1,2} = X_theta_freq(1,n);
 theta_array{n+1,3} = Y_theta_freq(1,n);
 theta_array{n+1,4} = X_theta_freq_stdev(1,n);
 theta_array{n+1,5} = Y_theta_freq_stdev(1,n);
 theta_array{n+1,6} = X_slope_frequency(1,n);
 theta_array{n+1,7} = Y_slope_frequency(1,n);
 theta_array{n+1,8} = X_intercept_frequency(1,n);
 theta_array{n+1,9} = Y_intercept_frequency(1,n);
 theta_array{n+1,10} = R_val_x_freq(1,n);
 theta_array{n+1,11} = R_val_y_freq(1,n);
 theta_array{n+1,12} = meanspeed_x;
 theta_array{n+1,13} = meanspeed_y;
 theta_array{n+1,14} = meanspeed_both;
end
end

% get overall standard deviation of the frequency estimates
 for z = 1:numCell
     N_x(z) = length(Theta_x_dir{1,z});
     SD_x(z) = X_theta_freq_stdev(1,z).^2;
     nominator_x = SD_x(z)*N_x(z);
     denominator_x = sum(N_x(z));
      N_y(z) = length(Theta_y_dir{1,z});
     SD_y(z) = Y_theta_freq_stdev(1,z).^2;
     nominator_y = SD_y(z)*N_y(z);
     denominator_y = sum(N_y(z));
 end
 
grand_stdev_x = (nominator_x/denominator_x).^0.5;
grand_stdev_y = (nominator_y/denominator_y).^0.5;
 
 % store the overall standard deviations in the cell array
 theta_array{1,12} = 'overall x deviation';
 theta_array{2,12} = grand_stdev_x;
 theta_array{1,13} = 'overall y deviation';
 theta_array{2,13} = grand_stdev_y;
 
 
 %save theta array as a speadsheet and the entire workspace as a .mat file
fileName = [spreadsheet,'_theta_analysis','.xls'];
xlswrite(fileName,theta_array);
matName = [spreadsheet,'.mat'];
save(matName);
end


