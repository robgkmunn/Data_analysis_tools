% speedModulatedThetaFrequency1_0('inputFile')
%
% This program calculates the relationship between EEG theta frequency and
% running speed. The method is described in the paper: "Grid Cells and
% Theta as Oscillatory Interference: Electrophysiological Data From Freely
% Moving Rats" by A. Jeewajee, C. Barry, J. O'Keefe, and N. Burgess. 
%
% But note that this program uses a different filter than the one described
% in the paper. The filter we use is an acuasal fft algorithm (i.e. no 
% phase shift) made by Ole Jensen, Brain Resarch Unit, Low Temperature 
% Laboratory, Helsinki University of Technology, 02015 HUT, Finland.
%
% This program was made by request from Lisa Giocomo.
% 
%  INPUT ARGUMENTS
%
% inputFile     Text file with information about the data to be analysed.
%               This must be an input file that you get from databaseMaker.
%               Look into that program for details. The name for the file
%               must be written in single quotes. Example: 'in1.txt'
%               The lines that must be included in the input file for each 
%               session are:
%               Session
%               Shape
%               Tetrode
%               Unit
%
%
% Version 1.0
% 31. Mar. 2011
%
% Version 1.1       Added regression line intercept value to the output.
% 06. Apr. 2011
%
% Version 1.2       Added acceleration modulated regression analysis of the
% 29. Apr. 2011     speed - theta frequency and speed - theta amplitude.
%                   Added calculation of the intrinsic frequency for the
%                   whole session (not speed modulated)
%
% Version 1.3       Added figure of power spectrum from the
% 04. May. 2011     autocorrelogram.
%                   Change the name on the figure files.
%
% Version 1.4       Added calculation of ISI.
% 05. May. 2011     Added cell name to the figures based on cell data.
%
% Version 1.5       Added option to set the amount of smoothing of the ISI
% 06. May. 2011     diagram.
%                   
%
% Created by Raymond Skjerpeng, KI/CBM, NTNU, 2011.
function speedModulatedThetaFrequency1_5(inputFile)



%__________________________________________________________________________
%
%                       Program parameters
%__________________________________________________________________________

% Bin width for the speed in theta anlysis. Theta frequency for each speed 
% bin is calculated.
p.binWidthSpeed = 5; % [cm/s]

% Bin width for the speed in the intrinsic frequency
p.intrinsicSpeedBinWidth = 10; % [cm/s]

% Minimum number of spikes in a speed bin. If the number of spikes is lower 
% than this the intrinsic frequency is set to NaN;
p.minimumNumSpikes = 100;

% Set the maximum speed bin value for the speed analysis. Threshold is used
% both for the theta frequency and the intrinsic frequency analysis.
p.maximumSpeedForPlot = 50; % [cm/s]

% Format for the images
% format = 1 -> bmp (24 bit)
% format = 2 -> png
% format = 3 -> eps
% format = 4 -> jpg
% format = 5 -> ai (Adobe Illustrator)
% format = 6 -> tiff (24 bit)
% format = 7 -> fig (Matlab figure)
p.imageFormat = 2;

% Low speed threshold. Segments of the path where the rat moves slower
% than this threshold will be removed. Set it to zero (0) to keep
% everything. Value in centimeters per second.
p.lowSpeedThreshold = 2.5; % [cm/s]

% High speed threshold. Segments of the path where the rat moves faster
% than this threshold will be removed. Set it to zero (0) to keep
% everything. Value in centimeters per second.
p.highSpeedThreshold = 100; % [cm/s]

% Folder to store the images created by this program. The folder will be
% created as a subfolder in the same directory as the input data is stored
% in.
p.imageFolder = 'thetaImages';

% Bandpass limits for the theta band filter.
% The passbands (Fp1 Fp2) and stop bands (Fs1 Fs2) are defined as
%                  ---------                      
%                /           \
%               /             \
%              /               \
%             /                 \
%   ----------                   ----------------- 
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
p.accelrationSampleOffset = 1;

% Set the amount of smoothing in the ISI. Must be an interger number. If
% set to zero no smoothing will be applied.
p.smoothingFactor = 3;

%__________________________________________________________________________


fprintf('%s%s\n','Start analysing at ', datestr(now));

% Set the operation system
if ispc
    % Windows
    p.computer = 0;
    % Directory delimiter
    p.delim = '\';
elseif isunix
    % Unix
    p.computer = 1;
    p.delim = '/';
else
    disp('ERROR: Sorry, this program can only be run on windows or Unix')
    return
end

% Check if the output directory for the images is present, if not make it.
dirInfo = dir(p.imageFolder);
if size(dirInfo,1) == 0
    mkdir(p.imageFolder);
end

disp('Reading and checking input file')
[status, sessionArray, unitArray, shapeArray] = inputFileReader(inputFile, p);

if status == 0
    return
end

% Check that the data listed exist
status = inputDataChecker(sessionArray, unitArray);
if status == 0
    return
end

% Output file name
sInd = strfind(inputFile,'.');
if ~isempty(sInd)
    outputFile = strcat('Results_',inputFile(1:sInd(end)),'xls');
else
    outputFile = strcat('Results_',inputFile,'.xls');
end

try
    fid2 = fopen(outputFile,'w');
    if fid2 == -1
        fprintf('%s%s\n','Error: Unable to open the output file: ', outputFile);
        return
    end
catch me
    fprintf('%s%s\n','Error: Unable to open the output file: ', outputFile)
    disp('Make sure it is not open in another program and that you have write permission to the destination folder')
    fprintf('%s%s\n','Matlab error message: ', me)
    return
end



% Write header to the output file
fprintf(fid2,'%s\t','Session');
fprintf(fid2,'%s\t','Cell');
fprintf(fid2,'%s\t','Theta Amplitude R-value');
fprintf(fid2,'%s\t','Theta Amplitude slope (uV/cm/s)');
fprintf(fid2,'%s\t','Theta Amplitude intersept');
fprintf(fid2,'%s\t','Theta Frequency R-value');
fprintf(fid2,'%s\t','Theta Frequency slop (Hz/cm/s)');
fprintf(fid2,'%s\t','Theta Frequency intersept');

fprintf(fid2,'%s\t','Acceleration Theta Amplitude R-value');
fprintf(fid2,'%s\t','Acceleration Theta Amplitude slope (uV/cm/s)');
fprintf(fid2,'%s\t','Acceleration Theta Amplitude intersept');
fprintf(fid2,'%s\t','Acceleration Theta Frequency R-value');
fprintf(fid2,'%s\t','Acceleration Theta Frequency slop (Hz/cm/s)');
fprintf(fid2,'%s\t','Acceleration Theta Frequency intersept');

fprintf(fid2,'%s\t','Deacceleration Theta Amplitude R-value');
fprintf(fid2,'%s\t','Deacceleration Theta Amplitude slope (uV/cm/s)');
fprintf(fid2,'%s\t','Deacceleration Theta Amplitude intersept');
fprintf(fid2,'%s\t','Deacceleration Theta Frequency R-value');
fprintf(fid2,'%s\t','Deacceleration Theta Frequency slop (Hz/cm/s)');
fprintf(fid2,'%s\t','Deacceleration Theta Frequency intersept');

numSpeedBins = ceil(p.maximumSpeedForPlot / p.intrinsicSpeedBinWidth);
fprintf(fid2,'%s\t','ISI');

startSpeed = 0;
stopSpeed = p.intrinsicSpeedBinWidth;
for ii = 1:numSpeedBins
    fprintf(fid2,'%s%3.1f%s%3.1f\t','ISI speed ',startSpeed,'-',stopSpeed);
    startSpeed = startSpeed + p.intrinsicSpeedBinWidth;
    stopSpeed = stopSpeed + p.intrinsicSpeedBinWidth;
end


fprintf(fid2,'%s\t','Intrinsic frequency');

startSpeed = 0;
stopSpeed = p.intrinsicSpeedBinWidth;
for ii = 1:numSpeedBins
    fprintf(fid2,'%s%3.1f%s%3.1f\t','Intrinsic frequency speed ',startSpeed,'-',stopSpeed);
    startSpeed = startSpeed + p.intrinsicSpeedBinWidth;
    stopSpeed = stopSpeed + p.intrinsicSpeedBinWidth;
end



fprintf(fid2,'\n');

% Number of session listed in the input file
numSessions = size(sessionArray,1);

if numSessions < 1
    disp('Error: No sessions listed')
    return
end


if ~strcmpi(p.imageFolder(end),'\')
    p.imageFolder = strcat(p.imageFolder,'\');
end



% Start loading data
for s = 1:numSessions
    fprintf('%s%s\n','Loading data for session ',sessionArray{s});
    
    sInd = strfind(sessionArray{s},'\');
    if ~isempty(sInd)
        if length(sInd) > 1
            subDir = sessionArray{s}(sInd(end-1)+1:sInd(end));
            dirPath = strcat(p.imageFolder,subDir);
            session = sessionArray{s}(sInd(end)+1:end);
        else
            dirPath = p.imageFolder;
            session = sessionArray{s}(sInd(end)+1:end);
        end
    else
        dirPath = p.imageFolder;
        session = sessionArray{s};
    end
    
    dirInfo = dir(dirPath);
    if size(dirInfo,1) == 0
        mkdir(dirPath);
    end
    
    dirPath = strcat(dirPath,session);
    
    
    
    
    % Load the position data
    posFile = strcat(sessionArray{s},'_pos.mat');
    load(posFile)
    
    % Make sure the correct variables were loaded from the file
    if ~exist('posx','var')
        disp('Error: The position file is missing the posx variable')
        return
    end
    if ~exist('posy','var')
        disp('Error: The position file is missing the posy variable')
        return
    end
    if ~exist('post','var')
        disp('Error: The position file is missing the post variable')
        return
    end
    if ~exist('recSystem','var')
        disp('Error: The position file is missing the recSystem variable')
        return
    end
    
    % Load the EEG data file
    eegFile = strcat(sessionArray{s,1},'_eeg.mat');
    load(eegFile);
    if ~exist('eegData','var')
        disp('Error: The EEG file is missing the eegData variable')
        return
    end
    
    if strcmpi(recSystem,'Axona')
        p.sampleTime = 0.02;
        p.videoSamplingRate = 50;
    else
        p.sampleTime = 0.04;
        p.videoSamplingRate = 25;
    end
    
    % Add the box shape and size to the parameter list
    p.boxType = shapeArray{s};
    
    % Scale the coordinates using the shape information
    minX = nanmin(posx);
    maxX = nanmax(posx);
    minY = nanmin(posy);
    maxY = nanmax(posy);
    xLength = maxX - minX;
    yLength = maxY - minY;
    sLength = max([xLength, yLength]);
    scale = shapeArray{s}(2) / sLength;
    posx = posx * scale;
    posy = posy * scale;
    
    % Calculate the speed of the animal
    speed = speed2D(posx,posy,post);
    
    % Calculate the sign of the acceleration
    accSign = accelrationSign(speed, p.accelrationSampleOffset);
    
    disp('Analysing the EEG data')

    
    % Band pass filter the EEG signal in the theta band from 5 to 12 Hz.
    theta = fftbandpass(eegData{1},eegData{2},3,5,11,12);
    
    %Calculate the phase of the filtered eeg signal and locate the peaks
    [~, peakInd] = thetaPhase(theta);
    
    % Calculate instantaneous frequency and amplitude
    thetaProperties = thetaFrequency(theta, peakInd, eegData{2});
    
    % Remove segments where the theta cycles are not present
    thetaProperties = thetaFiltering(thetaProperties);
    
    % Calculate the theta frequency and amplitude for each position sample
    [posTheta, posAmplitude]  = thetaFrequencyPosSample(post,thetaProperties);
    
    if p.lowSpeedThreshold > 0
        % Remove samples where the speed is to low
        ind = find(speed < p.lowSpeedThreshold);
        
        posTheta(ind) = NaN;
        posAmplitude(ind) = NaN;
    end
    
    if p.highSpeedThreshold > 0
        % Remove samples where the speed is to high
        ind = find(speed > p.highSpeedThreshold);
        
        posTheta(ind) = NaN;
        posAmplitude(ind) = NaN;
    end
    
    % Bin up the theta frequency and amplitude in speed bins
    [speedAxis, thetaFreqBin, thetaAmpBin] =  speedThetaBinning(speed, posTheta, posAmplitude, p);
    
    % Do the regression analysis for the speed vs theta amplitude
    [R2_amplitude, slope_amplitude, intercept_amplitude] = regressionAnalysis(speedAxis, thetaAmpBin);
    R_amplitude = sqrt(R2_amplitude);

    figure(1)
    bar(speedAxis, thetaAmpBin)
    xlabel('Speed [cm/s]');
    ylabel('Amplitude [microvolts]')
    title('Theta amplitude (EEG)')
    fName = strcat(dirPath,'_thetaAmplitude');
    imageStore(figure(1),p.imageFormat,fName,300);

    % Do the regression analysis for the speed vs theta frequency
    [R2_frequency, slope_frequency, intercept_frequency] = regressionAnalysis(speedAxis, thetaFreqBin);
    R_frequency = sqrt(R2_frequency);

    figure(2)
    bar(speedAxis, thetaFreqBin)
    xlabel('Speed [cm/s]');
    ylabel('Frequency [Hz]')
    title('Theta Frequency (EEG)')
    fName = strcat(dirPath,'_thetaFrequency');
    imageStore(figure(2),p.imageFormat,fName,300);
    
    % Find samples where the animal is accelerating
    accInd = find(accSign == 1);
    aSpeed = speed(accInd);
    aPosTheta = posTheta(accInd);
    aPosAmplitude = posAmplitude(accInd);
    
    % Bin up the theta frequency and amplitude in speed bins
    [speedAxis, thetaFreqBin, thetaAmpBin] =  speedThetaBinning(aSpeed, aPosTheta, aPosAmplitude, p);
    % Regression analysis
    [R2_amplitude_a, slope_amplitude_a, intercept_amplitude_a] = regressionAnalysis(speedAxis, thetaAmpBin);
    R_amplitude_a = sqrt(R2_amplitude_a);
    
    [R2_frequency_a, slope_frequency_a, intercept_frequency_a] = regressionAnalysis(speedAxis, thetaFreqBin);
    R_frequency_a = sqrt(R2_frequency_a);
    
    
    figure(1)
    bar(speedAxis, thetaAmpBin)
    xlabel('Speed [cm/s]');
    ylabel('Amplitude [microvolts]')
    title('Theta amplitude (EEG) when accelrating')
    fName = strcat(dirPath,'_thetaAmplitude_acceleration');
    imageStore(figure(1),p.imageFormat,fName,300);
    
    figure(2)
    bar(speedAxis, thetaFreqBin)
    xlabel('Speed [cm/s]');
    ylabel('Frequency [Hz]')
    title('Theta Frequency (EEG) when accelerating')
    fName = strcat(dirPath,'_thetaFrequency_acceleration');
    imageStore(figure(2),p.imageFormat,fName,300);
    
    
    % Find samples where the animal is deaccelration
    accInd = find(accSign == 0);
    
    dSpeed = speed(accInd);
    dPosTheta = posTheta(accInd);
    dPosAmplitude = posAmplitude(accInd);
    
    
    % Bin up the theta frequency and amplitude in speed bins
    [speedAxis, thetaFreqBin, thetaAmpBin] =  speedThetaBinning(dSpeed, dPosTheta, dPosAmplitude, p);
    [R2_amplitude_d, slope_amplitude_d, intercept_amplitude_d] = regressionAnalysis(speedAxis, thetaAmpBin);
    R_amplitude_d = sqrt(R2_amplitude_d);
    
    [R2_frequency_d, slope_frequency_d, intercept_frequency_d] = regressionAnalysis(speedAxis, thetaFreqBin);
    R_frequency_d = sqrt(R2_frequency_d);
    
    figure(1)
    bar(speedAxis, thetaAmpBin)
    xlabel('Speed [cm/s]');
    ylabel('Amplitude [microvolts]')
    title('Theta amplitude (EEG) when deaccelrating')
    fName = strcat(dirPath,'_thetaAmplitude_deacceleration');
    imageStore(figure(1),p.imageFormat,fName,300);
    
    figure(2)
    bar(speedAxis, thetaFreqBin)
    xlabel('Speed [cm/s]');
    ylabel('Frequency [Hz]')
    title('Theta Frequency (EEG) when deaccelerating')
    fName = strcat(dirPath,'_thetaFrequency_deacceleration');
    imageStore(figure(2),p.imageFormat,fName,300);
    
    
    if p.lowSpeedThreshold > 0 || p.highSpeedThreshold > 0
        % Speed filter the position data
        if p.lowSpeedThreshold > 0 && p.highSpeedThreshold > 0
            ind = find(speed < p.lowSpeedThreshold | speed > p.highSpeedThreshold);
        elseif p.lowSpeedThreshold > 0 && p.highSpeedThreshold == 0
            ind = find(speed < p.lowSpeedThreshold );
        else
            ind = find(speed > p.highSpeedThreshold );
        end
        
        % Remove the segments that have to high or to low speed
        posx(ind) = NaN;
        posy(ind) = NaN;
    end
    
    % Find cells part of this session
    cInd = find(unitArray(:,3) == s);
    
    % Number of cells in this session
    numCellsSession = length(cInd);
    
    for c = 1:numCellsSession
        cellId = sprintf('%s%u%s%u','T',unitArray(cInd(c),1),'C',unitArray(cInd(c),2));
        cellFileName = sprintf('%s%s%u%s%u%s',sessionArray{s},'_T',unitArray(cInd(c),1),'C',unitArray(cInd(c),2),'.mat');
        
        % Load the cell data
        load(cellFileName)
        
        % Make sure the correct variables were loaded from the file
        if ~exist('cellTS','var')
            disp('The cell file is missing the cellTS variable')
            return
        end
        
        fprintf(fid2,'%s\t', sessionArray{s});
        fprintf(fid2,'%s\t', cellId);
        fprintf(fid2,'%3.2f\t', R_amplitude);
        fprintf(fid2,'%3.2f\t', slope_amplitude);
        fprintf(fid2,'%3.2f\t', intercept_amplitude);
        fprintf(fid2,'%3.2f\t', R_frequency);
        fprintf(fid2,'%3.2f\t', slope_frequency);
        fprintf(fid2,'%3.2f\t', intercept_frequency);
        
        fprintf(fid2,'%3.2f\t', R_amplitude_a);
        fprintf(fid2,'%3.2f\t', slope_amplitude_a);
        fprintf(fid2,'%3.2f\t', intercept_amplitude_a);
        fprintf(fid2,'%3.2f\t', R_frequency_a);
        fprintf(fid2,'%3.2f\t', slope_frequency_a);
        fprintf(fid2,'%3.2f\t', intercept_frequency_a);
        
        fprintf(fid2,'%3.2f\t', R_amplitude_d);
        fprintf(fid2,'%3.2f\t', slope_amplitude_d);
        fprintf(fid2,'%3.2f\t', intercept_amplitude_d);
        fprintf(fid2,'%3.2f\t', R_frequency_d);
        fprintf(fid2,'%3.2f\t', slope_frequency_d);
        fprintf(fid2,'%3.2f\t', intercept_frequency_d);
        
        % Run the main function
        mainFunction(posx, post, speed, cellTS, cellId, p, dirPath, fid2)
        
        fprintf(fid2,'\n');
    end
end


fclose('all');
close all
fprintf('%s%s\n','Finished: ', datestr(now));
disp('====================================================================');


%__________________________________________________________________________
%
%                               Main
%__________________________________________________________________________


% Main function that runs all analyses for each cell
function mainFunction(posx, post, speed, ts, cellId, p, dirPath,fid2)

dirPath = sprintf('%s%s%s',dirPath,'_',cellId);

fprintf('%s%s\n','Analysing cell ', cellId);

if length(ts) < p.minimumNumSpikes
    disp('WARNING: Cell is skipped because of to few spikes');
    return
end

% Number of bins in the autocorrelogram
numBins = ceil(p.diagramLength/p.binWidth);


% Remove spikes that fall in removed sections of the path and get the
% position indices for the spikes
[ts, spkInd] = legalSpikes(posx, post, ts);

% Speed of animal at each spike time
spkSpeed = speed(spkInd);

% Find the spikes that fall into each speed bin
maxSpeed = min([nanmax(speed), p.maximumSpeedForPlot]);

numSpeedBins = ceil(maxSpeed / p.intrinsicSpeedBinWidth);

% Sampling frequency
Fs = 1/p.binWidth;

% Set fft length to 2^16
nfft = 2^16;


% Calculate the autocorrelogram for all the spikes
autoCorrGram = zeros(numBins,1);

N = length(ts);

% Find number of spikes in each bin of the correlogram
for ii = 1:N-1
    ref = ts(ii);
    for jj = 1:numBins
        low = (jj-1) * p.binWidth;
        high = jj * p.binWidth;
        autoCorrGram(jj) = autoCorrGram(jj) + length(find(ts(ii+1:end) >= (ref+low) & ts(ii+1:end) < (ref+high)));
    end
end

% set the axis for the autocorrelogram plot
autoCorrGramAxis = zeros(numBins,1);
low = 0;
high = p.binWidth;
for jj = 1:numBins
    autoCorrGramAxis(jj) = (low+high) / 2;
    low = low + p.binWidth;
    high = high + p.binWidth;
end

sLength = p.sampleTime * length(post);



% Convert from number of spikes to rate
autoCorrGram = autoCorrGram / sLength;



% Reduce the peak at zero lag to the highest values outside the peak
maxOutside = nanmax(autoCorrGram(2:end));
if autoCorrGram(1) > maxOutside
    autoCorrGram(1) = maxOutside;
end

% Smooth the ISI for plotting
if p.smoothingFactor > 0
    autoCorrGram = mapSmoothing(autoCorrGram, p);
end

% Remove DC component of the the ISI
autoCorrGramFFT = autoCorrGram - mean(autoCorrGram);

% Calculate the power spectrum of the autocorrelogram using fft
[fftACG,f] = fftPowerSpectrum(autoCorrGramFFT,Fs,nfft);

signalLength = length(fftACG);

% Calculate the peak of the power spectrum in the theta band
lowInd = round(p.thetaBand(2) * signalLength / (Fs/2)) + 1;
highInd = round(p.thetaBand(3) * signalLength / (Fs/2) + 1);
[~,peakInd] = nanmax(fftACG(lowInd:highInd));
peakIndFftACG = lowInd + peakInd - 1;
intrinsicFreq = f(peakIndFftACG);



% Make figure of the power spectrum
ind = find(f >= 20,1,'first');
figure(4)
plot(f(1:ind),fftACG(1:ind));
fName = sprintf('%s%s',dirPath,'_PowerSpectrum');
imageStore(figure(4),p.imageFormat,fName,300);




autoCorrGramAxis = autoCorrGramAxis * 1000;

figure(3)
bar(autoCorrGramAxis,autoCorrGram)
xlabel('Inter spike interval [milliseconds]')
ylabel('Rate [Hz]')
fName = sprintf('%s%s',dirPath,'_ISI');
imageStore(figure(3),p.imageFormat,fName,300);
    
% Calculate the isi using smoothed autocorrelogram 
isi = isiCalculation(autoCorrGram, autoCorrGramAxis, p);

intrinsicFreqSpeedBins = zeros(numSpeedBins,1);
isiFreqSpeedBins = zeros(numSpeedBins,1);

startSpeed = 0;
stopSpeed = p.intrinsicSpeedBinWidth;
for s = 1:numSpeedBins
    
    % Find spikes in this speed bin
    sTs = ts(spkSpeed >= startSpeed & spkSpeed < stopSpeed);
    
    % Number of position samples in this speed bin
    numSamples = length(find(speed>= startSpeed & speed < stopSpeed));
    sLength = p.sampleTime * numSamples;
    
    % -- Construct the autocorrelogram --
    autoCorrGram = zeros(numBins,1);
    

    % Number of spikes
    N = length(sTs);

    % Find number of spikes in each bin of the correlogram
    for ii = 1:N-1
        ref = sTs(ii);
        for jj = 1:numBins
            low = (jj-1) * p.binWidth;
            high = jj * p.binWidth;
            autoCorrGram(jj) = autoCorrGram(jj) + length(find(sTs(ii+1:end) >= (ref+low) & sTs(ii+1:end) < (ref+high)));
        end
    end
    
    % set the axis for the autocorrelogram plot
    autoCorrGramAxis = zeros(numBins,1);
    low = 0;
    high = p.binWidth;
    for jj = 1:numBins
        autoCorrGramAxis(jj) = (low+high) / 2;
        low = low + p.binWidth;
        high = high + p.binWidth;
    end
    
    % Convert from number of spikes to rate
    autoCorrGram = autoCorrGram / sLength;

    % Reduce the peak at zero lag to the highest values outside the peak
    maxOutside = nanmax(autoCorrGram(2:end));
    if autoCorrGram(1) > maxOutside
        autoCorrGram(1) = maxOutside;
    end
    
    % Smooth the ISI
    if p.smoothingFactor > 0
        autoCorrGram = mapSmoothing(autoCorrGram, p);
    end
    
    % Remove DC component of the the ISI
    autoCorrGramFFT = autoCorrGram - mean(autoCorrGram);
    
    % Calculate the power spectrum of the autocorrelogram using fft
    [fftACG,f] = fftPowerSpectrum(autoCorrGramFFT,Fs,nfft);
    
    signalLength = length(fftACG);

    % Calculate the peak of the power spectrum in the theta band
    lowInd = round(p.thetaBand(2) * signalLength / (Fs/2)) + 1;
    highInd = round(p.thetaBand(3) * signalLength / (Fs/2) + 1);
    [~,peakInd] = nanmax(fftACG(lowInd:highInd));
    peakIndFftACG = lowInd + peakInd - 1;
    intrinsicFreqSpeedBins(s) = f(peakIndFftACG);

    autoCorrGramAxis = autoCorrGramAxis * 1000;

    figure(3)
    bar(autoCorrGramAxis,autoCorrGram)
    xlabel('Inter spike interval [milliseconds]')
    ylabel('Rate [Hz]')
    fName = sprintf('%s%s%3.1f%s%3.1f',dirPath,'_ISI_Speed_',startSpeed,'-',stopSpeed);
    imageStore(figure(3),p.imageFormat,fName,300);
    
    % Calculate the isi using smoothed autocorrelogram
    isiFreqSpeedBins(s) = isiCalculation(autoCorrGram, autoCorrGramAxis, p);
    
    
    startSpeed = startSpeed + p.intrinsicSpeedBinWidth;
    stopSpeed = stopSpeed + p.intrinsicSpeedBinWidth;
end


% Write result to the output file
fprintf(fid2,'%1.1f\t',isi);
for ii = 1:numSpeedBins
    fprintf(fid2,'%1.1f\t',isiFreqSpeedBins(ii));
end


fprintf(fid2,'%1.1f\t',intrinsicFreq);
for ii = 1:numSpeedBins
    fprintf(fid2,'%1.1f\t',intrinsicFreqSpeedBins(ii));
end



function isi = isiCalculation(autoCorrGram, autoCorrGramAxis, p)



start = floor(0.090/p.binWidth);
stop = ceil(0.140/p.binWidth);

[~,peakInd] = max(autoCorrGram(start:stop));
peakInd = peakInd + start - 1;
isi = autoCorrGramAxis(peakInd);






% Removes spikes from the part of the recording where we are missing
% tracking
function [nTS, spkInd] = legalSpikes(posx, post, ts)

numSpikes = length(ts);
nTS = zeros(numSpikes, 1);
spkInd = zeros(numSpikes, 1);
spikeCounter = 0;

for ii = 1:numSpikes
    tdiff = (post-ts(ii)).^2;
    [~,ind] = min(tdiff);
    if ~isnan(posx(ind(1)))
        spikeCounter = spikeCounter + 1;
        nTS(spikeCounter) = ts(ii);
        spkInd(spikeCounter) = ind(1);
    end
end

nTS = nTS(1:spikeCounter);
spkInd = spkInd(1:spikeCounter);


%__________________________________________________________________________
%
%                          Statistic Functions
%__________________________________________________________________________

% Calculates the regression coefficient of determination (r^2)
function [R2, slope, intercept] = regressionAnalysis(x, y)

% Remove NaNs
ind = isfinite(x);
x = x(ind);
y = y(ind);

ind = isfinite(y);
x = x(ind);
y = y(ind);

N = length(x);

sumY = sum(y);
sumY2 = sum(y.^2);
sumXY = sum(x .* y);
sumX = sum(x);
sumX2 = sum(x.^2);

sum_xy = sumXY -  (sumX * sumY) / N;
sum_x2 = sumX2 - (sumX^2)/N;

% Total sum of squares
totSS = sumY2 - (sumY^2)/N;


% Regression sum of squares
regSS = sum_xy^2 / sum_x2;


R2 = regSS / totSS;


Lxx = sum_x2;
Lxy = sum_xy;

% Regression line slope
slope = Lxy / Lxx;

% Regression line y-axis intersect
intercept = (sumY - slope * sumX) / N;


% Calculates the power spectrum of the input signal using the FFT.
% The input signal is tapered with a Hamming window to reduce the spectral
% leakage.
%
% x     Input signal
% Fs    Sampling rate of the input signal
% nfft  Number of coeffisients of the fft (i.e the frequency resolution of
%       the fft)
%
% y     Power spectrum of the input signal
% f     Frequency array. 
%
% (c) Raymond Skjerpeng
function [y,f] = fftPowerSpectrum(x,Fs,nfft)

% Number of samples in the input signal
numSamp = length(x);

% Make a Hamming window with the same size as the input signal
w = hamming(numSamp);

if size(w,1) ~= size(x,1)
    w = w';
end

% Taper the input signal with the Hamming window to reduce spectral leakage
x = x .* w;


% Calculate the FFT (Automatically padding with zeros if necessary) and
% scale it to the length of the signal (Matlab does not do this
% automatically)
y = fft(x,nfft)/numSamp;


% FFT is symmetric, throw away second half
NumUniquePts = ceil((nfft+1)/2);
y = y(1:NumUniquePts);

% Take the magnitude of the FFT
y = abs(y);

% Take the square of the magnitude of the fft. 
y = y.^2;

% Since we dropped half the FFT, we multiply the power by 2 to keep the same energy.
% The DC component and Nyquist component, if it exists, are unique and should not
% be mulitplied by 2.
if rem(nfft, 2)
    % odd nfft excludes Nyquist point
    y(2:end) = y(2:end) * 2;
else
    y(2:end-1) = y(2:end-1) * 2;
end

% This is an evenly spaced frequency vector with NumUniquePts points.
f = (0:NumUniquePts-1) * Fs/nfft;


%__________________________________________________________________________
%
%                          Filter Functions
%__________________________________________________________________________




% Smooths the map with guassian smoothing
function sMap = mapSmoothing(map, p)

% Calculate the Gaussian smoothing weights
standardDev = max([1, p.smoothingFactor - 1]);
box = pdf('Normal',-p.smoothingFactor:p.smoothingFactor,0,standardDev);

boxLength = length(box);
offset = (boxLength+1) / 2;

numBins = length(map);
sMap = zeros(1,numBins);

for ii = 1:numBins
    for k = 1:boxLength
        % Bin shift
        sii = k-offset;
        % Bin index
        binInd = ii + sii;
        % Boundry check
        if binInd<1
            binInd = 1;
        end
        if binInd > numBins
            binInd = numBins;
        end
        
        sMap(ii) = sMap(ii) + map(binInd) * box(k);
    end
end





% Removes segments from the data where theta is not present and where the
% speed of the animal is outside the set range
function thetaProperties = thetaFiltering(thetaProperties)

% Set the cut off threshold to 2 standard deviations under the mean
meanAmp = nanmean(thetaProperties(:,4));
stdAmp = nanstd(thetaProperties(:,4));

threshold = meanAmp - 2 * stdAmp;

% Remove cycles where the amplitude is to low (No theta)
thetaProperties(thetaProperties(:,4) < threshold,3:4) = NaN;





function xf = fftbandpass(x,Fs,Fs1,Fp1,Fp2,Fs2)
% function XF = fftbandpass(X,FS,FS1,FP1,FP2,FS2)
%
% Bandpass filter for the signal X (time x trials). An acuasal fft 
% algorithm is applied (i.e. no phase shift). The filter functions is         
% constructed from a Hamming window. 
%
% Fs : sampling frequency
%
% The passbands (Fp1 Fp2) and stop bands (Fs1 Fs2) are defined as
%                 -----------                      
%                /           \
%               /             \
%              /               \
%             /                 \
%   ----------                   ----------------- 
%           Fs1  Fp1       Fp2  Fs2              
%
% If no output arguments are assigned the filter function H(f) and
% impulse response are plotted. 
%
% NOTE: for long data traces the filter is very slow.
%
%------------------------------------------------------------------------
% Ole Jensen, Brain Resarch Unit, Low Temperature Laboratory,
% Helsinki University of Technology, 02015 HUT, Finland,
% Report bugs to ojensen@neuro.hut.fi
%------------------------------------------------------------------------

%    Copyright (C) 2000 by Ole Jensen 
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You can find a copy of the GNU General Public License
%    along with this package (4DToolbox); if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

if size(x,1) == 1
    x = x';
end
% Make x even
Norig = size(x,1); 
if rem(Norig,2)
    x = [x' zeros(size(x,2),1)]';                
end

% Normalize frequencies  
Ns1 = Fs1/(Fs/2);
Ns2 = Fs2/(Fs/2);
Np1 = Fp1/(Fs/2);
Np2 = Fp2/(Fs/2);

% Construct the filter function H(f)
N = size(x,1);
Nh = N/2;

B = fir2(N-1,[0 Ns1 Np1 Np2 Ns2 1],[0 0 1 1 0 0]); 
H = abs(fft(B));  % Make zero-phase filter function
IPR = real(ifft(H));
if nargout == 0 
    subplot(2,1,1)
    f = Fs*(0:Nh-1)/(N);
    plot(f,H(1:Nh));
    xlim([0 2*Fs2])
    ylim([0 1]); 
    title('Filter function H(f)')
    xlabel('Frequency (Hz)')
    subplot(2,1,2)
    plot((1:Nh)/Fs,IPR(1:Nh))
    xlim([0 2/Fp1])
    xlabel('Time (sec)')
    ylim([min(IPR) max(IPR)])
    title('Impulse response')
end


if size(x,2) > 1
    for k=1:size(x,2)
        xf(:,k) = real(ifft(fft(x(:,k)) .* H'));
    end
    xf = xf(1:Norig,:);
else
    xf = real(ifft(fft(x') .* H));
    xf = xf(1:Norig);
end





%__________________________________________________________________________
%
%                               EEG
%__________________________________________________________________________




% Calculates mean theta amplitude for speed bins
function [speedAxis, posThetaBin, posAmplitudeBin] =  speedThetaBinning(speed, posTheta, posAmplitude, p)

% Maximum speed
maxSpeed = min([nanmax(speed), p.maximumSpeedForPlot]);

numSpeedBins = ceil(maxSpeed / p.binWidthSpeed);

speedAxis = zeros(numSpeedBins,1);
posThetaBin = zeros(numSpeedBins,1);
posAmplitudeBin = zeros(numSpeedBins,1);

start = 0;
stop = p.binWidthSpeed;
for ii = 1:numSpeedBins
    speedAxis(ii) = (start+stop) / 2;
    % Find samples in this speed bin
    ind = find(speed >= start & speed < stop);
    
    % Mean theta frequency in this speed bin
    posThetaBin(ii) = nanmean(posTheta(ind));
    % Mean theta amplitude in this speed bin
    posAmplitudeBin(ii) = nanmean(posAmplitude(ind));
    
    start = start + p.binWidthSpeed;
    stop = stop + p.binWidthSpeed;
end










% Calculates the theta frequency and theta amplitude for each position sample
function [posTheta, posAmplitude] = thetaFrequencyPosSample(post,thetaProperties)

% Number of theta cycles
N = size(thetaProperties,1);
% Number of position samples
numSamples = length(post);
posTheta = zeros(numSamples,1);
posAmplitude = zeros(numSamples,1);

posTheta(post < thetaProperties(1,2)) = thetaProperties(1,3);
posAmplitude(post < thetaProperties(1,2)) = thetaProperties(1,4);
for ii = 2:N
    posTheta(post >= thetaProperties(ii-1,2) & post < thetaProperties(ii,2)) = thetaProperties(ii,3);
    posAmplitude(post >= thetaProperties(ii-1,2) & post < thetaProperties(ii,2)) = thetaProperties(ii,4);
end

posTheta(posTheta == 0) = NaN;
posAmplitude(posAmplitude == 0) = NaN;




% Calculates instantinous theta filtered frequency
function thetaProperties = thetaFrequency(theta, peakInd, Fs)

numPeaks = length(peakInd);

% 1: Peak Index
% 2: Peak Time
% 3: Instantinous frequency
% 4: Instantinous amplitude
thetaProperties = zeros(numPeaks,4);
thetaProperties(:,1) = peakInd;
thetaProperties(:,2) = peakInd / Fs;

for ii = 2:numPeaks
    % Instantaneous frequency
    thetaProperties(ii,3) = 1 / (peakInd(ii)/Fs - peakInd(ii-1)/Fs);
    
    % Instantaneous amplitude
    thetaProperties(ii,4) = nanmax(theta(peakInd(ii-1):peakInd(ii))) - nanmin(theta(peakInd(ii-1):peakInd(ii)));
end




% Extract theta phase from the band pass filtered EEG signal by
% interpolating phases between peaks and troughs.
function [phase,peakPos] = thetaPhase(eeg)

% Number of samples in eeg signal
N = length(eeg);
phase = zeros(1,N);
peakPos = zeros(1,N);
peakCount = 0;

% Flag indicating if we are on the rising or falling part of the wave. 0 =
% falling part, 1 = rising part
rising = 0;
% Sample index to last known peak
lastPeak = 0;
% Sample index to last knows trough
lastTrough = 0;
% Index to the first sample after the first peak or trough
startIndex = 0;

% Find if first sample is on a rising or falling wave
first = eeg(1);
second = eeg(2);
if first < second
    rising = 1;
elseif first>second
    rising = 0;
else
    count = 2;
    while 1
        count = count + 1;
        if count == N
            break
        end
        if first < eeg(count)
            rising = 1;
            break
        end
        if first > eeg(count)
            rising = 0;
            break
        end
    end
end

% Locate the first peak or trough
for ii = 2:N
    if rising
        if eeg(ii) < eeg(ii-1)
            lastPeak = ii-1;
            rising = 0;
            startIndex = ii;
            % Add the peak to the peak array
            peakCount = peakCount + 1;
            peakPos(peakCount) = lastPeak;
            break;
        end
    else
        if eeg(ii) > eeg(ii-1)
            lastTrough = ii-1;
            rising = 1;
            startIndex = ii;
            break
        end
    end
end

% Set the phase before the first peak/trough to NaN, cause it is not
% possible to calculate the phase of this section with this method
 phase(1:startIndex-1) = NaN;

% Calculate the phase of the theta (bandpass filtered eeg) based on the
% location of peaks and troughs in the signal. Peaks are set as 0/360
% degrees and the troughs are set to 180 degrees. Any sample in between the
% peaks and troughs are interpolated between these two values.
for ii = startIndex:N
    if rising
        if eeg(ii) < eeg(ii-1)
            lastPeak = ii-1;
            phase(lastPeak) = 0;
            numSamp = lastPeak - lastTrough - 1;
            samps = 1:numSamp;
            phase(lastTrough+1:lastPeak-1) = 180 + samps*180/(numSamp+1);
            rising = 0;
            % Add the peak to the peak array
            peakCount = peakCount + 1;
            peakPos(peakCount) = lastPeak;
        end
    else
        if eeg(ii) > eeg(ii-1)
            lastTrough = ii-1;
            phase(lastTrough) = 180;
            numSamp = lastTrough - lastPeak - 1;
            samps = 1:numSamp;
            phase(lastPeak+1:lastTrough-1) = samps*180/(numSamp+1);
            rising = 1;
        end
    end
end

lastValue = max(lastPeak,lastTrough);
phase(lastValue+1:N) = NaN;

peakPos = peakPos(1:peakCount);






%__________________________________________________________________________
%
%                           Position processing
%__________________________________________________________________________






% Calculates the sign of the accelration
% 1 = Acceleration
% 0 = Deacceleration or no change in speed
function accSign = accelrationSign(s, offset)

N = length(s);
accSign = zeros(N,1);

for ii = 1+offset:N-offset
    if s(ii+offset) > s(ii-offset)
        accSign(ii) = 1;
    end
end


% Calculate the Speed of the rat in each position sample
%
% Version 1.0
% 3. Mar. 2008
% (c) Raymond Skjerpeng, CBM, NTNU, 2008.
function v = speed2D(x,y,t)

N = length(x);
v = zeros(N,1);

for ii = 2:N-1
    v(ii) = sqrt((x(ii+1)-x(ii-1))^2+(y(ii+1)-y(ii-1))^2)/(t(ii+1)-t(ii-1));
end
v(1) = v(2);
v(end) = v(end-1);






%__________________________________________________________________________
%
%                           Input data
%__________________________________________________________________________




function [status, sessionArray, unitArray, shapeArray] = inputFileReader(inputFile, p)


% Status = 0 -> Input file contain errors
% Status = 1 -> Input file is ok
status = 0;

% Number of sessions possible to have listed in the input file
N = 1000;

% Mean number of cell per session
M = 100;

% Session name array
% 1: Session name (whole path)
sessionArray    = cell(N, 1);

% Box shape
shapeArray      = cell(N, 1);

% Tetrode and cell number.
% 1: Tetrode.
% 2: Cell number.
% 3: Session number. Tell which session the cell belongs to.
unitArray       = zeros(M*N, 3);

% Open the input file for binary read access
fid = fopen(inputFile,'r');

if fid == -1
    msgbox('Could''n open the input file! Make sure the filname and path are correct.','File read error','error');
    disp('Input file could not be found.')
    disp('Failed')
    return
end

% Counts the number of sessions
sessionCounter = 0;
% Count the number of cells
unitCounter = 0;

% Keep track of the line number the programs reads from in the input file
currentLine = 0;


while ~feof(fid)
    % Read a line from the input file
    str = fgetl(fid);
    currentLine = currentLine + 1;
    
    % Remove space at end of line
    str = stripSpaceAtEndOfString(str);
    
    % Check that line is not empty
    if isempty(str)
        disp('Error: There can''t be any empty lines in the input file');
        fprintf('%s%u\n','Empty line was found in line number ',currentLine);
        return
    end
    
    % Check that the line is the "session" line
    if length(str)<7
        disp('Error: Expected keyword ''Session'' in the input file');
        fprintf('%s%u\n','Error on line ', currentLine);
        return
    end
    if ~strcmpi(str(1:7),'Session')
        disp('Error: Expected keyword ''Session'' in the input file');
        fprintf('%s%u\n','Error on line ', currentLine);
        return
    else
        sessionCounter = sessionCounter + 1;
        sessionArray{sessionCounter,1} = str(9:end);
        


        if strcmpi(p.delim,'/')
            sessionArray{sessionCounter,1} = strrep(sessionArray{sessionCounter,1},'\','/');
        else
            sessionArray{sessionCounter,1} = strrep(sessionArray{sessionCounter,1},'/','\');
        end

        
        % Read next line
        str = fgetl(fid);
        currentLine = currentLine + 1;
        
        % Remove space at end of line
        str = stripSpaceAtEndOfString(str);
        
        % Check that line is not empty
        if isempty(str)
            disp('Error: There can''t be any empty lines in the input file');
            fprintf('%s%u\n','Empty line was found in line number ',currentLine);
            return
        end
    end
    

    % Shape information should come next
    % 1 dim: shape. 1 = box, 2 = cylinder, 3 = linear track
    % 2 dim: Side length or diameter of the arena.
    shape = zeros(2,1);
    if length(str)<5 || ~strcmpi(str(1:5),'Shape')
        fprintf('%s%u\n','Error: Expected the ''Shape'' keyword at line ', currentLine)
        return
    else
        temp = str(7:end);
        if length(temp)>3 && strcmpi(temp(1:3),'Box')
            shape(1) = 1;
            shape(2) = str2double(temp(5:end));

        elseif length(temp)>5 && strcmpi(temp(1:5),'Track')
            shape(1) = 3;
            shape(2) = str2double(temp(7:end));
        elseif length(temp) > 6 && strcmpi(temp(1:6), 'Circle')
            shape(1) = 2;
            shape(2) = str2double(temp(8:end));
        elseif length(temp)>8 && strcmpi(temp(1:8),'Cylinder')
            shape(1) = 2;
            shape(2) = str2double(temp(10:end));
        else
            disp('Error: Missing shape information. Must be box, cylinder or Track');
            fprintf('%s%u\n','Error at line ', currentLine)
            return
        end


        % Add the shape information to the shape array
        shapeArray{sessionCounter} = shape;

        % Read next line
        str = fgetl(fid);
        currentLine = currentLine + 1;

        % Remove space at end of line
        str = stripSpaceAtEndOfString(str);

        % Check that line is not empty
        if isempty(str)
            disp('Error: There can''t be any empty lines in the input file');
            fprintf('%s%u\n','Empty line was found in line number ',currentLine);
            return
        end
    end
    
    
    while ~feof(fid)
        if strcmp(str,'---') % End of this block of data, start over.
            break
        end
        
        if length(str)>7
            if strcmpi(str(1:7),'Tetrode')
                tetrode = sscanf(str,'%*s %u');
                
                % Read next line
                str = fgetl(fid);
                currentLine = currentLine + 1;
                
                % Remove space at end of line
                str = stripSpaceAtEndOfString(str);

                % Check that line is not empty
                if isempty(str)
                    disp('Error: There can''t be any empty lines in the input file');
                    fprintf('%s%u\n','Empty line was found in line number ',currentLine);
                    return
                end
                
                while length(str) > 4 && strcmpi(str(1:4),'Unit')
                    unit = sscanf(str,'%*s %u');
                    unitCounter = unitCounter + 1;
                    
                    % Add tetrode and cell number to the cell array
                    unitArray(unitCounter,1) = tetrode;
                    unitArray(unitCounter,2) = unit;
                    unitArray(unitCounter,3) = sessionCounter;
                    
                    str = fgetl(fid);
                    currentLine = currentLine + 1;
                    
                    % Remove space at end of line
                    str = stripSpaceAtEndOfString(str);
                    
                    % Check that line is not empty
                    if isempty(str)
                        disp('Error: There can''t be any empty lines in the input file');
                        fprintf('%s%u\n','Empty line was found in line number ',currentLine);
                        return
                    end
                end
            else
                fprintf('%s%u\n','Error: Expected the Tetrode keyword at line ', currentLine);
                return
            end
        else
            fprintf('%s%u\n','Error: Expected the Tetrode keyword at line ', currentLine);
            return
        end
        
    end
    
end

sessionArray = sessionArray(1:sessionCounter,:);
shapeArray = shapeArray(1:sessionCounter,:);
unitArray = unitArray(1:unitCounter,:);

% Set status to success (1)
status = 1;


% Removes space at the end of the string input
function str = stripSpaceAtEndOfString(str)

if isempty(str)
    return
end

while ~isempty(str)
    if strcmp(str(end),' ')
        str = str(1:end-1);
    else
        break;
    end
end



function status = inputDataChecker(sessionArray, unitArray)

status = 0;

% Number of sessions in the input file
numSessions = size(sessionArray,1);

if numSessions == 0
    disp('Error: No sessions was listed')
    return
end

for ii = 1:numSessions

    % Look for position file
    videoFile = strcat(sessionArray{ii,1},'_pos.mat');
    d = dir(videoFile);
    if size(d,1) == 0
        fprintf('%s%s\n','Unable to find the position file: ',videoFile);
        disp('Please check your input file and data.')
        return
    end
    
    % Look for EEG files
    eegFile = strcat(sessionArray{ii,1},'_eeg.mat');
    d = dir(eegFile);
    if size(d,1) == 0
        fprintf('%s%s\n','Unable to find the EEG file: ',eegFile);
        disp('Please check your input file and data.')
        return
    end

    
    % Find cells listed for this session
    ind = find(unitArray(:,3) == ii);
    for c = 1:length(ind)
        cellFileName = sprintf('%s%s%u%s%u%s',sessionArray{ii},'_T',unitArray(ind(c),1),'C',unitArray(ind(c),2),'.mat');
        d = dir(cellFileName);
        if size(d,1) == 0
            fprintf('%s%s\n','Unable to find the cell file: ',cellFileName);
            disp('Please check your input file and data.')
            return
        end
    end
    
end

% Set status to success
status = 1;









% Function for storing figures to file
% figHanle  Figure handle (Ex: figure(1))
% format = 1 -> bmp (24 bit)
% format = 2 -> png
% format = 3 -> eps
% format = 4 -> jpg
% format = 5 -> ai (Adobe Illustrator)
% format = 6 -> tiff (24 bit)
% format = 7 -> fig (Matlab figure)
% figFile   Name (full path) for the file
% dpi       DPI setting for the image file
function imageStore(figHandle,format,figFile,dpi)

% Make the background of the figure white
set(figHandle,'color',[1 1 1]);
dpi = sprintf('%s%u','-r',dpi);

switch format
    case 1
        % Store the image as bmp (24 bit)
        figFile = strcat(figFile,'.bmp');
        print(figHandle, dpi, '-dbmp',figFile);
    case 2
        % Store image as png
        figFile = strcat(figFile,'.png');
        print(figHandle, dpi,'-dpng',figFile);
    case 3
        % Store image as eps (Vector format)
        figFile = strcat(figFile,'.eps');
        print(figHandle, dpi,'-depsc',figFile);
    case 4
        % Store image as jpg
        figFile = strcat(figFile,'.jpg');
        print(figHandle,dpi, '-djpeg',figFile);
    case 5
        % Store image as ai (Adobe Illustrator)
        figFile = strcat(figFile,'.ai');
        print(figHandle,dpi, '-dill',figFile);
    case 6
        % Store image as tiff (24 bit)
        figFile = strcat(figFile,'.tif');
        print(figHandle,dpi, '-dtiff',figFile);
    case 7
        % Store figure as Matlab figure
        figFile = strcat(figFile,'.fig');
        saveas(figHandle,figFile,'fig')
end