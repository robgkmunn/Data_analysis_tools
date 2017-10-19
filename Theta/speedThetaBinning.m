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