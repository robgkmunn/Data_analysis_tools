% Calculates instanteous theta filtered frequency
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

