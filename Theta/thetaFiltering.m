% Removes segments from the data where theta is not present and where the
% speed of the animal is outside the set range
function thetaProperties = thetaFiltering(thetaProperties)

% Set the cut off threshold to 2 standard deviations under the mean
meanAmp = nanmean(thetaProperties(:,4));
stdAmp = nanstd(thetaProperties(:,4));

threshold = meanAmp - 2 * stdAmp;

% Remove cycles where the amplitude is to low (No theta)
thetaProperties(thetaProperties(:,4) < threshold,3:4) = NaN;
