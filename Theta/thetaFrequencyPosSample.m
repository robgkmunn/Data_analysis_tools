%calclates theta frequency at given positional samples
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
end