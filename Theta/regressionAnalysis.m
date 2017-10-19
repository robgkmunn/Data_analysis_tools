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