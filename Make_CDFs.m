
A = uigetfile('*.mat','Choose your first group');
B = uigetfile('*.mat','Choose your second group');
A = load(A);
B = load(B);
% Compute the histogram of A and B.
[countsA, binsA] = hist(A);
[countsB, binsB] = hist(B);
% Compute the cumulative distribution function of A and B.
cdfA = cumsum(countsA) / sum(countsA);
cdfB = cumsum(countsB) / sum(countsB);
% Plot the probability distribution of A.
subplot(2,2, 1);
bar(binsA, countsA);
title('Histogram of A', '16', 16);
ylabel('Count A', '16', 16);
xlabel('Values of A', '16', 16);
grid on;
% Plot the probability distribution of B.
subplot(2,2, 2);
bar(binsB, countsB);
title('Histogram of B', '16', 16);
ylabel('Count B', '16', 16);
xlabel('Values of B', '16', 16);
grid on;
% Plot the cumulative distribution function of A.
subplot(2,2, 3);
bar(binsA, cdfA); % You can use plot() if you want to.
title('CDF of A', '16', 16);
grid on;
ylabel('Percentage A (/100)', '16', 16);
xlabel('Values of A', '16', 16);
% Plot the cumulative distribution function of B.
subplot(2,2, 4);
bar(binsB, cdfB); % You can use plot() if you want to.
title('CDF of B', '16', 16);
ylabel('Percentage B (/100)', '16', 16);
xlabel('Values of B', '16', 16);
grid on;
% Enlarge figure to full screen.
set(gcf, 'units','normalized','outerposition',[0 0 1 1]); % Maximize figure.
