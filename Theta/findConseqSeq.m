
function [x_consec] = findConseqSeq(x,max_con_seq)

% find the numbers in the sequence
x_num = find(isfinite(x));

x_cum_seq=cumsum([true;diff(x_num)~=1]);

xx = histc(x_cum_seq,1:x_cum_seq(end));
xx(:,2)=1:length(xx);
goodseqs = xx(xx(:,1 ) >= max_con_seq,2);

%lists indexes of "good" angle non-detects (that are in sets of <= maxconsec)
goodnons = x_num(ismember(x_cum_seq, goodseqs));

x_consec = nan(size(x));
x_consec(goodnons) = x(goodnons);

return