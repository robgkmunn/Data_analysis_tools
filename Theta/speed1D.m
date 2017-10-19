function v = speed1D_rob(dim,t)

N = length(dim);
M = length(t);

if N < M
    dim = dim(1:N);
    t = t(1:N);
else
    dim = dim(1:M);
    t = t(1:M);
end

v = zeros(min([N,M])-1,1);
n = length(dim);

for i = 1:n-1 
    dt(i) = t(i+1)-t(i);
end

for i = 1:n -1
v(i) = abs((dim(i+1) - dim(i)))/dt(i);
end
