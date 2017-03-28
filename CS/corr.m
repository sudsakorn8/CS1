function y = corr(x)
% sigma(x(n)*x(n-d))
% n=100;
n = length(x);
y = zeros(1,2*n+1);
idx = 1;
for num = n:-1:0
    xd1 = [zeros(1,num),x];
    xd2 = [x,zeros(1,num)];
    y(idx) = sum(xd1.*xd2);
    idx = idx+1;
    fprintf([num2str(idx),'\n']);
end
for num = 1: n
    xd1 = [x,zeros(1,num)];
    xd2 = [zeros(1,num),x];
    y(idx) = sum(xd1.*xd2);
    idx = idx+1;
    fprintf([num2str(idx),'\n']);
end
fprintf('yes \n');