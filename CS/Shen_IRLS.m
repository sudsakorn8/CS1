function [ xhat] = Shen_IRLS(y,Phi,eta,s)
%   input u(R*1) is the noisy sample vector
%   Based on algorithm proposed in "Iteratively Reweighted Least Squares Minimizaion for Sparse Recovery"
%   input Phi(R*W) is the measurement matrix 
%   input s is the sparsity level
%   output a(W*1) is the s-sparse approximation of the target signal(recovery signal)
%   y = Phi*x + e 

[R,W] = size(Phi); % Acquire the compressive parameters Phi(R*W)

%% initialize
Omega = ones(W,1);
Epsilon = 1;
D = zeros(W,W);
for num = 1:W
    D(num,num) = 1/Omega(num);%represents Omega
end
k=0;
while Epsilon>eta
    k = k+1;
    x = D*Phi'/(Phi*D*Phi')*y;
    x_sort=sort(abs(x),'descend');  
    Epsilon = min(Epsilon,x_sort(s+1)/W)% s should be K,but I don't know what is K?
    for num = 1:W
        D(num,num) = 1/(x(num)^2+Epsilon^2)^-0.5;
    end    
end
xhat = x;

