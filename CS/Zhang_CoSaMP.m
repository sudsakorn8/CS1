function [ a ] = Zhang_CoSaMP( u,Phi,s )  
%   CoSaOMP construct a signal approximation given the compressive samples.
%   Based on the CoSaMP algorithm proposed in "Iterative signal recovery from incomplete and inaccurate samples"
%   Detailed explanation goes here  
%   input u(R*1) is the noisy sample vector
%   input Phi(R*W) is the measurement matrix 
%   input s is the sparsity level
%   output a(W*1) is the s-sparse approximation of the target signal(recovery signal)
%   u = Phi*a + e  

if size(u,1) < size(u,2)  
    u = u'; % transpose u to a column vector  
end  
[R,W] = size(Phi); % Acquire the compressive parameters Phi(R*W)

%% initialize
a_k = zeros(W,1); % the approximation
v = u; % the current samples
T = []; % the support set

%% iteration
for k=1:s % the iteration time should be no more than the sparsity level 
    %(1) Identification  
    fprintf('CoSaMP iteration %d \n',k);
    y = Phi'*v; % form a proxy of the residual from the current samples
    [~,idx]=sort(abs(y),'descend');  
    Omega = idx(1:2*s); % identify the 2s largest components supp(y_2s)
    %(2) Support Merger  
    T = union(T,Omega); % merge support  
    %(3) Estimation  
    if length(T)<=R
        Phi_t = Phi(:,T);  
    else 
        if k == 1  
            a_k = 0; % initial a_k
        end  
        break;
    end 
    b_T = (Phi_t'*Phi_t)^(-1)*Phi_t'*u; % Signal estimation by least squares  
    %(4) Pruning  
    [~,idx]=sort(abs(b_T),'descend');
    T = T(idx(1:s));
    a_k = zeros(W,1);
    a_k(T) = b_T(idx(1:s)); %prune to obtain next approximation
    %(5) Sample Update
    v = u - Phi*a_k; % Update current samples  
    if norm(v)<1e-6 % halting criterion  
        break;   
    end  
end  
a = a_k;  
end


