function [ xhat] = Gao_RobustCSS(y1,A1,L,sigma,eta)
%% The function is used to recovery original signal from compressed measurements used for TSP paper...
%% ``Wideband Spectrum Sensing on Real-time Signals at Sub-Nyquist Sampling Rates in Single and Cooperative Multiple Nodes''


%%Inputs
%% y1 compressed measurements
%% A1: measurement matrix, AIC sampler is used in the paper
%% L? length of original signal
%% sigma: noise level
%% eta: recovery tolerance

%% Output
%% xhat: recovered signal

mm=length(y1); % number of compressed measurements

%% solve single node case
cvx_begin
    variable x(L) complex;
    minimize(norm(x,1));
    subject to
        norm(A1*x-y1,2)<=eta;
cvx_end

%             %% solve the multiple nodes case
%% J: number of SUs in multiple nodes case
%   cvx_begin
%             variable x(L,J);
%             minimize(norm_nuc(x));
%             subject to
%                 norm(y1-A1*vec(x),2) <=eta;
%         cvx_end

%% denoising
% lx1=find(x>sigma*sqrt(mm));
% A1T=A1(:,lx1);
% xhat=zeros(L,1);
% xhat(lx1)=A1T\y1;
xhat=x;
end

