clear,clc;
%% Generate the test signal with noise
% f(t) = sigma(a_w*exp(-2pi*j(w/W)t))
% y = Phi(s+e)
% Compression & Recovery
%% Parameters
W = 320; % Whole bandwidth
K = 32; %Sparsity level 
comp_ratio = K/W; % Compression ratio
% R = 160; % Sampling Rate R~2KlogW R divided W
R = 160; % Sampling Rate R~2KlogW R divided W
S = 32; % Recovery Sparsity

%%
% index matrix marks frequencies being occupied
idx = rand(W,1);
[temp,idx_sort] = sort(idx,'ascend');
idx_1 = idx_sort(1:K);
idx(:) = 0;
idx(idx_1) = 1;

%Create sparse signal s
s = zeros(W,1);
for w = 1:W
    s(w) = idx(w)*((exp(-2*pi*1i*(w-W/2)/W)-1)/(2*pi*1i*(w-W/2)));
end
s(W/2) = idx(W/2)/W; 
% w = 1:W;
% plot(w,abs(s));

% noise
% e = zeros(W,1);
SNR = 10; % SNR in dB
power_e = norm(s,2)^2/(10^(SNR/10)*W);% noise variance,namely average power spectral density
e = randn(W,1)+1i*randn(W,1);
e = sqrt(power_e)*e/norm(e,2);
s_n = s+e;

% Compressing through Random Demodulator
[y,Phi] = Zhang_RandomDemodulator(s_n,R);

% Recovery signal through CoSaMP
% a_cosamp = Zhang_CoSaMP(y,Phi,S);
a_css = Gao_RobustCSS(y,Phi,W,power_e,1e-6);

% error
% r_cosamp = norm(a_cosamp-s,2);
r_css = norm(a_css-s,2);

% if_singular = (At'*At)^(-1);
