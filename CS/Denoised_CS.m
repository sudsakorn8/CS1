function [ error_t,error_d,a_r,a_r_d ] = Denoised_CS( SNR,CS_ratio,Sparsity,Pf )
% Denoised_CS returns the error number of recovery signal (of the 40
% channels) in one trial using both traditional and proposed denoised CS algorithm
% input SNR is the signal to noise ratio
% input CS_ratio is the ratio of real sampling rate to Nyquist rate
% input Sparsity is the sparsity level of the signal (ratio)
% input Pf is the false alarm probability
% CS_ratio should be at least 2 times to the Sparsity level
% output error_t is the error channel number using traditional CS
% output error_r is the error channel number using denoised CS
%% Parameters
% White Space: 470MHz~790MHz
% Total Band Width: 8MHz/Channel * 40Channels = 320MHz
Channel_num = 40;
Channel_BW = 8*1e3; % kHz
BW = Channel_num * Channel_BW; % Whole band width
delta_f = 100; % kHz
N_OFDM = Channel_BW / delta_f; % length of a OFDM symbol
N = N_OFDM * Channel_num; % Original sampling rate
Channel_used = Channel_num * Sparsity; % number of used Channels
s = N_OFDM * Channel_used; % Sparsity level
R = N * CS_ratio;

%%
% index matrix marks frequencies being occupied
idx_Channel = rand(Channel_num,1);
[temp,idx_sort] = sort(idx_Channel,'ascend');
idx_1 = idx_sort(1:Channel_used);
idx_Channel(:) = 0;
idx_Channel(idx_1) = 1;
idx = zeros(N,1);
for i = 1:Channel_num
    idx((i-1)*N_OFDM+1:i*N_OFDM,1) = idx_Channel(i);
end

%Create sparse signal a_w
a_w = idx;

% noise
% e = zeros(N,1);
% SNR = 15; % SNR in dB
power_e = sum(conj(a_w).*a_w)/(10^(SNR/10)*N);
e = randn(N,1)+1i*randn(N,1);
e = sqrt(power_e)*e/sqrt(sum(conj(e).*e)/N);
a_w = a_w+e;

%Create amplitude vector u(with phase information)
u = zeros(N,1);
for k = 1:N
    u(k) = a_w(k)*((exp(-2*pi*1i*k/N)-1)/(2*pi*1i*k));
end


% Compressing through Random Demodulator
[y,Phi] = Zhang_RandomDemodulator(u,R);

% Recovery signal through CoSaMP
u_r = Zhang_CoSaMP(y,Phi,s);

% error
r = norm(u_r-u,2);

% Recover a
a_r = zeros(N,1);
for k = 1:N
    a_r(k) = u_r(k)/((exp(-2*pi*1i*k/N)-1)/(2*pi*1i*k));
end

%% Denoised
a_r_d = (a_r>power_e) .* a_r;

%% Spectrum Sensing: Energy Detection
sigma_n = power_e; % noise variance
% Pf = 0.01; % False alarm probability
lamda_d = sigma_n*(1+qfuncinv(Pf)/sqrt(N/2));
E_a = conj(a_r).*a_r;
E_a_d = conj(a_r_d).*a_r_d;
E = zeros(Channel_num,1);
E_d = zeros(Channel_num,1);
for i = 1:Channel_num
    E(i) = sum(E_a((i-1)*N_OFDM+1:i*N_OFDM,1))/N_OFDM;
    E_d(i) = sum(E_a_d((i-1)*N_OFDM+1:i*N_OFDM,1))/N_OFDM;
end
if_vacant = (E<lamda_d);
if_vacant_d = (E_d<lamda_d);

%%
detect = if_vacant+idx_Channel-ones(Channel_num,1);
error_t = sum(abs(detect));
detect_d = if_vacant_d+idx_Channel-ones(Channel_num,1);
error_d = sum(abs(detect_d));
end