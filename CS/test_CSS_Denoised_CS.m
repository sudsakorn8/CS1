function [ P_f,P_d,idx_Channel,if_vacant,E,E_ori,lamda_d,P_d_th,a_w,a_r ] = test_CSS_Denoised_CS( SNR,CS_ratio,Sparsity,Pf,lamdaRatio )
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
Channel_num = 10; % should divides 2 because of conjugate symmetry of spectrum
Channel_BW = 8*1e3; % kHz
BW = Channel_num * Channel_BW; % Whole band width
delta_f = 100; % kHz
N_OFDM = Channel_BW / delta_f; % length of a OFDM symbol
N = N_OFDM * Channel_num; % Original sampling rate
Channel_used = Channel_num * Sparsity; % number of used Channels
s = N_OFDM * Channel_used; % Sparsity level
R = round(N * CS_ratio);

%%
% index matrix marks frequencies being occupied
idx_Channel = rand(Channel_num,1);
[~,idx_sort] = sort(idx_Channel,'ascend');
idx_1 = idx_sort(1:Channel_used);
idx_Channel(:) = 0;
idx_Channel(idx_1) = 1;
idx = zeros(N,1);
for i = 1:Channel_num
    amp = 1;
    if idx_Channel(i)==1
        idx((i-1)*N_OFDM+1:i*N_OFDM,1) = ((rand(N_OFDM,1)>0.5)*2-1)*amp;%.*(1+rand(N_OFDM,1)); % if occupied: 1/-1
    else
        idx((i-1)*N_OFDM+1:i*N_OFDM,1) = 0; % if vacant: 0
    end
end

%Create sparse signal a_w
a_w = idx;
power_e = sum(conj(a_w).*a_w)/(10^(SNR/10)*N);
e = sqrt(0.5*power_e)*randn(N,1)+sqrt(0.5*power_e)*1i*randn(N,1);
a_w = a_w+e;
% a_w_o = a_w;
% a_w = awgn(a_w_o,SNR,'measured');
% power_e = var(a_w-a_w_o);

%Create amplitude vector u(with phase information)
u = zeros(N,1);
for k = 1:N
    u(k) = a_w(k)*((exp(-2*pi*1i*(k-N/2)/N)-1)/(2*pi*1i*(k-N/2)));
end
u(N/2) = a_w(N/2)/N; 
 
%% Compressing through Random Demodulator
% [ y2,y1,y,Phi,Omega_n ] = time_RandomDemodulator( a_w,R );
% figure
% stem(Omega_n,abs(a_w));
% title('signal in discrete frequency')
[y,Phi] = test_RandomDemodulator(a_w,R);
%% Recovery signal 
u_r = Gao_RobustCSS(y,Phi,N,sqrt(power_e),1e-2*norm(y,2)^2);
% u_r = zeros(N,1);
% u_r = Zhang_CoSaMP( y,Phi,s );
% u_r = Shen_IRLS(y,Phi,10e-8,s);
% u_r = cs_irls(y,Phi,s);
% u_r = cs_cosamp(y,Phi,s);
fprintf('Recovered \n');
% error
% r = norm(u_r-u,2);

% Recover a
% a_r = zeros(N,1);
% for k = 1:N
%     a_r(k) = u_r(k)/((exp(-2*pi*1i*(k-N/2)/N)-1)/(2*pi*1i*(k-N/2)));
% end
% a_r(N/2) = u_r(N/2)*N; 
a_r = u_r;
% a_r=0;
%% Denoised
%a_r_d = (a_r>power_e) .* a_r;

%% Spectrum Sensing: Energy Detection
sigma_n = power_e; % noise variance
% Pf = 0.01; % False alarm probability
lamda_d = lamdaRatio*sigma_n*(1+qfuncinv(Pf)/sqrt(N_OFDM));
% lamda_d =0.3;
E_a = conj(a_r).*a_r;
E_a_ori = conj(a_w).*a_w;
%E_a_d = conj(a_r_d).*a_r_d;
E = zeros(Channel_num,1);
E_ori = zeros(Channel_num,1);
%E_d = zeros(Channel_num,1);
for i = 1:Channel_num
    E(i) = sum(E_a((i-1)*N_OFDM+1:i*N_OFDM,1))/N_OFDM;
    E_ori(i) = sum(E_a_ori((i-1)*N_OFDM+1:i*N_OFDM,1))/N_OFDM;
    %E_d(i) = sum(E_a_d((i-1)*N_OFDM+1:i*N_OFDM,1))/N_OFDM;
end
if_vacant = (E<lamda_d);
%if_vacant_d = (E_d<lamda_d);
if_occupied=(E>lamda_d);

%%
% detect = if_vacant+idx_Channel-ones(Channel_num,1);
% error_t = sum(abs(detect));
% detect_d = if_vacant_d+idx_Channel-ones(Channel_num,1);
% error_d = sum(abs(detect_d));

detect_vacant = sum((ones(Channel_num,1)-idx_Channel).*if_vacant);
P_f = 1-detect_vacant/sum(ones(Channel_num,1)-idx_Channel);
detect_occupied = sum(idx_Channel.*if_occupied);
P_d = detect_occupied/sum(idx_Channel);
snr = 10^(SNR/10)/Sparsity;
P_d_th = 1-qfunc((1+snr-lamda_d/sigma_n)/sqrt((2*snr^2+1+2*snr)/N_OFDM));

% P_d_th = qfunc((qfuncinv(Pf)-10^(SNR/10)*sqrt(N_OFDM/2))/(1+10^(SNR/10)));
% P_d_th = qfunc(sqrt(N/2)*(lamda_d/power_e-(1+10^(SNR/10)))/(1+10^(SNR/10)));
%detect_d = sum((ones(Channel_num,1)-idx_Channel).*if_vacant_d);
%Pd_d = detect_d/sum(ones(Channel_num,1)-idx_Channel);



end