%% Parameters
% White Space: 470MHz~790MHz
% Total Band Width: 8MHz/Channel * 40Channels = 320MHz
Channel_num = 10; % should divides 2 because of conjugate symmetry of spectrum
Channel_BW = 8*1e3; % kHz
BW = Channel_num * Channel_BW; % Whole band width
delta_f = 100; % kHz
N_OFDM = Channel_BW / delta_f; % length of a OFDM symbol
N = N_OFDM * Channel_num; % Original sampling rate
Channel_used = Channel_num * 0.1; % number of used Channels
s = N_OFDM * Channel_used; % Sparsity level
R = round(N * 0.25);

%%
% index matrix marks frequencies being occupied
idx_Channel = rand(Channel_num,1);
[~,idx_sort] = sort(idx_Channel,'ascend');
idx_1 = idx_sort(1:Channel_used);
idx_1 = 5;
idx_Channel(:) = 0;
idx_Channel(idx_1) = 1;
idx = zeros(N,1);
for i = 1:Channel_num
    amp = 1;
    if idx_Channel(i)==1
        idx((i-1)*N_OFDM+1:i*N_OFDM,1) = ((rand(N_OFDM,1)>0.5)*2-1)*amp; % if occupied: 1/-1
    else
        idx((i-1)*N_OFDM+1:i*N_OFDM,1) = 0; % if vacant: 0
    end
end
f_f = zeros(1,4096);
f_f(1:410) = 1;
f_f(4096:3687) = 1;
f_t = ifft(f_f);
figure;
num_t = 1000*W;
delta_t = 1/num_t;
t = 0:delta_t:1-delta_t;
f_t = t*0;
y_n = f_t(1:L:num_t+1-L);
for i = 1:W
f_t = f_t+y_n(i)*sinc(W*t-i);
end

plot(abs(f_t))
%Create sparse signal a_w
a_w = idx;
% power_e = sum(conj(a_w).*a_w)/(10^(SNR/10)*N);
% e = sqrt(0.5*power_e)*randn(N,1)+sqrt(0.5*power_e)*1i*randn(N,1);
% a_w = a_w+e;
[ y,y1,Phi ] = Time_RandomDemodulator( a_w,R );