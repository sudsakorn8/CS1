function [ y,y1,y2,Phi,Omega_n ] = Time_RandomDemodulator( a_w,R )
% input a_w(W*1) is the original frequency cofficient
% R(1*1) is the sampling rate
% output y(R*1) is the compressive samples
W = length(a_w); % the whole bandwidth
%% Demodulation Matrix D(W*W)
chip_seq = rand(W,1);
chip_seq = ((chip_seq>0.5)-0.5)*2; %chipping sequence: take values +-1 with equal probability
D = diag(chip_seq);
%% Sampler Matrix H(R*W)
comp = W/R;
H = zeros(R,W);
for i = 1:R
    H(i,ceil(comp*(i-1))+1:floor(comp*i)) = 1;
    if i>1
        H(i,ceil(comp*(i-1)))=ceil(comp*(i-1))-comp*(i-1);
    end
    if i<R
        H(i,floor(comp*i)+1)=comp*(i)-floor(comp*i);
    end
end
M = H*D; 
% DFT Matrix F
n = 1:W;
w = 1-W/2:W/2;
F = exp(-2*pi*1i*n'*w/W);%1/sqrt(W)*(

% Overall action Matrix on s Phi(R*W)
Phi = M*F;% the action of the hardware system Matrix on x M(R*W)

%Create amplitude vector u(with phase information)
u = zeros(W,1);
N=W;
for k = 1:N
    u(k) = a_w(k)*((exp(-2*pi*1i*(k-N/2)/N)-1)/(2*pi*1i*(k-N/2)));
end
u(N/2) = a_w(N/2)/N; 
y1 = Phi*u;%for comparation
%% time_random_modulator
num_t = 100*W;
delta_t = 1/num_t;
t = 0:delta_t:1-delta_t;
f_t = t*0;
% generate input time-sequence
Omega_n = zeros(1,W);
for w_n = -W/2+1:W/2 % integer-valued frequncies
%     f_t = f_t+a_w(w_n+W/2)*(exp(-2*pi*1i*(w_n)*t));
    Omega_n(w_n+W/2)=w_n+((rand(1,1)>0.5)-0.5)*rand(1,1);
    f_t = f_t+a_w(w_n+W/2)*(exp(-2*pi*1i*(Omega_n(w_n+W/2))*t));
end
% f_t = 1*exp(-2*pi*1i*33.5*t)+2*exp(-2*pi*1i*88.8*t)+3*exp(-2*pi*1i*156*t)+3*exp(2*pi*1i*255.3*t);

% generate random time-sequence 
L = num_t/W;
% DFT
% y_n = f_t(1:L:num_t+1-L);
% fft_y = fft(y_n,4096);
% figure;
% plot(abs(fft_y));
% title('DFT of the original signal in Nyquist sampling' );

pc_t = t*0;
for num = 1:W
    pc_t((num-1)*L+1:num*L) = chip_seq(num);
end
% mix
f_mix = f_t.*pc_t;
%integral
L1=num_t/R;
y = (1:R)'*0;
for i = 1:R
    y(i) = sum(f_mix(round((i-1)*L1)+1:round(i*L1)))*delta_t;
end

%% for continues frequency signals
delta_w = 0.1;
w = -399:delta_w:400;
num_w = length(w);
F_w = w*0;
F_w(399*num_w/800:450*num_w/800) = 1;
delta_t1 = 1/num_t;
f_t1 = (1:1/delta_t1)*0;
for t_num = 1:1/delta_t1;
    f_t1(t_num) = sum(F_w.*exp(1j*2*pi*w*t_num*delta_t1))*delta_w;  
end
figure;
% plot(abs(fft(y_n)));
plot(abs(fft(f_t1(1:L:num_t+1-L))));
%% interpolation
% f_t1 = f_t*0;
% % y_n = f_t(1:L:num_t+1-L);
% for i = 1:W
% f_t1 = f_t1+y_n(i)*sinc(W*t-i);
% end
% f_t1 = f_t;

% figure;
% plot(abs(f_t));
% figure;
% plot(abs(f_t1));
% mix
f_mix1 = f_t1.*pc_t;
%% integral
L1=num_t/R;
y2 = (1:R)'*0;
for i = 1:R
    y2(i) = sum(f_mix1(round((i-1)*L1)+1:round(i*L1)))*delta_t;
end

aa=1;
% figure;
% plot(abs(y));
% figure;
% plot(abs(y1));
% figure;
% plot(abs(y./y1));