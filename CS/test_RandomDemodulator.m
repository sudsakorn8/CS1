function [ y,Phi ] = test_RandomDemodulator( s,R )
% RANDOMDEMODULATOR Compress the signal s using Random Demodulator
% Based on algorithm proposed in "Beyond Nyquist: Efficient Sampling of Sparse Bandlimited Signals"
% input s(W*1) is the original signal in frequency domain
% input R(1*1) is the sampling rate. Specially R divided W
% input e(R*1) is the noise signal
% output y(R*1) is the compressive samples
% output Phi(R*W) is the measurement matrix
% y = Phi*s+e = HD*x+e = HDF*s+e

W = length(s); % the whole bandwidth

% Demodulation Matrix D(W*W)
chip_seq = rand(W,1);
chip_seq = ((chip_seq>0.5)-0.5)*2; %chipping sequence: take values +-1 with equal probability
D = diag(chip_seq);

% Sampler Matrix H(R*W)
comp = W/R;
H = zeros(R,W);
% rand_W = rand(W,1);
% [~,idx_sort] = sort(rand_W,'ascend');
% [idx_ss,~] = sort(idx_sort(1:R),'ascend');
for i = 1:R
    H(i,i*comp) = 1;
end
% for i = 1:R

%     H(i,ceil(comp*(i-1))+1:floor(comp*i)) = 1;
%     if i>1
%         H(i,ceil(comp*(i-1)))=ceil(comp*(i-1))-comp*(i-1);
%     end
%     if i<R
%         H(i,floor(comp*i)+1)=comp*(i)-floor(comp*i);
%     end
% end


M = H;%*D; % the action of the hardware system Matrix on x M(R*W)

% DFT Matrix F
n = 1:W;
w = 1-W/2:W/2;
F = (exp(-2*pi*1i*n'*w/W));%1/sqrt(W)*

% Overall action Matrix on s Phi(R*W)
Phi = M*F;

% Vector of samples with noise y(R*1)
y = Phi*s;
end

