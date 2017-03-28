N = 1024;
y = zeros(1,2047);
ftx = zeros(1,N);
psd2 = zeros(1,N);
trial_num=10000;
for num =1:trial_num
 xn = randn(1,N); 
 y = y+xcorr(xn,'biased')*sum(xn.*xn); 
 ftx = fft(xn);
 psd2 = psd2+ftx.*conj(ftx);
end
y=y/trial_num;
psd2 = psd2/trial_num;
 subplot(3,1,1)
 stem(y);
 fty = fft(y);
 psd1 = abs(fty);
 subplot(3,1,2)
 plot(psd1);
%  axis([0 201 0 10e4]);
 subplot(3,1,3)
 plot(psd2);
%  
%  
%  figure(2) 
%  subplot(3,1,1) 
%  plot(xn);
%  axis([0 N -5 5])
%  [r,lag]=xcorr(xn,100,'biased'); 
% subplot(3,1,2) 
%   hndl=stem(lag,r);
%  set(hndl,'Marker','.') 
%  set(hndl,'MarkerSize',2);
%  ylabel('¡Á??¨¤??????R(m)');
%  tf=fft(r);
%  mx_x=abs(tf);
%  subplot(3,1,3) 
%   plot(mx_x);
%   axis([0 201 -5 5])
%   xlabel('?????¡Á') 


% close all;
% x=ones(1,65536);
% f=zeros(1,65536);
% for i=1:1
% y = awgn(x,0) ;
% e = y-x;
% % e = randn(1,4096);
% f = f+abs(fft(e));
% 
% end
% figure;
% plot(e);
% % figure;
% % plot(abs(f));
% figure;
% hist(e,100);
% figure;
% plot(f/100);
% % hist(abs(f),100);