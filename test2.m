N=120;
P_d_th = [];
for SNR = -20:2:0
power_e = 1/(10^((SNR+10)/10)*N*10);
lamda_d = power_e*(1+qfuncinv(0.01)/sqrt(N/2));
P_d_th = [P_d_th;qfunc(sqrt(N/2)*(lamda_d/power_e-(1+10^((SNR+10)/10)))/(1+10^((SNR+10)/10)))];
end
figure;
plot(-20:2:0,P_d_th);

N=80;
P_d_th = [];
for SNR = -20:2:0
power_e = 1/(10^((SNR)/10)*N*40);
lamda_d = power_e*(1+qfuncinv(0.01)/sqrt(N/2));
P_d_th = [P_d_th;qfunc(sqrt(N/2)*(lamda_d/power_e-(1+10^((SNR)/10)))/(1+10^((SNR)/10)))];
end
figure;
plot(-20:2:0,P_d_th);