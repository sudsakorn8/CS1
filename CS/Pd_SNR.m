% draw the detection probability Pd to SNR curve
clear,clc;
tic
Channel_num = 10;% should divides 2 because of conjugate symmetry of spectrum
Sparsity = 0.1; 
Pf = 0.01;
SNR = -14:2:0;
%SNR=-6:1:-5;
SNR_num = length(SNR);
% CS_ratio = [0.2;0.25];
CS_ratio = 0.39;
CS_num = length(CS_ratio);
trial_num = 20;
%trial_num = 1; % trial number for calculating each Pd
Pd = zeros(SNR_num,CS_num); % detection probability with different SNR and Compression ratio
% Pd_denoised = zeros(SNR_num,CS_num);
Pd_occupied = zeros(SNR_num,CS_num);
error = zeros(trial_num,2);
pd_trial = zeros(trial_num,2);
idx_Channel = zeros(Channel_num,SNR_num);
if_vacant = zeros(Channel_num,SNR_num);
E = zeros(Channel_num,SNR_num);
E_ori = zeros(Channel_num,SNR_num);
a_r = zeros(80*Channel_num,SNR_num);
lamda_d = zeros(1,SNR_num); 
fid1=['log','.txt'];
cRec=fopen(fid1,'w');
for n = 1:CS_num
    for i = 1:SNR_num
        for j = 1:trial_num
            fprintf('CS_ratio=%f,SNR=%f,trial_num=%d  \n',CS_ratio(n),SNR(i),j);
            [pd_trial(j,1),pd_trial(j,2),idx_Channel(:,i),if_vacant(:,i),E(:,i),E_ori(:,i),lamda_d(:,i),a_r(:,i)] = CSS_Denoised_CS( SNR(i),CS_ratio(n),Sparsity,Pf);
%             [error(j,1),error(j,2),a_r,a_r_d] = Denoised_CS( SNR(i),CS_ratio(n),Sparsity,Pf);
%             [error(j,1),error(j,2)] = CSS_Denoised_CS( SNR(i),CS_ratio(n),Sparsity,Pf);
        end
%         Pd(i,n) =  1 - sum(error(:,1))/(Channel_num*trial_num);
%         Pd_denoised(i,n) =  1 - sum(error(:,2))/(Channel_num*trial_num);
        Pd(i,n) = sum(pd_trial(:,1))/(trial_num);
        Pd_o(i,n) = sum(pd_trial(:,2))/(trial_num);
        %Pd_denoised(i,n) = sum(pd_trial(:,2))/(trial_num);
        fprintf(cRec,'CS_ratio=%f,SNR=%f,trial_num=%d  \n',CS_ratio(n),SNR(i),j);
    end
end
toc
fclose(cRec);
figure;
plot(SNR,1-Pd(:,1),'o-b'); hold on
plot(SNR,Pd_o(:,1),'square-b'); hold on
% plot(SNR,Pd(:,2),'square-b'); hold on
%plot(SNR,Pd_denoised(:,1),'*-g'); hold on
% plot(SNR,Pd_denoised(:,2),'square-g'); hold on
xlabel('SNR(dB)'); ylabel('P');
% legend('Traditional CS based SS 20%','Traditional CS based SS 25%','Denoised CS based SS 20%','Denoised CS based SS 25%','Location','SouthEast');
legend('P_f','P_d','Location','SouthEast');
title('Traditional CS based on SS 40%,recovered by CoSaMP' );

