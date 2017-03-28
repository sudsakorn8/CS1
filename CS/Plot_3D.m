load('data.mat')
CS_ratio = [50,40,30,25,20,15,10,5];
SNR = -20:2:0;
Spa = 0.1:0.1:0.5;
CS_ratio_mesh = 0.1:0.1:0.5;
[x,y] = meshgrid(CS_ratio_mesh,Spa);
z(:,1)=SNR_minus2_spa_10.P_d;
z(:,2)=SNR_minus2_spa_20.P_d;
z(:,3)=SNR_minus2_spa_30.P_d;
z(:,4)=SNR_minus2_spa_40.P_d;
z(:,5)=SNR_minus2_spa_50.P_d;
figure;
% for CS_ratio_num = 1:length(CS_ratio)
%     var_str = ['struct_data_',num2str(CS_ratio(CS_ratio_num))];
%     eval(['temp=',var_str]);
%     SNR_num = length(temp.SNR);
% %     for SNR_num = 1:length(temp.SNR)
% %     scatter3(temp.SNR,CS_ratio(CS_ratio_num)*ones(1,SNR_num),temp.P_d);hold on;
%     
% end
surf(x,y,z);
title('Pd-Sparsity-CSratio based on Traditional CSS' );
x1=xlabel('Spa');
x2=ylabel('CS ratio');
zlabel('P_d');
set(x1,'Rotation',30);
set(x2,'Rotation',-30);
    
         
         