load('data.mat')
CS_ratio = [50,40,30,25,20,15,10,5];
SNR = -20:2:0;
CS_ratio_mesh = 10:10:50;
[x,y] = meshgrid(SNR,CS_ratio_mesh);
z(1,:)=struct_data_10.P_d(1:2:21);
z(2,:)=struct_data_20.P_d(1:2:21);
z(3,:)=struct_data_30.P_d;
z(4,:)=struct_data_40.P_d;
z(5,:)=struct_data_50.P_d;
% z(1,:)=struct_data_10.P_f(1:2:21);
% z(2,:)=struct_data_20.P_f(1:2:21);
% z(3,:)=struct_data_30.P_f;
% z(4,:)=struct_data_40.P_f;
% z(5,:)=struct_data_50.P_f;
figure;
for CS_ratio_num = 1:length(CS_ratio)
    var_str = ['struct_data_',num2str(CS_ratio(CS_ratio_num))];
    eval(['temp=',var_str]);
    SNR_num = length(temp.SNR);
%     for SNR_num = 1:length(temp.SNR)
%     scatter3(temp.SNR,CS_ratio(CS_ratio_num)*ones(1,SNR_num),temp.P_d);hold on;
    surf(x,y,z);
end
title('Pd-SNR-CSratio based on Traditional CSS' );
x1=xlabel('SNR(dB)');
x2=ylabel('CS ratio');
zlabel('P_d');
set(x1,'Rotation',30);
set(x2,'Rotation',-30);
    
         
         