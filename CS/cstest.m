clc;clear
%%  1. ʱ������ź�����
K=7;      %  ϡ���(��FFT���Կ�����)
N=256;    %  �źų���
M=128;     %  ������
f1=5e6;    %  �ź�Ƶ��1
f2=10e6;   %  �ź�Ƶ��2
f3=14e6;   %  �ź�Ƶ��3
f4=20e6;   %  �ź�Ƶ��4
fs=40e6;   %  ����Ƶ��
ts=1/fs;  %  �������
Ts=1:N;   %  ��������
snr_avgdB =10;%����ȣ�dBֵ��
snr_avg=10^(snr_avgdB/10);%�����
Dn=1;         %�������ʣ���˹�����ķ���
Ds=snr_avg*Dn;%�źŹ���
x=sqrt(Ds)*cos(2*pi*f1*Ts*ts)+sqrt(Ds)*cos(2*pi*f2*Ts*ts)+sqrt(Ds)*cos(2*pi*f3*Ts*ts)+sqrt(Ds)*cos(2*pi*f4*Ts*ts);  %  �����ź�
noi1 =sqrt(Dn)* randn(1,N);
n=0:511;f=40000000*n/512;
y = noi1 + x;
X=fft(x,512)
figure(1)
plot(f,abs(X))
title('�źŵ�Ƶ��')
%%  2.  ʱ���ź�ѹ������
Phi=randn(M,N);                                   %  ��������(��˹�ֲ�������)
s=Phi*y.'                                       %  ������Բ���
P=fft(s,128)
% h=0:127;
% f=40000000*h/128;
%figure(2);
%plot(f,abs(P))
%%  3.  ����ƥ��׷�ٷ��ع��ź�(��������L_1�������Ż�����)
%���ź����ҵ�һ��������r_n����ص�С��1�������T���ڻ����ֵ����������ȥ����С��������ӡ�ǣ��ػ�в����r_n��
%�����ظ�֪����������С����ǡ����͡��ռ�����data
Psi=fft(eye(N,N))/sqrt(N);                        %  ����Ҷ���任����
T=Phi*Psi';                                       %  �ָ�����(��������*�������任����)

hat_y=zeros(1,N);                                 %  ���ع�������(�任��)����                     
Aug_t=[];                                         %  ֧�ż�
r_n=s;                                            %  �в�ֵ

for times=1:K;                                    %  ��������(�������������,�õ�������ΪK)
    for col=1:N;                                  %  �ָ����������������
        product(col)=abs(T(:,col)'*r_n);          %  �ָ�������������Ͳв��ͶӰϵ��(�ڻ�ֵ)
    end
    [val,pos]=max(product);                       %  ���ͶӰϵ����Ӧ��λ�ã��������ֵ��ͬʱ������һ���±�����
    Aug_t=[Aug_t,T(:,pos)];                       %  ����֧�ż� ,�����������С��  
    T(:,pos)=zeros(M,1);                          %  ���ָ����������0����ȥ��ӡ��
    aug_y=(Aug_t'*Aug_t)^(-1)*Aug_t'*s;           %  ��С����,ʹ�в���С
    r_n=s-Aug_t*aug_y;                            %  �в�
    pos_array(times)=pos;                         %  ��¼���������С����λ�ã����ڻָ�
end
hat_y(pos_array)=aug_y                           %  �ع�����������
hat_x=real(Psi'*hat_y.')                        %  ���渵��Ҷ�任�ع��õ�ʱ���źţ�realΪȡ������ʵ��������ÿ��Ԫ�ص�ʵ��
%%  4.  �ָ��źź�ԭʼ�źŶԱ�
figure(3);
hold on;
figure;
f_r = fs*(1:N)/N;
plot(f_r,abs(hat_y))
title('�ָ��źŵ�Ƶ��')
plot(hat_x,'*')                                 %  �ؽ��ź�
plot(y,'k-o')                                       %  ԭʼ�ź�
legend('Recovery','Original')
xlabel('�ŵ��ź�')
ylabel('�źŴ�С')
norm(hat_x.'-x)/norm(x)                           %  �ع����
%%  Detection Process   ������
%��֪ʱ��
T=0.0000025;
% �źŴ���W
W=4*10^7;
m =2*T*W;
Sim_Times=10000;%�������
Pf=0:0.05:1-0.05;
h=length(Pf);
L=zeros(1,length(Pf));
Pd=zeros(1,length(Pf));

for i=1:length(Pf)

    for kk = 1:Sim_Times

        power=sum(abs(hat_x).^2);
        Threshold=qfuncinv(Pf(i))*sqrt(2*m)+m;%������pf���Ƶ�Th
        if power>Threshold

            L(i)= L(i) +1;
        end
    end
    Pd(i)=L (i)/Sim_Times;
end
Threshold_matrix=qfuncinv(Pf)*sqrt(2*m)+m;
Pd_theory=qfunc((Threshold_matrix-m-snr_avg)/(sqrt(2*m+snr_avg)));
figure(4)
plot(Pf,Pd,'-*r')
hold on
plot(Pf,Pd_theory,'-*b')
hold off
xlabel('�龯����Pf');
ylabel('������Pd');
legend('����ֵ','����ֵ');
title('���û��������ROC����ͼ')
