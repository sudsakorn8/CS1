clc;clear
%%  1. 时域测试信号生成
K=7;      %  稀疏度(做FFT可以看出来)
N=256;    %  信号长度
M=128;     %  测量数
f1=5e6;    %  信号频率1
f2=10e6;   %  信号频率2
f3=14e6;   %  信号频率3
f4=20e6;   %  信号频率4
fs=40e6;   %  采样频率
ts=1/fs;  %  采样间隔
Ts=1:N;   %  采样序列
snr_avgdB =10;%信噪比（dB值）
snr_avg=10^(snr_avgdB/10);%信噪比
Dn=1;         %噪声功率，高斯噪声的方差
Ds=snr_avg*Dn;%信号功率
x=sqrt(Ds)*cos(2*pi*f1*Ts*ts)+sqrt(Ds)*cos(2*pi*f2*Ts*ts)+sqrt(Ds)*cos(2*pi*f3*Ts*ts)+sqrt(Ds)*cos(2*pi*f4*Ts*ts);  %  完整信号
noi1 =sqrt(Dn)* randn(1,N);
n=0:511;f=40000000*n/512;
y = noi1 + x;
X=fft(x,512)
figure(1)
plot(f,abs(X))
title('信号的频谱')
%%  2.  时域信号压缩传感
Phi=randn(M,N);                                   %  测量矩阵(高斯分布白噪声)
s=Phi*y.'                                       %  获得线性测量
P=fft(s,128)
% h=0:127;
% f=40000000*h/128;
%figure(2);
%plot(f,abs(P))
%%  3.  正交匹配追踪法重构信号(本质上是L_1范数最优化问题)
%在信号中找到一个与数据r_n最相关的小波1与测量阵T的内积最大值，在数据中去掉该小波的所有印记，重获残差矩阵r_n。
%不断重复知道我们能用小波标记“解释”收集到的data
Psi=fft(eye(N,N))/sqrt(N);                        %  傅里叶正变换矩阵
T=Phi*Psi';                                       %  恢复矩阵(测量矩阵*正交反变换矩阵)

hat_y=zeros(1,N);                                 %  待重构的谱域(变换域)向量                     
Aug_t=[];                                         %  支撑集
r_n=s;                                            %  残差值

for times=1:K;                                    %  迭代次数(有噪声的情况下,该迭代次数为K)
    for col=1:N;                                  %  恢复矩阵的所有列向量
        product(col)=abs(T(:,col)'*r_n);          %  恢复矩阵的列向量和残差的投影系数(内积值)
    end
    [val,pos]=max(product);                       %  最大投影系数对应的位置，返回最大值的同时，返回一个下标向量
    Aug_t=[Aug_t,T(:,pos)];                       %  更新支撑集 ,加入关联最大的小波  
    T(:,pos)=zeros(M,1);                          %  将恢复矩阵该列置0，即去除印迹
    aug_y=(Aug_t'*Aug_t)^(-1)*Aug_t'*s;           %  最小二乘,使残差最小
    r_n=s-Aug_t*aug_y;                            %  残差
    pos_array(times)=pos;                         %  纪录各次最相关小波的位置，用于恢复
end
hat_y(pos_array)=aug_y                           %  重构的谱域向量
hat_x=real(Psi'*hat_y.')                        %  做逆傅里叶变换重构得到时域信号，real为取负数的实部，返回每个元素的实部
%%  4.  恢复信号和原始信号对比
figure(3);
hold on;
figure;
f_r = fs*(1:N)/N;
plot(f_r,abs(hat_y))
title('恢复信号的频谱')
plot(hat_x,'*')                                 %  重建信号
plot(y,'k-o')                                       %  原始信号
legend('Recovery','Original')
xlabel('信道信号')
ylabel('信号大小')
norm(hat_x.'-x)/norm(x)                           %  重构误差
%%  Detection Process   检测过程
%感知时间
T=0.0000025;
% 信号带宽W
W=4*10^7;
m =2*T*W;
Sim_Times=10000;%仿真次数
Pf=0:0.05:1-0.05;
h=length(Pf);
L=zeros(1,length(Pf));
Pd=zeros(1,length(Pf));

for i=1:length(Pf)

    for kk = 1:Sim_Times

        power=sum(abs(hat_x).^2);
        Threshold=qfuncinv(Pf(i))*sqrt(2*m)+m;%给定的pf反推的Th
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
xlabel('虚警概率Pf');
ylabel('检测概率Pd');
legend('仿真值','理论值');
title('单用户能量检测ROC曲线图')
