%% 控制体计算函数

function [t_f2, DNBR, t_cs, t_ci,t_u, t_o, H_f] = CTRLvolume(t_f1, Phi, N_cv, q_max, q_lmax)

% 常数定义
t_fin = 279.4;  %冷却剂入口平均温度
t_fout = 315.6; %冷却剂出口平均温度
Pres = 15.5e3;  %函数内定义运行压力(kPa)
Wt = 14314;     %冷却剂流量(kg/s)
Bf = 0.059;     %旁流系数(Bypass factor)
Hv = 4.2672;    %堆芯高度(m)
d_cs = 9.5e-3;  %燃料元件包壳外径(m)
d_ci = 8.36e-3; %包壳内表面直径
d_u = 8.19e-3;  %燃料芯块直径
GP = 12.6e-3;   %燃料棒直径/栅距(m)
Kg = 5678;      %包壳和芯块间气隙等效放热系数
delH = Hv/N_cv; %控制体单元步长
m = 157;        %燃料组件数
n = 264;        %单个燃料组件燃料棒数
n0 = 17;        %燃料组件为17x17规格

% 平均管情况
Af = m*n*(GP^2-0.25*pi*d_cs^2)+4*n0*GP*(GP-d_cs)/2;     %总流通截面积
tf_ = 0.5*(t_fin+t_fout);                               %冷却剂平均温度
rhof_ = refpropm('D','T',tf_+273.15,'P',Pres,'water');  %冷却剂平均密度
v_ = Wt*(1-Bf)/(Af*rhof_);                              %冷却剂平均流速
Au = GP^2-0.25*pi*d_cs^2;                               %单元流通截面积
Wu = Wt*(1-Bf)*Au/Af;                                   %单元截面流量
De = 4*(GP^2-pi/4*d_cs^2)/(pi*d_cs);                    %单元通道当量直径

% 二氧化铀积分热导率参考表拟合公式
UO2table = readtable('二氧化铀积分热导率参考表.xlsx');
UO2matrix = table2array(UO2table);

xData = UO2matrix(:,1);
yData = UO2matrix(:,2);
ft = fittype('poly5');
[UO2fit, ~] = fit(xData, yData, ft);
[UO2fit_v, ~] = fit(yData, xData, ft);

% 求出该控制体出口处的温度
t_f2 = 380;     %设置初始控制体出口温度
e_cv1 = 0.1;    %设置出口误差限
while  e_cv1>0.001
    ti = 0.5*(t_f1+t_f2);
    Cpi = refpropm('C','T',ti+273.15,'P',Pres,'WATER');
    tmpi = t_f1+q_max*Phi*pi*d_cs*delH/(Wu*Cpi);
    e_cv1 = (tmpi-t_f2)/t_f2;
    t_f2 = tmpi;
end

% 求解控制体焓升

H_f = refpropm('H','T',t_f2+273.15,'P',Pres,'water');        %热管出口处冷却剂焓值

% 求解临界热流量与烧毁比
Hfin = refpropm('H','T',t_f1+273.15,'P',Pres,'water'); %控制体入口冷却剂的比焓
Hfs = refpropm('H','P',Pres,'Q',0,'water');            %运行压力下的饱和水比焓
Hgs = refpropm('H','P',Pres,'Q',1,'water');            %运行压力下的饱和蒸汽比焓
H = refpropm('H','T',t_f2+273.15,'P',Pres,'water');    %控制体出口处冷却剂比焓
x_e = (H-Hfs)/(Hgs-Hfs);                               %该控制体处含汽量
G = rhof_*v_;    %冷却剂质量流速(kg/m^2*s)
p = Pres*1000;   %换单位
% 根据W-3公式计算出临界热流量
qDNB = 3.154e6*((2.022-6.238e-8*p)+(0.1722-1.43e-8*p)*exp((18.177- 5.987e-7*p)*x_e))*((0.1484-1.596*x_e+0.1729*x_e*abs(x_e))*(737.64*G/10e6)+1.037)*(1.157-0.869*x_e)*(0.2664+0.8357*exp(-124*De))*(0.8258+0.341e-6*(Hfs-Hfin));
% 计算烧毁比
DNBR = qDNB/(q_max*Phi); 

% 查找该压力温度下冷却剂的热物性
vis = refpropm('$','T',t_f2+273.15,'P',Pres,'water');   %运动粘度
lamda = refpropm('L','T',t_f2+273.15,'P',Pres,'water'); %热导率
Pr = refpropm('^','T',t_f2+273.15,'P',Pres,'water');    %普朗特数
t_s = refpropm('T','P',Pres,'Q',0,'water')-273.15;      %冷却剂饱和温度(℃)

Re = v_*De/(0.0001*vis);                                %计算该处的雷诺数
h = 0.023*(Re^0.8)*(Pr^0.4)*lamda/De;                   %该处的对流换热系数

% 单相强迫对流放热公式算得的温压
delT1 = q_max*Phi/h;
% 采用詹斯-洛特斯传热方程算得的过冷沸腾膜温压
delT2 = 25*((q_max*Phi/10^6)^0.25)*exp(-Pres/6.2)+t_s-t_f2;
% 膜温压取两个中较小值，算得包壳外表面温度
if delT1<delT2
    t_cs = t_f2+delT1;
else 
    t_cs = t_f2+delT2;
end

% 计算包壳内表面温度
t_ci = 349;     %设定初始包壳内表面温度
e_cv2 = 0.1;    %设定初始误差限
while e_cv2>0.001
    tc_ = 0.5*(t_ci+t_cs);
    Kci = 0.0547*(1.8*tc_+32)+13.8;
    tmptci = t_cs+q_lmax*Phi*log(d_cs/d_ci)/(2*pi*Kci);
    e_cv2 = (tmptci-t_ci)/tmptci;
    t_ci = tmptci;                          %采用迭代算法求得包壳内表面温度
end

% 计算燃料芯块中心温度
t_u = t_ci+q_lmax*Phi*2/(pi*(d_ci+d_u)*Kg); %燃料芯块表面温度
I_tu = UO2fit(t_u);                         %燃料芯块表面积分热导率
dI_tou = q_lmax*Phi/(4*pi);
I_to = I_tu+dI_tou;
t_o = UO2fit_v(I_to);                       %根据积分热导率拟合公式反函数求解芯块中心温度

end