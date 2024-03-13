clc
clear
%% 堆芯冷却剂出口温度计算

t_fin = 279.4;  %堆芯冷却剂入口温度(℃)
Fa = 0.974;     %燃料元件发热占总发热的份额
Nt = 2700e+6;   %堆芯输出功率(W)
Pres = 15.5e3;  %反应堆运行压力(kPa)
Wt = 14314;     %冷却剂流量(kg/s)
Bf = 0.059;     %旁流系数(Bypass factor)

% 冷却剂出口温度求解
t_fout = 323;   %设定初始冷却剂出口温度(℃)
e_fout = 1;     %设定初始误差限(error)
while  abs(e_fout)>0.001
    tm = 0.5*(t_fout+t_fin);
    Cp = refpropm('C','T',tm+273.15,'P',Pres,'WATER');
    tmp_f = t_fin+Fa*Nt/(Wt*(1-Bf)*Cp);
    e_fout = (t_fout-tmp_f)/t_fout;
    t_fout=tmp_f; 
end

%% 热流密度计算

m = 157;        %燃料组件数
n = 264;        %单个燃料组件燃料棒数
N_rods = m*n;   %燃料棒总数
d_cs = 0.0095;  %燃料元件包壳外径(m)
GP = 12.6e-3;   %燃料棒直径/栅距(m)
Hv = 4.2672;    %堆芯高度(m)

F_qN = 2.524;   %热流量核热点因子
F_qE = 1.03;    %热流量工程热点因子
F_dHE = 1.085;  %焓升核热管因子

q_ = Nt*Fa/(pi*d_cs*Hv*N_rods);  %平均热流密度(W/m2)
q_max = q_*F_qN*F_qE;            %最大热流密度(W/m2)
ql_ = q_*pi*d_cs;                %平均线功率(W/m)
ql_max = ql_*F_qN*F_qE;          %最大线功率(W/m)

%% 控制体计算

N_cv = 20;      %控制体数目

% 修正功率归一化因子读取
mPNFtable = readtable('修正的功率归一化因子.xlsx');
mPNFmatrix = table2array(mPNFtable);
Phi = mPNFmatrix(:,2);

ANStable = zeros(8,N_cv);  %定义结果表
ANStable(2,1) = t_fin;     %入口温度赋值(℃)

for i = 1:1:N_cv
    [t_f2, DNBR, t_cs, t_ci, t_u, t_o, H_f] = CTRLvolume(ANStable(2,i), Phi(i,1), N_cv, q_max, ql_max);
    ANStable(1,i) = i;
    ANStable(2,i) = t_f2;
    ANStable(3,i) = DNBR;
    ANStable(4,i) = t_cs;
    ANStable(5,i) = t_ci;
    ANStable(6,i) = t_u;
    ANStable(7,i) = t_o;
    ANStable(8,i) = H_f;
    if i<N_cv
        ANStable(2,i+1) = t_f2;
    else
        ANStable(2,i) = t_f2;
    end
end

%% 热管最大焓升

H_in = refpropm('H','T',t_fin+273.15,'P',Pres,'water');            %热管冷却剂入口处焓值(J/kg)
H_out = refpropm('H','T',ANStable(2,20)+273.15,'P',Pres,'water');  %热管出口处冷却剂焓值(J/kg)

delH_fmax = H_out-H_in;  %热管最大焓升(J/kg)

%% 压降计算

g = 9.8;        %重力加速度(N/kg)
Kin = 0.75;     %入口局部阻力系数
Kout = 1.0;     %出口局部阻力系数
Kgr = 1.05;     %定位格架局部阻力系数

% 摩擦压降
tf_ = 0.5*(t_fin+t_fout);                               %冷却剂平均温度(℃)
Af = m*n*(GP^2-0.25*pi*d_cs^2);                         %总流通截面积(m^2)
De = 4*(GP^2-pi/4*d_cs^2)/(pi*d_cs);                    %单元通道当量直径(m)
rhof_ = refpropm('D','T',tf_+273.15,'P',Pres,'water');  %冷却剂平均密度(kg/m^3)
vf_ = Wt*(1-Bf)/(Af*rhof_);                             %冷却剂平均流速(m/s)
vis = refpropm('$','T',tf_+273.15,'P',Pres,'water');    %运动粘度(m^2/s)
Re = vf_*De/(0.0001*vis);                               %计算该处的雷诺数

f = 0.3164/Re^0.25;
dPf = f*Hv*rhof_*vf_^2/(2*De);
% 单相流体提升压降计算
dPel = rhof_*g*Hv;
% 进口局部压降计算
rho_in = refpropm('D','T',t_fin+273.15,'P',Pres,'water');   %冷却剂平均密度(kg/m^3)
vf_in = Wt*(1-Bf)/(Af*rho_in);                              %冷却剂平均流速(m/s)
dPin = 0.5*Kin*rhof_*vf_in^2;
% 出口局部压降计算
rho_out = refpropm('D','T',t_fout+273.15,'P',Pres,'water'); %冷却剂平均密度(kg/m^3)
vf_out = Wt*(1-Bf)/(Af*rho_out);                            %冷却剂平均流速(m/s)
dPout = 0.5*Kout*rho_out*vf_out^2;
% 定位搁架出口压降计算
dPgr = 0.5*Kgr*rhof_*vf_^2;
% 总的压降计算
dP = dPf+dPel+dPin+dPout+dPgr;

%% 打印结果

disp('控制体计算结果如下表所示：');
disp(' ');
Nam = {'控制体编号';'节点出口温度(℃)';'DNBR';'包壳外表面温度(℃)';'包壳内表面温度(℃)';'燃料芯块表面温度(℃)';'燃料芯块中心温度(℃)';'控制体焓(J/kg)'};
DISPLAY = [Nam, num2cell(ANStable)]';
disp(DISPLAY);
disp(' ');
disp('堆芯计算结果如下所示：');
disp(' ');
disp(['堆芯冷却剂出口的温度为：',num2str(t_fout),'℃.']);
disp(['燃料棒表面平均热流密度为：',num2str(q_),'W/m2.']);
disp(['燃料棒表面最大热流密度为：',num2str(q_max),'W/m2.']);
disp(['燃料棒表面平均线功率为：',num2str(ql_),'W/m.']);
disp(['燃料棒表面最大线功率为：',num2str(ql_max),'W/m.']);
disp(['堆芯热管的最大焓升为：',num2str(delH_fmax),'J/K.']);
disp(['包壳外表面最高温度为：',num2str(max(ANStable(4,:))),'℃.']);
disp(['包壳内表面最高温度为：',num2str(max(ANStable(5,:))),'℃.']);
disp(['燃料芯块表面最高温度为：',num2str(max(ANStable(6,:))),'℃.']);
disp(['燃料芯块中心最高温度为：',num2str(max(ANStable(7,:))),'℃.']);
disp(['堆芯的总压降为：',num2str(dP),'Pa.']);

figure(1)

plot(ANStable(1,:),ANStable(5,:),'r^-');
xlabel('控制体编号');
ylabel('包壳内表面温度（℃）');
legend('包壳内表面温度');

figure(2)
plot(ANStable(1,:),ANStable(4,:),'r^-');
xlabel('控制体编号');
ylabel('包壳外表面温度（℃）');
legend('包壳外表面温度');

figure(3)
plot(ANStable(1,:),ANStable(7,:),'r^-');
xlabel('控制体编号');
ylabel('芯块中心温度（℃）');
legend('芯块中心温度');

figure(4)
plot(ANStable(1,:),ANStable(3,:),'b^-');
xlabel('控制体编号');
ylabel('DNBR');
legend('DNBR');

figure(5)
plot(ANStable(1,:),ANStable(8,:),'r^-');
xlabel('控制体编号');
ylabel('控制体出口焓值(J/kg)');
legend('控制体出口焓');
