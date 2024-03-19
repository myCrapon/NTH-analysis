clc
clear
% 轴向功率归一化因子拟合(PNF)
Hv = 4.2672;    %堆芯高度(m)
PNFtable = readtable('轴向功率归一化因子分布.csv');
PNFmatrix = table2array(PNFtable);

px = PNFmatrix(:, 1);
py = PNFmatrix(:, 2);

[xData, yData] = prepareCurveData( px, py );

% 设置 fittype 和选项
ft = fittype( 'gauss7' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0];
opts.StartPoint = [1.186026198 0.557951246 0.0538033643079901 1.17273860505418 0.680498322 0.0556634616251109 1.1501407267418 0.402730884 0.0450506312682384 1.10652945899263 0.807075401 0.0743599571543443 1.08603781851748 0.267897939 0.0347157816290735 0.998732647008271 0.190933211 0.043679007815396 0.944393552165349 0.332633184 0.0426507016475911];

% 对数据进行模型拟合
[PNFfit, PNFgof] = fit( xData, yData, ft, opts );

% 生成功率归一化因子表
N_cv = 20;     %控制体数目(control volume number)
deltaH = Hv/N_cv;
Hi = 0;
sum = 0;
for i = 1:1:N_cv

    Hi = i*deltaH-0.5*deltaH;
    hi = Hi/Hv;
    tmpPhi = PNFfit(hi);
    Phi(i) = tmpPhi;
    sum = Phi(i)+sum;
end
plot(Phi)

%% 随后根据计算数值手动微调，使得功率归一化因子的和为N_cv.
