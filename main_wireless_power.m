clear all; close all; clc;

fprintf('=== 无线能量传输系统仿真 ===\n');
fprintf('开始初始化系统...\n');


wpt = WirelessPowerSystem();

% 显示系统参数
wpt.displaySystemParameters();

% 绘制频率响应
fprintf('\n绘制频率响应曲线...\n');
figure(1);
wpt.plotFrequencyResponse(80e3, 120e3);
title('无线能量传输系统 - 频率特性');

% 绘制时域波形
fprintf('绘制时域波形...\n');
figure(2);
wpt.plotTimeDomainWaveforms();
sgtitle('无线能量传输系统 - 时域波形');

% 参数扫描：不同耦合系数的影响
fprintf('分析耦合系数影响...\n');
plotCouplingEffect(wpt);

fprintf('\n=== 仿真完成 ===\n');
fprintf('所有图形已生成，请查看弹出的图形窗口。\n');

% 耦合系数影响分析函数
function plotCouplingEffect(wpt)
    k_range = linspace(0.1, 0.9, 20);
    efficiency_k = zeros(1, length(k_range));
    P_out_k = zeros(1, length(k_range));
    
    original_M = wpt.M;
    
    for i = 1:length(k_range)
        k = k_range(i);
        wpt.M = k * sqrt(wpt.L1 * wpt.L2);
        wpt.calculateCompensationCapacitors();
        [~, P_out, efficiency] = wpt.calculatePower();
        efficiency_k(i) = efficiency;
        P_out_k(i) = P_out;
    end
    
    % 恢复原始M值
    wpt.M = original_M;
    wpt.calculateCompensationCapacitors();
    
    figure('Position', [100, 100, 800, 600]);
    yyaxis left
    plot(k_range, efficiency_k, 'b-o', 'LineWidth', 2, 'MarkerSize', 4);
    ylabel('效率 (%)');
    yyaxis right
    plot(k_range, P_out_k, 'r-s', 'LineWidth', 2, 'MarkerSize', 4);
    ylabel('输出功率 (W)');
    xlabel('耦合系数 k');
    title('耦合系数对系统性能的影响');
    grid on;
    legend('效率', '输出功率', 'Location', 'best');
end
