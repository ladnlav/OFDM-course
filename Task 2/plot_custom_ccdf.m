function plot_custom_ccdf(PAPRs)

%Plot CDF
numBins = 100;
binEdges = linspace(floor(min(PAPRs)), ceil(max(PAPRs)), numBins + 1);

figure();
[Values,edges]=hist(PAPRs,binEdges);
prob_func=Values/numel(PAPRs);
CCDF_sig=1-cumsum(prob_func);
semilogy(edges, CCDF_sig,'-o','DisplayName','Обычный сигнал', 'LineWidth', 2);
hold on;

% binEdges = linspace(floor(min(sc_PAPRs)), ceil(max(sc_PAPRs)), numBins + 1);
% [Values,edges]=hist(sc_PAPRs,binEdges);
% prob_func=Values/numel(sc_PAPRs);
% sc_CCDF=1-cumsum(prob_func);
% semilogy(edges, sc_CCDF,'-o','DisplayName', 'Сигнал с рандомизатором', 'LineWidth', 2);

hold on;
grid on;
xlabel('PAPR (dB)');
ylabel('CCDF');
set(gca, 'Fontsize', 20)
legend('Location','southeast');
title('CCDF(PAPR) в режиме скользящего окна');

end