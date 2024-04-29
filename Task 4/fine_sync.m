function sync_signal = fine_sync(rx_signal, pilotCarriers, pilotValues,time_desync, freq_desync)
tx_pilotValues = pilotValues;

rx_pilotValues = rx_signal(pilotCarriers,:);

deltak = pilotCarriers(2)-pilotCarriers(1);

taus = zeros(1,length(tx_pilotValues(:,1)));

% for i=2:length(tx_pilotValues(:,2))
% 
%     q_k_1 = tx_pilotValues(i-1,2)*conj(rx_pilotValues(i-1,2));
%     q_k = tx_pilotValues(i,2)*conj(rx_pilotValues(i,2));
%     
%     taus(i-1) = angle(q_k*conj(q_k_1))/(2*pi*deltak);
% end
% 
% figure(); plot(pilotCarriers,(taus),'LineWidth', 1.5);
% title("\tau - временная задержка при SNR = 25 dB для эффективной полосы одного ofdm-символа")
% xlabel("k - номер поднесущей")
% ylabel("\tau - значение временной задержки")
% grid on;
% set(gca, 'Fontsize', 20);

for i=2:length(tx_pilotValues(:))
    q_k_1 = tx_pilotValues(i-1)*conj(rx_pilotValues(i-1));
    q_k = tx_pilotValues(i)*conj(rx_pilotValues(i));
    
    taus(i-1) = angle(q_k*conj(q_k_1))/(2*pi*deltak);
end

diffs = diff(taus);
mask = [false, (abs(diffs) < 1e-3 & abs(diffs)~=0)];
taus_result = taus(mask);
tau = mean(taus_result(length(pilotCarriers)+1:end)); 

% Тонкая временная синхронизация
if time_desync
    nn = 0:1024-1;
    nn_exp = exp(-2j*pi*tau*nn);
    for i=1:size(rx_signal,2)
        rx_signal(:,i)=rx_signal(:,i).*nn_exp';
    end
end

% Тонкая фазовая синхронизация (компенсация постоянного фазового сдвига)
rx_pilotValues = rx_signal(pilotCarriers,:);
qks = zeros(1,length(tx_pilotValues(:)));
for i=1:length(tx_pilotValues(:))
    qks(i) = angle(tx_pilotValues(i)*conj(rx_pilotValues(i)));
end
phase_shift = mean(qks(abs(qks)>1e-3));

if freq_desync
    sync_signal = rx_signal.*exp(1j*phase_shift);
else
    sync_signal = rx_signal;
end

end

