% Рассчитать оконный PAPR
function PAPRs = calculate_window_PAPR(Tx_OFDM_Signal,Nfft)    
    % Сколько всего будет значений PAPR для данного сигнала
    n_PAPRs=length(Tx_OFDM_Signal) - Nfft + 1;
    PAPRs = zeros(1, n_PAPRs);

    % Цикл расчёта PAPR скользящим окном
    for i = 1:n_PAPRs
            % Извлечь отсчёты в окне 
            window = Tx_OFDM_Signal(i:i+Nfft-1);
            
            % Рассчитать PAPR для текущего окна
            PAPRs(i) = calculatePAPR(window);
    end
end