% Расчитать PAPR
function PAPR = calculatePAPR(OFDM_signal)    
    % Поиск квадрата максимальной амплитуды
    peak_power = max(abs(OFDM_signal))^2;
    
    % Нахождение средней мощности сигнала
    average_power = mean(abs(OFDM_signal).^2);
    
    % Расчёт PAPR
    PAPR = 10*log10(peak_power / average_power);
end