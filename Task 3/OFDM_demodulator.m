% Переход OFDM-symbol в частотную область
function TX_IQ=OFDM_demodulator(OFDM_time_guarded, T_guard) 
    
    % Удалть cp 
    OFDM_time = OFDM_time_guarded(T_guard+1:end,:); % take the last samples without CP samples

    % FFT
    TX_IQ = fft(OFDM_time,[],1);

end