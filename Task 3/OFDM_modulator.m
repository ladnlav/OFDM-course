% Создать OFDM-symbol во времени
function OFDM_time_guarded=OFDM_modulator(OFDM_symbols, T_guard) 
    
    % IFFT
    OFDM_time = ifft(OFDM_symbols,[],1);
    
    % Добавить cp
    cp = OFDM_time(end-T_guard+1:end,:);  % take the last CP samples
    OFDM_time_guarded = [cp; OFDM_time];

end