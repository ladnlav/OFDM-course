% Извлечь данные из OFDM-символов
function RX_IQ=get_payload(RX_OFDM_symbols, dataCarriers) 
    RX_IQ=RX_OFDM_symbols(dataCarriers,:);
end