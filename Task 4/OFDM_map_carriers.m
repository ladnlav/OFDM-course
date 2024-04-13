% Распределить данные по поднесущим 
function mapped_carriers=OFDM_map_carriers(QAM_payload, N_symb, Nfft, dataCarriers,pilotCarriers, pilotValues) 
    mapped_carriers = zeros(Nfft, N_symb);
    data =  reshape(QAM_payload.',[length(dataCarriers), N_symb]);
    
    mapped_carriers(dataCarriers,:)=data;
    
    mapped_carriers(pilotCarriers,:)=pilotValues;
end