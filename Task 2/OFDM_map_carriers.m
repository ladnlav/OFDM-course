% Распределить данные по поднесущим 
function mapped_carriers=OFDM_map_carriers(QAM_payload, N_symb, Nfft, dataCarriers, pilotCarriers,amp_pilots) 
    mapped_carriers = zeros(Nfft, N_symb);
    data =  reshape(QAM_payload.',[length(dataCarriers), N_symb]);
    
    mapped_carriers(dataCarriers,:)=data;
    pilotValues =zeros(1,length(pilotCarriers));
    pilotValues(1:2:end) = amp_pilots*exp(1i*0);
    pilotValues(2:2:end) = amp_pilots*exp(1i*pi);
    
    pilotValues = repmat(pilotValues', 1, 50);
    mapped_carriers(pilotCarriers,:)=pilotValues;
end