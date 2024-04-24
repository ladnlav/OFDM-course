function equalized_Hest = equalize_signal(OFDM_demod, Hest,N_carrier)
    % Perform equalization by element-wise division
    equalized_Hest = zeros(size(OFDM_demod));

    for i=1:size(OFDM_demod,2)
        equalized_Hest(1:N_carrier,i) = OFDM_demod(1:N_carrier,i) ./ Hest(1:N_carrier).';
    end
end
