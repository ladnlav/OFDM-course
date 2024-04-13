function equalized_Hest = equalize_signal(OFDM_demod, Hest)
    % Perform equalization by element-wise division
    equalized_Hest = zeros(size(OFDM_demod));

    for i=1:size(OFDM_demod,2)
        equalized_Hest(:,i) = OFDM_demod(:,i) ./ Hest.';
    end
end
