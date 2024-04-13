function [IQ_RX, N_var] = Noise(SNR, IQ_TX)
    % calculate power of the signal
    IQ_TXPower = mean(abs(IQ_TX).^2);
    % calculate power of the noise
    NoisePower = IQ_TXPower / (10^(SNR/10));
    % generating white noise with the certain SNR
    Noise = sqrt(NoisePower/2) * normrnd(0, 1, size(IQ_TX)) + ...
            1i * sqrt(NoisePower/2) * normrnd(0, 1, size(IQ_TX));
    % adding noise to the signal
    IQ_RX = IQ_TX + Noise;
    N_var=sqrt(NoisePower);
end