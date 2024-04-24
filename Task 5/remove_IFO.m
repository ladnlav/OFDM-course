function [fixed_rx_signal, IFO] = remove_IFO(rx_signal,Nfft)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

spectrum = ((abs(fft(rx_signal(Nfft+1:2*Nfft),[],1))));
inds = find(spectrum>0.77);

IFO = (inds(1)-1);
fixed_rx_signal = add_CFO(rx_signal,-IFO,Nfft);

end

