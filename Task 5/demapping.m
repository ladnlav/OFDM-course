function [de_bits] = demapping(pad, IQ, Constellation)
% Make the different dictionary for BPSK, QPSK, 8PSK, 16QAM constellations
% calculate the Bit_depth for each contellation
[Dictionary, Bit_depth_Dict] = constellation_func(Constellation);

% Find the closest constellation point for each IQ value
dist = zeros(2^Bit_depth_Dict,length(IQ));
for idx = 1:(2^Bit_depth_Dict)
    dist(idx,:) = (real(IQ) - (real(Dictionary(idx)))).^2 + (imag(IQ) - (imag(Dictionary(idx)))).^2;
end

[~, idx] = min(dist, [], 1);

% Convert the indices to binary representation
de_bits = int2bit(idx-1, Bit_depth_Dict);

% Reshape the bits into a row vector
de_bits = reshape(de_bits, 1, []);


if pad~=-1
    de_bits=de_bits(1:end-pad);
end

end

