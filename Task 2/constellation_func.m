%% constellation_func
% Make the different dictionary for BPSK, QPSK, 8PSK, 16QAM constellations
% calculate the Bit_depth for each contellation
function [Dictionary, Bit_depth_Dict] = constellation_func(Constellation)
    switch Constellation
        case 'BPSK'
            Dictionary = [-1+0i 1+0i];
            Bit_depth_Dict = 1;
        case 'QPSK'
            Dictionary = [-1-1i -1+1i 1-1i 1+1i];
            Bit_depth_Dict = 2;
        case '8PSK'
            gray_map = [5 4 2 3 6 7 1 0];
            Dictionary = exp(1i*(gray_map*2*pi/8));
            Bit_depth_Dict = log2(8);
        case '16QAM'
            Dictionary = [-3+3i, -3+1i, -3-3i, -3-1i, -1+3i, -1+1i, -1-3i, ... 
                -1-1i, 3+3i, 3+1i, 3-3i, 3-1i, 1+3i, 1+1i, 1-3i, 1-1i];
            Bit_depth_Dict = 4;
    end
 % Normalise the constellation.
    % Mean power of every constellation must be equel 1.
    % Make the function to calculate the norm, 
    % which can be applied for every constellation

% Make the function to calculate the norm, which can be applied for every constellation
 N=size(Dictionary,2);
 norm = sqrt(sum(Dictionary .* conj(Dictionary))/N);
 Dictionary = Dictionary./norm;

%Mean power P
 % Calculate P
 %P = 1/N * sum(Dictionary .* conj(Dictionary));
 %disp(P);
end