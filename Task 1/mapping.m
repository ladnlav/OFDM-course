function [IQ,pad] = mapping(bits, constellation)
    % Make the different dictionary for BPSK, QPSK, 8PSK, 16QAM constellations
    % calculate the Bit_depth for each contellation

    [dictionary, bit_depth] = constellation_func(constellation);

    pad = -1;
    remainder=mod(numel(bits),bit_depth);
    if remainder~=0
        pad=bit_depth-remainder;
        bits=vertcat(bits,zeros(pad,1));
    end

    % Reshape the bits into a matrix with bit_depth rows
    symbols_index = reshape(bits, bit_depth, []).';
    
    % Convert the rows of the matrix to decimal values
    symbols_index = bi2de(symbols_index, 'left-msb');
    
    % Use the decimal values as indices into the dictionary
    symbols = dictionary(symbols_index+1);

    % write  the function of mapping from bit vector to IQ vector
    IQ=symbols;
end