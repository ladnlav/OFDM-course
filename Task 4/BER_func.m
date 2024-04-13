function [BER] = BER_func(Bit_Tx, Bit_Rx)
    % Calculate the number of errors in the stream
    numErrors = sum(Bit_Tx ~= Bit_Rx);
    
    % Calculate the probability of an error (Bit Error Rate)
    BER = numErrors / length(Bit_Tx);
end

