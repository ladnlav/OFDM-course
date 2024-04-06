function [dsc_sequence, Register] = DeScrambler(Register, sequence)
    % Инициализация параметров
    %sequence_length = 2^12 - 1; % длина псевдослучайной последовательности
    dsc_sequence=zeros(size(sequence));

    % Процесс Дескрэмблинга
    for i=1:length(sequence)
        feedback = sequence(i); 
        symb_reg=array_xor(Register); % выходной бит
        dsc_sequence(i) = xor(symb_reg,sequence(i));
    
        Register = circshift(Register,1); % сдвиг регистра
        Register(1) = feedback; % обновление первого бита
    end

end

function xor_result = array_xor(Register)
    
    coeff_mask=[1 0 0 0 0 0 0 0 0 0 0 0 0 1 1]; % was 77141 (oct) in bits
    xor_in=Register(logical(coeff_mask(2:end)));

    xor_result = xor_in(1);
    for i = 2:length(xor_in)
        xor_result=xor(xor_result,xor_in(i));
    end

end