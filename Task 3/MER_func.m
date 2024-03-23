function [MER] = MER_func(IQ_RX, Constellation)

[Dictionary, ~] = constellation_func(Constellation);

%finding "ideal" point for each rx point
IQ_RX_ideal = zeros(size(IQ_RX));
for i = 1:length(IQ_RX)
    min_distance = inf;
    min_index = 0;
    for j = 1:length(Dictionary)
        distance = abs(IQ_RX(i) - Dictionary(j));
        if distance < min_distance
            min_distance = distance;
            min_index = j;
        end
    end
    IQ_RX_ideal(i) = Dictionary(min_index);
end

%finding Power of the "Ideal" IQ-points
sum1=sum(real(IQ_RX_ideal).^2+imag(IQ_RX_ideal).^2);
%finding Noise Power 
sum2=sum(real(IQ_RX_ideal-IQ_RX).^2+imag(IQ_RX_ideal-IQ_RX).^2);
%MER itself
MER = 10 * log10(sum1 / sum2);
end

