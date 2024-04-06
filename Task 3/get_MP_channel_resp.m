% Рассчитать канальную характеристику при многолучевом распространении
function [impulse_response,frequency_response] = get_MP_channel_resp(channel_taps, Nfft)
    % Вычислить максимальную задержку среди всех лучей
    max_delay = max(channel_taps(:, 1));
    impulse_response_length = max_delay + 1;

    % Инициализация массива с импульсной характеристикой канала
    impulse_response = zeros(1, impulse_response_length);

    % Создание массива с импульсной характеристикой канала
    for i = 1:size(channel_taps, 1)
        delay = channel_taps(i, 1);
        amplitude = channel_taps(i, 2);
        impulse_response(delay + 1) = amplitude;
    end

    % Вычислим частотную характеристику канала для будущих сравнений
    frequency_response = fft(impulse_response, Nfft);
end
