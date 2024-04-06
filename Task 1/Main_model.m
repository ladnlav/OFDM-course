clc; % чистка командного окна
close all; % закрыть дополнительные окна 
clear all; % очистить память
rng(1); % фиксирование начального состояния генератора случайных чисел Матлаба
%% 0. Параметры OFDM
Nfft = 1024;
N_carrier = 400;
T_Guard = Nfft /8;
Amount_OFDM_Frames = 10;
Amount_ODFM_SpF = 5;
N_symb=Amount_OFDM_Frames*Amount_ODFM_SpF;

% Параметры полосы
Percent_pilot = 25;

allCarriers = linspace(1,Nfft,Nfft); % индексы поднесущих
amount_pilots = round(Percent_pilot/100*N_carrier);   % кол-во пилотов на всех Nfft отсчётах
pilot_step = floor(N_carrier/amount_pilots);
pilotCarriers = allCarriers(1:pilot_step:N_carrier-2);
pilotCarriers = [pilotCarriers,allCarriers(N_carrier)];
amount_pilots = length(pilotCarriers);

dataCarriers = allCarriers(~ismember(allCarriers(1:N_carrier), pilotCarriers));
amount_data_carriers = length(dataCarriers);



File = 'eagle.tiff';
Constellation = "16QAM";
[dict,bps] = constellation_func(Constellation); % bps - bits per symbol - битов на один символ созвездия

%% 1. Чтение файла
Size_Buffer = Amount_ODFM_SpF*Amount_OFDM_Frames*amount_data_carriers*bps;
input_bits = file_reader(File, Size_Buffer);

%% 2. Mapping
[TX_IQ,pad] = mapping(input_bits,Constellation);

%% 3. Формирование полосы
max_amplitude = max(abs(dict));
amp_pilots = 2*max_amplitude ;        % пилотное значение

OFDM_mapped_carriers = OFDM_map_carriers(TX_IQ, N_symb, Nfft, dataCarriers, pilotCarriers,amp_pilots);

%% 4. OFDM-модуляция: переход во временную область + добавление CP
Tx_OFDM_Signal_matrix = OFDM_modulator(OFDM_mapped_carriers, T_Guard); % lab 2

Tx_OFDM_Signal=Tx_OFDM_Signal_matrix(:);

%% 5.0 Временная рассинхронизация
nSTO = 0;
% nSTO = 1;                % рассинхронизация на 1 отсчёт
% nSTO = T_Guard/2;      % рассинхронизация на T_Guard/2 отсчётов
% nSTO = T_Guard;        % рассинхронизация на T_Guard отсчётов
% Tx_OFDM_Signal=add_STO(Tx_OFDM_Signal_matrix(:),nSTO);

%% 5. Канал
%SNR_dB = 25;
%[Rx_OFDM_Signal,~] = Noise(SNR_dB, Tx_OFDM_Signal);
figure();
plot(abs(fft(Tx_OFDM_Signal_matrix(1:Nfft+T_Guard,1),[],1)), 'LineWidth', 2);
xlabel('Номер отсчёта');
ylabel('Амплитуда')
grid on
set(gca, 'Fontsize', 20)
title('АЧХ OFDM сигнала')

Rx_OFDM_Signal = Tx_OFDM_Signal;

%% 6. Разделение сигнала на OFDM-символы
Rx_OFDM_Signal =  reshape(Rx_OFDM_Signal,[Nfft+T_Guard, N_symb]);

%% 7. OFDM-Демодуляция: удаление CP + переход в частотную область - получение точек созвездия
RX_OFDM_mapped_carriers = OFDM_demodulator(Rx_OFDM_Signal, T_Guard);

figure();
plot(abs(RX_OFDM_mapped_carriers(:,1)), 'LineWidth', 2);
xlabel('Номер отсчёта');
ylabel('Амплитуда')
xlim([1,Nfft])
grid on
set(gca, 'Fontsize', 20)
title_name = "АЧХ OFDM сигнала на этапе получения точек созвездия (рассинхронизация на " + num2str(nSTO) +" отсчётов)";
title(title_name)

RX_IQ = get_payload(RX_OFDM_mapped_carriers,dataCarriers);

RX_IQ = RX_IQ(:).';

title_name = "Созвездие на этапе получения IQ точек (рассинхронизация на " + num2str(nSTO) +" отсчётов)";
scatterplot(RX_IQ(1:400));
title(title_name)
%% 8. Demapping
output_bits = demapping(pad, RX_IQ, Constellation);


%% 9. Проверка полученных данных

if (output_bits==input_bits)
    disp("Проверка пройдена!");
    BER = BER_func(input_bits, output_bits);
    disp("BER="+num2str(BER));
    figure();
    display_pic(output_bits) % отображение полученной части картинки
else
    disp("Проверка НЕ пройдена!");
    BER = BER_func(input_bits, output_bits);
    disp("BER="+ num2str(BER));
end