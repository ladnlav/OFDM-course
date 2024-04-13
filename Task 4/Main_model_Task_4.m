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
Percent_pilot = 15;

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

% Генерируем массив пилотов
max_amplitude = max(abs(dict));
amp_pilots = 4/3*max_amplitude ;        % пилотное значение
pilotValues = zeros(1,length(pilotCarriers));
pilotValues(1:2:end) = amp_pilots*exp(1i*0);
pilotValues(2:2:end) = amp_pilots*exp(1i*pi);
pilotValues = repmat(pilotValues', 1, N_symb);
%% 1. Чтение файла
Size_Buffer = Amount_ODFM_SpF*Amount_OFDM_Frames*amount_data_carriers*bps;
input_bits = file_reader(File, Size_Buffer);

%% 2. Scrambling
% Scrambling
Register = [1 0 0 1 0 1 0 1 0 0 0 0 0 0 0];
sc_bits_matrix=zeros(Amount_ODFM_SpF*amount_data_carriers*bps,Amount_OFDM_Frames);

% Для каждого нового OFDM кадра сбрасывать состояние регистра РСЛОС до
% начального
for i = 1:Amount_OFDM_Frames
    start_index = (i-1)*Amount_ODFM_SpF*amount_data_carriers*bps+1;
    end_index = i*Amount_ODFM_SpF*amount_data_carriers*bps;
    
    if i == Amount_OFDM_Frames
        sc_bits_matrix(:,i) = Scrambler(Register, input_bits(start_index:end));
    else
        sc_bits_matrix(:,i) = Scrambler(Register, input_bits(start_index:end_index));
    end
end
sc_bits = sc_bits_matrix(:).';
%% 2. Mapping
[TX_IQ,pad] = mapping(sc_bits,Constellation);

%% 3. Формирование полосы
OFDM_mapped_carriers = OFDM_map_carriers(TX_IQ, N_symb, Nfft, dataCarriers,pilotCarriers, pilotValues);

%% 4. OFDM-модуляция: переход во временную область + добавление CP
Tx_OFDM_Signal_matrix = OFDM_modulator(OFDM_mapped_carriers, T_Guard);

Tx_OFDM_Signal=Tx_OFDM_Signal_matrix(:);
%% 5. Расчёт PAPR

% PAPR всего сигнала
PAPR_sig  = calculatePAPR(Tx_OFDM_Signal);
disp(['PAPR всего сигнала равен: ',int2str(PAPR_sig),' dB']);

% Оконный PAPR
% PAPRs = calculate_window_PAPR(Tx_OFDM_Signal,Nfft);
% [PAPR_ccdf, CCDF] = calculateCCDF(PAPRs);

%% 6. Канал
% Чистый канал
Rx_OFDM_Signal = Tx_OFDM_Signal;

% Добавление АБГШ
SNR_dB = 25;
% [Rx_OFDM_Signal,~] = Noise(SNR_dB, Tx_OFDM_Signal);

%% 6.1 Временная рассинхронизация на Time_Delay
% Time_Delay выбирается случайным образом и принимает целые 
% числа в диапазоне от 0 до Nfft+T_guard
% Time_Delay = randi([0, Nfft + T_Guard], 1);
% Time_Delay = 20;
% Rx_OFDM_Signal=add_STO(Rx_OFDM_Signal,Time_Delay); % рассинхронизация на Time_Delay отсчётов
% disp(strcat("Текущая временная задержка: ",num2str(Time_Delay)))

%% 6.2 Частотная рассинхронизация на Freq_Shift
%Freq_Shift = randi([-30, 30], 1);
% Freq_Shift=25.24;
% Rx_OFDM_Signal = add_CFO(Rx_OFDM_Signal,Freq_Shift,Nfft);
% disp(strcat("Текущий частотный сдвиг в канале: ",num2str(Freq_Shift)))
%% Влияние шума на оценку частотного сдвига
% SNRs=(0:1:60);
% MSEs=zeros(length(SNRs));
% Estimated_freq_offsets = zeros(1,length(SNRs));
% for i=1:length(SNRs)
%     SNR_dB = SNRs(i);
%     [Rx_OFDM_Signal2,~] = Noise(SNR_dB, Tx_OFDM_Signal);
%     Time_Delay = 150;
%     Rx_OFDM_Signal2=add_STO(Rx_OFDM_Signal2,Time_Delay); % рассинхронизация на Time_Delay отсчётов
%     Freq_Shift=25.24;
%     Rx_OFDM_Signal2 = add_CFO(Rx_OFDM_Signal2,Freq_Shift,Nfft);
% 
%     [~, ~, FreqOffset] = AutoCorrFunction(Rx_OFDM_Signal2, T_Guard/2, Nfft);
%     [~, e_IFO] = remove_IFO(Rx_OFDM_Signal2,Nfft);
% 
%     Estimated_freq_offsets(i) = FreqOffset+e_IFO;
% end
% figure();
% plot(SNRs, abs(Estimated_freq_offsets-Freq_Shift), 'LineWidth', 2);
% grid on;
% xlabel('SNR (dB)');
% ylabel('|Estimated \delta f - \delta f|');
% set(gca, 'Fontsize', 20)
% title('Влияние шума на оценку частотного сдвига');
%% 6.3 Многолучевое распространение
channel_taps = [
    0, 1;   % (delay, amplitude)
    4, 0.6;
    10, 0.3;
    16, 0.7
];

% Получаем импульсную и частотную характеристики для заданной конфигурации
% лучей
[H_tau, H_freq] = get_MP_channel_resp(channel_taps, Nfft);

% Прохождение через многолучевой канал
Rx_OFDM_Signal = conv(Rx_OFDM_Signal,H_tau.',"full");
Rx_OFDM_Signal = Rx_OFDM_Signal(length(H_tau):end);
%% 6.4 АЧХ на приёмнике
% figure();
% plot((abs(fft(Rx_OFDM_Signal(1:Nfft+T_Guard),[],1))), 'LineWidth', 2);
% xlabel('Номер отсчёта');
% ylabel('Амплитуда')
% grid on
% set(gca, 'Fontsize', 20)
% title('АЧХ OFDM сигнала')

%% 7. Разделение сигнала на OFDM-символы
% 7.1 Грубая временная синхронизация
% [AutoCorr, TgPosition, FreqOffset] = AutoCorrFunction(Rx_OFDM_Signal, T_Guard, Nfft);
% disp(strcat("Оценённая позиция защитного интервала: ", num2str(TgPosition)))
% disp(strcat("Оценённый дробный частотный сдвиг: ",num2str(FreqOffset)))
% 
% figure()
% plot(abs(AutoCorr),'LineWidth', 1.5);
% grid on;
% xlabel('Номер отсчёта');
% ylabel('|ACF|');
% title('ACF принятого сигнала');

% Компенсация когда мы точно знаем позицию начала первого ofdm-символа
% Rx_OFDM_Signal = add_STO(Rx_OFDM_Signal,-(Nfft+T_Guard-(TgPosition-1)));

% Компенсация, когда мы знаем только точную позицию T_Guard в первом
% ofdm-символе
% Rx_OFDM_Signal = add_STO(Rx_OFDM_Signal,TgPosition);
% Первый ofdm-символ, который мы "прослушали" заменим нулями
% Rx_OFDM_Signal = add_STO(Rx_OFDM_Signal,-(Nfft+T_Guard));


% 7.2 Частотная синхронизация
% Устранение дробного сдвига частотного смещения
% Rx_OFDM_Signal = add_CFO(Rx_OFDM_Signal,-FreqOffset,Nfft);
% Устранение целочисленного сдвига частотного смещения
% [Rx_OFDM_Signal, e_IFO] = remove_IFO(Rx_OFDM_Signal,Nfft);
% disp(strcat("Оценённый целочисленный частотный сдвиг: ",num2str(e_IFO)));

% plot(abs(fft(Rx_OFDM_Signal(1:Nfft+T_Guard),[],1)));
Rx_OFDM_Signal =  reshape(Rx_OFDM_Signal,[Nfft+T_Guard, N_symb]);
%% 8. OFDM-Демодуляция: удаление CP + переход в частотную область - получение точек созвездия
RX_OFDM_mapped_carriers = OFDM_demodulator(Rx_OFDM_Signal, T_Guard);

%% 8.5 Тонкая временная синхронизация
% RX_OFDM_mapped_carriers=fine_time_sync(RX_OFDM_mapped_carriers,pilotCarriers,pilotValues);

%% 8.7 Эквализация
[H_est, Hest_at_pilots] = estimate_channel(RX_OFDM_mapped_carriers,allCarriers,pilotCarriers,pilotValues);

figure;
plot(allCarriers, abs(H_freq), 'LineWidth', 2, 'DisplayName', 'Correct Channel');
hold on;
stem(pilotCarriers, abs(Hest_at_pilots), 'filled', 'DisplayName', 'Pilot estimates');
plot(allCarriers, abs(H_est), '--', 'LineWidth', 1.5, 'DisplayName', 'Estimated channel via interpolation');
grid on;
xlabel('Carrier index');
ylabel('|H(f)|');
legend('FontSize', 10);
ylim([0, 2]);
title('Channel Estimation');
hold off;

RX_OFDM_mapped_carriers = equalize_signal(RX_OFDM_mapped_carriers, H_est);
%% 8.9 Получение точек созвездия
RX_IQ = get_payload(RX_OFDM_mapped_carriers,dataCarriers);
RX_IQ = RX_IQ(:).';

scatterplot(RX_IQ(length(dataCarriers)+1:length(dataCarriers)+332));
% здесь видно, что "нулевые" поднесущие берут на себя часть мощности АБГШ
% scatterplot(RX_OFDM_mapped_carriers(1:400)); 
%% 9. Demapping
output_bits = demapping(pad, RX_IQ, Constellation);

%% 10. DeScrambler
dsc_bits_matrix=zeros(Amount_ODFM_SpF*amount_data_carriers*bps,Amount_OFDM_Frames);

% Для каждого нового OFDM кадра сбрасывать состояние регистра РСЛОС до
% начального
for i = 1:Amount_OFDM_Frames
    start_index = (i-1)*Amount_ODFM_SpF*amount_data_carriers*bps+1;
    end_index = i*Amount_ODFM_SpF*amount_data_carriers*bps;
    
    if i == Amount_OFDM_Frames
        dsc_bits_matrix(:,i) = DeScrambler(Register, output_bits(start_index:end));
    else
        dsc_bits_matrix(:,i) = DeScrambler(Register, output_bits(start_index:end_index));
    end
end
dsc_bits = dsc_bits_matrix(:).';
%% 11. Проверка полученных данных
BER = BER_func(input_bits, dsc_bits);
if BER<0.2
    figure();
    disp("Проверка пройдена!");
    display_pic(dsc_bits) % отображение полученной части картинки
else
    disp("Проверка НЕ пройдена!");
end
MER = MER_func(RX_IQ,Constellation);
% Разница в MER и SNR вызвана тем, что шум накладывается на сигнал с CP, а
% MER рассчитывается по сигналу уже без CP
disp("При значениии SNR=" +num2str(SNR_dB)+" dB; "+" MER="+num2str(MER)+" dB; "+ "BER="+ num2str(BER));