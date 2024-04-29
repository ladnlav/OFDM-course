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

%% 6. Канал. Параметры рассинхронизации
% АБГШ
noise_desync = 0;                % 1 - on, 0- off
% Временная
time_desync = 0;               
% Частотная
freq_desync = 0;
% Многолучевое распространение
mp_desync = 0;

% Чистый канал
Rx_OFDM_Signal = Tx_OFDM_Signal;

%% 6.0 Добавление АБГШ
SNR_dB = 25;
if noise_desync
    [Rx_OFDM_Signal,~] = Noise(SNR_dB, Tx_OFDM_Signal);
end
%% 6.1 Временная рассинхронизация на Time_Delay
% Time_Delay выбирается случайным образом и принимает целые 
% числа в диапазоне от 0 до Nfft+T_guard
if time_desync
    Time_Delay = randi([0, Nfft + T_Guard], 1);
%     Time_Delay = 20;
    Rx_OFDM_Signal=add_STO(Rx_OFDM_Signal,Time_Delay); % рассинхронизация на Time_Delay отсчётов
    disp(strcat("Текущая временная задержка: ",num2str(Time_Delay)))
end
%% 6.2 Частотная рассинхронизация на Freq_Shift
if freq_desync
    Freq_Shift = randi([0, 30], 1)+ (rand(1)-0.5);
%     Freq_Shift=25.24;
    Rx_OFDM_Signal = add_CFO(Rx_OFDM_Signal,Freq_Shift,Nfft);
    disp(strcat("Текущий частотный сдвиг в канале: ",num2str(Freq_Shift)))
end
%% Влияние шума на оценку частотного сдвига
% SNRs=(0:1:60);
% Estimated_freq_offsets = zeros(1,length(SNRs));
% for i=1:length(SNRs)
%     SNR_dB = SNRs(i);
%     [Rx_OFDM_Signal2,~] = Noise(SNR_dB, Tx_OFDM_Signal);
%     Time_Delay = 150;
%     Rx_OFDM_Signal2=add_STO(Rx_OFDM_Signal2,Time_Delay); % рассинхронизация на Time_Delay отсчётов
%     Freq_Shift=0.24;
%     Rx_OFDM_Signal2 = add_CFO(Rx_OFDM_Signal2,Freq_Shift,Nfft);
% 
%     [~, ~, FreqOffset] = AutoCorrFunction(Rx_OFDM_Signal2, T_Guard, Nfft);
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

%% Влияние шума на тонкую временную синхронизацию
% SNRs=(0:1:40);
% MERs = zeros(1,length(SNRs));
% gifFile = 'myAnimation.gif';
% 
% % Create a cell array to store scatter plot handles
% scatterplots = cell(1, length(SNRs));
% 
% for i = 1:length(SNRs)
%     SNR_dB = SNRs(i);
%     [Rx_OFDM_Signal3, ~] = Noise(SNR_dB, Tx_OFDM_Signal);
%     Time_Delay = 12;
%     Rx_OFDM_Signal3 = add_STO(Rx_OFDM_Signal3, Time_Delay); % рассинхронизация на Time_Delay отсчётов
% 
%     [AutoCorr, TgPosition, FreqOffset] = AutoCorrFunction(Rx_OFDM_Signal3, T_Guard, Nfft);
%     Rx_OFDM_Signal3 = add_STO(Rx_OFDM_Signal3, TgPosition);
%     % Первый ofdm-символ, который мы "прослушали" заменим нулями
%     Rx_OFDM_Signal3 = add_STO(Rx_OFDM_Signal3, -(Nfft + T_Guard));
% 
%     Rx_OFDM_Signal3 = reshape(Rx_OFDM_Signal3, [Nfft + T_Guard, N_symb]);
%     RX_OFDM_mapped_carriers3 = OFDM_demodulator(Rx_OFDM_Signal3, T_Guard);
%     RX_OFDM_mapped_carriers3 = fine_sync(RX_OFDM_mapped_carriers3, pilotCarriers, pilotValues, 1, freq_desync);
% 
%     RX_IQ = get_payload(RX_OFDM_mapped_carriers3, dataCarriers);
%     RX_IQ = RX_IQ(:).';
%     
%     MERs(i) = MER_func(RX_IQ(Nfft+T_Guard+1:end),Constellation);
%     
%     scatterplots{i} = scatterplot(RX_IQ(length(dataCarriers) + 1:2 * length(dataCarriers)));
%     hold on;
%     title1 =  sprintf('SNR = %d dB', SNR_dB); % Generate title
%     title(title1); % Add title to scatter plot
%     hold off;
% end
% 
% % Save scatter plots as frames of the GIF
% for i = 1:length(SNRs)
%     frame = scatterplots{i};
%     frame_filename = sprintf('frame_%d.png', i);
%     saveas(frame, frame_filename);
%     % Close the figure to avoid accumulation in memory
%     close(frame);
% end
% 
% % Create GIF from saved frames
% frames = cell(1, length(SNRs));
% for i = 1:length(SNRs)
%     frames{i} = rgb2gray(imread(sprintf('frame_%d.png', i))); % Convert to grayscale
% end
% delete frame_*.png; % Clean up the frame images
% 
% % Write GIF file
% gif_speed = 1; % Adjust speed as needed
% for i = 1:length(SNRs)
%     if i == 1
%         imwrite(frames{i}, gifFile, 'Loop', Inf, 'DelayTime', gif_speed);
%     else
%         imwrite(frames{i}, gifFile, 'WriteMode', 'append', 'DelayTime', gif_speed);
%     end
% end
% figure();
% plot(SNRs, abs(MERs-SNRs), 'LineWidth', 2);
% grid on;
% xlabel('SNR (dB)');
% ylabel('(MER-SNR)');
% set(gca, 'Fontsize', 20)
% title('Влияние шума на тонкую временную синхронизацию');

%% Влияние шума на оценку канальной характеристики
% SNRs=(0:0.5:30);
% MSEs=zeros(1,length(SNRs));
% for i=1:length(SNRs)
%     SNR_dB = SNRs(i);
%     % АБГШ
%     [Rx_OFDM_Signal,~] = Noise(SNR_dB, Tx_OFDM_Signal);
%     
%     channel_taps = [
%         0, 1;   % (delay, amplitude)
%         4, 0.6;
%         10, 0.3
%     ];
%     
%     % Получаем импульсную и частотную характеристики для заданной конфигурации
%     % лучей
%     [H_tau, H_freq] = get_MP_channel_resp(channel_taps, Nfft);
% 
%     % Прохождение через многолучевой канал
%     Rx_OFDM_Signal = conv(Rx_OFDM_Signal,H_tau.',"full");
%     Rx_OFDM_Signal = Rx_OFDM_Signal(1:end-length(H_tau)+1);
% 
%     Rx_OFDM_Signal =  reshape(Rx_OFDM_Signal,[Nfft+T_Guard, N_symb]);
%     RX_OFDM_mapped_carriers = OFDM_demodulator(Rx_OFDM_Signal, T_Guard);
%    
%     [H_est, ~] = estimate_channel(RX_OFDM_mapped_carriers,allCarriers,pilotCarriers,pilotValues);
% 
%     
%     H_freq_useful = H_freq(1:N_carrier);
%     H_est_useful = H_est(1:N_carrier);
%     
%     MSEs(i) = (H_freq_useful-H_est_useful)*(H_freq_useful-H_est_useful)'/N_carrier;
% 
% end
% figure();
% plot(SNRs, MSEs, 'LineWidth', 2);
% grid on;
% xlabel('SNR (dB)');
% ylabel('NMSE');
% set(gca, 'Fontsize', 20)
% title(strcat('Влияние шума на оценку канала: NMSE(SNR)',  '\newline',...
%     ' Nfft=', num2str(Nfft), ...
%     ' N\_carrier=',num2str(N_carrier), ...
%     ' pilot\_step=',num2str(pilot_step)));

%% 6.3 Многолучевое распространение
if mp_desync
    channel_taps = [
        0, 1;   % (delay, amplitude)
        4, 0.6;
        10, 0.3
    ];
    
    % Получаем импульсную и частотную характеристики для заданной конфигурации
    % лучей
    [H_tau, H_freq] = get_MP_channel_resp(channel_taps, Nfft);
    
    % Прохождение через многолучевой канал
    Rx_OFDM_Signal = conv(Rx_OFDM_Signal,H_tau.',"full");
    Rx_OFDM_Signal = Rx_OFDM_Signal(1:end-length(H_tau)+1);
end
%% 6.4 АЧХ на приёмнике
figure();
plot((abs(fft(Rx_OFDM_Signal(Nfft+2*T_Guard+1:2*(Nfft+T_Guard)),[],1))), 'LineWidth', 2);
xlabel('Номер отсчёта');
ylabel('Амплитуда')
grid on
set(gca, 'Fontsize', 20)
title('АЧХ полезной части OFDM-символа')

%% 7. Разделение сигнала на OFDM-символы
% 7.1 Грубая временная синхронизация
if time_desync || freq_desync
    [AutoCorr, TgPosition, FreqOffset] = AutoCorrFunction(Rx_OFDM_Signal, T_Guard, Nfft);
    disp(strcat("Оценённая позиция защитного интервала: ", num2str(TgPosition)))
    disp(strcat("Оценённый дробный частотный сдвиг: ",num2str(FreqOffset)))
    
    figure()
    plot(abs(AutoCorr),'LineWidth', 1.5);
    grid on;
    xlabel('Номер отсчёта');
    ylabel('|ACF|');
    title('ACF принятого сигнала');
    
    if time_desync
        % Компенсация, когда мы знаем только точную позицию T_Guard в первом
        % ofdm-символе
        Rx_OFDM_Signal = add_STO(Rx_OFDM_Signal,TgPosition);
        % Первый ofdm-символ, который мы "прослушали" заменим нулями
        Rx_OFDM_Signal = add_STO(Rx_OFDM_Signal,-(Nfft+T_Guard));
    end
end

% 7.2 Частотная синхронизация
if freq_desync
%     Устранение дробного сдвига частотного смещения
    Rx_OFDM_Signal = add_CFO(Rx_OFDM_Signal,-FreqOffset,Nfft);
%     Устранение целочисленного сдвига частотного смещения
    [Rx_OFDM_Signal, e_IFO] = remove_IFO(Rx_OFDM_Signal,Nfft);
    disp(strcat("Оценённый целочисленный частотный сдвиг: ",num2str(e_IFO)));
end

% plot(abs(fft(Rx_OFDM_Signal(1:Nfft+T_Guard),[],1)));
Rx_OFDM_Signal =  reshape(Rx_OFDM_Signal,[Nfft+T_Guard, N_symb]);
%% 8. OFDM-Демодуляция: удаление CP + переход в частотную область - получение точек созвездия
RX_OFDM_mapped_carriers = OFDM_demodulator(Rx_OFDM_Signal, T_Guard);

%% 8.5 Тонкая временная синхронизация
if time_desync || freq_desync
    RX_OFDM_mapped_carriers=fine_sync(RX_OFDM_mapped_carriers,pilotCarriers,pilotValues,time_desync, freq_desync);
end
%% 8.7 Эквализация
if mp_desync
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
    ylim([0, max(abs(H_freq))]);
    xlim([0,N_carrier]);
    title('Channel Estimation');
    hold off;
    
    RX_OFDM_mapped_carriers = equalize_signal(RX_OFDM_mapped_carriers, H_est, N_carrier);
end
%% 8.9 Получение точек созвездия
scatterplot(RX_OFDM_mapped_carriers(:,2));
title("Сигнальное созвездие с пилотами")

RX_IQ = get_payload(RX_OFDM_mapped_carriers,dataCarriers);
RX_IQ = RX_IQ(:).';

scatterplot(RX_IQ(length(dataCarriers)+1:2*length(dataCarriers)));
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