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
Constellation = "16QAM"; % 'BPSK' / 'QPSK' / '8PSK' / '16QAM'
[dict,bps] = constellation_func(Constellation); % bps - bits per symbol - битов на один символ созвездия

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
max_amplitude = max(abs(dict));
amp_pilots = 4/3*max_amplitude ;        % пилотное значение

OFDM_mapped_carriers = OFDM_map_carriers(TX_IQ, N_symb, Nfft, dataCarriers, pilotCarriers,amp_pilots);

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

% % Построить график CCDF кривой
% figure();
% semilogy(PAPR_ccdf, CCDF, 'LineWidth', 2);
% grid on;
% xlabel('PAPR (dB)');
% ylabel('CCDF');
% set(gca, 'Fontsize', 20)
% legend('Обычный сигнал','Location','southeast');
% title('CCDF(PAPR) в режиме скользящего окна');

%% 6. Канал. Параметры рассинхронизации
% АБГШ
noise_desync = 0;                % 1 - on, 0- off
% Временная
time_desync = 0;               
% Частотная
freq_desync = 0;
% Многолучевое распространение
mp_desync = 1;

% Чистый канал
Rx_OFDM_Signal = Tx_OFDM_Signal;

%% 6.0 Добавление АБГШ
SNR_dB = 25;
if noise_desync
    [Rx_OFDM_Signal,~] = Noise(SNR_dB, Tx_OFDM_Signal);
end

%% 6.1 Временная рассинхронизация на Time_Delay
if time_desync
%     Time_Delay = randi([0, Nfft + T_Guard], 1);
    Time_Delay = 37;
    Rx_OFDM_Signal=add_STO(Rx_OFDM_Signal,Time_Delay); % рассинхронизация на Time_Delay отсчётов
    disp(strcat("Текущая временная задержка: ",num2str(Time_Delay)))
end

%% 6.2 Частотная рассинхронизация на Freq_Shift
if freq_desync
%     Freq_Shift = randi([-30, 30], 1)+ (rand(1)-0.5);
    Freq_Shift=100;
    Rx_OFDM_Signal = add_CFO(Rx_OFDM_Signal,Freq_Shift,Nfft);
    disp(strcat("Текущий частотный сдвиг в канале: ",num2str(Freq_Shift)))
end

%% 6.3 Многолучевое распространение
if mp_desync
    channel_taps = [
        0, 1;   % (delay, amplitude)
        2, 0.4;
        4, 0.01
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
plot((abs(fft(Rx_OFDM_Signal(T_Guard+1:Nfft+T_Guard),[],1))), 'LineWidth', 2);
xlabel('Номер отсчёта');
ylabel('Амплитуда')
grid on
set(gca, 'Fontsize', 20)
title('АЧХ полезной части OFDM-символа')

%% 7. Разделение сигнала на OFDM-символы
Rx_OFDM_Signal =  reshape(Rx_OFDM_Signal,[Nfft+T_Guard, N_symb]);

%% 8. OFDM-Демодуляция: удаление CP + переход в частотную область - получение точек созвездия
RX_OFDM_mapped_carriers = OFDM_demodulator(Rx_OFDM_Signal, T_Guard);

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
if input_bits==dsc_bits
    figure();
    disp("Проверка пройдена!");
    display_pic(dsc_bits) % отображение полученной части картинки
else
    disp("Проверка НЕ пройдена!");
end

BER = BER_func(input_bits, dsc_bits);
MER = MER_func(RX_IQ,Constellation);
% Разница в MER и SNR вызвана тем, что шум накладывается на сигнал с CP, а
% MER рассчитывается по сигналу уже без CP
disp("При значениии SNR=" +num2str(SNR_dB)+" dB; "+" MER="+num2str(MER)+" dB; "+ "BER="+ num2str(BER));

%% BER (SNR)
SNRs=(0:0.5:30);
BERs=zeros(4,length(SNRs));
Constellations = ["BPSK"  "QPSK"  "8PSK"  "16QAM"];
for const=1:length(Constellations)
    Constellation = Constellations(const);
    [dict,bps] = constellation_func(Constellation); % bps - bits per symbol - битов на один символ созвездия

    Size_Buffer = Amount_ODFM_SpF*Amount_OFDM_Frames*amount_data_carriers*bps;
    input_bits = file_reader(File, Size_Buffer);

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
    
    % mapping
    [TX_IQ,pad] = mapping(sc_bits,Constellation);

    % Формирование полосы
    max_amplitude = max(abs(dict));
    amp_pilots = 4/3*max_amplitude ;        % пилотное значение

    OFDM_mapped_carriers = OFDM_map_carriers(TX_IQ, N_symb, Nfft, dataCarriers, pilotCarriers,amp_pilots);
    
    % OFDM-модуляция: переход во временную область + добавление CP
    Tx_OFDM_Signal_matrix = OFDM_modulator(OFDM_mapped_carriers, T_Guard);
    
    Tx_OFDM_Signal=Tx_OFDM_Signal_matrix(:);
    
    PAPR_sig  = calculatePAPR(Tx_OFDM_Signal);
    disp(['PAPR всего сигнала равен: ',int2str(PAPR_sig),' dB']);

    for i=1:length(SNRs)
        SNR_dB = SNRs(i);
        [Rx_OFDM_Signal,~] = Noise(SNR_dB, Tx_OFDM_Signal);
        Rx_OFDM_Signal =  reshape(Rx_OFDM_Signal,[Nfft+T_Guard, N_symb]);
    
        RX_OFDM_symbols = OFDM_demodulator(Rx_OFDM_Signal, T_Guard);
    
        RX_IQ = get_payload(RX_OFDM_symbols,dataCarriers);
        RX_IQ = RX_IQ(:).';
        
        output_bits = demapping(pad, RX_IQ, Constellation);
    
        dsc_bits_matrix=zeros(Amount_ODFM_SpF*amount_data_carriers*bps,Amount_OFDM_Frames);
    
        % Для каждого нового OFDM кадра сбрасывать состояние регистра РСЛОС до
        % начального
        for j = 1:Amount_OFDM_Frames
            start_index = (j-1)*Amount_ODFM_SpF*amount_data_carriers*bps+1;
            end_index = j*Amount_ODFM_SpF*amount_data_carriers*bps;
            
            if j == Amount_OFDM_Frames
                dsc_bits_matrix(:,j) = DeScrambler(Register, output_bits(start_index:end));
            else
                dsc_bits_matrix(:,j) = DeScrambler(Register, output_bits(start_index:end_index));
            end
        end
        dsc_bits = dsc_bits_matrix(:).';
    
    
        BERs(const,i) = BER_func(input_bits, dsc_bits);
    end
end
figure();
for i=1:length(Constellations)
    semilogy(SNRs, BERs(i,:), 'LineWidth', 2);
    hold on;
end
grid on;
xlabel('SNR (dB)');
ylabel('BER');
set(gca, 'Fontsize', 20)
legend("BPSK", "QPSK", "8PSK", "16QAM", 'Location','southwest');
title('BER(SNR) for different constellations');
