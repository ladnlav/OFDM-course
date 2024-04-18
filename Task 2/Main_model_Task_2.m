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
Percent_pilot = 1;

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
[TX_IQ,~] = mapping(input_bits,Constellation);
[sc_TX_IQ,pad] = mapping(sc_bits,Constellation);

%% 3. Формирование полосы
max_amplitude = max(abs(dict));
amp_pilots = 2*max_amplitude ;        % пилотное значение

OFDM_mapped_carriers = OFDM_map_carriers(TX_IQ, N_symb, Nfft, dataCarriers, pilotCarriers,amp_pilots);
sc_OFDM_mapped_carriers = OFDM_map_carriers(sc_TX_IQ, N_symb, Nfft, dataCarriers, pilotCarriers,amp_pilots);

%% 4. OFDM-модуляция: переход во временную область + добавление CP
Tx_OFDM_Signal_matrix = OFDM_modulator(OFDM_mapped_carriers, T_Guard);
sc_Tx_OFDM_Signal_matrix = OFDM_modulator(sc_OFDM_mapped_carriers, T_Guard);

Tx_OFDM_Signal=Tx_OFDM_Signal_matrix(:);
sc_Tx_OFDM_Signal=sc_Tx_OFDM_Signal_matrix(:);
%% 5. Расчёт PAPR

% PAPR всего сигнала
PAPR_sig  = calculatePAPR(Tx_OFDM_Signal);
PAPR_rand = calculatePAPR(sc_Tx_OFDM_Signal);
disp(['PAPR обычного сигнала равен: ',int2str(PAPR_sig),' dB']);
disp(['PAPR рандомизированного сигнала равен: ',int2str(PAPR_rand),' dB']);

% Оконный PAPR
PAPRs = calculate_window_PAPR(Tx_OFDM_Signal,Nfft);
[PAPR_ccdf, CCDF] = calculateCCDF(PAPRs);

sc_PAPRs = calculate_window_PAPR(sc_Tx_OFDM_Signal,Nfft);
[sc_PAPR_ccdf, sc_CCDF] = calculateCCDF(sc_PAPRs);

% Построить график CCDF кривой
figure();
semilogy(PAPR_ccdf, CCDF, 'LineWidth', 2);
hold on;
semilogy(sc_PAPR_ccdf, sc_CCDF, 'LineWidth', 2);
grid on;
xlabel('PAPR (dB)');
ylabel('CCDF');
set(gca, 'Fontsize', 20)
legend('Обычный сигнал','Сигнал с рандомизатором','Location','southeast');
title('CCDF(PAPR) в режиме скользящего окна');
% plot_custom_ccdf(PAPRs); % самостоятельно написанная функция построения
% ccdf

%% 6. Канал
%SNR_dB = 25;
%[Rx_OFDM_Signal,~] = Noise(SNR_dB, Tx_OFDM_Signal);

Rx_OFDM_Signal = Tx_OFDM_Signal;
sc_Rx_OFDM_Signal=sc_Tx_OFDM_Signal;

%% 7. Разделение сигнала на OFDM-символы
Rx_OFDM_Signal =  reshape(Rx_OFDM_Signal,[Nfft+T_Guard, N_symb]);
sc_Rx_OFDM_Signal =  reshape(sc_Rx_OFDM_Signal,[Nfft+T_Guard, N_symb]);

%% 8. OFDM-Демодуляция: удаление CP + переход в частотную область - получение точек созвездия
RX_OFDM_mapped_carriers = OFDM_demodulator(Rx_OFDM_Signal, T_Guard);

RX_IQ = get_payload(RX_OFDM_mapped_carriers,dataCarriers);
RX_IQ = RX_IQ(:).';

sc_OFDM_mapped_carriers = OFDM_demodulator(sc_Rx_OFDM_Signal, T_Guard);

sc_RX_IQ = get_payload(sc_OFDM_mapped_carriers,dataCarriers);
sc_RX_IQ = sc_RX_IQ(:).';

%% 9. Demapping
output_bits = demapping(pad, RX_IQ, Constellation);
sc_output_bits = demapping(pad, sc_RX_IQ, Constellation);

%% 10. DeScrambler
dsc_bits_matrix=zeros(Amount_ODFM_SpF*amount_data_carriers*bps,Amount_OFDM_Frames);

% Для каждого нового OFDM кадра сбрасывать состояние регистра РСЛОС до
% начального
for i = 1:Amount_OFDM_Frames
    start_index = (i-1)*Amount_ODFM_SpF*amount_data_carriers*bps+1;
    end_index = i*Amount_ODFM_SpF*amount_data_carriers*bps;
    
    if i == Amount_OFDM_Frames
        dsc_bits_matrix(:,i) = DeScrambler(Register, sc_output_bits(start_index:end));
    else
        dsc_bits_matrix(:,i) = DeScrambler(Register, sc_output_bits(start_index:end_index));
    end
end
dsc_bits = dsc_bits_matrix(:).';
%% 11. Проверка полученных данных
if output_bits==input_bits
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

if dsc_bits==input_bits
    figure();
    disp("Проверка с рандомизатором пройдена!");
    display_pic(dsc_bits) % отображение полученной части картинки
    BER = BER_func(input_bits, output_bits);
    disp("BER="+num2str(BER));
else
    disp("Проверка с рандомизатором НЕ пройдена!");
    BER = BER_func(input_bits, output_bits);
    disp("BER="+ num2str(BER));
end