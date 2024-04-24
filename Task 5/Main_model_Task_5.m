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
pilot_step = 2; % comb

allCarriers = linspace(1,Nfft,Nfft); % индексы поднесущих
amount_pilots = floor(N_carrier/pilot_step);   % кол-во пилотов на всех Nfft отсчётах
% pilot_step = floor(N_carrier/amount_pilots);
pilotCarriers = allCarriers(1:pilot_step:N_carrier);

if pilot_step == 1
    Percent_pilot = 100;
    
    allCarriers = linspace(1,Nfft,Nfft); % индексы поднесущих
    amount_pilots = round(Percent_pilot/100*N_carrier);   % кол-во пилотов на всех Nfft отсчётах
    pilot_step = floor(N_carrier/amount_pilots);
    pilotCarriers = allCarriers(1:pilot_step:N_carrier-2);
    pilotCarriers = [pilotCarriers,allCarriers(N_carrier)];
    amount_pilots = length(pilotCarriers);
end
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
noise_desync = 1;                % 1 - on, 0- off
% Многолучевое распространение
mp_desync = 1;

% Чистый канал
Rx_OFDM_Signal = Tx_OFDM_Signal;

%% 6.0 Добавление АБГШ
SNR_dB = 20;
if noise_desync
    [Rx_OFDM_Signal,~] = Noise(SNR_dB, Tx_OFDM_Signal);
end
%% 6.3 Многолучевое распространение
if mp_desync
    channel_taps = [
        0, 1;   % (delay, amplitude)
        4, 0.8;
        10, 0.6;
        15, 0.4;
        21, 0.2;
        25, 0.1;
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
title('АЧХ OFDM сигнала')

%% 7. Разделение сигнала на OFDM-символы
% plot(abs(fft(Rx_OFDM_Signal(1:Nfft+T_Guard),[],1)));
Rx_OFDM_Signal =  reshape(Rx_OFDM_Signal,[Nfft+T_Guard, N_symb]);
%% 8. OFDM-Демодуляция: удаление CP + переход в частотную область - получение точек созвездия
RX_OFDM_mapped_carriers = OFDM_demodulator(Rx_OFDM_Signal, T_Guard);
%% Test
% H_f_test = H_freq;
% H_f_test_pilots = H_f_test(pilotCarriers);
% H_t_ls_test = LS_CE(H_f_test_pilots,pilotValues,pilotCarriers,Nfft);
% 
% figure();
% plot(abs(H_tau),'b'); 
% hold on;
% plot(abs(H_t_ls_test),'r:+'); 
% legend('True Channel','ls');
% title('Impulse Response');
% 
% figure();
% plot(abs(H_freq),'b');
% hold on;
% plot(abs(fft(H_t_ls_test)),'r:+');
% legend('True Channel','ls');
% title('Frequency Response');
%% 8.7 Эквализация
if mp_desync
    [H_est, Hest_at_pilots] = estimate_channel(RX_OFDM_mapped_carriers,allCarriers,pilotCarriers,pilotValues);
%     figure;
%     plot(allCarriers, abs(H_freq), 'LineWidth', 2, 'DisplayName', 'Correct Channel');
%     hold on;
%     stem(pilotCarriers, abs(Hest_at_pilots), 'filled', 'DisplayName', 'Pilot estimates');
%     plot(allCarriers, abs(H_est), '--', 'LineWidth', 1.5, 'DisplayName', 'Estimated channel via interpolation');
%     grid on;
%     xlabel('Carrier index');
%     ylabel('|H(f)|');
%     legend('FontSize', 10);
%     ylim([0, max(abs(H_freq))]);
%     xlim([0,N_carrier]);
%     title('Channel Estimation');
%     hold off;
    
    H_est_LS_l = LS_CE(RX_OFDM_mapped_carriers,pilotValues,pilotCarriers,Nfft);
    H_est_MMSE = MMSE_CE(RX_OFDM_mapped_carriers,pilotValues,pilotCarriers,Nfft,H_tau,SNR_dB);
    [H_est_MP, h_t_MP] =  MP_estimate(RX_OFDM_mapped_carriers, pilotCarriers, Nfft,N_carrier, length(channel_taps));
    [H_est_OMP,h_t_OMP,kk1] =  OMP_estimate(RX_OFDM_mapped_carriers, pilotCarriers, Nfft,N_carrier, length(channel_taps));
%     [x,kk2] =  OMP2(RX_OFDM_mapped_carriers, pilotCarriers, Nfft,N_carrier, length(channel_taps));
    %     [h_t_s, h_f_s] = OMPForEach(RX_OFDM_mapped_carriers, Nfft, pilotCarriers, length(channel_taps));

    H_freq_useful = H_freq(1:N_carrier);
    H_est_LS_l_useful = H_est_LS_l(1:N_carrier);
    H_est_MMSE_useful = H_est_MMSE(1:N_carrier);
    H_est_MP_useful = H_est_MP(1:N_carrier);
    H_est_OMP_useful = H_est_OMP(1:N_carrier);
    
    MSE_l  = (H_freq_useful-H_est_LS_l_useful)*(H_freq_useful-H_est_LS_l_useful)'/N_carrier;
    MSE_mmse = (H_freq_useful-H_est_MMSE_useful)*(H_freq_useful-H_est_MMSE_useful)'/N_carrier;
    MSE_mp = (H_freq_useful-H_est_MP_useful)*(H_freq_useful-H_est_MP_useful)'/N_carrier;
    MSE_omp = (H_freq_useful-H_est_OMP_useful)*(H_freq_useful-H_est_OMP_useful)'/N_carrier;

    figure();
    subplot(4,1,1);
    plot(abs(H_freq_useful),'b'); 
    hold on;
    plot(abs(H_est_LS_l_useful),'r:+'); 
    legend('True Channel','linear');
    title(strcat("MSE=",num2str(MSE_l)));
    subplot(4,1,2);
    plot(abs(H_freq_useful),'b'); 
    hold on;
    plot(abs(H_est_MMSE_useful),'r:+');
    legend('True Channel','MMSE');
    title(strcat("MSE=",num2str(MSE_mmse)));
    subplot(4,1,3);
    plot(abs(H_freq_useful),'b'); 
    hold on;
    plot(abs(H_est_MP_useful),'r:+');
    legend('True Channel','MP');
    title(strcat("MSE=",num2str(MSE_mp)));
    subplot(4,1,4);
    plot(abs(H_freq_useful),'b'); 
    hold on;
    plot(abs(H_est_OMP_useful),'r:+');
    legend('True Channel','OMP');
    title(strcat("MSE=",num2str(MSE_omp)));
    hold off;
    
%     figure()
%     plot(abs(H_freq))
%     hold on;
%     plot(abs(H_est_OMP))
    figure()
    plot(abs(H_tau))
    hold on;
    plot(abs(h_t_MP))
    plot(abs(h_t_OMP))
    legend('True Channel','MP','OMP');
    RX_OFDM_mapped_carriers = equalize_signal(RX_OFDM_mapped_carriers, H_est, N_carrier);
end
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