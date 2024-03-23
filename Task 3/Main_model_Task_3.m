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

File = 'eagle.tiff';
Constellation = "16QAM";
[~,bps] = constellation_func(Constellation); % bps - bits per symbol - битов на один символ созвездия

%% 1. Чтение файла
Size_Buffer = Amount_ODFM_SpF*Amount_OFDM_Frames*N_carrier*bps;
input_bits = file_reader(File, Size_Buffer);

%% 2. Scrambling
% Scrambling
Register = [1 0 0 1 0 1 0 1 0 0 0 0 0 0 0];
sc_bits_matrix=zeros(Amount_ODFM_SpF*N_carrier*bps,Amount_OFDM_Frames);

% Для каждого нового OFDM кадра сбрасывать состояние регистра РСЛОС до
% начального
for i = 1:Amount_OFDM_Frames
    start_index = (i-1)*Amount_ODFM_SpF*N_carrier*bps+1;
    end_index = i*Amount_ODFM_SpF*N_carrier*bps;
    
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
allCarriers = linspace(1,Nfft,Nfft); % индексы поднесущих
P = 32;                   % кол-во пилотов
pilotValue = 3+3j ;        % пилотное значение

pilotCarriers = allCarriers(1:round(Nfft/(P)):end);
pilotCarriers = [pilotCarriers,allCarriers(end)]; % поднесущие пилотных сигналов
P=P+1;

dataCarriers = allCarriers(~ismember(allCarriers, pilotCarriers));
OFDM_symbols = OFDM_symbol(TX_IQ, N_carrier, N_symb, Nfft, dataCarriers, pilotCarriers);

%% 4. OFDM-модуляция: переход во временную область + добавление CP
Tx_OFDM_Signal_matrix = OFDM_modulator(OFDM_symbols, T_Guard);

Tx_OFDM_Signal=Tx_OFDM_Signal_matrix(:);
%% 5. Расчёт PAPR

% PAPR всего сигнала
PAPR_sig  = calculatePAPR(Tx_OFDM_Signal);
disp(['PAPR всего сигнала равен: ',int2str(PAPR_sig),' dB']);

% Оконный PAPR
PAPRs = calculate_window_PAPR(Tx_OFDM_Signal,Nfft);
[PAPR_ccdf, CCDF] = calculateCCDF(PAPRs);

% % Построить график CCDF кривой
% figure();
% semilogy(PAPR_ccdf, CCDF, 'LineWidth', 2);
% grid on;
% xlabel('PAPR (dB)');
% ylabel('CCDF');
% set(gca, 'Fontsize', 20)
% legend('Обычный сигнал','Location','southeast');
% title('CCDF(PAPR) в режиме скользящего окна');

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
% Time_Delay = 3;
% Rx_OFDM_Signal=[Rx_OFDM_Signal(Time_Delay+1:end).',zeros(1,Time_Delay)]; % рассинхронизация на Time_Delay отсчётов

%% 6.2 Частотная рассинхронизация на Freq_Shift
%Freq_Shift = randi([-30, 30], 1);
Freq_Shift=0.01;
nn=0:length(Rx_OFDM_Signal)-1; 
Rx_OFDM_Signal = Rx_OFDM_Signal.*exp(j*2*pi*Freq_Shift*nn'/Nfft);

%% 6.3 Многолучевое распространение
channel_taps = [
    0, 1;   % (delay, amplitude)
    4, 0.6;
    10, 0.3
];

% Получаем импульсную и частотную характеристики для заданной конфигурации
% лучей
[H_tau, H_freq] = get_MP_channel_resp(channel_taps, Nfft);

% Прохождение через многолучевой канал
% Rx_OFDM_Signal = conv(Rx_OFDM_Signal,H_tau.',"full");
% Rx_OFDM_Signal = Rx_OFDM_Signal(length(H_tau):end);
%% 6.4 АЧХ на приёмнике
figure();
plot((abs(fft(Rx_OFDM_Signal(1:Nfft+T_Guard),[],1))), 'LineWidth', 2);
xlabel('Номер отсчёта');
ylabel('Амплитуда')
grid on
set(gca, 'Fontsize', 20)
title('АЧХ OFDM сигнала')

%% 7. Разделение сигнала на OFDM-символы
Rx_OFDM_Signal =  reshape(Rx_OFDM_Signal,[Nfft+T_Guard, N_symb]);

%% 8. OFDM-Демодуляция: удаление CP + переход в частотную область - получение точек созвездия
RX_OFDM_symbols = OFDM_demodulator(Rx_OFDM_Signal, T_Guard);

RX_IQ = get_payload(RX_OFDM_symbols,dataCarriers,N_carrier);
RX_IQ = RX_IQ(:).';

scatterplot(RX_IQ(1:400));
% здесь видно, что "нулевые" поднесущие берут на себя часть мощности АБГШ
% scatterplot(RX_OFDM_symbols(1:400)); 
%% 9. Demapping
output_bits = demapping(pad, RX_IQ, Constellation);

%% 10. DeScrambler
dsc_bits_matrix=zeros(Amount_ODFM_SpF*N_carrier*bps,Amount_OFDM_Frames);

% Для каждого нового OFDM кадра сбрасывать состояние регистра РСЛОС до
% начального
for i = 1:Amount_OFDM_Frames
    start_index = (i-1)*Amount_ODFM_SpF*N_carrier*bps+1;
    end_index = i*Amount_ODFM_SpF*N_carrier*bps;
    
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
BERs=zeros(length(SNRs));
for i=1:length(SNRs)
    SNR_dB = SNRs(i);
    [Rx_OFDM_Signal,~] = Noise(SNR_dB, Tx_OFDM_Signal);
    Rx_OFDM_Signal =  reshape(Rx_OFDM_Signal,[Nfft+T_Guard, N_symb]);

    RX_OFDM_symbols = OFDM_demodulator(Rx_OFDM_Signal, T_Guard);

    RX_IQ = get_payload(RX_OFDM_symbols,dataCarriers,N_carrier);
    RX_IQ = RX_IQ(:).';
    
    output_bits = demapping(pad, RX_IQ, Constellation);

    dsc_bits_matrix=zeros(Amount_ODFM_SpF*N_carrier*bps,Amount_OFDM_Frames);

    % Для каждого нового OFDM кадра сбрасывать состояние регистра РСЛОС до
    % начального
    for j = 1:Amount_OFDM_Frames
        start_index = (j-1)*Amount_ODFM_SpF*N_carrier*bps+1;
        end_index = j*Amount_ODFM_SpF*N_carrier*bps;
        
        if j == Amount_OFDM_Frames
            dsc_bits_matrix(:,j) = DeScrambler(Register, output_bits(start_index:end));
        else
            dsc_bits_matrix(:,j) = DeScrambler(Register, output_bits(start_index:end_index));
        end
    end
    dsc_bits = dsc_bits_matrix(:).';


    BERs(i) = BER_func(input_bits, dsc_bits);
end
figure();
semilogy(SNRs, BERs, 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('BER');
set(gca, 'Fontsize', 20)
%legend('Обычный сигнал','Location','southeast');
title('BER(SNR)');
%% Функции

% Прочитать файл
function input_bits = file_reader(File, Size_Buffer)
    % Прочитать черно-белую картинку
    grayImage=imread(File);
    
    % Перевести картинку в битовый формат
    binaryImage = imbinarize(grayImage);
    input_bits = double(binaryImage); % из логического массива в массив 0 и 1

    % Взять нужное количество бит из всей картинки
    input_bits = input_bits(1:Size_Buffer);
    
end
% Отобразить картинку
function display_pic(binaryImage)
    
    % Добиваем до 360*360 битов нулями
    padAmount = 360*360-size(binaryImage,2);
    binaryImage=padarray(binaryImage, [0,padAmount], 0,'post');

    % Делаем из колбасы матрицу - картинку
    binaryImage=reshape(binaryImage,360,360);

    % Перевести картинку в обычный формат
    grayImage2 = uint8(binaryImage) * 255;
    
    % отобразить картинку
    imshow(grayImage2);
end
% Распределить данные по поднесущим 
function matrix=OFDM_symbol(QAM_payload, N_carrier, N_symb, Nfft, dataCarriers, pilotCarriers) 
    symbols = zeros(Nfft, N_symb);
    data =  reshape(QAM_payload.',[N_carrier, N_symb]);
    
    symbols(dataCarriers(1:N_carrier),:)=data;
    symbols(pilotCarriers,:)=0+0j;

    matrix = symbols;
end
% Создать OFDM-symbol во времени
function OFDM_time_guarded=OFDM_modulator(OFDM_symbols, T_guard) 
    
    % IFFT
    OFDM_time = ifft(OFDM_symbols,[],1);
    
    % Добавить cp
    cp = OFDM_time(end-T_guard+1:end,:);  % take the last CP samples
    OFDM_time_guarded = [cp; OFDM_time];

end
% Переход OFDM-symbol в частотную область
function TX_IQ=OFDM_demodulator(OFDM_time_guarded, T_guard) 
    
    % Удалть cp 
    OFDM_time = OFDM_time_guarded(T_guard+1:end,:); % take the last samples without CP samples

    % FFT
    TX_IQ = fft(OFDM_time,[],1);

end
% Извлечь данные из OFDM-символов
function RX_IQ=get_payload(RX_OFDM_symbols, dataCarriers, N_carrier) 
    RX_IQ=RX_OFDM_symbols(dataCarriers(1:N_carrier),:);
end

% Расчитать PAPR
function PAPR = calculatePAPR(OFDM_signal)    
    % Поиск квадрата максимальной амплитуды
    peak_power = max(abs(OFDM_signal))^2;
    
    % Нахождение средней мощности сигнала
    average_power = mean(abs(OFDM_signal).^2);
    
    % Расчёт PAPR
    PAPR = 10*log10(peak_power / average_power);
end

% Рассчитать оконный PAPR
function PAPRs = calculate_window_PAPR(Tx_OFDM_Signal,Nfft)    
    % Сколько всего будет значений PAPR для данного сигнала
    n_PAPRs=length(Tx_OFDM_Signal) - Nfft + 1;
    PAPRs = zeros(1, n_PAPRs);

    % Цикл расчёта PAPR скользящим окном
    for i = 1:n_PAPRs
            % Извлечь отсчёты в окне 
            window = Tx_OFDM_Signal(i:i+Nfft-1);
            
            % Рассчитать PAPR для текущего окна
            PAPRs(i) = calculatePAPR(window);
    end
end
% Рассчитать CCDF
function [PAPR_ccdf, CCDF] = calculateCCDF(PAPR_values)
    
    [CCDF,PAPR_ccdf] = ecdf(PAPR_values);
    CCDF = 1 - CCDF;
end
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
