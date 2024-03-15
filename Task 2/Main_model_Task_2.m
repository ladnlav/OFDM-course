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
[TX_IQ,~] = mapping(input_bits,Constellation);
[sc_TX_IQ,pad] = mapping(sc_bits,Constellation);

%% 3. Формирование полосы
allCarriers = linspace(1,Nfft,Nfft); % индексы поднесущих
P = 32;                   % кол-во пилотов
pilotValue = 3+3j ;        % пилотное значение

pilotCarriers = allCarriers(1:round(Nfft/(P)):end);
pilotCarriers = [pilotCarriers,allCarriers(end)]; % поднесущие пилотных сигналов
P=P+1;

dataCarriers = allCarriers(~ismember(allCarriers, pilotCarriers));
OFDM_symbols = OFDM_symbol(TX_IQ, N_carrier, N_symb, Nfft, dataCarriers, pilotCarriers);
sc_OFDM_symbols = OFDM_symbol(sc_TX_IQ, N_carrier, N_symb, Nfft, dataCarriers, pilotCarriers);

%% 4. OFDM-модуляция: переход во временную область + добавление CP
Tx_OFDM_Signal_matrix = OFDM_modulator(OFDM_symbols, T_Guard);
sc_Tx_OFDM_Signal_matrix = OFDM_modulator(sc_OFDM_symbols, T_Guard);

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
RX_OFDM_symbols = OFDM_demodulator(Rx_OFDM_Signal, T_Guard);

RX_IQ = get_payload(RX_OFDM_symbols,dataCarriers,N_carrier);
RX_IQ = RX_IQ(:).';

sc_RX_OFDM_symbols = OFDM_demodulator(sc_Rx_OFDM_Signal, T_Guard);

sc_RX_IQ = get_payload(sc_RX_OFDM_symbols,dataCarriers,N_carrier);
sc_RX_IQ = sc_RX_IQ(:).';

%% 9. Demapping
output_bits = demapping(pad, RX_IQ, Constellation);
sc_output_bits = demapping(pad, sc_RX_IQ, Constellation);

%% 10. DeScrambler
dsc_bits_matrix=zeros(Amount_ODFM_SpF*N_carrier*bps,Amount_OFDM_Frames);

% Для каждого нового OFDM кадра сбрасывать состояние регистра РСЛОС до
% начального
for i = 1:Amount_OFDM_Frames
    start_index = (i-1)*Amount_ODFM_SpF*N_carrier*bps+1;
    end_index = i*Amount_ODFM_SpF*N_carrier*bps;
    
    if i == Amount_OFDM_Frames
        dsc_bits_matrix(:,i) = DeScrambler(Register, sc_output_bits(start_index:end));
    else
        dsc_bits_matrix(:,i) = DeScrambler(Register, sc_output_bits(start_index:end_index));
    end
end
dsc_bits = dsc_bits_matrix(:).';
%% 11. Проверка полученных данных
if output_bits==input_bits
    figure();
    disp("Проверка пройдена!");
    display_pic(output_bits) % отображение полученной части картинки
else
    disp("Проверка НЕ пройдена!");
end

if dsc_bits==input_bits
    figure();
    disp("Проверка с рандомизатором пройдена!");
    display_pic(dsc_bits) % отображение полученной части картинки
else
    disp("Проверка с рандомизатором НЕ пройдена!");
end
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