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

%% 2. Mapping
[TX_IQ,pad] = mapping(input_bits,Constellation);

%% 3. Формирование полосы
allCarriers = linspace(1,Nfft,Nfft); % индексы поднесущих
P = 32;                   % кол-во пилотов на всех Nfft отсчётах
pilotValue = 3+3j ;        % пилотное значение

pilotCarriers = allCarriers(1:round(Nfft/(P)):end);
pilotCarriers = [pilotCarriers,allCarriers(end)]; % поднесущие пилотных сигналов
P=P+1;

dataCarriers = allCarriers(~ismember(allCarriers, pilotCarriers));
OFDM_symbols = OFDM_symbol(TX_IQ, N_carrier, N_symb, Nfft, dataCarriers, pilotCarriers);

%% 4. OFDM-модуляция: переход во временную область + добавление CP
Tx_OFDM_Signal_matrix = OFDM_modulator(OFDM_symbols, T_Guard); % lab 2

% Tx_OFDM_Signal=Tx_OFDM_Signal_matrix(:);

%% 5.0 Временная рассинхронизация
% nSTO = 0;
% nSTO = 1;                % рассинхронизация на 1 отсчёт
% nSTO = T_Guard/2;      % рассинхронизация на T_Guard/2 отсчётов
nSTO = -T_Guard+1;        % рассинхронизация на T_Guard отсчётов
Tx_OFDM_Signal=add_STO(Tx_OFDM_Signal_matrix(:),nSTO);

%% 5. Канал
%SNR_dB = 25;
%[Rx_OFDM_Signal,~] = Noise(SNR_dB, Tx_OFDM_Signal);
% figure();
% plot(mean(abs(fft(Tx_OFDM_Signal_matrix,Nfft,1)),2), 'LineWidth', 2);
% xlabel('Номер отсчёта');
% ylabel('Амплитуда')
% grid on
% set(gca, 'Fontsize', 20)
% title('АЧХ OFDM сигнала')

Rx_OFDM_Signal = Tx_OFDM_Signal;

%% 6. Разделение сигнала на OFDM-символы
Rx_OFDM_Signal =  reshape(Rx_OFDM_Signal,[Nfft+T_Guard, N_symb]);

%% 7. OFDM-Демодуляция: удаление CP + переход в частотную область - получение точек созвездия
RX_OFDM_symbols = OFDM_demodulator(Rx_OFDM_Signal, T_Guard);

figure();
% plot(mean(abs(RX_OFDM_symbols),2), 'LineWidth', 2);
plot(abs(RX_OFDM_symbols(:,1)), 'LineWidth', 2);
xlabel('Номер отсчёта');
ylabel('Амплитуда')
xlim([1,Nfft])
grid on
set(gca, 'Fontsize', 20)
title_name = "АЧХ OFDM сигнала на этапе получения точек созвездия (рассинхронизация на " + num2str(nSTO) +" отсчётов)";
title(title_name)

RX_IQ = get_payload(RX_OFDM_symbols,dataCarriers,N_carrier);

RX_IQ = RX_IQ(:).';

title_name = "Созвездие на этапе получения IQ точек (рассинхронизация на " + num2str(nSTO) +" отсчётов)";
scatterplot(RX_IQ(1:400));
title(title_name)
%% 8. Demapping
output_bits = demapping(pad, RX_IQ, Constellation);


%% 9. Проверка полученных данных

if (output_bits==input_bits)
    disp("Проверка пройдена!");
    figure();
    display_pic(output_bits) % отображение полученной части картинки
else
    disp("Проверка НЕ пройдена!");
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

%     figure();
%     data_plot=abs(fft(OFDM_time,[],1));
%     plot(data_plot(:,2), 'LineWidth', 2);
%     xlabel('Номер отсчёта');
%     ylabel('Амплитуда')
%     grid on
%     set(gca, 'Fontsize', 20)
%     title('АЧХ одного OFDM символа до добавления T_Guard')

    % Adding cp
    cp = OFDM_time(end-T_guard+1:end,:);  % take the last CP samples
    OFDM_time_guarded = [cp; OFDM_time];
    
%     figure();
%     data_plot=abs(fft(OFDM_time_guarded,[],1));
%     plot(data_plot(:,2), 'LineWidth', 2);
%     xlabel('Номер отсчёта');
%     ylabel('Амплитуда')
%     grid on
%     set(gca, 'Fontsize', 20)
%     title('АЧХ одного OFDM символа ПОСЛЕ добавления T_Guard')
end
% Переход OFDM-symbol в частотную область
function TX_IQ=OFDM_demodulator(OFDM_time_guarded, T_guard) 
    
    % Remove cp 
    OFDM_time = OFDM_time_guarded(T_guard+1:end,:); % take the last samples without CP samples

    % FFT
    TX_IQ = fft(OFDM_time,[],1);

end
% Извлечь данные из OFDM-символов
function RX_IQ=get_payload(RX_OFDM_symbols, dataCarriers, N_carrier) 
    RX_IQ=RX_OFDM_symbols(dataCarriers(1:N_carrier),:);
end

function y_STO=add_STO(y, nSTO)
    % add STO (symbol time offset)
    % y : Received signal
    % nSTO : Number of samples corresponding to STO
    if nSTO>=0
        y_STO=[y(nSTO+1:end); zeros(1,nSTO)']; % advance
    else 
        y_STO=[zeros(1,-nSTO)'; y(1:end+nSTO)]; % delay
    end
end