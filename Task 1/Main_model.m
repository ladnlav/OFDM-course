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
function mapped_carriers=OFDM_map_carriers(QAM_payload, N_symb, Nfft, dataCarriers, pilotCarriers,amp_pilots) 
    mapped_carriers = zeros(Nfft, N_symb);
    data =  reshape(QAM_payload.',[length(dataCarriers), N_symb]);
    
    mapped_carriers(dataCarriers,:)=data;
    pilotValues =zeros(1,length(pilotCarriers));
    pilotValues(1:2:end) = amp_pilots*exp(1i*0);
    pilotValues(2:2:end) = amp_pilots*exp(1i*pi);
    
    pilotValues = repmat(pilotValues', 1, 50);
    mapped_carriers(pilotCarriers,:)=pilotValues;
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
function RX_IQ=get_payload(RX_OFDM_symbols, dataCarriers) 
    RX_IQ=RX_OFDM_symbols(dataCarriers,:);
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