clc; % чистка командного окна
close all; % закрыть дополнительные окна 
clear all; % очистить память
%% MP/OMP на случайных пилотах (часть-2)
Nfft = 4096;
N_carrier = 1024;
T_Guard = Nfft /8;
Amount_OFDM_Frames = 2;
Amount_ODFM_SpF = 7;
N_symb=Amount_OFDM_Frames*Amount_ODFM_SpF;
verbose = 1; % 1 - display plots, 0 - silence mode
reg_pilot = 1; % 1 - regular pilot mask, 0 - random pilot mask
combs = 4:1:256;
amounts_pilots =  floor(N_carrier./combs);
[~,ia,~]=unique(amounts_pilots);
combs = combs(sort(ia));
amounts_pilots =  floor(N_carrier./combs);

monteCarloRuns = 1e2;
SNR_dB = 20;
Nps=amounts_pilots;

rng(5);
channel_Seeds = randi(2^16,length(amounts_pilots),monteCarloRuns);


channel.Seed = 0;
channel.NRxAnts = 1;
channel.DelayProfile= 'EPA'; %  EPA/EVA/ETU
channel.DopplerFreq = 0;
channel.MIMOCorrelation = 'Low';
channel.SamplingRate = 4e7;
channel.InitTime = 0;
channel.InitPhase = "Random";

if reg_pilot
    scenario_loop = length(combs);
    NMSEs=zeros(4,length(combs));
    BERs=zeros(4,length(combs));
else
    scenario_loop = length(Nps);
    NMSEs=zeros(4,length(Nps));
    BERs=zeros(4,length(Nps));
end

for kk = 1:scenario_loop

    %% Параметры полосы
    allCarriers = linspace(1,Nfft,Nfft); % индексы поднесущих
    if reg_pilot
        comb = combs(kk);
        pilot_step = comb; % comb
        %Регулярные пилоты
        amount_pilots = floor(N_carrier/pilot_step);   % кол-во пилотов на всех Nfft отсчётах
        pilotCarriers = allCarriers(1:pilot_step:N_carrier);
        amount_pilots = length(pilotCarriers);
    else 
        % Случайные пилоты
        amount_pilots = Nps(kk);
        num_values = amount_pilots;
        range_start = 1;
        range_end = N_carrier;
        pilotCarriers = sort(randperm(range_end - range_start + 1, num_values));
        pilot_step = pilotCarriers(3)-pilotCarriers(2); % comb
    end

    if pilot_step == 1
        Percent_pilot = 100;
        
        allCarriers = linspace(1,Nfft,Nfft); % индексы поднесущих
        amount_pilots = round(Percent_pilot/100*N_carrier);   % кол-во пилотов на всех Nfft отсчётах
        pilot_step = floor(N_carrier/amount_pilots);
        pilotCarriers = allCarriers(1:pilot_step:N_carrier-1);
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
    amp_pilots = 2*max_amplitude ;        % пилотное значение
    pilotValues = zeros(1,length(pilotCarriers));
%     pilotValues(1:end) = amp_pilots*exp(1i*0);
    pilotValues(1:2:end) = amp_pilots*exp(1i*0);
    pilotValues(2:2:end) = amp_pilots*exp(1i*pi);
    pilotValues = repmat(pilotValues', 1, N_symb);
    
    %% 1. Чтение файла
    if pilot_step ~= 1
    Size_Buffer = Amount_ODFM_SpF*Amount_OFDM_Frames*amount_data_carriers*bps;
    input_bits = file_reader(File, Size_Buffer);
    
    %% 2. Scrambling
%     % Scrambling
%     Register = [1 0 0 1 0 1 0 1 0 0 0 0 0 0 0];
%     sc_bits_matrix=zeros(Amount_ODFM_SpF*amount_data_carriers*bps,Amount_OFDM_Frames);
%     
%     % Для каждого нового OFDM кадра сбрасывать состояние регистра РСЛОС до
%     % начального
%     for i = 1:Amount_OFDM_Frames
%         start_index = (i-1)*Amount_ODFM_SpF*amount_data_carriers*bps+1;
%         end_index = i*Amount_ODFM_SpF*amount_data_carriers*bps;
%         
%         if i == Amount_OFDM_Frames
%             sc_bits_matrix(:,i) = Scrambler(Register, input_bits(start_index:end));
%         else
%             sc_bits_matrix(:,i) = Scrambler(Register, input_bits(start_index:end_index));
%         end
%     end
%     sc_bits = sc_bits_matrix(:).';


    %% 2. Mapping
    [TX_IQ,pad] = mapping(input_bits,Constellation);
   

    %% 3. Формирование полосы
    OFDM_mapped_carriers = OFDM_map_carriers(TX_IQ, N_symb, Nfft, dataCarriers,pilotCarriers, pilotValues);

    else
        OFDM_mapped_carriers = zeros(Nfft,N_symb);
        OFDM_mapped_carriers(pilotCarriers,:) = pilotValues;
    end
    %% 4. OFDM-модуляция: переход во временную область + добавление CP
    Tx_OFDM_Signal_matrix = OFDM_modulator(OFDM_mapped_carriers, T_Guard);
    Tx_OFDM_Signal=Tx_OFDM_Signal_matrix(:);

    %% АБГШ
    [Tx_OFDM_Signal_noised,~] = Noise(SNR_dB, Tx_OFDM_Signal);

    %% 5. Многолучевое распространение
    NMSEs_jj=zeros(4,length(monteCarloRuns));
    BERs_jj=zeros(4,length(monteCarloRuns));

    % создаём временные переменные для parfor
    NMSEs_jj_ls=zeros(1,monteCarloRuns);
    NMSEs_jj_mmse=zeros(1,monteCarloRuns);
    NMSEs_jj_mp=zeros(1,monteCarloRuns);
    NMSEs_jj_omp=zeros(1,monteCarloRuns);

    % можно переключить в parfor, только нужно закомментировать все, что
    % под verbose
    for jj = 1:monteCarloRuns     
        local_channel = channel;
%         channel_Seed = randi(2^16);
        local_channel.Seed = channel_Seeds(kk,jj);
        [rx_signal, info] = lteFadingChannel(local_channel,Tx_OFDM_Signal_noised);

        h_t = lteFadingChannel(local_channel, [1; zeros(Nfft-1,1)]);
        H_f = fft(h_t);
        %% 8. Оценка канала
        % Power delay profile of the channel
        if verbose 
            figure()
            stem(info.PathSampleDelays, abs(info.PathGains(1,:)).^2, 'LineWidth',2);
            xlabel("Delays in samples")
            ylabel("CIR power")
            title(strcat(channel.DelayProfile, ' seed= ', num2str(channel.Seed)))
            set(gca, 'Fontsize', 20);
            grid on;
        end
        
        %% 6. Разделение сигнала на OFDM-символы
        rx_signal =  reshape(rx_signal,[Nfft+T_Guard, N_symb]);
%         random_rx_signal =  reshape(random_rx_signal,[Nfft+T_Guard, N_symb]);
         %% 7. OFDM-Демодуляция: удаление CP + переход в частотную область - получение точек созвездия
        RX_OFDM_mapped_carriers = OFDM_demodulator(rx_signal, T_Guard);

        H_est_LS_l = LS_CE(RX_OFDM_mapped_carriers,pilotValues,pilotCarriers,N_carrier);
%         h_t_mmse = ifft(H_est_LS_l);
        h_t_mmse=h_t(1:N_carrier).';
        H_est_MMSE = MMSE_CE(RX_OFDM_mapped_carriers,pilotValues,pilotCarriers,Nfft,N_carrier,h_t_mmse,SNR_dB);
        h_t_est = ifft(H_est_LS_l(1:N_carrier));
        
        %Данные для MP/OMP
        F = dftmtx(Nfft);
        if reg_pilot
            F = F(:,1:ceil(Nfft/comb));
        end
        P = zeros(amount_pilots, Nfft);
        for i=1:amount_pilots
            P(i,pilotCarriers(i)) = 1;
        end
        sensing_matrix = P*F;
        Y = RX_OFDM_mapped_carriers(pilotCarriers,1)./pilotValues(:,1);
        
        [H_est_MP, h_t_MP] = MP_estimate(Y, sensing_matrix, Nfft, length(info.PathSampleDelays));
        [H_est_OMP,h_t_OMP,kk1] =  OMP_estimate(Y, sensing_matrix, Nfft, length(info.PathSampleDelays),SNR_dB);
        

        H_freq_useful = H_f(1:N_carrier).';
        H_est_LS_l_useful = H_est_LS_l(1:N_carrier);
        H_est_OMP_useful = H_est_OMP(1:N_carrier);
        H_est_MP_useful = H_est_MP(1:N_carrier);
        H_est_MMSE_useful = H_est_MMSE(1:N_carrier);

        NMSEs_jj_ls(jj) = (H_freq_useful-H_est_LS_l_useful)*(H_freq_useful-H_est_LS_l_useful)'/N_carrier; %LS
        NMSEs_jj_mmse(jj)= (H_freq_useful-H_est_MMSE_useful)*(H_freq_useful-H_est_MMSE_useful)'/N_carrier; %MMSE
        NMSEs_jj_mp(jj) = (H_freq_useful-H_est_MP_useful)*(H_freq_useful-H_est_MP_useful)'/N_carrier; %MP
        NMSEs_jj_omp(jj) = (H_freq_useful-H_est_OMP_useful)*(H_freq_useful-H_est_OMP_useful)'/N_carrier; %OMP
        
        % Закомментировать при запуске parfor
        if verbose 
            close all; % закрыть дополнительные окна 
            figure()
            subplot(2,2,1)
            plot(abs(H_f(1:N_carrier)), 'LineWidth',2)
            hold on;
            plot(abs(H_est_LS_l(1:N_carrier)), 'LineWidth',2)
            legend('True Channel','LS');
            xlabel("frequency sample")
            ylabel("|CFR|")
            grid on;
            set(gca, 'Fontsize', 20);
            if reg_pilot
                title(strcat("MSE=",num2str(mean(NMSEs_jj_ls)),' comb=',num2str(comb)));
            else
                title(strcat("MSE=",num2str(mean(NMSEs_jj_ls)),' Np=',num2str(amount_pilots)));
            end

            subplot(2,2,2)
            plot(abs(H_f(1:N_carrier)), 'LineWidth',2)
            hold on;
            plot(abs(H_est_MMSE(1:N_carrier)), 'LineWidth',2)
            legend('True Channel','MMSE');
            xlabel("frequency sample")
            ylabel("|CFR|")
            grid on;
            set(gca, 'Fontsize', 20);
            title(strcat("MSE=",num2str(mean(NMSEs_jj_mmse))));

            subplot(2,2,3)
            plot(abs(H_f(1:N_carrier)), 'LineWidth',2)
            hold on;
            plot(abs(H_est_MP(1:N_carrier)), 'LineWidth',2)
            xlabel("frequency sample")
            ylabel("|CFR|")
            legend('True Channel','MP');
            grid on;
            set(gca, 'Fontsize', 20);
            title(strcat("MSE=",num2str(mean(NMSEs_jj_mp))));

            subplot(2,2,4)
            plot(abs(H_f(1:N_carrier)), 'LineWidth',2)
            hold on;
            plot(abs(H_est_OMP(1:N_carrier)), 'LineWidth',2)
            xlabel("frequency sample")
            ylabel("|CFR|")
            legend('True Channel','OMP');
            set(gca, 'Fontsize', 20);
            grid on;
            title(strcat("MSE=",num2str(mean(NMSEs_jj_omp))));
            
            figure()
            stem(info.PathSampleDelays, abs(info.PathGains(1,:))/max(abs(info.PathGains(1,:))), 'LineWidth',2)
            hold on
            stem(abs(h_t_est)/max(abs(h_t_est)), 'LineWidth',2)
            legend("Ground Truth CIR", "CIR LS estimation")
            set(gca, 'Fontsize', 20);
            grid on;
            title(strcat("F\_s =", num2str(channel.SamplingRate ), " Hz"))
        end
        %% 9. Эквализация
            ls_RX_OFDM_mapped_carriers = equalize_signal(RX_OFDM_mapped_carriers, H_est_LS_l, N_carrier);
            mmse_RX_OFDM_mapped_carriers = equalize_signal(RX_OFDM_mapped_carriers, H_est_MMSE, N_carrier);
            mp_RX_OFDM_mapped_carriers = equalize_signal(RX_OFDM_mapped_carriers, H_est_MP, N_carrier);
            omp_RX_OFDM_mapped_carriers = equalize_signal(RX_OFDM_mapped_carriers, H_est_OMP, N_carrier);

            equalized = [[ls_RX_OFDM_mapped_carriers]; [mmse_RX_OFDM_mapped_carriers]; [mp_RX_OFDM_mapped_carriers]; [omp_RX_OFDM_mapped_carriers]];
            if pilot_step ~= 1
            %% 10. Получение точек созвездия
                for ll=1:4
                    step = size(ls_RX_OFDM_mapped_carriers,1);
                    RX_IQ = get_payload(equalized((ll-1)*step+1:ll*step,:),dataCarriers);
                    RX_IQ = RX_IQ(:).';
                    
                    %% 11. Demapping
                    output_bits = demapping(pad, RX_IQ, Constellation);
%                     
                    %% 10. DeScrambler
%                     dsc_bits_matrix=zeros(Amount_ODFM_SpF*amount_data_carriers*bps,Amount_OFDM_Frames);
%                     
%                     % Для каждого нового OFDM кадра сбрасывать состояние регистра РСЛОС до
%                     % начального
%                     for i = 1:Amount_OFDM_Frames
%                         start_index = (i-1)*Amount_ODFM_SpF*amount_data_carriers*bps+1;
%                         end_index = i*Amount_ODFM_SpF*amount_data_carriers*bps;
%                         
%                         if i == Amount_OFDM_Frames
%                             dsc_bits_matrix(:,i) = DeScrambler(Register, output_bits(start_index:end));
%                         else
%                             dsc_bits_matrix(:,i) = DeScrambler(Register, output_bits(start_index:end_index));
%                         end
%                     end
%                     dsc_bits = dsc_bits_matrix(:).';
%     
                    %% 13. Проверка полученных данных
                    BERs_jj(ll,jj) = BER_func(input_bits, output_bits);
                end
            end
    end


    BERs(1,kk) = mean(BERs_jj(1,:));
    BERs(2,kk) = mean(BERs_jj(2,:));
    BERs(3,kk) = mean(BERs_jj(3,:));
    BERs(4,kk) = mean(BERs_jj(4,:));


    NMSEs(1,kk) = mean(NMSEs_jj_ls); %LS
    NMSEs(2,kk) = mean(NMSEs_jj_mmse); %MMSE
    NMSEs(3,kk) = mean(NMSEs_jj_mp); %MP
    NMSEs(4,kk) = mean(NMSEs_jj_omp); %OMP

    disp(num2str(kk)+"### NMSE (LS)="+ NMSEs(1,kk)+"; NMSE (MMSE)="+ NMSEs(2,kk)+"; NMSE (MP)="+ NMSEs(3,kk)+"; NMSE (OMP)="+ NMSEs(4,kk));
end
%% Результаты

if reg_pilot
    figure();
    for i=1:size(NMSEs,1)
        plot(amounts_pilots, NMSEs(i,:), 'LineWidth', 2);
        hold on;
    end
    grid on;
    xlabel("amount of pilots")
    ylabel("NMSE")
    set(gca, 'Fontsize', 20)
    legend("LS", "MMSE","MP","OMP", 'Location','northwest');
%     legend("LS", "MMSE",'Location','northwest');
    title(strcat('NMSE(Np) for different channel estimation method (comb-like pilots)', ...
        'channel F\_s= ', num2str(channel.SamplingRate), '\newline',...
        ' Nfft=', num2str(Nfft), ...
        ' N\_carrier=',num2str(N_carrier), ...
        ' channel\_model: ',num2str(channel.DelayProfile)));
    set(gca,'XDir','reverse')

    figure();
    for i=1:size(BERs,1)
        plot(amounts_pilots, BERs(i,:), 'LineWidth', 2);
        hold on;
    end
    xlabel("amount of pilots")
    ylabel("BER")
    set(gca, 'Fontsize', 20);
    grid on;
    legend("LS", "MMSE","MP","OMP", 'Location','northwest');
%     legend("LS", "MMSE",'Location','northwest');
    title(strcat("BER(Np) for channel estimation (comb-like pilots) ", ...
        'channel F\_s= ', num2str(channel.SamplingRate), '\newline',...
        ' Nfft=', num2str(Nfft), ...
        ' N\_carrier=',num2str(N_carrier),...
        ' channel\_model: ',num2str(channel.DelayProfile)));
    set(gca,'XDir','reverse')

else
    figure();
    for i=1:size(NMSEs,1)
        plot(Nps, NMSEs(i,:), 'LineWidth', 2);
        hold on;
    end
    grid on;
    xlabel("amount of pilots")
    ylabel("NMSE")
    set(gca, 'Fontsize', 20)
    legend("LS", "MMSE","MP","OMP", 'Location','northwest');
%     legend("MP","OMP", 'Location','northwest');
    title(strcat('NMSE(Np) for different channel estimation method ', ...
        'channel F\_s= ', num2str(channel.SamplingRate),  '\newline',...
        ' Nfft=', num2str(Nfft), ...
        ' N\_carrier=',num2str(N_carrier), ...
        ' channel\_model: ',num2str(channel.DelayProfile)));
    set(gca,'XDir','reverse')
    
    figure();
    for i=1:size(BERs,1)
        plot(Nps, BERs(i,:), 'LineWidth', 2);
        hold on;
    end
    xlabel("amount of pilots")
    ylabel("BER")
    set(gca, 'Fontsize', 20);
    grid on;
    legend("LS", "MMSE","MP","OMP", 'Location','northwest');
%     legend("MP","OMP", 'Location','northwest');
    title(strcat("BER(Np) for channel estimation ", ...
        'channel F\_s= ', num2str(channel.SamplingRate), '\newline',...
        ' Nfft=', num2str(Nfft), ...
        ' N\_carrier=',num2str(N_carrier),...
        ' channel\_model: ',num2str(channel.DelayProfile)));
    set(gca,'XDir','reverse')
end