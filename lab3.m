clc;
clear;
close all;
%% формирование сетки

Nfft = 2048;
Gi = 1/16;

all_sc = -Nfft/2 : Nfft/2 - 1; 
all_sc_shifted = all_sc+Nfft/2+1;
zero_sc = Nfft/2+1;

% ИСПРАВЛЕНО: приведение к int32 убирает Warning об операторе двоеточия
Gi_left = 1:int32(Nfft*Gi/2+1);
Gi_right = int32(Nfft-Nfft*Gi/2+1):Nfft;

pilot_sc = max(Gi_left):16:min(Gi_right)-1;

exc_sc = [Gi_left,Gi_right,zero_sc,pilot_sc];
Y_inf = all_sc_shifted(~ismember(all_sc_shifted,exc_sc));
%% Формирование парметров битов
M_inf = 16;
k_inf = log2(M_inf);
n_inf_bits = length(Y_inf);
M_pilot = 4;
k_pilot = log2(M_pilot);
n_pilot_bits = length(pilot_sc);


%% Константы кадра
CP_len = Nfft/4;
T_symb = 100e-6;
Fs = Nfft/T_symb;
N_symbs = 5;

time_shift = 600;
t_delay = zeros(time_shift,1);

punkt = 0; % выбор пункта лабы
delta = 20;

%% канал 
SNR_vector = 1:20; 
BER_vector = zeros(1, length(SNR_vector));

N_frames = 100;   

for idx_ber = 1:length(SNR_vector)
    SNR = SNR_vector(idx_ber);
    
    total_errors = 0;
    total_bits = 0;

    for frame = 1:N_frames
        %% формирование пилотов
        pilot_bits = randi([0,1],n_pilot_bits*k_pilot,N_symbs);
        qampilot = qammod(pilot_bits,M_pilot,'InputType', 'bit','UnitAveragePower', true);

        bits = randi([0,1],n_inf_bits*k_inf,N_symbs);
        qamsymb = qammod(bits,M_inf,'InputType', 'bit','UnitAveragePower', true);

        % Формируем матрицу переданных бит в нужном формате
        bits_dem = qamdemod(qamsymb.', 16,'OutputType', 'bit','UnitAveragePower', true);
        bits_dem_cut = bits_dem(1:k_inf*(N_symbs-1),:);
        % формирование OFDM символа
        tx_fd_sig = zeros(N_symbs,Nfft);
        tx_fd_sig(:,pilot_sc) = qampilot.';
        tx_fd_sig(:,Y_inf) = qamsymb.';

        % Сигнал во времени
        OFDM_symb = ifft(ifftshift(tx_fd_sig,2),[],2);

        % добавление ЦП
        OFDM_symb_cp = [OFDM_symb(:,end-CP_len+1:end),OFDM_symb];

        %преобразование в один сигнал для передачи
        OFDM_vector = OFDM_symb_cp.';
        OFDM_vector = OFDM_vector(:);

        % Добавляем временной сдвиг
        OFDM_vector_shift = [t_delay; OFDM_vector];

        % АБГШ       
        Eb =  mean(abs(qammod(0:M_inf-1, M_inf)).^2) / log2(M_inf) / Nfft;
        N0 = Eb / 10^(SNR/10);
        sigma = sqrt(N0/2);
        noise = sigma*randn(size(OFDM_vector_shift)) + 1j*sigma*randn(size(OFDM_vector_shift));

        OFDM_noise = OFDM_vector_shift+noise;
        if punkt == 1
            true_delay = time_shift+1+delta;
        else
            T_frame = Nfft + CP_len;
            N_avg = N_symbs - 1;
            N_search = Nfft - 1;

            x = OFDM_noise;

            shifts = (0:N_avg-1)' * T_frame + (1:N_search);

            z = x(1:end-Nfft) .* conj(x(Nfft+1:end));

            Sz = [0; cumsum(z)];

            R = Sz(shifts + CP_len) - Sz(shifts);

            % Когерентное накопление
            if N_symbs == 2
                ACF = abs(R);
            else
                ACF = abs(sum(R, 1));
            end
            % Приведение к строке
            ACF = reshape(ACF, 1, []);

            [~, true_delay] = max(ACF);
        end
                
        % компенсация задержки
        OFDM_compens = OFDM_noise(true_delay:(true_delay+(Nfft+CP_len)*(N_symbs-1)-1));
        % приемник
        OFDM_clean = reshape(OFDM_compens,size(OFDM_symb_cp,2),size(OFDM_symb_cp,1)-1).';

        % удаление цп
        OFDM_rx = OFDM_clean(:,CP_len+1:end);

        % переход в частоты
        OFDM_rx_fd = fftshift(fft(OFDM_rx,[],2),2);
        
        % Извлечение данных
        rx_symbs = OFDM_rx_fd(:,Y_inf);
        rx_pilots = OFDM_rx_fd(:,pilot_sc);
        H = rx_pilots./tx_fd_sig(1:N_symbs-1,pilot_sc);
        if N_symbs == 2
            H_all = interp1(double(pilot_sc),H,1:Nfft,"linear","extrap");
        else
            H_all = interp1(double(pilot_sc),mean(H),1:Nfft,"linear","extrap");
        end
        rx_data_equalized = OFDM_rx_fd./H_all;
        rx_data_equalized = rx_data_equalized(:,Y_inf);
        rx_data_demap = qamdemod(rx_data_equalized, 16,'OutputType', 'bit','UnitAveragePower', true);
        
        % Считаем ошибки текущего фрейма
        [numErrors, ~] = biterr(bits_dem_cut(:), rx_data_demap(:));

        total_errors = total_errors + numErrors;
        total_bits = total_bits + numel(bits_dem_cut);
    end
    BER_vector(idx_ber) = total_errors / total_bits;
end

%% АКФ при фикс снр
figure;
plot(linspace(1,length(ACF),length(ACF)),abs(ACF));
xlabel("time_shift");
ylabel("ACF");
%% Теоретическая кривая
P_const = mean(abs(qammod(0:M_inf-1, M_inf)).^2);   
EbN0_vector = SNR_vector - 10*log10(P_const);

% Теоретический BER для 16-QAM в канале АБГШ
BER_theory = berawgn(EbN0_vector, 'qam', M_inf);

% Построение графика
figure;
semilogy(SNR_vector, BER_vector, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
hold on;
semilogy(SNR_vector, BER_theory, 'r:', 'LineWidth', 1.5);
hold off;
grid on;

xlim([min(SNR_vector), max(SNR_vector)]);
ylim([1e-4, 1]); 

xlabel('SNR, дБ');
ylabel('BER');
title('Кривая помехоустойчивости OFDM-системы');
legend('Экспериментальный BER', 'Теоретический BER', 'Location', 'SouthWest');
