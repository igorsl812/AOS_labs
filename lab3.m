clc;
clear;
close all;
%% формирование сетки

Nfft = 2048;
Gi = 1/16;

all_sc = -Nfft/2 : Nfft/2 - 1; 
all_sc_shifted = all_sc+Nfft/2+1;
zero_sc = Nfft/2+1;

Gi_left = 1:Nfft*Gi/2+1;
Gi_right = (Nfft-Nfft*Gi/2+1):Nfft;

exc_sc = [Gi_left,Gi_right,zero_sc];
Y_inf = all_sc_shifted(~ismember(all_sc_shifted,exc_sc));
%% Формирование парметров битов
M_inf = 16;
k_inf = log2(M_inf);
n_inf_bits = length(Y_inf);

%% Константы кадра
CP_len = 16;%Nfft/4;
T_symb = 100e-6;
Fs = Nfft/T_symb;
N_symbs = 10;

time_shift = 15;
t_delay = zeros(time_shift,1);

%% канал 
SNR_vector = 1:30; 
BER_vector = zeros(1, length(SNR_vector));

N_frames = 100;   

for idx_ber = 1:length(SNR_vector)
    SNR = SNR_vector(idx_ber);
%     SNR = 30;
    
    total_errors = 0;
    total_bits = 0;

    for frame = 1:N_frames
        
        bits = randi([0,1],n_inf_bits*k_inf,N_symbs);
        
        qamsymb = qammod(bits,M_inf,'InputType', 'bit','UnitAveragePower', true);
        % Формируем матрицу переданных бит в нужном формате
        bits_dem = qamdemod(qamsymb.', 16,'OutputType', 'bit','UnitAveragePower', true);
        bits_dem_cut = bits_dem(1:k_inf*(N_symbs-1),:);
        % формирование OFDM символа
        tx_fd_sig = zeros(N_symbs,Nfft);
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
        ACF = zeros(Nfft,1);
        % АКФ
        for j_time_shift = 1:Nfft-1
            left = OFDM_noise((j_time_shift-1) + (1:CP_len));
            right = OFDM_noise((j_time_shift-1) +(1:CP_len) + Nfft);
            R_identical = corrcoef(left, right);
            % Выделяем конкретное значение корреляции между ними
            ACF(j_time_shift) = R_identical(1,2);
        end
        [~,true_delay] = max(abs(ACF));
        % figure;
        % plot(linspace(1,Nfft,Nfft),abs(ACF));
        
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
%         scatterplot(rx_symbs(1,:));
        rx_data_demap = qamdemod(rx_symbs, 16,'OutputType', 'bit','UnitAveragePower', true);
        % Считаем ошибки текущего фрейма
        [numErrors, ~] = biterr(bits_dem_cut(:), rx_data_demap(:));

        total_errors = total_errors + numErrors;
        total_bits = total_bits + numel(bits);
%         scatterplot(rx_data_equalized(1,:));
    end
    BER_vector(idx_ber) = total_errors / total_bits;
end

%% Теоретическая кривая

P_const = mean(abs(qammod(0:M_inf-1, M_inf)).^2);   % = 10
EbN0_vector = SNR_vector - 10*log10(P_const);

% Теоретический BER для 16-QAM в канале АБГШ
BER_theory = berawgn(EbN0_vector, 'qam', M_inf);

BER_theory_Rayleigh = berfading(EbN0_vector, 'qam', M_inf,1);
% Построение графика

figure;

% 1. текущий эксперимент 
semilogy(SNR_vector, BER_vector, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
hold on;

% 3. Теория для канала АБГШ (Красный точечный) — ваша проверенная идеальная прямая
semilogy(SNR_vector, BER_theory, 'r:', 'LineWidth', 1.5);

hold off;
grid on;

% Настройка масштаба осей
xlim([min(SNR_vector), max(SNR_vector)]);
ylim([1e-4, 1]); 

xlabel('SNR, дБ');
ylabel('BER');
title('Кривая помехоустойчивости OFDM-системы ');

legend('Экспериментальный BER ', ...
       'Теоретический BER ', ...
       'Location', 'SouthWest');