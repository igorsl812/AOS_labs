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
CP_len = 128;%Nfft/4;
T_symb = 100e-6;
Fs = Nfft/T_symb;
t_delay = zeros(10,1);
num_symbs = 10;
%% канал 
SNR_vector = 30:30; 
BER_vector = zeros(1, length(SNR_vector));
N_frames = 100;   

for idx_ber = 1:length(SNR_vector)
    SNR = SNR_vector(idx_ber);
%     SNR = 30;
    
    total_errors = 0;
    total_bits = 0;

    for frame = 1:100 
        
        bits = randi([0,1],n_inf_bits*k_inf,num_symbs);

        qamsymb = qammod(bits,M_inf,'InputType', 'bit','UnitAveragePower', true);
        % Формируем матрицу переданных бит в нужном формате
        bits_dem = qamdemod(qamsymb.', 16,'OutputType', 'bit','UnitAveragePower', true);

        % формирование OFDM символа
        tx_fd_sig = zeros(num_symbs,Nfft);
        tx_fd_sig(:,Y_inf) = qamsymb.';
        OFDM_symb = ifft(ifftshift(tx_fd_sig,2),[],2);

        % добавление ЦП
        OFDM_symb_cp = [OFDM_symb(:,end-CP_len+1:end),OFDM_symb];%zeros(num_symbs,CP_len) или OFDM_symb(:,end-CP_len+1:end)
        OFDM_vector = OFDM_symb_cp.';
        OFDM_vector = OFDM_vector(:);
        OFDM_W_t_delay = [t_delay;OFDM_vector]; %t_delay;
        % АБГШ       
        SNR_linear = 10^(SNR / 10);
        
        % Измеряем точную среднюю мощность чистого сигнала до добавления ЦП
        P_sig = mean(abs(OFDM_symb).^2,'all'); 
        
        % Рассчитываем необходимую мощность шума
        P_noise = P_sig / SNR_linear;
        
        % Амплитуда комплексного шума
        sigma = sqrt(P_noise / 2);
        
        % Генерируем комплексный гауссовский шум
        noise = sigma * (randn(size(OFDM_W_t_delay)) + 1i * randn(size(OFDM_W_t_delay)));
        %преобразование в один сигнал для передачи
        OFDM_noise = OFDM_W_t_delay+noise;
       
        % синхронизация
        AllDelay = 0:15;
        for j_time_shift = 1:length(AllDelay)+1
            left = OFDM_noise(j_time_shift-1 + (1:CP_len));
            right = OFDM_noise(j_time_shift-1 +(1:CP_len) + Nfft);
            ACF(j_time_shift) = conv(left, fliplr(conj(right)), "valid");
        end
        [~,true_delay] = max(abs(ACF));
        OFDM_clean = reshape(OFDM_noise,size(OFDM_symb_cp,2),size(OFDM_symb_cp,1)).';
        % приемник

        % удаление цп
        OFDM_rx = OFDM_clean(:,CP_len+1:end);

        % переход в частоты
        OFDM_rx_fd = fftshift(fft(OFDM_rx,[],2),2);
        
        % Извлечение данных
        rx_symbs = OFDM_rx_fd(:,Y_inf);
%         scatterplot(rx_symbs(1,:));
        
        rx_data_demap = qamdemod(rx_symbs, 16,'OutputType', 'bit','UnitAveragePower', true);
        % Считаем ошибки текущего фрейма
        [numErrors, ~] = biterr(bits_dem(:), rx_data_demap(:));

        total_errors = total_errors + numErrors;
        total_bits = total_bits + numel(bits);
%         scatterplot(rx_data_equalized(1,:));
    end
%     Сохраняем усредненный BER 
    BER_vector(idx_ber) = total_errors / total_bits;
end



%% Теоретическая кривая

% Коэффициент частотной избыточности
L_data = length(Y_inf) / Nfft; 

EbN0_vector = SNR_vector-10*log10(k_inf)-10*log10(L_data);
EbN0_linear = 10.^(EbN0_vector/10);

% Теоретический BER для 16-QAM в канале АБГШ
BER_theory = berawgn(EbN0_vector, 'qam', M_inf);

% Построение графика

figure;

% 1. текущий эксперимент в канале 
semilogy(SNR_vector, BER_vector, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
hold on;

% 2. Теория для канала  
semilogy(SNR_vector, BER_theory, 'g--', 'LineWidth', 2);


hold off;
grid on;

% Настройка масштаба осей
xlim([min(SNR_vector), max(SNR_vector)]);
ylim([1e-4, 1]); 

xlabel('SNR, дБ');
ylabel('BER');
title('Кривая помехоустойчивости OFDM-системы в канале с АБГШ');

legend('Экспериментальный BER ', ...
       'Теоретический BER ', ...
       'Location', 'SouthWest');

