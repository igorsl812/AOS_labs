clc;
clear;
% close all;
%% формирование сетки

Nfft = 2048;
Gi = 1/8;

all_sc = -Nfft/2 : Nfft/2 - 1; 
all_sc_shifted = all_sc+Nfft/2+1;
zero_sc = Nfft/2+1;

Gi_left = 1:(Nfft*Gi/2+1);
Gi_right = (Nfft-Nfft*Gi/2+1):Nfft;

pilot_sc = max(Gi_left):2048:min(Gi_right)-1;

exc_sc = [Gi_left,Gi_right,zero_sc,pilot_sc];
Y_inf = all_sc_shifted(~ismember(all_sc_shifted,exc_sc));
%% Формирование парметров битов
M_inf = 16;
k_inf = log2(M_inf);
n_inf_bits = length(Y_inf);
M_pilot = 4;
k_pilot = log2(M_pilot);
n_pilot_bits = length(pilot_sc);


%% параметры кадра
CP_len = Nfft/4;
T_symb = 100e-6;
Fs = Nfft/T_symb;
N_symbs = 1;

alg = 'pts'; % 'none' 'clp' 'slm' 'pts'
max_clip_iters = 3;
Nslm = 100;
Nblk = 256;
%% канал 
SNR_vector = 1:1; 
BER_vector = zeros(1, length(SNR_vector));

N_frames = 100;   
rng(1)
for idx_ber = 1:length(SNR_vector)
    SNR = SNR_vector(idx_ber);
    
    total_errors = 0;
    total_bits = 0;
    PAPR = [];

    for frame = 1:N_frames

        bits = randi([0,1],n_inf_bits*k_inf,N_symbs);
        qamsymb = qammod(bits,M_inf,'InputType', 'bit','UnitAveragePower', true);

        % Формируем матрицу переданных бит в нужном формате
        bits_dem = qamdemod(qamsymb.', 16,'OutputType', 'bit','UnitAveragePower', true);
        bits_dem_cut = bits_dem(1:k_inf*(N_symbs),:);
        % формирование OFDM символа
        tx_fd_sig = zeros(N_symbs,Nfft);
        tx_fd_sig(:,Y_inf) = qamsymb.';

        

        % подавление PAPR
        switch (alg)
            case 'none'
                OFDM_symb = ifft(ifftshift(tx_fd_sig,2),[],2);
            case 'clp'
                OFDM_symb = ifft(ifftshift(tx_fd_sig,2),[],2);
                for i = 1:max_clip_iters
                    OFDM_to_clip = OFDM_symb(:);
                    clip_lvl = 0.03; %mean(abs(OFDM_to_clip))
    
                    idx = abs(OFDM_symb) > clip_lvl;
                    OFDM_symb(idx) = clip_lvl * (OFDM_symb(idx) ./ abs(OFDM_symb(idx)));
    
                    OFDM_clip = fftshift(fft(OFDM_symb,[],2),2);
    
                    OFDM_clip(:,Gi_right) =  0;
                    OFDM_clip(:,Gi_left) =  0;
                    OFDM_clip(:,zero_sc) =  0;
    
                    OFDM_symb = ifft(ifftshift(OFDM_clip,2),[],2);   
                end
            case 'slm'
                best_papr = inf;
                for j = 1:Nslm
                    slm_bits = 2 * randi([0,1],1,Nfft) - 1;
                    OFDM_slm = tx_fd_sig.*slm_bits;

                    OFDM_time_slm = ifft(ifftshift(OFDM_slm,2),[],2);
                    OFDM_time_slm_vec = OFDM_time_slm(:);

                    current_papr = 10*log10(max(abs(OFDM_time_slm_vec).^2)/mean(abs(OFDM_time_slm_vec).^2));
                    if current_papr < best_papr
                        best_papr = current_papr;
                        OFDM_symb = OFDM_time_slm;
                    end
                end
            case 'pts'
                OFDM_symb = zeros(1, Nfft);
                M = zeros(Nblk,Nfft);
                for n = 1:Nblk
                    M(n,(Nfft/Nblk*(n-1))+1:(Nfft/Nblk*(n))) = tx_fd_sig(1, (Nfft/Nblk*(n-1))+1:(Nfft/Nblk*(n)));
                    M_time = ifft(ifftshift(M,2),[],2);
               
                    cand_1 = OFDM_symb + M_time(n,:).*1;
                    cand_2 = OFDM_symb + M_time(n,:).*-1;
                    PAPR1 = 10*log10(max(abs(cand_1).^2)/mean(abs(cand_1).^2));
                    PAPR2 = 10*log10(max(abs(cand_2).^2)/mean(abs(cand_2).^2));

                    if PAPR1 < PAPR2
                        OFDM_symb = cand_1;
                    else
                        OFDM_symb = cand_2;
                    end
                end
                
        end

        % добавление ЦП
        OFDM_symb_cp = [OFDM_symb(:,end-CP_len+1:end),OFDM_symb];

        %преобразование в один сигнал для передачи
        OFDM_vector = OFDM_symb_cp.';
        OFDM_vector = OFDM_vector(:);

        %PAPR
        PAPR_current = (abs(OFDM_vector).^2) / mean(abs(OFDM_vector).^2);
        PAPR = [PAPR; 10*log10(PAPR_current)];

        % АБГШ       
        Eb =  mean(abs(qammod(0:M_inf-1, M_inf)).^2) / log2(M_inf) / Nfft;
        N0 = Eb / 10^(SNR/10);
        sigma = sqrt(N0/2);
        noise = sigma*randn(size(OFDM_vector)) + 1j*sigma*randn(size(OFDM_vector));

        OFDM_noise = OFDM_vector+noise;

        % приемник
        OFDM_clean = reshape(OFDM_noise,size(OFDM_symb_cp,2),size(OFDM_symb_cp,1)).';

        % удаление цп
        OFDM_rx = OFDM_clean(:,CP_len+1:end);

        % переход в частоты
        OFDM_rx_fd = fftshift(fft(OFDM_rx,[],2),2);
        
        % Извлечение данных
        rx_symbs = OFDM_rx_fd(:,Y_inf);
        rx_pilots = OFDM_rx_fd(:,pilot_sc);
        rx_data_demap = qamdemod(rx_symbs, 16,'OutputType', 'bit','UnitAveragePower', true);
        
        % Считаем ошибки текущего фрейма
        [numErrors, ~] = biterr(bits_dem_cut(:), rx_data_demap(:));

        total_errors = total_errors + numErrors;
        total_bits = total_bits + numel(bits_dem_cut);
    end
    BER_vector(idx_ber) = total_errors / total_bits;
end


%% Рассчет PAPR и CCDF
X_db = 0:1:12; 
CCDF = zeros(length(X_db),1);
for k=X_db
    CCDF(k+1) = mean(PAPR>k);
    CCDF(CCDF==0) = 1e-5;
end
figure
semilogy(X_db,CCDF)
xlim([0 12])
ylim([1e-4 1])
xlabel('X, дБ');
ylabel('P(PAPR > X)');
%% спеткр
% figure
% frequencyLimits = [-1 1]*pi;
% pspectrum(OFDM_vector,'FrequencyLimits',frequencyLimits);
% title('Спектр мощности OFDM-сигнала');
%% Теоретическая кривая
% P_const = mean(abs(qammod(0:M_inf-1, M_inf)).^2);   
% EbN0_vector = SNR_vector - 10*log10(P_const);
% 
% % Теоретический BER для 16-QAM в канале АБГШ
% BER_theory = berawgn(EbN0_vector, 'qam', M_inf);
% 
% % Построение графика
% figure;
% semilogy(SNR_vector, BER_vector, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
% hold on;
% semilogy(SNR_vector, BER_theory, 'r:', 'LineWidth', 1.5);
% hold off;
% grid on;
% 
% xlim([min(SNR_vector), max(SNR_vector)]);
% ylim([1e-4, 1]); 
% 
% xlabel('SNR, дБ');
% ylabel('BER');
% title('Кривая помехоустойчивости OFDM-системы');
% legend('Экспериментальный BER', 'Теоретический BER', 'Location', 'SouthWest');
