clc;
clear;
close all;
%% формирование сетки

Nfft = 2048;
Gi = 1/16;

all_sc = -Nfft/2 : Nfft/2 - 1; 
all_sc_shifted = all_sc+Nfft/2+1;
zero_sc = Nfft/2+1;

Gi_left = 1:(Nfft*Gi/2+1);
Gi_right = (Nfft-Nfft*Gi/2+1):Nfft;

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


%% параметры кадра
CP_len = Nfft/4;
T_symb = 100e-6;
Fs = Nfft/T_symb;
N_symbs = 10;

time_shift = 600;
t_delay = zeros(time_shift,1);
f_n = 0;
alg = 2; %1 - мин 2 - ЦП
punkt = 1; % 1 - без синхронизации
delta = -20;
%% пилоты
rng(1);
pilot_bits = randi([0,1],n_pilot_bits*k_pilot,N_symbs);
qampilot = qammod(pilot_bits,M_pilot,'InputType', 'bit','UnitAveragePower', true);

%% канал 
SNR_vector = 1:30; 
BER_vector = zeros(1, length(SNR_vector));

N_frames = 100;   

for idx_ber = 1:length(SNR_vector)
    SNR = SNR_vector(idx_ber);
    
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
        tx_fd_sig(:,pilot_sc) = qampilot.';
        tx_fd_sig(:,Y_inf) = qamsymb.';

        % Сигнал во времени
        OFDM_symb = ifft(ifftshift(tx_fd_sig,2),[],2);

        % добавление ЦП
        if alg == 1
            A = OFDM_symb(1,10:9+Nfft/4);
            sync_meen = [A,A,-A,-A];
            OFDM_symb_cp = OFDM_symb;
        else
            OFDM_symb_cp = [OFDM_symb(:,end-CP_len+1:end),OFDM_symb];
        end

        %преобразование в один сигнал для передачи
        OFDM_vector = OFDM_symb_cp.';
        OFDM_vector = OFDM_vector(:);

        % Добавляем временной сдвиг
        if alg == 1
            OFDM_vector_shift = [t_delay;sync_meen(:); OFDM_vector];
        else
            OFDM_vector_shift = [t_delay;OFDM_vector];
        end

        % частотная расстройка
        
        delta_f = f_n * (Fs / Nfft);
        t = (0 : length(OFDM_vector_shift) - 1).' / Fs;
        OFDM_vector_shift = OFDM_vector_shift .* exp(1j * 2 * pi * delta_f * t);

        % АБГШ       
        Eb =  mean(abs(qammod(0:M_inf-1, M_inf)).^2) / log2(M_inf) / Nfft;
        N0 = Eb / 10^(SNR/10);
        sigma = sqrt(N0/2);
        noise = sigma*randn(size(OFDM_vector_shift)) + 1j*sigma*randn(size(OFDM_vector_shift));

        OFDM_noise = OFDM_vector_shift+noise;

        if alg == 1

            x = OFDM_noise;
            L = Nfft/4;

            N_search = length(x) - Nfft + 1;
            shifts = 1:N_search;

            z = x(1:end-L) .* conj(x(L+1:end));

            Sz = [0; cumsum(z)];

            P1 = Sz(shifts+L) - Sz(shifts);

            P2 = Sz(shifts+3*L) - Sz(shifts+2*L);

            P = P1 + P2;

            E = [0; cumsum(abs(x).^2)];

            E1 = E(shifts+L)   - E(shifts);
            E2 = E(shifts+2*L) - E(shifts+L);
            E3 = E(shifts+3*L) - E(shifts+2*L);
            E4 = E(shifts+4*L) - E(shifts+3*L);

            M = abs(P).^2 ./ ((E1+E3).*(E2+E4) + eps);

            [~, sync_position] = max(M);
            % --- БЛОК ЧАСТОТНОЙ СИНХРОНИЗАЦИИ ---
            % Извлекаем комплексное значение пика
            CF_max = P(sync_position);
            
            % Находим фазу пика (угол поворота)
            theta = angle(CF_max);
            
            % Оценка нормированного частотного сдвига по формуле из методички
            f_n_est = -4*theta / (2 * pi);
            
            % Переводим оценку обратно в Герцы
            delta_f_est = f_n_est * (Fs / Nfft);
            
            % Домножаем ВЕСЬ принятый сигнал во времени на компенсирующую экспоненту
            t_rx = (0 : length(OFDM_noise) - 1).' / Fs;
            OFDM_noise = OFDM_noise .* exp(-1j * 2 * pi * delta_f_est * t_rx);

        else
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
                    ACF = R;
                else
                    ACF = sum(R, 1);
                end
                % Приведение к строке
                ACF = reshape(ACF, 1, []);

                [~, true_delay] = max(abs(ACF));
                % --- БЛОК ЧАСТОТНОЙ СИНХРОНИЗАЦИИ ---
                % Извлекаем комплексное значение пика
                CF_max = ACF(true_delay);
                
                % Находим фазу пика (угол поворота)
                theta = angle(CF_max);
                
                % Оценка нормированного частотного сдвига по формуле из методички
                f_n_est = -theta / (2 * pi);
                
                % Переводим оценку обратно в Герцы
                delta_f_est = f_n_est * (Fs / Nfft);
                
                % Домножаем ВЕСЬ принятый сигнал во времени на компенсирующую экспоненту
                t_rx = (0 : length(OFDM_noise) - 1).' / Fs;
                OFDM_noise = OFDM_noise .* exp(-1j * 2 * pi * delta_f_est * t_rx);
        
            end
        end       
        % компенсация задержки
        if alg == 2
            OFDM_compens = OFDM_noise(true_delay:(true_delay+(Nfft+CP_len)*(N_symbs-1)-1));
        else
            OFDM_compens = OFDM_noise(sync_position+Nfft:(sync_position+Nfft+(Nfft)*(N_symbs-1)-1));
        end
        % приемник
        OFDM_clean = reshape(OFDM_compens,size(OFDM_symb_cp,2),size(OFDM_symb_cp,1)-1).';

        % удаление цп
        if alg == 1
            OFDM_rx = OFDM_clean;
        else
            OFDM_rx = OFDM_clean(:,CP_len+1:end);
        end
        

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
        rx_data_demap = qamdemod(rx_data_equalized, 16,'OutputType', 'bit','UnitAveragePower', true);% если без экв то = rx_symbs, с экв rx_data_equalized
        
        % Считаем ошибки текущего фрейма
        [numErrors, ~] = biterr(bits_dem_cut(:), rx_data_demap(:));

        total_errors = total_errors + numErrors;
        total_bits = total_bits + numel(bits_dem_cut);
    end
    BER_vector(idx_ber) = total_errors / total_bits;
end

%% АКФ 
% if alg == 2
%     figure;
%     plot(linspace(1,length(ACF),length(ACF)),abs(ACF));
%     xlabel("time_shift");
%     ylabel("ACF");
% else
%     figure;
%     plot(linspace(1,length(M),length(M)),abs(M));
%     xlabel("time_shift");
%     ylabel("Metric_meen");
% end
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
