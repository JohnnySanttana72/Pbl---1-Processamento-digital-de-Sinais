clc;
clear all;
pkg load signal;

Fs = 1e3; % Fequência de amostragem de 1MHz
Ts = 1/Fs; % Período de amostragem de 1 microsegundos
f0 = 10; % Frequência da primeira senoide em Hz
f1 = 50; % Frequência da segunda senoide em Hz
f2 = 520; % Frequência da segunda senoide em Hz
fc = 500; % Frequência de corte que vai de encontro à regra de Nyquist

t = [0:Ts:0.1]; % valores para o eixo que representa o tempo

x0 = 3*sin(2*pi*f0*t); % primeira senoide
x1 = 5*sin(2*pi*f1*t + pi/4); % segunda senoide
x2 = 5*sin(2*pi*fc*t); % terceira senoide

sinal = x0 + x1 + x2;

% Janelamento com a janela de Hamming
janela_hamming = hamming(length(t))';
sinal_janelado = sinal .* janela_hamming;

ordem = 5;
lowcut = fc*2/Fs; % frequência de corte

[b,a] = butter(ordem, lowcut, 'low'); % filtro Butterworth de ordem 5

sinal_lowfilter = filtfilt(b, a, sinal_janelado);

N_filter = length(sinal_lowfilter); % quantidade total de amostras do sinal filtrado
f_filter = Fs*(-N_filter/2:N_filter/2-1)/N_filter; % Frequencia normalizada
Y = fftshift(fft(sinal_lowfilter/N_filter)); % Transformada de Fourier do sinal filtrado com a amplitude normalizada

N_sampling = length(sinal_lowfilter(1:100)); % quantidade total do sinal que será amostrado
f_sampling = Fs*(-N_sampling/2:N_sampling/2-1)/N_sampling; % Frequencia normalizada
Y_sampling = fftshift(fft(sinal_lowfilter(1:100)/N_sampling)); % Transformada de Fourier do sinal amostrado

% Plotagem das senoides individuais
figure;
subplot(2,2,1);
plot(t, x0, 'b');
title('Senoide 1 - Tempo');
xlabel('Tempo');
ylabel('Amplitude');

subplot(2,2,2);
plot(t, x1, 'r');
title('Senoide 2 - Tempo');
xlabel('Tempo');
ylabel('Amplitude');

subplot(2,2,3);
plot(t, x2, 'g');
title('Senoide 3 - Tempo');
xlabel('Tempo');
ylabel('Amplitude');

% Plotagem da soma das senoides
subplot(2,2,4);
plot(t, sinal, 'm');
title('Soma das Senóides - Tempo');
xlabel('Tempo');
ylabel('Amplitude');

% Plotagem do sinal filtrado e sinal PAM
figure;
subplot(2,1,1);
plot(t, sinal_lowfilter, 'k', 'LineWidth', 2);
title('Sinal Senoidal x(t) filtrado - Tempo');
xlabel('Tempo');
ylabel('Amplitude');

% Plotagem do sinal PAM amostrado
subplot(2,1,2);
stem(t(1:100), sinal_lowfilter(1:100), 'b', 'Marker', 'o');
title('Sinal PAM Amostrado');
xlabel('Tempo');
ylabel('Amplitude');

% Plotagem do espectro do sinal amostrado
figure;
subplot(2,1,1);
stem(f_filter, abs(Y), 'b', 'Marker', 'o');
title('Espectro do Sinal Amostrado - PAM');
xlabel('Frequência (Hz)');
ylabel('Amplitude FFT');

subplot(2,1,2);
stem(f_sampling, abs(Y_sampling), 'r', 'Marker', 'o');
title('Espectro do Sinal Amostrado Janelado');
xlabel('Frequência (Hz)');
ylabel('Amplitude FFT');

% Plotagem do janelamento
figure;
subplot(2,1,1);
plot(t, sinal, 'b', 'LineWidth', 1.5);
hold on;
plot(t, janela_hamming, 'r--', 'LineWidth', 1.5);
title('Sinal Senoidal x(t) e Janelamento de Hamming');
xlabel('Tempo');
ylabel('Amplitude');
legend('Sinal x(t)', 'Janela de Hamming');

subplot(2,1,2);
plot(t, sinal_janelado, 'm', 'LineWidth', 1.5);
title('Sinal Senoidal x(t) Janelado com Hamming');
xlabel('Tempo');
ylabel('Amplitude');

