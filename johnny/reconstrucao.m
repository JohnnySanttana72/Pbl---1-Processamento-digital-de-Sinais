clc;
clear all;
pkg load signal;

Fs = 1e3; % Frequência de amostragem de 1kHz
Ts = 1/Fs; % Período de amostragem de 1 milissegundo
f0 = 10; % Frequência da primeira senoide em Hz
f1 = 50; % Frequência da segunda senoide em Hz
f2 = 520; % Frequência da terceira senoide em Hz
fc = 500; % Frequência de corte que atende à regra de Nyquist

t = [0:Ts:0.1]; % Valores para o eixo que representa o tempo

x0 = 3*sin(2*pi*f0*t); % primeira senoide
x1 = 5*sin(2*pi*f1*t + pi/4); % segunda senoide
x2 = 5*sin(2*pi*f2*t); % terceira senoide

sinal_original = x0 + x1 + x2; % Soma das senoides

% Aplicação do Janelamento de Hamming
janela_hamming = hamming(length(t))';
sinal_janelado = sinal_original .* janela_hamming;

ordem = 5; % Aumento da ordem do filtro
lowcut = fc*2/Fs; % Frequência de corte

[b,a] = butter(ordem, lowcut, 'low'); % Filtro Butterworth de ordem 10

sinal_lowfilter = filtfilt(b, a, sinal_janelado);

% Equalização para realçar a frequência f1 = 50 Hz
frequencia_desejada = f1 / (Fs/2); % Frequência normalizada
ganho_desejado = 10; % Ganho desejado para a frequência

% Projeto do filtro de equalização
[b_eq, a_eq] = butter(ordem, frequencia_desejada, 'low');
[h_eq, frequencies_eq] = freqz(b_eq, a_eq, length(sinal_lowfilter), Fs);

% Aplicação do filtro de equalização
sinal_equalizado = filter(b_eq, a_eq, sinal_lowfilter);

% Amostragem PAM natural
amplitude_pam = 1;
ganho = amplitude_pam / max(abs(sinal_equalizado));
sinal_pam = ganho * sinal_equalizado;

% Aplicação da FFT
sinal_fft = fft(sinal_equalizado);
N = length(sinal_fft);
frequencies = (0:N-1)*(Fs/N);

% Reconstrução do sinal original
sinal_reconstruido = ifft(sinal_fft);

% Plotagem de Todos os Sinais
figure;

subplot(3,3,1);
plot(t, x0, 'b');
title('Senoide x0 - Original');
xlabel('Tempo');
ylabel('Amplitude');

subplot(3,3,2);
plot(t, x1, 'r');
title('Senoide x1 - Original');
xlabel('Tempo');
ylabel('Amplitude');

subplot(3,3,3);
plot(t, x2, 'g');
title('Senoide x2 - Original');
xlabel('Tempo');
ylabel('Amplitude');

subplot(3,3,4);
plot(t, sinal_original, 'm');
title('Soma das Senoides Originais');
xlabel('Tempo');
ylabel('Amplitude');

subplot(3,3,5);
plot(t, sinal_janelado, 'c');
title('Sinal Janelado (Após Janelamento)');
xlabel('Tempo');
ylabel('Amplitude');

subplot(3,3,6);
plot(t, sinal_lowfilter, 'k');
title('Sinal Filtrado');
xlabel('Tempo');
ylabel('Amplitude');

subplot(3,3,7);
plot(t, sinal_equalizado, 'y');
title('Sinal Equalizado');
xlabel('Tempo');
ylabel('Amplitude');

subplot(3,3,8);
stem(t, sinal_pam, '.');
title('Sinal Amostrado (PAM)');
xlabel('Tempo');
ylabel('Amplitude');

% Plotagem do sinal no Domínio da Frequência (FFT)
figure;
stem(frequencies, abs(sinal_fft), '.');
title('Sinal no Domínio da Frequência (FFT)');
xlabel('Frequência (Hz)');
ylabel('Amplitude FFT');

% Plotagem do sinal reconstruído
figure;
plot(t, real(sinal_reconstruido), 'g');
title('Sinal Reconstruído');
xlabel('Tempo');
ylabel('Amplitude');

