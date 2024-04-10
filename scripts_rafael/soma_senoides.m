clc;
clear all;
pkg load signal;

Fs = 1e6; % Fequência de amostragem de 1MHz
Ts = 1/Fs; % Período de amostragem de 1 microsegundos

f0 = 20; % Frequência da primeira senoide em Hz
f1 = 60; % Frequência da segunda senoide em Hz
f2 = 501; % Frequência da terceira senoide em Hz
f3 = 560; % Frequência da quarta senoide em Hz

t = [0:Ts:0.1]; % valores para o eixo que representa o tempo

y = [1:1:length(t)];


x0 = 3*sin(2*pi*f0*t);
x1 = 5*sin(2*pi*f1*t + pi/4);
x2 = 7*cos(2*pi*f2*t);
x3 = 2*cos(2*pi*f3*t + pi/4);

x = x0 + x1 + x2 + x3;

ordem = 5;
lowcut = 300*2/Fs; % frequência de corte

[b,a] = butter(ordem, lowcut, 'low'); % filtro Butterworth de ordem 5

lowfilter = filtfilt(b, a, x);

subplot(2,2,1);
plot(t,x);
title('Sinal Senoidal x(t)');
xlabel('t');
ylabel('x(t)');

subplot(2,2,2);
plot(t,lowfilter, 'r');
title('Sinal Senoidal y(t)');
xlabel('t');
ylabel('y(t)');

N = length(x)
F = Fs*(-N/2:N/2-1)/N; % Frequencia normalizada
X = fftshift(fft(x/N)); % Transformada de Fourier deslocado com a amplitude normalizada

N_filter = length(lowfilter)
f_filter = Fs*(-N_filter/2:N_filter/2-1)/N_filter; % Frequencia normalizada
Y = fftshift(fft(lowfilter/N_filter)); % Transformada de Fourier do sinal filtrado com a amplitude normalizada

subplot(2,2,3);
plot(F, abs(X));
xlim ([-1000 1000]); % limites de visualização do espectro de -1000 à 1000 hertz
title('Sinal x(t) Domínio da Frequência (FFT)');
xlabel('F(Hz)');
ylabel('Amplitude (dB)');

subplot(2,2,4);
plot(f_filter, abs(Y), 'r');
xlim ([-1000 1000]); % limites de visualização do espectro de -1000 à 1000 hertz
title('Sinal y(t) Domínio da Frequência (FFT) após o filtro');
xlabel('F(Hz)');
ylabel('Amplitude (dB)');

