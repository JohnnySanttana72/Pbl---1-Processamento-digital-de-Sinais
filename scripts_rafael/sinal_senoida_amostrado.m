clc;
clear all;
pkg load signal;

Fs = 1e3; % Fequência de amostragem de 1MHz
Ts = 1/Fs; % Período de amostragem de 1 microsegundos
f0 = 10; % Frequência da primeira senoide em Hz
f1 = 50; % Frequência da segunda senoide em Hz
f2 = 520; % Frequência da segunda senoide em Hz
fc = 500; % Frequência dde corte que vai de encontro a regra de Nyquist

t = [0:Ts:0.1]; % valores para o eixo que representa o tempo


x0 = 3*sin(2*pi*f0*t); % primeira senoide
x1 = 5*sin(2*pi*f1*t + pi/4); % segunda senoide
x2 = 5*sin(2*pi*fc*t); % terceira senoide

sinal = x0 + x1 + x2;

ordem = 5;
lowcut = fc*2/Fs % frequência de corte calculo igual a 1

[b,a] = butter(ordem, lowcut, 'low'); % filtro Butterworth de ordem 5
##[c, d] = cheby1(ordem, 20, 500*2/Fs); % filtro Chebyshev 1 de ordem 5

sinal_lowfilter = filtfilt(b, a, sinal);

##sinal_lowfilter2 = filtfilt(c, d, sinal);

subplot(2,2,1);
plot(t, sinal_lowfilter);
title('Sinal Senoidal x(t) filtrado');
xlabel('t');
ylabel('x(t)');

subplot(2,2,2);
stem(t(1:100),sinal_lowfilter(1:100), 'r'); % plotar o sinal amostrado usando stem
hold on;
title('Sinal Senoidal x(t) filtrado amostrado');
xlabel('n');
ylabel('y(t)');


N_filter = length(sinal_lowfilter); % quantidade total de amostras do sinal filtrado
f_filter = Fs*(-N_filter/2:N_filter/2-1)/N_filter; % Frequencia normalizada
Y = fftshift(fft(sinal_lowfilter/N_filter)); % Transformada de Fourier do sinal filtrado com a amplitude normalizada

N_sampling = length(sinal_lowfilter(1:100)); % quantidade total do sinal que será amostrado
f_sampling = Fs*(-N_sampling/2:N_sampling/2-1)/N_sampling; % Frequencia normalizada
Y_sampling = fftshift(fft(sinal_lowfilter(1:100)/N_sampling)); % Transformada de Fourier do sinal amostrado


subplot(2,2,3);
stem(f_filter, abs(Y)); % plotar stem do sinal no domínio da frequência
xlim ([-1000 1000]); % limites de visualização do espectro de -1000 à 1000 hertz
title('Sinal x(t) filtrado contínuo Domínio da Frequência (FFT) amostrado');
xlabel('F(Hz)');
ylabel('Amplitude');

subplot(2,2,4);
stem(f_sampling, abs(Y_sampling), 'r'); % plotar stem do sinal no domínio da frequência
xlim ([-1000 1000]); % limites de visualização do espectro de -1000 à 1000 hertz
title('Sinal x(t) filtrado discretizado no Domínio da Frequência (FFT) também amostrado');
xlabel('F(Hz)');
ylabel('Amplitude');

