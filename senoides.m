% Carrega o pacote de processamento de sinal
pkg load signal

% Definição dos parâmetros das senoides
A1 = 1;     % Amplitude da primeira senoide
f1 = 10;    % Frequência da primeira senoide

A2 = 0.5;   % Amplitude da segunda senoide
f2 = 3;     % Frequência da segunda senoide

A3 = 0.8;   % Amplitude da terceira senoide
f3 = 10;    % Frequência da terceira senoide

% Definição do vetor de tempo
t = 0:0.01:2*pi;  % De 0 a 2*pi, com incremento de 0.01

% Cálculo das senoides
senoide1 = A1 * sin(2*pi*f1*t);
senoide2 = A2 * sin(2*pi*f2*t);
senoide3 = A3 * sin(2*pi*f3*t);

% Soma das senoides
soma_senoides = senoide1 + senoide2 + senoide3;

% Definição do filtro passa-baixas Butterworth
ordem = 4; % Ordem do filtro
fc = 5;    % Frequência de corte
fs = 100;  % Frequência de amostragem

% Normalização da frequência de corte
fc_norm = fc / (fs/2);

% Coeficientes do filtro passa-baixas Butterworth
[b, a] = butter(ordem, fc_norm, 'low');

% Aplicação do filtro na soma das senoides
soma_filtrada = filter(b, a, soma_senoides);

% Calcula a FFT da soma filtrada para o domínio da frequência
soma_fft = fft(soma_filtrada);
N = length(soma_fft);
frequencies = (0:N-1)*(fs/N);

% Plotagem das senoides individuais, soma e resultado filtrado
figure;
subplot(2,2,1);
plot(t, senoide1, 'b');
title('Senóide 1 - Tempo');
xlabel('Tempo');
ylabel('Amplitude');

subplot(2,2,2);
plot(t, senoide2, 'r');
title('Senóide 2 - Tempo');
xlabel('Tempo');
ylabel('Amplitude');

subplot(2,2,3);
plot(t, senoide3, 'g');
title('Senóide 3 - Tempo');
xlabel('Tempo');
ylabel('Amplitude');

subplot(2,2,4);
plot(t, soma_senoides, 'm');
hold on;
plot(t, soma_filtrada, 'k', 'LineWidth', 2);
hold off;
title('Soma das Senóides e Resultado Filtrado - Tempo');
xlabel('Tempo');
ylabel('Amplitude');
legend('Soma das Senóides', 'Resultado Filtrado');

% Plotagem do resultado da FFT separado
figure;
subplot(2,1,1);
plot(frequencies, abs(soma_fft));
title('FFT da Soma Filtrada - Frequência');
xlabel('Frequência (Hz)');
ylabel('Amplitude FFT');

% Plotagem da soma filtrada no domínio do tempo
subplot(2,1,2);
plot(t, soma_filtrada, 'k', 'LineWidth', 2);
title('Soma Filtrada no Domínio do Tempo');
xlabel('Tempo');
ylabel('Amplitude');

