% Carrega o pacote de processamento de sinal
pkg load signal

% Definição dos parâmetros das senoides
A1 = 1;     % Amplitude da primeira senoide
f1 = 10;    % Frequência da primeira senoide

A2 = 0.5;   % Amplitude da segunda senoide
f2 = 3;     % Frequência da segunda senoide

A3 = 0.8;   % Amplitude da terceira senoide
f3 = 30;    % Frequência da terceira senoide

% Definição do vetor de tempo
t = 0:0.005:4*2*pi;  % De 0 a 4 períodos completos, com incremento de 0.005

% Cálculo das senoides
senoide1 = A1 * sin(2*pi*f1*t);
senoide2 = A2 * sin(2*pi*f2*t);
senoide3 = A3 * sin(2*pi*f3*t);

% Soma das senoides
soma_senoides = senoide1 + senoide2 + senoide3;

% Definição do filtro passa-baixas Butterworth
ordem = 6;  % Ordem do filtro
fc = 20;    % Frequência de corte
fs = 40;    % Frequência de amostragem (não obedecendo ao Teorema de Nyquist)

% Normalização da frequência de corte
fc_norm = fc / (fs/2);

% Coeficientes do filtro passa-baixas Butterworth
[b, a] = butter(ordem, fc_norm, 'low');

% Aplicação do filtro na soma das senoides
soma_filtrada = filter(b, a, soma_senoides);

% Modulação PAM natural no sinal filtrado
amplitude_pam = 1;
ganho = amplitude_pam / max(abs(soma_filtrada));
sinal_pam = ganho * soma_filtrada;

% Calcula a FFT da soma filtrada para o domínio da frequência
soma_fft = fft(sinal_pam);
N = length(soma_fft);
frequencies = (0:N-1)*(fs/N);

% Plotagem das senoides individuais
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

% Plotagem da soma das senoides
subplot(2,2,4);
plot(t, soma_senoides, 'm');
title('Soma das Senóides - Tempo');
xlabel('Tempo');
ylabel('Amplitude');

% Plotagem do sinal filtrado e sinal PAM
figure;
subplot(2,1,1);
plot(t, soma_filtrada, 'k', 'LineWidth', 2);
title('Soma Filtrada - Tempo');
xlabel('Tempo');
ylabel('Amplitude');

% Plotagem do sinal PAM amostrado
subplot(2,1,2);
stem(t, sinal_pam, 'b', 'Marker', 'o');
title('Sinal PAM Amostrado');
xlabel('Tempo');
ylabel('Amplitude');

% Plotagem do espectro do sinal amostrado
figure;
stem(frequencies, abs(soma_fft), 'b', 'Marker', 'o');
title('Espectro do Sinal Amostrado - PAM');
xlabel('Frequência (Hz)');
ylabel('Amplitude FFT');

