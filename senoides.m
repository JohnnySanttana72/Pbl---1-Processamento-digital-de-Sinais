% Carrega o pacote de processamento de sinal  pacote fornece várias funções úteis para análise, filtragem e processamento de sinais.

pkg load signal

% Definição dos parâmetros das senoides
A1 = 1;     % Amplitude da primeira senoide
f1 = 10;    % Frequência da primeira senoide

A2 = 0.5;   % Amplitude da segunda senoide
f2 = 3;     % Frequência da segunda senoide

A3 = 0.8;   % Amplitude da terceira senoide
f3 = 10;    % Frequência da terceira senoide

% Definição do vetor de tempo
%Criamos um vetor de tempo t que varia de 0 a 2*pi (um período completo de uma senoide) com um intervalo de 0.01.
%Este vetor de tempo será usado para calcular os valores das senoides ao longo do tempo.
t = 0:0.01:2*pi;  % De 0 a 2*pi, com incremento de 0.01

% Cálculo das senoides
senoide1 = A1 * sin(2*pi*f1*t);
senoide2 = A2 * sin(2*pi*f2*t);
senoide3 = A3 * sin(2*pi*f3*t);

% Soma das senoides
soma_senoides = senoide1 + senoide2 + senoide3;

% Definição do filtro passa-baixas Butterworth

%A ordem do filtro indica o quão rápido o filtro pode atenuar as frequências fora da banda de passagem.
%Uma ordem maior geralmente resulta em uma resposta mais nítida na frequência de corte,
%mas também pode aumentar a complexidade computacional do filtro.
ordem = 4; % Ordem do filtro

%A frequência de corte (fc) é a frequência em Hertz na qual o filtro começa a atenuar as frequências
%do sinal de entrada. Qualquer componente de frequência acima dessa frequência de corte é atenuado pelo filtro.
fc = 5;    % Frequência de corte

%A frequência de amostragem (fs) é a taxa na qual o sinal foi originalmente amostrado. É importante para a normalização da
%frequência de corte e para garantir que o filtro esteja funcionando corretamente com base nas propriedades do sinal original.
fs = 100;  % Frequência de amostragem

% Normalização da frequência de corte
%Este passo é necessário porque as funções de filtro em muitos softwares,
%requerem a frequência de corte normalizada em relação à frequência de Nyquist.
fc_norm = fc / (fs/2);


% Coeficientes do filtro passa-baixas Butterworth
%Com todos esses parâmetros definidos, calculamos os coeficientes do filtro passa-baixas Butterworth usando a função butter.
[b, a] = butter(ordem, fc_norm, 'low');
%b são os coeficientes dos termos de realimentação do filtro.
%a são os coeficientes dos termos de entrada do filtro.


% Aplicação do filtro na soma das senoides
%A função "filter" é usada para convolucionar o sinal de entrada com os coeficientes do filtro.
soma_filtrada = filter(b, a, soma_senoides);

% Calcula a FFT da soma filtrada para o domínio da frequência
soma_fft = fft(soma_filtrada);
N = length(soma_fft); %N é o comprimento da FFT, que é igual ao número de pontos

%é um vetor que contém as frequências correspondentes aos coeficientes da FFT.
%Estas frequências são calculadas a partir da frequência de amostragem fs e são distribuídas uniformemente de 0 até a frequência de Nyquist (fs/2)
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

