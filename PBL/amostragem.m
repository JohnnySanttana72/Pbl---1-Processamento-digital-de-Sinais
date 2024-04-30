clear all; close all; clc

%% dados do sinal
[y_audio, f] = audioread('audio_entrada.wav');
xt = y_audio(:,1)';
t = linspace(0,length(xt)/f,length(xt));

%% Etapa 1 => Filtro de Entrada
% filtro = designfilt('lowpassiir', 'PassbandFrequency', 1000, 'StopbandFrequency', 1500, 'PassbandRipple', 1, 'StopbandAttenuation', 60, 'SampleRate', f);
wp = (2*pi*200)/f; % Frequência da banda passante normalizada
ws = (2*pi*1500)/f; % Frequência da banda de parada normalizada

[N, wc] = buttord(wp, ws, 3, 60); % ripple de banda passante 3 db e atenuação na banda de parada 60 Hz
[b, a] = butter(N, wc, 'low'); % filtro butterworth de ordem N com frequência de corte wc
sinal_filtrado = filter(b, a, xt);

[h, w] = freqz(b, a, (0:0.01:pi)); % recuperação parâmetros

% Converção de radianos/amostra para Hz/plot
freq = w * f / (2*pi);

figure(1)
plot(freq, abs(h), [200, 200], [0, 1], [1500, 1500], [0, 0.08357])
grid on;
xlabel('$f$(Hz)','Interpreter','LaTex');
ylabel('Magnitude','Interpreter','LaTex');
legend('', 'w_{p}','w_{s}')
title('Filtro Anti-Alising (Passa-Baixa)');


%% Etapa 2 => Amostragem PAM
fs = 4000; % Frequência de amostragem em Hz
% ffc = f/fs; % Pega o quanto fs é menor que f

% Ts = 1/fs;
% N = length(t)/ffc; % A quantidade de amostras é 1/fcc menor do que o tamanho de amostras iniciais
% n = 0 : 1 : N-1;
% t_amostrado = 0 : Ts : n(N)*Ts;
% sinal_amostrado1 = zeros(1);
% k = 0;
% 
% for i = 1 : 1 : N
%    sinal_amostrado1(i) = sinal_filtrado(k + 1);
%    k = k + ffc;
% end

sinal_retangular = square(2*pi*fs*t);
n = length(sinal_retangular);

for i = 1 : n
   if (sinal_retangular(i) <= 0)
       sinal_retangular(i) = 0;
   else
       sinal_retangular(i) = 1;
   end
end

sinal_amostrado = sinal_filtrado .* sinal_retangular;

%% Faz o plot dos resultados

% Figura 1: x(t) filtrado
figure(2)
subplot(3,1,1)
plot(t, sinal_filtrado);
xlabel('$t$','Interpreter','LaTex');
ylabel('$x(t)$','Interpreter','LaTex');
title('Sinal x(t) filtrado');

% Figura 1: x(nt) amostrado
subplot(3,1,2)
% plot(n, sinal_amostrado,'.'); 
plot(t, sinal_amostrado); 
xlabel('$nT$','Interpreter','LaTex')
ylabel('$x(nT)$','Interpreter','LaTex')
title('Sinal x(nT) amostrado');

% Figura 1: x(t) e x[t]
subplot(3,1,3)
plot(t, sinal_filtrado);
hold;
% plot(t_amostrado, sinal_amostrado, 'o');
plot(t, sinal_amostrado);
xlabel('$t$','Interpreter','LaTex')
ylabel('$x[nT_s],x(t)$','Interpreter','LaTex')
title('Sinal x(t) e x(nT)');

% Figura 2: Espectro de x(t) sem o filtro
figure(3)
subplot(3,1,1);

y=fft(xt); grid on;
yaux=fliplr(y(1,2:end)); % inverte da esquerda para a direita os elementos da posição 2 até 12800
X=[yaux y];
X(1,1:length(X)/4)=0; % adição de zeros da posição 1 até 32000
X(1,3*length(X)/4:end)=0; % adição de zeros da posição 9600 até 12800 (último termo)
length(X);
omega=0:f/length(y):f-(f/length(y));
waux=-fliplr(omega(1,2:end));
w=[waux omega];
length(w);

stem(w,abs(2*X/length(t)), '.');
xlim([-6000, 6000]); % limita a exibição do gráfico de -6000 à 6000 Hz
xlabel('$f$(Hz)','interpreter','latex');
ylabel('Magnitude');
title('Espectro do sinal X(jw) sem o filtro');

% Figura 2: Espectro de x(t) com o filtro
figure(3)
subplot(3,1,2);

y=fft(sinal_filtrado); grid on;
yaux=fliplr(y(1,2:end)); % inverte da esquerda para a direita os elementos da posição 2 até 12800
X=[yaux y];
X(1,1:length(X)/4)=0; % adição de zeros da posição 1 até 32000
X(1,3*length(X)/4:end)=0; % adição de zeros da posição 9600 até 12800 (ultimo termo)
length(X);
omega=0:f/length(y):f-(f/length(y));
waux=-fliplr(omega(1,2:end));
w=[waux omega];
length(w);

stem(w,abs(2*X/length(t)), '.');
xlim([-6000, 6000]); % limita a exibição do gráfico de -6000 à 6000 Hz
xlabel('$f$(Hz)','interpreter','latex');
ylabel('Magnitude');
title('Espectro do sinal X(jw) filtrado');

% Figura 2: Espectro de x(nT)
figure(3)
subplot(3,1,3);

y=fft(sinal_amostrado); grid on;
yaux=fliplr(y(1,2:end)); % % inverte da esquerda para a direita os elementos da posição 2 até 12800
y(1,2:4)
fliplr(y(1,2:4))
X=[yaux y];
X(1,1:length(X)/4)=0; % adição de zeros da posição 1 até 32000
X(1,3*length(X)/4:end)=0; % adição de zeros da posição 9600 até 12800 (ultimo termo)
length(X);
omega=0:fs/length(y):fs-(fs/length(y));
waux=-fliplr(omega(1,2:end));
w=[waux omega];
length(w);

stem(w,abs(2*X/length(t)), '.');
xlim([-6000, 6000]); % limita a exibição do gráfico de -6000 à 6000 Hz
xlabel('$f$(Hz)','interpreter','latex');
ylabel('Magnitude');
title('Espectro do sinal X(jw) filtrado e amostrado');

%% Mudança da taxa de amostragem (acho que não vamos precisar do processo de subamostragem do sinal)

% M = 2
% sinal_subamostrado = downsample(sinal_amostrado,M);
% t_subamostrado = downsample(t_amostrado,M);
% Td = M*Ts;
% fd = 1/Td;
% Nd = N/M;
% nd = 0 : 1 : Nd-1;
% figure(3)
% plot(t_amostrado, sinal_amostrado, t_subamostrado, sinal_subamostrado);
% xlabel('$t_d$','Interpreter','LaTex','FontSize',14);
% ylabel('$y(nT), y(nT_d)$','Interpreter','LaTex','FontSize',14);


%% Reconstrução
% sinal_sinc = sinc(fd*(ones(length(nd),1)*t-(nd*Td)'*ones(1, length(t))));

[N, wc] = buttord(wp, ws, 3, 60); % ripple de banda passante 3 db e atenuação na banda de parada 60 Hz
[b, a] = butter(N, wc, 'low'); % filtro butterworth de ordem N com frequência de corte wc

% sinal_reconstruido = sinal_subamostrado*sinal_sinc;

% figure(4)
% subplot(3,1,1);
% plot(t, sinal_filtrado);
% xlabel('$t$','Interpreter','LaTex','FontSize',14);
% ylabel('$x_c(t)$','Interpreter','LaTex','FontSize',14);
% 
% subplot(3,1,2);
% plot(t, sinal_reconstruido);
% xlabel('$t$','Interpreter','LaTex');
% ylabel('$x_r(t)$','Interpreter','LaTex');
% 
% subplot(3,1,3);
% plot(t, sinal_filtrado);
% hold;
% plot(t, sinal_reconstruido);
% xlabel('$t$','Interpreter','LaTex')
% ylabel('$x_r(t),x_c(t)$','Interpreter','LaTex')
% 
% plot(t, sinal_reconstruido);

% audiowrite('audio_reconstruido.wav', sinal_reconstruido, f);
% sound(y_audio, f);
% pause(length(xt)/f);
% pause(2);
% sound(sinal_reconstruido, f);
