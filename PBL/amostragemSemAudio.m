clear all; close all; clc

%% dados do sinal
f = 8000; % Frequência entrada em Hz
[xt,fc,phi,t] = sinal(2,10,3000);

%% Etapa 1 => Filtro de Entrada
% filtro = designfilt('lowpassiir', 'PassbandFrequency', 1000, 'StopbandFrequency', 1500, 'PassbandRipple', 1, 'StopbandAttenuation', 60, 'SampleRate', f);
wp = 2/f*1000; % Frequência da banda passante normalizada
ws = 2/f*1500; % Frequência da banda de parada normalizada

[N, wc] = buttord(wp, ws, 1, 60); % ripple de banda passante 3 db e atenuação na banda de parada 60 Hz
[b, a] = butter(N, wc, 'low'); % filtro butterworth de ordem N com frequência de corte wc
sinal_filtrado = filter(b,a, xt);
%fvtool(filtro);

[h, w] = freqz(b, a, (0:0.01:pi)); % recuperação parâmetros
%
% Converção de radianos/amostra para Hz/plot
freq = w * f / (2*pi);

% Figura 2: Gáfico do filtro passa-baixa
figure(1)
plot(freq, abs(h), [1000, 1000], [0, 1.2], [1500, 1500], [0, 1.2])
grid on;
xlabel('$f$(Hz)','Interpreter','LaTex');
ylabel('Magnitude','Interpreter','LaTex');
legend('', 'w_{p}','w_{s}')
title('Filtro Anti-Alising (Passa-Baixa)');


%% Etapa 2 => Amostragem
fs = 10000; % Frequência de amostragem em Hz
% ffc = f/fs; % Pega o quanto fs é menor que f
%
% Ts = 1/fs;
% N = length(t)/ffc; % A quantidade de amostras é 1/fcc menor do que o tamanho de amostras iniciais
% n = [0 : 1 : N-1];
% t_amostrado = [0 : Ts : n(N)*Ts];
% sinal_amostrado = zeros(1);
% k = 0;
%
% for i = 1 : 1 : N
%   sinal_amostrado(i) = sinal_filtrado(k + 1);
%   k = k + ffc;
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

% Faz o plot dos resultados
% Figura 2:  x(t) sem filtro no domínio do tempo e seu espectro emfrequência
figure(2);
subplot(2,1,1);
plot(t, xt);
xlabel('$t$(s)','Interpreter','LaTex')
ylabel('Amplitude')
title('$x(t) = x_{f}(t) + x_{r}(t)$ ','Interpreter','LaTex','FontSize',14);
axis([0 0.05 -inf inf]);

% Aplicar a FFT ao sinal com ruído, espelhar o espectro e manter o espectro central
subplot(2,1,2);
y=fft(xt); grid on;
yaux=fliplr(y(1,2:end)); % inverte da esquerda para a direita os elementos da posição 2 até 31999
X=[yaux y];
X(1,1:length(X)/4)=0; % adição de zeros da posição 1 até 800
X(1,3*length(X)/4:end)=0; % adição de zeros da posição 2400 até 32000
length(X)
omega=0:f/length(y):f-(f/length(y));
waux=-fliplr(omega(1,2:end));
w=[waux omega];
length(w);
stem(w,abs(2*X/length(t)), '.');
xlabel('$f$(Hz)','interpreter','latex');
ylabel('Magnitude');
title('$|X(j\omega)|$ ','Interpreter','LaTex','FontSize',14);
axis([-8000 8000 -inf inf]);

% Figura 2: x(t) filtrado e seu espectro em frequência
figure(3);
subplot(2,1,1);
plot(t, sinal_filtrado);
xlabel('$t$(s)','Interpreter','LaTex')
ylabel('Amplitude')
title('$x_{cc}(t)$','Interpreter','LaTex','FontSize',14);
axis([0 0.05 -inf inf]);

% Aplicar a FFT ao sinal filtrado, espelhar o espectro e manter o espectro central
subplot(2,1,2);
y=fft(sinal_filtrado); grid on;
yaux=fliplr(y(1,2:end)); % inverte da esquerda para a direita os elementos da posição 2 até 31999
X=[yaux y];
X(1,1:length(X)/4)=0; % adição de zeros da posição 1 até 800
X(1,3*length(X)/4:end)=0; % adição de zeros da posição 2400 até 32000
length(X);
omega=0:f/length(y):f-(f/length(y));
waux=-fliplr(omega(1,2:end));
w=[waux omega];
length(w);
stem(w,abs(2*X/length(t)), '.');
xlabel('$f$(Hz)','interpreter','latex');
ylabel('Magnitude');
title('$|X_{cc}(j\omega)|$','Interpreter','LaTex','FontSize',14);
axis([-8000 8000 -inf inf]);

% Plot sinal amostrado e discretizado
figure(4);
subplot(2,1,1);
plot(t, sinal_amostrado);
xlabel('$n$','Interpreter','LaTex')
ylabel('Amplitude')
title('$x_{cc}(t) = x_{cc}(nT) = x[n]$ amostrado e discretizado', 'Interpreter','LaTex','FontSize',14);
% axis([0 50 -inf inf]);
axis([0 0.05 -inf inf]);

% Aplicar a FFT ao sinal amostrado, espelhar o espectro e manter o espectro central
subplot(2,1,2);
y=fft(sinal_amostrado); grid on;
yaux=fliplr(y(1,2:end)); % inverte da esquerda para a direita os elementos da posição 2 até 31999
X=[yaux y];
X(1,1:length(X)/4)=0; % adição de zeros da posição 1 até 800
X(1,3*length(X)/4:end)=0; % adição de zeros da posição 2400 até 32000
length(X);
omega=0:fs/length(y):fs-(fs/length(y));
waux=-fliplr(omega(1,2:end));
w=[waux omega];
length(w);
stem(w,abs(2*X/length(t)), '.');
xlabel('$f$(Hz)','interpreter','latex');
ylabel('Magnitude');
title('$|X(e^{j\omega})|$','Interpreter','LaTex','FontSize',14);
axis([-8000 8000 -inf inf])

% Plot sinal filtrado e sinal amostrado
figure(5);
plot(t, sinal_filtrado);
hold;
plot(t, sinal_amostrado);
xlabel('$t$(s)','Interpreter','LaTex')
ylabel('Amplitude')
title('$x_{cc}(t)$ e $x_{amostrado}(t)$','Interpreter','LaTex','FontSize',14);
axis([0 0.05 -inf inf]);
legend('x_{cc}','x_{amostrado}');

%% Mudança da taxa de amostragem

%M = 2;
%sinal_subamostrado = downsample(sinal_amostrado,M);
%t_subamostrado = downsample(t_amostrado,M);
%Td = M*Ts;
%fd = 1/Td;
%Nd = N/M;
%nd = 0 : 1 : Nd-1;
%figure(5)
%plot(t_amostrado, sinal_amostrado, t_subamostrado, sinal_subamostrado);
%xlabel('$t$(s)','Interpreter','LaTex')
%ylabel('Amplitude')
%title('$y(nT)$, $y(nT_d)$','Interpreter','LaTex','FontSize',14);
%axis([0 0.05 -inf inf]);
%legend('y(nT)','y(nT_d)');

% Reconstrução
wp = 2/f*1000; % Frequência da banda passante normalizada
ws = 2/f*1500; % Frequência da banda de parada normalizada

[N, wc] = buttord(wp, ws, 1, 60); % ripple de banda passante 3 db e atenuação na banda de parada 60 Hz
[b, a] = butter(N, wc, 'low'); % filtro butterworth de ordem N com frequência de corte wc

% pw = pulsewidth(sinal_filtrado, t);
% bw = bandwidth(sinal_filtrado);

% ganho = db2mag(max(sinal_amostrado)/max(sinal_filtrado));
sinal_reconstruido = (1.55823)*filter(b,a, sinal_amostrado);

% sinal_reconstruido = sinal_amostrado*sinc(fd*(ones(length(nd),1)*t-(nd*Td)'*ones(1, length(t))));

% Plot sinal reconstruído
figure(6);
subplot(2,1,1);
plot(t, sinal_reconstruido);
xlabel('$t$(s)','Interpreter','LaTex')
ylabel('Amplitude')
title('$y_{cc}(t)$','Interpreter','LaTex','FontSize',14);
axis([0 0.05 -inf inf]);

% Aplicar a FFT ao sinal reconstruído, espelhar o espectro e manter o espectro central
subplot(2,1,2);
y=fft(sinal_reconstruido); grid on;
yaux=fliplr(y(1,2:end)); % inverte da esquerda para a direita os elementos da posição 2 até 31999
X=[yaux y];
X(1,1:length(X)/4)=0; % adição de zeros da posição 1 até 800
X(1,3*length(X)/4:end)=0; % adição de zeros da posição 2400 até 32000
length(X);
omega=0:fs/length(y):fs-(fs/length(y));
waux=-fliplr(omega(1,2:end));
w=[waux omega];
length(w);
stem(w,abs(2*X/length(t)), '.');
xlabel('$f$(Hz)','interpreter','latex');
ylabel('Magnitude');
title('$|Y_{cc}(j\omega)|$','Interpreter','LaTex','FontSize',14);
axis([-8000 8000 -inf inf])

% Plot sinal reconstruído e filtrado juntos
figure(7);
plot(t, sinal_filtrado);
hold;
plot(t, sinal_reconstruido);
xlabel('$t$(s)','Interpreter','LaTex')
ylabel('Amplitude')
title('$y_{cc}(t)$ e $x_{cc}(t)$ ','Interpreter','LaTex','FontSize',14);
axis([0 0.05 -inf inf]);
legend('x_{cc}','y_{cc}');

figure(8);
subplot(2,1,1);
plot(t, sinal_reconstruido);
hold;
plot(t, xt);
xlabel('$t$(s)','Interpreter','LaTex')
ylabel('Amplitude')
title('$y_{cc}(t)$ e $x(t)$ ','Interpreter','LaTex','FontSize',14);
axis([0 0.05 -inf inf]);
legend('y','x');

subplot(2,1,2);
stem(w,abs(2*X/length(t)), '.');
hold;
y=fft(xt); grid on;
yaux=fliplr(y(1,2:end)); % inverte da esquerda para a direita os elementos da posição 2 até 31999
X=[yaux y];
X(1,1:length(X)/4)=0; % adição de zeros da posição 1 até 800
X(1,3*length(X)/4:end)=0; % adição de zeros da posição 2400 até 32000
length(X);
omega=0:f/length(y):f-(f/length(y));
waux=-fliplr(omega(1,2:end));
w=[waux omega];
length(w);
stem(w,abs(2*X/length(t)), '.');
xlabel('$f$(Hz)','interpreter','latex');
ylabel('Magnitude');
title('$|Y(j\omega)|$ e $|X(j\omega)|$','Interpreter','LaTex','FontSize',14);
axis([-8000 8000 -inf inf])
legend('Y','X');

sound(sinal_reconstruido, f);

