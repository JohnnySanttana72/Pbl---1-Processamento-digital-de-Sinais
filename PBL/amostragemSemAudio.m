clear all; close all; clc
%% dados do sinal
f = 8000; % Frequência entrada em Hz
[xt,fc,phi,t] = sinal(2,10,3000);

%% Etapa 1 => Filtro de Entrada
filtro = designfilt('lowpassiir', 'PassbandFrequency', 1000, 'StopbandFrequency', 1500, 'PassbandRipple', 1, 'StopbandAttenuation', 60, 'SampleRate', f);
sinal_filtrado = filter(filtro, xt);
%fvtool(filtro);
%% Etapa 2 => Amostragem
fs = 4000; % Frequência de amostragem em Hz
ffc = f/fs; % Pega o quanto fs é menor que f

Ts = 1/fs;
N = length(t)/ffc; % A quantidade de amostras é 1/fcc menor do que o tamanho de amostras iniciais
n = [0 : 1 : N-1];
t_amostrado = [0 : Ts : n(N)*Ts];
sinal_amostrado = zeros(1);
k = 0;

for i = 1 : 1 : N
   sinal_amostrado(i) = sinal_filtrado(k + 1);
   k = k + ffc;
end

%% Faz o plot dos resultados
% Figura 1: x(t) sem filtro no domínio do tempo e seu espectro em
% frequência
figure(1);
subplot(2,1,1);
plot(t, xt);
xlabel('$t$(s)','Interpreter','LaTex')
ylabel('Amplitude')
title('$x(t) = x_{f}(t) + x_{r}(t)$ ','Interpreter','LaTex','FontSize',14);
axis([0 0.05 -inf inf]);

subplot(2,1,2);
y=fft(xt); grid on;
yaux=fliplr(y(1,2:end));
X=[yaux y];
X(1,1:length(X)/4)=0;
X(1,3*length(X)/4:end)=0;
length(X);
omega=0:f/length(y):f-(f/length(y));
waux=-fliplr(omega(1,2:end));
w=[waux omega];
length(w);
plot(w,abs(2*X/length(t)));
xlabel('$f$(Hz)','interpreter','latex');
ylabel('Magnitude');
title('$|X(j\omega)|$ ','Interpreter','LaTex','FontSize',14);
axis([-8000 8000 -inf inf]);

% Figura 2: x(t) filtrado e seu espectro em frequência
figure(2);
subplot(2,1,1);
plot(t, sinal_filtrado);
xlabel('$t$(s)','Interpreter','LaTex')
ylabel('Amplitude')
title('$x_{cc}(t)$','Interpreter','LaTex','FontSize',14);
axis([0 0.05 -inf inf]);

subplot(2,1,2);
y=fft(sinal_filtrado); grid on;
yaux=fliplr(y(1,2:end));
X=[yaux y];
X(1,1:length(X)/4)=0;
X(1,3*length(X)/4:end)=0;
length(X);
omega=0:f/length(y):f-(f/length(y));
waux=-fliplr(omega(1,2:end));
w=[waux omega];
length(w);
plot(w,abs(2*X/length(t)));
xlabel('$f$(Hz)','interpreter','latex');
ylabel('Magnitude');
title('$|X_{cc}(j\omega)|$','Interpreter','LaTex','FontSize',14);
axis([-8000 8000 -inf inf]);

% Plot sinal amostrado e discretizado
figure(3);
subplot(2,1,1);
stem(n, sinal_amostrado,'.'); 
xlabel('$n$','Interpreter','LaTex')
ylabel('Amplitude')
title('$x_{cc}(t) = x_{cc}(nT) = x[n]$ amostrado e discretizado', 'Interpreter','LaTex','FontSize',14);
axis([0 50 -inf inf]);


subplot(2,1,2);
y=fft(sinal_amostrado); grid on;
yaux=fliplr(y(1,2:end));
X=[yaux y];
X(1,1:length(X)/4)=0;
X(1,3*length(X)/4:end)=0;
length(X);
omega=0:fs/length(y):fs-(fs/length(y));
waux=-fliplr(omega(1,2:end));
w=[waux omega];
length(w);
plot(w,abs(2*X/length(t)));
xlabel('$f$(Hz)','interpreter','latex');
ylabel('Magnitude');
title('$|X(e^{j\omega})|$','Interpreter','LaTex','FontSize',14);
axis([-8000 8000 -inf inf])

% Plot sinal filtrado e sinal amostrado
figure(4);
plot(t, sinal_filtrado);
hold;
plot(t_amostrado, sinal_amostrado, '.');
xlabel('$t$(s)','Interpreter','LaTex')
ylabel('Amplitude')
title('$x_{cc}(t)$ e $x_{amostrado}(t)$','Interpreter','LaTex','FontSize',14);
axis([0 0.05 -inf inf]);
legend('x_{cc}','x_{amostrado}');

%% Mudança da taxa de amostragem

M = 2;
sinal_subamostrado = downsample(sinal_amostrado,M);
t_subamostrado = downsample(t_amostrado,M);
Td = M*Ts;
fd = 1/Td;
Nd = N/M;
nd = 0 : 1 : Nd-1;
figure(5)
plot(t_amostrado, sinal_amostrado, t_subamostrado, sinal_subamostrado);
xlabel('$t$(s)','Interpreter','LaTex')
ylabel('Amplitude')
title('$y(nT)$, $y(nT_d)$','Interpreter','LaTex','FontSize',14);
axis([0 0.05 -inf inf]);
legend('y(nT)','y(nT_d)');

%% Reconstrução

sinal_reconstruido = sinal_subamostrado*sinc(fd*(ones(length(nd),1)*t-(nd*Td)'*ones(1, length(t))));
% Plot sinal reconstruído
figure(6);
subplot(2,1,1);
plot(t, sinal_reconstruido);
xlabel('$t$(s)','Interpreter','LaTex')
ylabel('Amplitude')
title('$y_{cc}(t)$','Interpreter','LaTex','FontSize',14);
axis([0 0.05 -inf inf]);

subplot(2,1,2);
y=fft(sinal_reconstruido); grid on;
yaux=fliplr(y(1,2:end));
X=[yaux y];
X(1,1:length(X)/4)=0;
X(1,3*length(X)/4:end)=0;
length(X);
omega=0:fs/length(y):fs-(fs/length(y));
waux=-fliplr(omega(1,2:end));
w=[waux omega];
length(w);
plot(w,abs(2*X/length(t)));
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
plot(w,abs(2*X/length(t)));
hold;
y=fft(xt); grid on;
yaux=fliplr(y(1,2:end));
X=[yaux y];
X(1,1:length(X)/4)=0;
X(1,3*length(X)/4:end)=0;
length(X);
omega=0:f/length(y):f-(f/length(y));
waux=-fliplr(omega(1,2:end));
w=[waux omega];
length(w);
plot(w,abs(2*X/length(t)));
xlabel('$f$(Hz)','interpreter','latex');
ylabel('Magnitude');
title('$|Y(j\omega)|$ e $|X(j\omega)|$','Interpreter','LaTex','FontSize',14);
axis([-8000 8000 -inf inf])
legend('y','x');

%sound(sinal_reconstruido, f);

