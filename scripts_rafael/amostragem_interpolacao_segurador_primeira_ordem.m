clc;
clear all;
pkg load signal;

n = 1000;
t = 1/n*(0:1:999); % vetor do tempo com passo de 0,001
f_par1 = 50; % Frequência da primeira senoide em Hz
f_par2 = 60; % Frequência da segunda senoide em Hz

sinal = 2*sin(2*pi*f_par1*t) + 3 *sin(2*pi*f_par2*t); % sinal de entrada

% função para criar sinal retangular
function y=square(t)
  t = t*(1/(pi));
  y = ones(size(t));
  y(find(bitand(abs(floor(t)),1))) = -1;
endfunction

fs = 200; % Frequencia de amostragem do sinal retangular

sq = square(2*pi*fs*t); % sinal retângular

% remove parte negativa do sinal retangular
N_square = length(sq);
for k = 1:N_square
  if sq(k) <= 0
    sq(k) = 0;
   else
    sq(k) = 1;
  endif
endfor


f_sinal = (-n/2:n/2-1); % Frequencia do sinal de entrada
X = fftshift(fft(sinal)); % Transformada de Fourier do sinal de entrada

f_sq = (-n/2:n/2-1); % Frequencia do sinal retangular
S = fftshift(fft(sq)); % Transformada de Fourier do sinal retangular

figure;
subplot(3,2,1);
plot(t, sinal);
title('Sinal Senoidal x(t) de entrada');
xlabel('F(Hz)');
ylabel('Magnitude');

subplot(3,2,2);
plot(f_sinal, abs(X),'r'); % O espectro estará centrado em 60 e 50 Hz simetricamente
title('Sinal Senoidal X(jw) no domínio da Frequência');
xlabel('F(Hz)');
ylabel('Magnitude');

subplot(3,2,3);
plot(t, sq); % O espectro etará centrado em zero e nos multiplos inteiros de 200 Hz simetricamente
title('Sinal retangular s(t)');
xlabel('Magnitude');
ylabel('s(t)');

subplot(3,2,4);
plot(f_sq, abs(S), 'r'); % Os espectro estarão centrados na origem e nos multiplos inteiros de 200 Hz simetricamente
title('Sinal retangular S(jw) domínio da Frequência');
xlabel('F(Hz)');
ylabel('Magnitude');

## AMOSTRAGEM DO SINAL

sinal_pam = sinal .* sq; % amostragem do sinal de entrada por meio de pulsos retangulares

f = (-n/2:n/2-1); % Frequencia normalizada
sinal_sampling = fftshift(fft(sinal_pam/n)); % Transformada de Fourier do sinal amostrado

subplot(3,2,5);
plot(t, sinal_pam);  % Os espectro estarão em 60 e 50 Hz e serão duplicados nos multiplos inteiros de 400 Hz  e centrados em 400 Hz (polo positivo e negativo)
title('Sinal Senoidal x(nT)');
xlabel('x(nT)');
ylabel('t');

subplot(3,2,6);
plot(f, abs(sinal_sampling), 'r');
title('Sinal Senoidal Amostrado no domínio da Frequência');
xlabel('F(Hz)');
ylabel('Magnitude');


% Para o filtro passa-baixa a frequência de corte será Ωc = Ωs - ΩM (nesse caso Ωc = 400 - 60)
% Para evitar o Alising (sobreposição das réplicas):  Ωs - ΩM >= ΩM , logo Ωs >= 2ΩM (teorema de Nyquist)
% A frequência de de amostragem tem que ser no mínimo 2 vezes 60 ou maior do que isso
% O efeito de alising só se manisfesta quando o sinal é recuperado quando Ωs < 2ΩM, ou seja, um sinal de frequência maior irá se comportar como um sinal de frequência menor
% Normalmente se amostra com uma taxa bem maior do que o teorema de nyquist e antes da amostragem um filtro anti-alising é aplicado para retirar componentes de alta frequência


## RECONSTRUÇÃO DO SINAL
## Escolher Ωc (frequência de corte) e dar um ganho de T, pois a amplitude será reduzida para 1/T
## O Ωc terá que ser maior que ΩM e menor que Ωs, ou seja, Ωc = Ωs/2

f = (-n/2:n/2-1); % vetor de frequência
t2 = (-0.5:1/n:0.499); % vetor do tempo do pulso triangular

T_max = 1/(fs/2); % largura máxima do pulso retangular
retangulo = rectpuls(t2, T_max); % construção do pulso retangular com largura de 0,007

% filtro de primeira ordem
pulso_triangular = conv(retangulo, retangulo, "same"); % convolução de dois pulsos retângulares para obter o pulso triangular

sinc_triangulo = fftshift(fft(pulso_triangular)); % Transformada de Fourier do pulso triangular

sinal_interpolado = fftshift(sinc_triangulo).* sinal_sampling; % realizar interpolação para reconstrução no domínio da frequência

sinal_reconstruido = ifft(sinal_interpolado); % Transformada Inversa de Fourier aplicada ao sinal reconstruído

figure;
subplot(2,2,1);
plot(t2, pulso_triangular/(T_max*length(t2)));  % Os espectro estarão em 60 e 50 Hz e serão duplicados nos multiplos inteiros de 400 Hz  e centrados em 400 Hz (polo positivo e negativo)
title('Pulso Triangular');
ylabel('t');

subplot(2,2,2);
plot(f, abs(sinc_triangulo), 'r', [fs/2,fs/2], [0,25], 'b', [-fs/2,-fs/2], [0,25], 'b');  % sinal sinc no domínio da frequência com frequência de corte igual a 100 Hz
title('Pulso Triangular no domínio da Frequência');
xlabel('F(Hz)');
ylabel('Magnitude');

subplot(2,2,3);
plot(t, sinal_reconstruido); % O sinal x(t)reconstruído não ficou exatamente igual
title('Sinal x(t) reconstruído');
xlabel('x(t)');
ylabel('t');

subplot(2,2,4);
plot(f, abs(fftshift(sinal_interpolado)), 'r');  % Os espectro ficou centrado em 50 e 40 Hz simetricamente e com ruídos
title('Sinal x(t) reconstruído no domínio da Frequência');
xlabel('F(Hz)');
ylabel('Magnitude');
