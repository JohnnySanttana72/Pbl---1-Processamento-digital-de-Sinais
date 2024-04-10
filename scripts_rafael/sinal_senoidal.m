f0 = 1500; % Frequência do sinal
fs = 100000; % Frequência de amostragem
w0 = 2*pi*f0;
t = [0:(1/fs):1];
n = length(t)

x = 3*sin(w0*t);

plot(t,x);
title('Sinal Senoidal');
xlabel('t');
ylabel('x(t)');


y_fft = abs(fft(x)/n); % espectro com a amplitude eixo X (somente o módulo e não a parte complexa)

f = fs*(0:(n /2-1))/n; % eixo Y lei de nixtre para evitar aliasing não vai passar de 100 (50*2)

plot(f, y_fft(1:(n/2)));
title('Transformada de Fourier do Sinal Senoidal');
xlabel('f(Hz)');
ylabel('Amplitude');
% no sinal de audio não é apresentado um único espectro de frequência, mas vários espectros de amplitude da frequencia
