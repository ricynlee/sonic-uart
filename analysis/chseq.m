clc;
clear;

snr = 0;

fs = 48e3;
fc = 18e3;
n = 2048;
t = (0:n-1)/fs;
u0 = 1200/((n-1)/fs); % bandwidth = 1200Hz

bbt = cos(2*pi*((-600)*t+u0*t.^2/2)).';
bbr = exp(1j*2*pi*((-600+120*rand()-60)*t+u0*t.^2/2)+1j*pi*(2*rand()-1)).' + ... 
      10^(-snr/20)*randn(n, 1)./sqrt(2) + 1j*10^(-snr/20)*randn(n, 1)./sqrt(2);

[m, r] = filter(flipud(bbt), 1, bbr, zeros(n-1, 1));

subplot(1, 2, 1);
% plot(t, real(sig));
spectrum = abs(fft(bbr,n));
plot((-n/2:n/2-1)*fs/n, [spectrum(n/2+1:n); spectrum(1:n/2)]);
grid on;
title('Spectrum of noisy, freq shifted and phase shifted RX chirp');

subplot(1, 2, 2);
plot(abs([m; r]), '.-');
title('Match filtered');
grid on;
