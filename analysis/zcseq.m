R = 1; % this is not any integer. 1 is legal
N = 64;
t = transpose(0:N-1);

zc = exp(-1j*pi*R*t.*(t+1)/(N-1));

plot(real(zc), imag(zc), 'o', 'LineSmooth', 'on');
axis([-1.2 1.2 -1.2 1.2]);
grid on;
axis square;

title('Press any key to continue');
pause;

zc_phase_shifted = zc.*exp(-1j*2*pi*rand());

ac = zeros(1, N);
for i=1:N
    sh = i-1+floor(-N/2);
    ac(i) = abs(circshift(zc, sh)'*zc);
end

plot(ac, 'LineSmooth', 'on');
grid on;
