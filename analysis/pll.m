clc;
clear;

%% stimulus params
FS = 12e3; % sampling rate
N = 65536;
F = rand()*200-100; % freq offset
PHI = 2*pi*rand(); % phase
SNR =  -30;

%% pll params
LB = 128; % loop bandwidth in hz
DF = sqrt(2)/2; % damping factor
G = 1024+128; % gain

init_beta = 2*DF*2*pi*LB/FS - (2*pi*LB)^2 / (FS^2*G);
beta = init_beta; % proportional
alpha = (2*pi*LB)^2 / (FS*G); % integral

%% sim
t = (0:N-1).'/FS;
x = exp(1j*2*pi*F.*t+PHI) + wgn(N, 1, -SNR);% awgn(exp(1j*2*pi*F.*t+PHI), SNR, 'measured');
y = ones(N, 1);
f = zeros(N, 1);
phi = zeros(N, 1);
lock = false(N, 1);
phierror = zeros(N, 1);

mav = 1024*ones(32, 1);

for i = 1:N
    delta = atan2(imag(x(i)/y(i)), real(x(i)/y(i))); % phase detector
    phierror(i) = delta;
    if i<N
        f(i+1) = f(i) + beta*delta;
        phi(i+1) = phi(i) + alpha*delta + 2*pi*f(i+1)/FS;
        while phi(i+1)>pi
            phi(i+1) = phi(i+1)-2*pi;
        end
        while phi(i+1)<-pi
            phi(i+1) = phi(i+1)+2*pi;
        end
    end
    
    if mod(i-1, 128)==0
        mav = [delta; mav(1:31)];
        if mav(32)~=1024
            phistdd = std(mav);
            if phistdd<pi/2
                beta = beta*0.995;
                f(i+1) = f(i) + beta*delta;
            else
                beta = beta*1.005;
                if beta>init_beta; beta=init_beta; end
                f(i+1) = f(i) + beta*delta;
            end
            if phistdd<0.5
                lock(i) = true;
            end
        end
    else
        lock(i) = lock(i-1);
    end
    
    y(i+1) = exp(1j*phi(i));
end

%%
clf;
subplot(311);
hold on;
plot([0, t(N)], [F, F]);
plot(t, f, 'r');
plot(t, lock*200-100, 'g');
hold off;
ylim([-120,120]);
grid on;

subplot(312);
plot(t, f-F);
grid on;
title(sprintf('\\Delta f=%.3f',abs(f(end)-F)));

subplot(313);
plot(t, phierror, '.');
grid on;
