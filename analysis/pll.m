clc;
clear;

%% stimulus params
FS = 12e3; % sampling rate
N = 65536-4096;
F = rand()*6-3; % freq offset
PHI = 2*pi*rand(); % phase
SNR = 10;

%% pll params
LB = 1; % loop bandwidth in hz
DF = sqrt(2)/2; % damping factor
G = 3; % gain

beta = 2*DF*2*pi*LB/FS - (2*pi*LB)^2 / (FS^2*G); % proportional
alpha = (2*pi*LB)^2/(FS*G); % integral

%% sim
t = (0:N-1).'/FS;
x = awgn(exp(1j*2*pi*F.*t+PHI), SNR, 'measured');
y = ones(N, 1);
f = zeros(N, 1);
phi = zeros(N, 1);
lock = false(N, 1);
stable = false(N, 1);
delta = zeros(N, 1);

mav0 = 1024*ones(32, 1);
mav1 = 1024*ones(32, 1);

for i = 1:N
    delta(i) = atan2(imag(x(i)/y(i)), real(x(i)/y(i))); % phase detector
    
    if mod(i-1, 256)==0
        mav0 = [delta(i); mav0(1:31)];
        if mav0(32)~=1024
            phistdd = std(mav0);
            if phistdd<pi/9
                lock(i) = true;
            end
        end
    elseif i>1
        lock(i) = lock(i-1);
    end

    if mod(i-1, 256)==127
        mav1 = [f(i); mav1(1:31)];
        if mav1(32)~=1024
            fstdd = std(mav1);
            if fstdd<0.04
                stable(i) = true;
                final = mean(mav1);
            end
        end
    elseif i>1
        stable(i) = stable(i-1);
    end

    if i<N
        f(i+1) = f(i) + beta*delta(i);
        phi(i+1) = phi(i) + alpha*delta(i) + 2*pi*f(i+1)/FS;
        while phi(i+1)>pi
            phi(i+1) = phi(i+1)-2*pi;
        end
        while phi(i+1)<-pi
            phi(i+1) = phi(i+1)+2*pi;
        end
        y(i+1) = exp(1j*phi(i+1));
    end
end

%%
clf;
subplot(311);
plot([0, t(N)], [F, F], ':', 'linewidth', 2);
hold on;
plot(t, f, 'r');
plot(t, lock*8-4, 'g', 'displayname', 'Phase lock');
plot(t, stable*8-4, 'c', 'displayname', 'Frq stable');
hold off;
ylim([-5,5]);
grid on;
title(sprintf('final=%.3fHz \\Delta=%.3fHz',final, abs(final-F)));
labelled = findobj(gca, 'Type', 'line', '-not', 'DisplayName', '');
legend(labelled, 'location', 'best');

subplot(312);
plot(t, delta, '.', 'markersize', 1);
grid on;
title('Phase error');

subplot(313);
plot(t, real(x), '.', 'markersize', 1);
hold on;
plot(t, imag(x), 'r.', 'markersize', 1);
hold off;
grid on;
title('Input signal');
legend('I', 'Q', 'location', 'southeast');
