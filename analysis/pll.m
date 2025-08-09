clc;

%% stimulus params
FS = 12e3; % sampling rate
N = (65536-3072)/4;
F = rand()*4-2; % freq offset
PHI = 2*pi*rand(); % phase
SNR = -15;

%% pll params
LB = 2; % loop bandwidth in hz
DF = sqrt(2)/2; % damping factor
G = 6; % gain

CC_UB = 96;
CC_LB = 32;

beta = 2*DF*2*pi*LB/FS - (2*pi*LB)^2 / (FS^2*G); % proportional
alpha = (2*pi*LB)^2/(FS*G); % integral

%% sim
t = (0:N-1).'/FS;
if exist('rx', 'var')
    x = rx(1:4:end,4) + 1j*rx(1:4:end,5);
else
    x = awgn(exp(1j*2*pi*F.*t+PHI), SNR, 'measured');
end
y = ones(N, 1);
f = zeros(N, 1);
phi = zeros(N, 1);
lock = false(N, 1);
stable = false(N, 1);
delta = zeros(N, 1);

mav0 = 1024*ones(32, 1);
mav1 = 1024*ones(32, 1);

cc = CC_UB; % confidence coef, smaller is better (indicating higher SNR)
sc = 0; % stable count, larger is better (indicating higher SNR)
final = 0; % final freq offset in hz

for i = 1:N
    delta(i) = atan2(imag(x(i)/y(i)), real(x(i)/y(i))); % phase detector

    if mod(i-1, 64)==0
        mav0 = [delta(i); mav0(1:31)];
        if mav0(32)~=1024
            phistdd = std(mav0);
            if phistdd<pi/9
                lock(i) = true;
                cc = cc - 1;
                if cc<CC_LB; cc=CC_LB; end
            else
                cc = cc + 1;
                if cc>CC_UB; cc=CC_UB; end
            end
        end
    elseif i>1
        lock(i) = lock(i-1);
    end

    if mod(i-1, 64)==31
        mav1 = [f(i); mav1(1:31)];
        if mav1(32)~=1024
            fstdd = std(mav1);
            if fstdd<0.025
                stable(i) = true;
                final = mean(mav1);
                sc = sc + 1;
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
plot(t, f, 'b');
if exist('rx', 'var')
    plot(t, rx(1:4:end,7), 'r.', 'markersize', 1);
end
% plot(t, lock*4-2, 'g', 'displayname', 'Phase lock');
% plot(t, stable*4-2, 'c', 'displayname', 'Frq stable');
hold off;
axis([min(t) max(t) -2.2, 2.2]);
grid on;
if F~=0
    if abs(final-F)<=abs(F)
        improvement = 1-abs(final-F)/abs(F);
        improvement = sprintf('%.1f%% better', improvement*100);
    else
        improvement = abs(final-F)/abs(F);
        improvement = sprintf('%.2fx worse', improvement);
    end
else
    improvement = abs(final);
    improvement = sprintf('%.3fhz deviated', improvement);
end
fast_exp = @(x) typecast(int32(single(x)*single(12102203)) + int32(1065353216), 'single');
title(sprintf('final=%.3fHz confidence=%.3f\n\\Delta=%.3fHz %s',final, 1-fast_exp(-sc/cc), abs(final-F), improvement));
labelled = findobj(gca, 'Type', 'line', '-not', 'DisplayName', '');
if ~isempty(labelled); legend(labelled, 'location', 'best'); end

subplot(312);
plot(t, delta, '.', 'markersize', 1);
if exist('rx', 'var')
    hold on;
    plot(t, rx(1:4:end,6), 'r.', 'markersize', 1);
    hold off;
end
grid on;
axis([min(t) max(t) -pi pi]);
title('Phase error detected');

subplot(313);
plot(t, real(x), '.', 'markersize', 1);
hold on;
plot(t, imag(x), 'r.', 'markersize', 1);
hold off;
grid on;
xlim([min(t) max(t)]);
title('Input signal');
