%% Matching filtering
clc;
rx = load('log');
% spectrogram(rx(:,1),4096);

t=(1:8:65536)'/48000;
ch = cos(2*pi*(-600*t+(1200. * 48000 / 65536 / 2)*t.*t));

t=(1:length(rx))'/48000;
sig = rx(1:8:length(rx), 2)+1j*rx(1:8:length(rx), 3);
x = abs(filter(ch, 1, sig));

peak = zeros(65,1);
psum = zeros(fix(length(rx)/8),1);
pulse = zeros(fix(length(rx)/8),1);
for i=2:(length(rx)/8-1)
    if ((x(i)>x(i-1) && x(i)>=x(i+1)) || (x(i)>=x(i-1) && x(i)>x(i+1)))
        psum(i) = psum(i-1)-peak(1)+peak(65);
        peak = [peak(2:65); x(i)];
        if peak(1)>0 && psum(i)*.21<x(i)
            pulse(i)=1;
        end
    else
        psum(i) = psum(i-1);
    end
end

clf;
hold on;
plot(x, '.-');
plot(pulse*20, 'r.-');
plot(psum(5:end)*.21, 'g.-');
hold off;
grid on;

%% Carrier Sync
% clc;
% rx = load('log');
% plot(rx(:,2));