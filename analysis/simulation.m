%%
clear;
clc;
clf;

vldt = @(v) assert(abs(v) < 2^15);

%% Parameters
fs = 48e3; % sampling rate
fc = 18e3; % central freq

t = 0:1/fs:(128*33-1)/fs; % time axis

%% Based-band signal (DSSS)
MSEQ = [1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0];

BIT = 0;

BB = [zeros(1, 128), reshape(repmat(xor(MSEQ, BIT), 128, 1), 1, 128*length(MSEQ)), zeros(1, 128)];

subplot(4,2,1);
plot(t,BB,'LineSmoothing','On');
grid on;
axis([min(t), max(t), -0.5, 1.5]);
title(sprintf('1.Base band signal (bit): %d', BIT));

clear BIT;

%% Generate S16.10 Fixed Pointed Signal
A = 0.5;
X = A*cos(2*pi*fc*t).*BB;

subplot(4,2,3);
plot(t,X,'LineSmoothing','On');
grid on;
axis([min(t), max(t), -1, 1]);
title('2.Audio Signal');

clear BB fc;

%% Add Noise
SNR=0;
X = awgn(X, SNR-db(0.5*A^2));
X = round(X*2048);

subplot(4,2,5);
plot(t,X,'LineSmoothing','On');
grid on;
axis([min(t), max(t), -5000, 5000]);
title(sprintf('3.Audio Signal w/ Noise, SNR=%ddB', SNR));

clear SNR;

%% S12.9 Band-pass IIR Filter Coefficients
B=[0.00176620483398438,0,-0.00529837608337402,0,0.00529837608337402,0,-0.00176620483398438;];
A=[1,3.90625000000000,7.58105468750000,8.73242187500000,6.36035156250000,2.74853515625000,0.590820312500000;];

ORDER=length(B)-1;

KB=2^22;
KA=2^11;

B=round(B*KB);
A=round(A*KA);

vldt(max(abs(B)));
vldt(max(abs(A)));

N=1024;
subplot(4,2,7);
plot(0:fs/N/2:(fs/2-fs/N/2), db(freqz(B/KB,A/KA,N)),'LineSmoothing','On');
hold on;
grid on;
axis([0 fs/2 -120 10]);
title('4.IIR Filter');

clear N fs;

%% Filter
Y=zeros(1,length(X));
Z=zeros(1,ORDER);

for i=1:length(X)
    for n=1:ORDER+1
        if n==1
            Y(i)=Z(n)+floor(X(i)*B(n)/KB);
            vldt( X(i)*B(n)/KB );
            vldt( Y(i) );
        elseif n==ORDER+1
            Z(n-1)=floor(X(i)*B(n)/KB)-floor(Y(i)*A(n)/KA);
            vldt( Y(i)*A(n)/KA );
            vldt( X(i)*B(n)/KB );
            vldt( Z(n-1) );
        else
            Z(n-1)=Z(n)+floor(X(i)*B(n)/KB)-floor(Y(i)*A(n)/KA);
            vldt( Y(i)*A(n)/KA );
            vldt( X(i)*B(n)/KB );
            vldt( Z(n-1) );
        end
    end
end

subplot(4,2,2);
plot(t,Y,'LineSmoothing','On', 'Color',[0.6,0.6,1]);
grid on;
axis([min(t), max(t), -3000, 3000]);

clear KA KB ORDER Z n vldt;

%% Rectification
for i=1:length(X)
    Y(i) = abs(Y(i));
end
clear Yq;

% subplot(4,2,2);
hold on;
plot(t,Y,'LineSmoothing','On');
grid on;
% axis([min(t), max(t), -1000, 2500]);
title('5.Filtered & Rectified');

%% Envelope detection
linear_envelope_detection_thresh = 128;
for i=1:length(X)
    if i==1
        prev=0;
    else
        prev=Y(i-1);
    end

    Y(i) = floor(mean([prev,Y(i)]));

    if Y(i)>prev+linear_envelope_detection_thresh
        Y(i) = prev+linear_envelope_detection_thresh;
    elseif Y(i)<prev-linear_envelope_detection_thresh
        Y(i) = prev-linear_envelope_detection_thresh;
    end
end

subplot(4,2,4);
plot(t,Y,'LineSmoothing','On','Color',[0.6,0.6,1]);
grid on;
axis([min(t), max(t), -1000, 1500]);
title('7.Envelope Detection');

clear X prev;

%% Accumulation
ACC=16;
for i=1:floor(length(Y)/ACC)
    Y(i) = sum(Y((i*ACC-ACC+1):(i*ACC)))/ACC;
    t(i) = t(i*ACC-ACC+1);
end

Y = Y(1:floor(length(Y)/ACC));
t = t(1:floor(length(t)/ACC));

subplot(4,2,4);
hold on;
plot(t, Y, 'LineSmoothing','On','Color', 'b');
hold off;
title('6.Envelope Detection & Accumulation');

clear ACC;

%% Despread
MSEQ = reshape(repmat(MSEQ, 8, 1), 1, 8*31)*2-1;
Y_match = zeros(1, length(Y));
for i=1:length(Y)
    for j=1:length(MSEQ)
        if i+j-124>0 && i+j-124<=length(Y)
            Y_match(i) = Y_match(i)+MSEQ(j)*Y(i+j-124);
        end
    end
end

Y = Y_match;

subplot(4,2,6);
plot(t, Y, '-', 'LineSmooth', 'on');
axis([min(t), max(t), min(Y)-20000, max(Y)+20000]);
grid on;
title('7.Despread');

clear Y_match MSEQ linear_envelope_detection_thresh;

%% Decisioned
Y = Y./max(abs(Y));

decision_thresh = 0.8;

subplot(4,2,8);
plot(t,Y,'LineSmoothing','On');
hold on;
plot([min(t), max(t)], [decision_thresh, decision_thresh],'--','LineSmoothing','On','Color',[1,0.4,0.4]);
plot([min(t), max(t)], [-decision_thresh, -decision_thresh],'--','LineSmoothing','On','Color',[1,0.4,0.4]);
hold off;
grid on;
axis([min(t), max(t), -1.2, 1.2]);
title('8.Decisioned');

clear Y decision_thresh i t;
