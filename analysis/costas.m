% point-by-point loop
clc;
clear;

fs = 48e3; % sampling freq

snr = 0;

n = 27000;

alpha = 0.03; % phi ctrl
beta = 0.12; % frq ctrl

lpf = [0.00323248474800766,0.00341121901342408,0.00354043287435739,0.00361486147941774,0.00362990165094776,0.00358172662821808,0.00346739301787512,0.00328493784523272,0.00303346369827736,0.00271321009190362,0.00232560935148027,0.00187332552027135,0.00136027503176679,0.000791628152273114,0.000173790487284370,-0.000485635847193166,-0.00117791146487034,-0.00189323910564847,-0.00262086259714419,-0.00334918445573330,-0.00406590087447166,-0.00475815253104020,-0.00541268935309986,-0.00601604710646551,-0.00655473342833562,-0.00701542071806652,-0.00738514312577998,-0.00765149474800123,-0.00780282605245685,-0.00782843551334176,-0.00771875344526418,-0.00746551507939154,-0.00706192002893901,-0.00650277544214257,-0.00578462033751411,-0.00490582885597019,-0.00386669044409368,-0.00266946529835247,-0.00131841374694390,0.000180201381138810,0.00181813995689834,0.00358523747143389,0.00546948540925826,0.00745713839187177,0.00953284699712183,0.0116798147603736,0.0138799774825452,0.0161142026138886,0.0183625061569499,0.0206042842428824,0.0228185562871214,0.0249842164276578,0.0270802897953071,0.0290861900632833,0.0309819746749972,0.0327485941554144,0.0343681319727539,0.0358240315330925,0.0371013070590076,0.0381867353223180,0.0390690254670064,0.0397389644675092,0.0401895361149986,0.0404160118046930,0.0404160118046930,0.0401895361149986,0.0397389644675092,0.0390690254670064,0.0381867353223180,0.0371013070590076,0.0358240315330925,0.0343681319727539,0.0327485941554144,0.0309819746749972,0.0290861900632833,0.0270802897953071,0.0249842164276578,0.0228185562871214,0.0206042842428824,0.0183625061569499,0.0161142026138886,0.0138799774825452,0.0116798147603736,0.00953284699712183,0.00745713839187177,0.00546948540925826,0.00358523747143389,0.00181813995689834,0.000180201381138810,-0.00131841374694390,-0.00266946529835247,-0.00386669044409368,-0.00490582885597019,-0.00578462033751411,-0.00650277544214257,-0.00706192002893901,-0.00746551507939154,-0.00771875344526418,-0.00782843551334176,-0.00780282605245685,-0.00765149474800123,-0.00738514312577998,-0.00701542071806652,-0.00655473342833562,-0.00601604710646551,-0.00541268935309986,-0.00475815253104020,-0.00406590087447166,-0.00334918445573330,-0.00262086259714419,-0.00189323910564847,-0.00117791146487034,-0.000485635847193166,0.000173790487284370,0.000791628152273114,0.00136027503176679,0.00187332552027135,0.00232560935148027,0.00271321009190362,0.00303346369827736,0.00328493784523272,0.00346739301787512,0.00358172662821808,0.00362990165094776,0.00361486147941774,0.00354043287435739,0.00341121901342408,0.00323248474800766];
ord = length(lpf) - 1;

t = (0:n-1)/fs;

frqt = 18000 + 120*rand()-60;
[pol, tmp] = filter(lpf, 1, mod(fix((0:n-1)/1024), 2)*2-1, zeros(1, ord));
pol = [pol(round(ord/2):n), tmp(1:round(ord/2)-1)'];
sig = pol.*cos(2*pi*frqt.*t + rand()*2*pi) + 10^(-snr/20)*randn(1, n)./sqrt(2);

phir = zeros(1, n);
frqr = [18000, zeros(1, n-1)];
loi = zeros(1, n);
loq = zeros(1, n);
bbi = zeros(1, n);
bbq = zeros(1, n);

mixi = zeros(1, ord+1);
mixq = zeros(1, ord+1);

frqr_sample = zeros(1, 12);
lock = zeros(1, n);

phie = zeros(1, n);
for i=1:n
    loi(i) = cos(phir(i));
    loq(i) = -sin(phir(i));

    mixi = [loi(i)*sig(i), mixi(1:ord)];
    mixq = [loq(i)*sig(i), mixq(1:ord)];

    bbi(i) = sum(lpf.*mixi);
    bbq(i) = sum(lpf.*mixq);
    
    phie(i) = bbi(i)*bbq(i)/sqrt(bbi(i)^2+bbq(i)^2);
    
    if mod(i, 300)==0
        frqr_sample = [frqr(i), frqr_sample(1:length(frqr_sample)-1)];
        lock(i) = (max(frqr_sample) - min(frqr_sample) < 4);
    elseif i>1
        lock(i) = lock(i-1);
    end
    
    if i~=n
        frqr(i+1) = frqr(i) + beta*phie(i);
        phir(i+1) = phir(i) + 2*pi*frqr(i)/fs + alpha*phie(i);
    end
end

clf;
set(gcf, 'renderer', 'painters');

subplot(1,3,1);
title('Freq following');
hold on;
plot(t, frqr, 'g', 'linesmooth', 'on');
plot(t, frqt+lock*20, 'm-');
plot([t(1), t(n)], [frqt, frqt]);
hold off;
ylim([18000-100, 18000+100]);
xlim([t(1), t(n)]);
grid on;
set(gca, 'yticklabel', get(gca, 'ytick')); % cannnot get rid of yaxis exponent if opengl is used

subplot(1,3,2);
title('Phase error');
hold on;
plot(t, phie, 'g', 'linesmooth', 'on');
plot([t(1), t(n)], [0, 0], 'linesmooth', 'on');
hold off;
ylim([-.5, .5]);
xlim([t(1), t(n)]);
grid on;

subplot(1,3,3);
title('Waveform');
hold on;
plot(t(n-16:n), abs(loi(n-16:n)), 'g', 'linesmooth', 'on');
plot(t(n-16:n), abs(sig(n-16:n)), 'b.-', 'linesmooth', 'on');
hold off;
ylim([-.2, 1.6]);
xlim([t(n-16), t(n)]);
grid on;
