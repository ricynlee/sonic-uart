% point-by-point loop
% 48ksps -> 12ksps filtered down sample (extra latency introduced)
clc;
clear;

fs = 48e3; % sampling freq

snr = -15;

n = 32768; 

init_beta = 0.0052;
init_bblim = 0.116;
adjcoef_beta = 0.958;
adjcoef_bblim = 0.998;

alpha = 0.0026; % phi ctrl
beta = init_beta; % frq ctrl
bblim = init_bblim; % phie fluctuation lim

% kaiser len=16, fs=24000, fpass=750, fstop=1250, apass=0.1, astop=80
hbf = [-1.00565929571123e-05,-0.000321623413532287,0.00199470588834347,0.00750258131670068,-0.0214172449746208,-0.0521987304211924,0.123744402545264,0.440705965651995,0.440705965651995,0.123744402545264,-0.0521987304211924,-0.0214172449746208,0.00750258131670068,0.00199470588834347,-0.000321623413532287,-1.00565929571123e-05];
% kaiser len=128, fs=12000, fpass=750, fstop=1250, apass=0.1, astop=83.5
lpf12k = [9.02066018368156e-06,1.47146669044470e-05,1.61263137389339e-05,8.36971724868580e-06,-1.14243566349565e-05,-4.14028243603962e-05,-7.33568493648115e-05,-9.34228974243534e-05,-8.57779880072831e-05,-3.88699872491052e-05,4.75818004373267e-05,0.000157561627528463,0.000258585874502042,0.000308226459248327,0.000267039417592598,0.000114928782541558,-0.000134334583601472,-0.000426637736610233,-0.000674070310034469,-0.000775990612521527,-0.000651111180603454,-0.000272055513271638,0.000309383871499698,0.000957812776710881,0.00147769250649465,0.00166367467102975,0.00136715170910193,0.000560190518785614,-0.000625486697300812,-0.00190341734857556,-0.00288956489440250,-0.00320442537962824,-0.00259627872245236,-0.00104985191066245,0.00115787661807726,0.00348351469031493,0.00523287158194523,0.00574734583134558,0.00461602924949001,0.00185201785725053,-0.00202858461573343,-0.00606729811070622,-0.00907033021918942,-0.00992536811557056,-0.00795199199938506,-0.00318686416700240,0.00349195854790376,0.0104653748736718,0.0157069302183590,0.0172931437219755,0.0139756395029551,0.00566684550283199,-0.00630534642857666,-0.0192748016756516,-0.0296710625174893,-0.0337440854734465,-0.0284317144017534,-0.0121708268783935,0.0145507107854594,0.0490564027534398,0.0868227419137602,0.122253252530521,0.149732253247598,0.164732783595156,0.164732783595156,0.149732253247598,0.122253252530521,0.0868227419137602,0.0490564027534398,0.0145507107854594,-0.0121708268783935,-0.0284317144017534,-0.0337440854734465,-0.0296710625174893,-0.0192748016756516,-0.00630534642857666,0.00566684550283199,0.0139756395029551,0.0172931437219755,0.0157069302183590,0.0104653748736718,0.00349195854790376,-0.00318686416700240,-0.00795199199938506,-0.00992536811557056,-0.00907033021918942,-0.00606729811070622,-0.00202858461573343,0.00185201785725053,0.00461602924949001,0.00574734583134558,0.00523287158194523,0.00348351469031493,0.00115787661807726,-0.00104985191066245,-0.00259627872245236,-0.00320442537962824,-0.00288956489440250,-0.00190341734857556,-0.000625486697300812,0.000560190518785614,0.00136715170910193,0.00166367467102975,0.00147769250649465,0.000957812776710881,0.000309383871499698,-0.000272055513271638,-0.000651111180603454,-0.000775990612521527,-0.000674070310034469,-0.000426637736610233,-0.000134334583601472,0.000114928782541558,0.000267039417592598,0.000308226459248327,0.000258585874502042,0.000157561627528463,4.75818004373267e-05,-3.88699872491052e-05,-8.57779880072831e-05,-9.34228974243534e-05,-7.33568493648115e-05,-4.14028243603962e-05,-1.14243566349565e-05,8.36971724868580e-06,1.61263137389339e-05,1.47146669044470e-05,9.02066018368156e-06];
% cheb ii, direct ii, cascade, 2 sections, fs=12k, fpass=12, fstop=110, apass=0.1, astop=42
nbf12k = [1,-1.99557936191559,1,1,-1.98248481750488,0.982806622982025;1,1,0,1,-0.981456756591797,0];
g = [0.0727965161204338;0.00927162170410156;1];
for j = 1:length(g)-1
    nbf12k(j, 1:3) = nbf12k(j, 1:3).*g(j);
end

t = (0:n-1)/fs;

frqt = 18000 + 6*rand()-3;
amp = 10^(2*rand()-1);
sig = cos(2*pi*frqt.*t + rand()*2*pi) + 10^(-snr/20)*randn(1, n)./sqrt(2);
sig = amp*sig;

phir = zeros(1, n);
frqr = [18000, zeros(1, n-1)];
loi = zeros(1, n);
loq = zeros(1, n);

z0 = zeros(2, length(hbf));
z1 = zeros(2, length(hbf));
z2 = zeros(2, length(lpf12k));
z3 = zeros(2, size(nbf12k, 1), 6);

bbi = zeros(1, n);
bbq = zeros(1, n);

bb_sample = ones(2, 32);
lock = zeros(1, n);

phie = zeros(1, n);
for i=1:n
    % 48k
    loi(i) = cos(phir(i));
    loq(i) = -sin(phir(i));
    
    if i<n
        phir(i+1) = phir(i) + 2*pi*frqr(i)/fs;
        frqr(i+1) = frqr(i);
    end
    if i>1
        lock(i) = lock(i-1);
    end

    z0 = [[loi(i)*sig(i); loq(i)*sig(i)], z0(:, 1:length(hbf)-1)];
    filtered = sum([hbf;hbf].*z0, 2);

    if mod(i-1, 2)==0
        % 24k
        z1 = [filtered, z1(:, 1:length(hbf)-1)];
        filtered = sum([hbf;hbf].*z1, 2);

        if mod(i-1, 4)==0
            % 12k

            % z2 = [filtered, z2(:, 1:length(lpf12k)-1)];
            % filtered = sum([lpf12k;lpf12k].*z2, 2);
            
            for iq=1:2
                for sec=1:size(nbf12k, 1)
                    z3(iq, sec, 1) = filtered(iq) - z3(iq, sec, 2)*nbf12k(sec, 5) - z3(iq, sec, 3)*nbf12k(sec, 6);
                    filtered(iq) = z3(iq, sec, 1)*nbf12k(sec, 1) + z3(iq, sec, 2)*nbf12k(sec, 2) + z3(iq, sec, 3)*nbf12k(sec, 3);
                    z3(iq, sec, 3) = z3(iq, sec, 2);
                    z3(iq, sec, 2) = z3(iq, sec, 1);
                end
            end

            bbi(i) = filtered(1);
            bbq(i) = filtered(2);

            phie(i) = atan(bbq(i)/bbi(i));
            if i<n
                frqr(i+1) = frqr(i+1) + beta*phie(i);
                phir(i+1) = phir(i+1) + alpha*phie(i);
            end

            if mod(i-1, 360)==0
                bb_sample = [[bbi(i);bbq(i)], bb_sample(:, 1:length(bb_sample)-1)];
                lock(i) = ( ...
                    abs( ...
                        mean(bb_sample(2,:)) / mean(bb_sample(1,:)) ...
                    ) < bblim ...
                );
                if lock(i)
                    beta = beta * adjcoef_beta;
                    bblim = bblim * adjcoef_bblim;
                else
                    beta = beta / adjcoef_beta;
                    bblim = bblim / adjcoef_bblim;
                    if beta > init_beta, beta = init_beta; end
                    if bblim > init_bblim, bblim = init_bblim; end
                end
            end
        end
    end
end

clf;
set(gcf, 'renderer', 'painters');

subplot(1,3,1);
title(sprintf('Freq following\n\\Delta(1)=%.2f \\Delta(n)=%.2f', frqr(1)-frqt, frqr(n)-frqt), 'interpreter', 'tex');
hold on;
plot(t, frqr, 'g-', 'linesmooth', 'on');
plot(t, frqt+lock-.5, 'c');
plot([t(1), t(n)], [frqt, frqt]);
hold off;
ylim([18000-10, 18000+10]);
xlim([t(1), t(n)]);
grid on;
set(gca, 'yticklabel', get(gca, 'ytick')); % cannnot get rid of yaxis exponent if opengl is used

subplot(1,3,2);
title('Measured phase error');
hold on;
plot(t, phie, 'g', 'linesmooth', 'on');
plot([t(1), t(n)], [0, 0], 'linesmooth', 'on');
hold off;
ylim([-pi/2, pi/2]);
xlim([t(1), t(n)]);
grid on;

subplot(1,3,3);
hold on;
plot(t, bbi, 'g-');
plot(t, bbq, 'r-');
hold off;

% title(sprintf('Normalized waveform\nAmplitude=%.2f', amp));
% hold on;
% plot(t(n-16:n), abs(loi(n-16:n)), 'g', 'linesmooth', 'on');
% plot(t(n-16:n), abs(sig(n-16:n)/amp), 'b.-', 'linesmooth', 'on');
% hold off;
% ylim([-.2, 1.6]);
% xlim([t(n-16), t(n)]);
% grid on;
