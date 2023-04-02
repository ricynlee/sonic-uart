clc;

ORDER = 6; % 6 or 8

lfsr=[1, zeros(1,ORDER-1)];

state = zeros(2^ORDER-1, 1);
m = zeros(2^ORDER-1, 1);

for i=1:(2^ORDER-1)
    state(i)=sum(lfsr.*(2.^(0:ORDER-1)));
    m(i)=lfsr(1);
    if ORDER==6
        % f=xor(lfsr(1), lfsr(2)); % 000011
        % f=xor(xor(lfsr(1), lfsr(2)), xor(lfsr(3), lfsr(6))); % 100111
        f=xor(xor(lfsr(1), lfsr(3)), xor(lfsr(4), lfsr(6))); % 101101        
    elseif ORDER==8
        f=xor(xor(lfsr(1), lfsr(2)), xor(lfsr(6), lfsr(7)));
    end
    lfsr=[lfsr(2:ORDER), f];
end

if length(unique(state))==2^ORDER-1
    fprintf('%d,', m);
    fprintf('\n');
else
    fprintf('NOT MSEQ!\n');
end

m = 2*m-1;

ac = zeros(1, 2^ORDER-1);
for i=1:2^ORDER-1
    sh = i-1+floor(-(2^ORDER-1)/2);
    ac(i) = abs(circshift(m, sh)'*m);
end

plot(ac, 'LineSmooth', 'on');
grid on;
