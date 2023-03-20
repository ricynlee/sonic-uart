clc;

ORDER = 6; % 6 or 8

a=[1, zeros(1,ORDER-1)];

state = zeros(1, 2^ORDER);
m = zeros(1, 2^ORDER);

for i=1:(2^ORDER-1)
    state(i+1)=sum(a.*(2.^(0:ORDER-1)));
    m(i+1)=a(1);
    if ORDER==6
        f=xor(a(1), a(2));
    elseif ORDER==8
        f=xor(xor(a(1), a(2)), xor(a(6), a(7)));
    end
    a=[a(2:ORDER), f];
end

if length(unique(state))==2^ORDER
    fprintf('m seq:\n');
    fprintf('%d,', m);
    fprintf('\n');
else
    fprintf('NOT MSEQ!\n');
end

m = reshape(repmat(m, 16, 1), 1, length(m)*16);

s = zeros(1, length(m));
for i=0:(length(m)-1)
    for j=1:length(m)
        s(i+1) = s(i+1) + m(j)*m(mod(j+i, length(m))+1);
    end
end

plot(s, '.');
