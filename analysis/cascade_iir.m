clear;
clc;

n=32768;
t=(0:n-1)/3000;

s = cos(2*pi*26*t);
filtered = s;

lpf3k = [1,2,1,1,-1.98580607763713,0.987430879565983;1,2,1,1,-1.96234827644227,0.963953885016314;1,2,1,1,-1.94282576708918,0.944415402194706;1,2,1,1,-1.92888371699001,0.930461944602056;1,2,1,1,-1.92163389881911,0.923206194573720];
g = [0.000406200482212439;0.000401402143510437;0.000397408776382334;0.000394556903011068;0.000393073938652093;1];

for j = 1:length(g)-1
    lpf3k(j, 1:3) = lpf3k(j, 1:3).*g(j);
end

clf;
hold on;
plot(s, 'r');
plot(sosfilt(lpf3k,s));
hold off;

%%
z3 = zeros(size(lpf3k));

for i=1:n
    for sec=1:size(lpf3k, 1)
        z3(sec, 1) = filtered(i) - z3(sec, 2)*lpf3k(sec, 5) - z3(sec, 3)*lpf3k(sec, 6);
        filtered(i) = z3(sec, 1)*lpf3k(sec, 1) + z3(sec, 2)*lpf3k(sec, 2) + z3(sec, 3)*lpf3k(sec, 3);
        z3(sec, 3) = z3(sec, 2);
        z3(sec, 2) = z3(sec, 1);
    end
end

clf;
hold on;
plot(t, s, 'b');
plot(t, filtered, 'r');
hold off;
