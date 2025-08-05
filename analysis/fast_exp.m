clc;
clear;
clf;

fexp = @(x) typecast(int32(single(x)*single(12102203)) + int32(1065353216), 'single');

x = 0:1:255;
k = single(-1/96);

hold on;
plot(x, fexp(x*k), 'r-', 'linesmooth', 'on');
plot(x, exp(x*k), 'b-', 'linesmooth', 'on');
hold off;
axis([0 255 0 1]);
grid on;
legend('Fast exp', 'Std exp');
