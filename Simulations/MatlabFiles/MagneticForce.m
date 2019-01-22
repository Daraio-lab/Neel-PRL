%% Magnetic Force
clc
clear all
close all

A = 6*6.25*10^-6;
n = -2.73;

x1 = 0;
x2 = 8;
l = x2 - x1; % x2 - x1

u1 = 0;
u2 = -6.5:0.01:1;

y1 = x1 + u1;
y2 = x2 + u2;

F2 = A*((y2-y1)*0.01).^n - A*(l*0.01).^n;

figure
plot(u2,F2)
xlabel('displacement (cm)')
ylabel('force (10 N)')

