%% Force displacement
close all
clear all
clc

Displacement3 = [1.7 2.2 2.7 3.2 3.7 4.2 4.7 5.2 5.7 6.2];
Force3 = [2.170 0.990 0.510 0.290 0.165 0.104 0.057 0.032 0.017 0.006];

logD3 = log(Displacement3);
logF3 = log(Force3);

coeff = polyfit(logD3,logF3,1);
fit = polyval(coeff,logD3);

figure
plot(logD3,logF3,'o',logD3,fit,'r')
xlabel('Log(D)')
ylabel('Log(F)')

p3 = coeff(1);
A3 = exp(coeff(2));

d = 1:0.001:7.5;
F3 = A3*(d).^p3;

figure
h=plot(Displacement3,Force3,'o',d,F3,'r');
set(h,'LineWidth',2)
set(gca,'FontSize',16)
ylabel('Force [10 N]')
xlabel('Displacement [cm]')
legend('experiments','numerical fit')