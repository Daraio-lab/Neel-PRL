%% forcedata
close all
clear all
clc

file3 = importdata('MagneticForce.txt');
zwickdata = file3.data;

n = 9500;

Displacement = 0.1*(142-zwickdata(:,3));
Force = 0.1*(zwickdata(:,4));

logD = log(Displacement(n:11457));
logF = log(Force(n:11457));

coeff = polyfit(logD,logF,1);
fit = polyval(coeff,logD);

% figure
% plot(logD,logF,logD,fit,'r')
% xlabel('Log(D)')
% ylabel('Log(F)')

p = coeff(1);
A = exp(coeff(2));

d = 1:0.001:12.5;
F = A*(d).^p;

figure
hold on
plot(Displacement,10*Force,'sr','markers',4)
h=plot(d,10*F,'--k');
set(h,'LineWidth',2.5)
set(gca,'FontSize',26)
axis([1, 8, 0, 30])
xlabel('Displacement (cm)')
ylabel('Force (N)')
[legh,objh,outh,outm] = legend('experiments','theory');
set(objh,'linewidth',2.5);
box on