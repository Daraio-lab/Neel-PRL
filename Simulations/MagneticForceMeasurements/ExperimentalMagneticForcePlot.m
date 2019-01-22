%% Magnetic force measurements
clc
close all
clear all

% file1 = importdata('MagForceMeasurementsData.mat');
file2 = importdata('MagForceMeasurementsCorrectData.mat');
file = importdata('MagForceMeasurementsCorrectDataMay2015.mat');
Force = 10^-1*file(2:9,3); % W in mN
Displacement = 0.1*file(2:9,1); % D in mm 

Force2 = 10^-3*file2.Mass; % W in mN
Displacement2 = 0.1*file2.Distance; % D in mm 

file3 = importdata('MagneticForce.txt');
zwickdata = file3.data;

logD = log(Displacement);
logF = log(Force);

logD2 = log(Displacement2(2:12));
logF2 = log(Force2(2:12));

coeff = polyfit(logD,logF,1);
fit = polyval(coeff,logD);

coeff2 = polyfit(logD2,logF2,1);
fit2 = polyval(coeff2,logD2);

figure
plot(logD,logF,'o',logD,fit,'r')
xlabel('Log(D)')
ylabel('Log(F)')

p = coeff(1);
A = exp(coeff(2));

p2 = coeff2(1);
A2 = exp(coeff2(2));

d = 1:0.001:7.5;
F = A*(d).^p;
F2 = A2*(d).^p2;

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

F3 = A3*(d).^p3;

figure
h=plot(Displacement,Force,'o',d,F,Displacement2(2:12),Force2(2:12),'+',d,F2,Displacement3,Force3,'s',d,F3,0.1*(141-zwickdata(:,3)),0.1*(zwickdata(:,4)+0.001),'y');
set(h,'LineWidth',2)
set(gca,'FontSize',16)
ylabel('Force [10 N]')
xlabel('Displacement [cm]')
legend('experiments Andres1','numerical fit Andres1','experiments Andres2','numerical fit Andres2','experiments Neel','numerical fit Neel')

figure
plot(Displacement2(2:12),Force2(2:12),'o',0.1*(141-zwickdata(:,3)),0.1*(zwickdata(:,4)+0.001),'r');

