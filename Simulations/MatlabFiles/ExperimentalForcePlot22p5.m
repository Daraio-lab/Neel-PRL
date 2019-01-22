% Experimental Force Plot: lattice distance = 22.5
clear all
close all
clc

data5 = xlsread('../Experiments/Quasi-StaticTest225mmSpec5Run2.xlsx');

p = polyfit(0.1*data5(:,3),data5(:,1)*10^-2,8);

a = roots(polyder(p));
eqmpoints = sort(a(a==real(a)));

u = -0.4:0.001:4.5;

energyfit = polyval(p,u);
figure
h = plot(0.1*data5(:,3),data5(:,1),'x',u,energyfit*10^2,'--r');
set(h,'LineWidth',2)
set(gca,'FontSize',16)
xlabel('Out-of-plane displacement (cm)')
ylabel('Energy (mJ)')
legend('experimental result','numerical fit')


q = polyder(p);
forcefit = polyval(q,u);

figure
plot(0.1*data5(:,3),0.1*data5(:,4),u,forcefit)
xlabel('displacement (cm)')
ylabel('force (10 N)')

r = polyder(q);
stiffnessfit = polyval(r,u);

figure
plot(u,stiffnessfit)
xlabel('displacement (cm)')
ylabel('stiffness (N/mm)')

