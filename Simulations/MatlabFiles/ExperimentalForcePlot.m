%% Experimental Bistable Force Plot
% close all
clear all
clc

file = load('../Experiments/Quasi-staicTestsRealGeometrySpec5.mat');
data1 = file.DispContrSnapPreStress175mmSpec5Run33FL19;
data2 = file.DispContrSnapPreStress175mmSpec5Run35FL20;
data3 = file.DispContrSnapPreStress180mmSpec5Run30FL17;
data4 = file.DispContrSnapPreStress180mmSpec5Run29FL16;
data5 = file.DispContrSnapPreStress185mmSpec5Run22FL12;
data6 = file.DispContrSnapPreStress180mmSpec5Run32FL18;


for i = 2:length(data1(:,1))
    energy(i) = trapz(data1(1:i,3),data1(1:i,4));
end

plot(data1(:,3),energy,data1(:,3),data1(:,1));

figure
plot(data1(:,3),10*data1(:,1),data2(:,3),10*data2(:,1),data3(:,3),10*data3(:,1),data4(:,3),10*data4(:,1),data5(:,3),10*data5(:,1),data6(:,3),10*data6(:,1))
xlabel('displacement (cm)')
ylabel('energy (10^-1 J)')
legend('sample175_1','sample175_2','sample180_1','sample180_2','sample185_1','sample180_3')

p = polyfit(0.1*data3(:,3),data3(:,1)*10^-2,8);

a = roots(polyder(p));
eqmpoints = sort(a(a==real(a)));

u = -0.4:0.001:5.0;

energyfit = polyval(p,u);
figure
h = plot(0.1*data3(:,3),data3(:,1)*10^-2,'x',u,energyfit);
set(h,'LineWidth',2)
set(gca,'FontSize',16)
xlabel('Out-of-plane displacement (cm)')
ylabel('Energy (10^{-1} J)')
legend('Experimental result','Numerical fit')

q = polyder(p);
forcefit = polyval(q,u);

figure
plot(0.1*data3(:,3),0.1*data3(:,4),u,forcefit)
xlabel('displacement (cm)')
ylabel('force (10 N)')