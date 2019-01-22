clc
close all
clear all

displacements1 = importdata('Stiffness10Nonlinearity1Gamma1.txt');
displacements2 = importdata('Stiffness150Nonlinearitym3Gamma1.txt');
displacements3 = importdata('Stiffness10Nonlinearity5Gamma1.txt');
displacements4 = importdata('Stiffness10Nonlinearity7Gamma1.txt');
x = 0:length(displacements1(:,1))-1;
[E1, Eth] = findEnergy(displacements1);
E2 = findEnergy(displacements2);
E3 = findEnergy(displacements3);
E4 = findEnergy(displacements4);

figure
plot(x,E1,'*',x,E2,'*',x,E3,'*',x,E4,'*',x,Eth)
xlabel('Nodal Point [n]')
axis([0 200 1 4])
ylabel('Energy at each nodal point [E(n) = \int_{-\infty}^{\infty} u_{n,t}(t)^2 dt]')
legend('nonlinearity = 1','nonlinearity = 3','nonlinearity = 5','nonlinearity = 7','theoretical energy')
% legend('stiffness = 10','stiffness = 50','stiffness = 100','stiffness = 150','theoretical energy')
timestep = 0.1;
time = timestep*(0:length(displacements2(2,:))-1);

figure
for i = 1:length(time)/10;
   plot(displacements2(:,i*10))
   axis([0 200 0 2])
   pause(0.1)  
end