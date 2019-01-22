% Displacement check
clear all
close all
clc

displacements = importdata('../Cases/HyperelasticLatticeDiffusive.txt');

figure
plot(displacements(:,4000),'linewidth',2)
axis([0,600,-0.1,2.1])
set(gca,'fontsize', 24);
xlabel('Nodal position')
ylabel('Displacements')

mass = 1;
timestep = 0.01;

[E,v] = EnergyVelocity(displacements,timestep,1000,1500);

sum = zeros(length(displacements(1,:)),1);

for j = 1:length(sum)
    for i = 1:length(disp)-1
        sum(j) = sum(j) + (displacements(i+1,j)-displacements(i,j)).^2;
    end
end
