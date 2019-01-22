%% Force displacement curve of simulations

clc
clear all
close all

energy = load('../SimLatticeDistance8Lattice22p5/EnergyDisplacement.txt');
force = load('../SimLatticeDistance8Lattice22p5/ForceDisplacement.txt');

figure
plot(energy(:,1),energy(:,2))
xlabel('displacement (cm)')
ylabel('energy (10^{-1} J)')

figure
plot(force(:,1),force(:,2))
xlabel('displacement (cm)')
ylabel('force (10 N)')