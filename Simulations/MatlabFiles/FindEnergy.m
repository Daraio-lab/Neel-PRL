%% Open displacements
clc
close all
clear all

% file = importdata('../SimLatticeDistance6Lattice22p5/RailDistance22.5Lattice6test14000fps.csv');
% data = file.data;
% displacements = [data(:,7), data(:,17), data(:,27), data(:,37)];
% x = 0:length(displacements(1,:))-1;
% timestep = 0.25; % 0.25 ms
% time = timestep*(0:length(displacements(:,1))-1);
% step = 25;

data = importdata('../SimLatticeDistance6Lattice22p5/ExperimentalRandomizedBistableMagneticChainDisplacementsLattice6Rail22p5.txt');
displacements = data(1:length(data(:,1)),:);
x = 0:length(displacements(:,1))-1;
timestep = 0.25; % 0.25 ms
time = timestep*(0:length(displacements(1,:))-1);
step = 25;

for j = 1:length(x)
    for i = 2:length(time)-1
        particlev(j,i) = (displacements(j,i+1)-displacements(j,i-1))/(2*timestep);
    end
end

figure
hold on
for i = 6:length(x)-4
   plot(0.5*particlev(i,:).^2);   
end
hold off

figure
hold on
for i = 6:length(x)
   plot(i,sum(0.5*particlev(i,:).^2)*timestep,'o');   
end
hold off