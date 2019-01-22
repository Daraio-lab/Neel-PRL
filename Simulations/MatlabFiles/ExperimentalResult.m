%% Experimental wave propagation

clc
clear all
close all

file = importdata('../Experiments/RailDistance225Lattice6test14000fps.csv');

data = file.data;
time = data(:,1)/4000;
h = plot(time,0.1*data(:,7),time,0.1*data(:,17),time,0.1*data(:,27),time,0.1*data(:,37));
set(h,'LineWidth',2)
set(gca,'FontSize',16)
xlabel('Time (s)')
ylabel('Out-of-plane displacement (cm)')
legend('Element 1','Element 2','Element 3','Element 4')