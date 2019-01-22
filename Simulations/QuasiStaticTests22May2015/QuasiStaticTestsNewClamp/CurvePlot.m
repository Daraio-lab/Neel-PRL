%% Quasistatic results
clear all
close all
clc

file = importdata('21p5NormalTestSpecimenm4.txt');

A = 17.5180;
p = -3.2744;

dinit = 11;

data = file.data;
[b,a] = butter(6,0.1);
force = filter(b,a,0.1*data(:,4));
force = force+A*dinit^p;%-force(1);
displacement = filter(b,a,0.1*data(:,3));
displacement = displacement-displacement(1);
plot(displacement,force)

displacement2 = (force/A).^(1/p)+displacement-dinit;
forcenew = force-A*dinit^p;

[M, idx] = max(forcenew);

figure
plot(displacement2,forcenew)
xlabel('displacement (cm)')
ylabel('force (N)')