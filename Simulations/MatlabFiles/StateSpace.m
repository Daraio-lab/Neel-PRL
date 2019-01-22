clc
close all
clear all

latticed = 8;
raild = 21.5;

if (rem(raild,2)==0)    
    data = importdata(sprintf('../SimulationsWithoutExperiments/ExperimentalSplinedBistableMagneticChainDisplacementsLattice%d.0Rail%d.txt',latticed,raild));
else    
    data = importdata(sprintf('../SimulationsWithoutExperiments/ExperimentalSplinedBistableMagneticChainDisplacementsLattice%d.0Rail%dp5.txt',latticed,raild-rem(raild,1)));
end

displacements1 = data(1:length(data(:,1)),:);
disp = 10^-1*displacements1(10,:);
vel = 400*diff(displacements1(10,:));
acc = 400*diff(vel);

figure
plot(disp)
hold on
plot(vel,'r')
hold on
plot(acc,'g')


vel = [vel 0];
acc = [acc 0 0];

figure
plot(disp,vel)
xlabel('u')
ylabel('u_\xi')

figure
plot(vel,acc)
xlabel('u_\xi')
ylabel('u_{\xi\xi}')

figure
plot(disp,acc)
xlabel('u')
ylabel('u_{\xi\xi}')


figure
n=50;
disp_pT = [zeros(1,n-1) disp(1:10000-n+1)];
plot(disp_pT)
hold on
plot(disp,'r')

figure
plot(disp,disp_pT)

q = 6;
n = 133;
disp_pT = [zeros(1,n-1) displacements1(q,1:10000-n+1)];

figure
hold on
plot(displacements1(q+1,:),'b')
plot(displacements1(17,:),'r')
plot(disp_pT,'r')