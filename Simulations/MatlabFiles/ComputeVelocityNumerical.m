%% Velocity of numerical simulation
clear all
close all
clc

latticed = 6;
raild = 21.5;
test = 3;

if (rem(raild,2)==0)    
    file = importdata(sprintf('../SimLatticeDistance%dLattice%d/RailDistance%dLattice%dtest%d4000fps.csv',latticed,raild,raild,latticed,test));
    numdata = importdata(sprintf('../SimLatticeDistance%dLattice%d/ExperimentalSplinedBistableMagneticChainDisplacementsLattice%d.00Rail%d.txt',latticed,raild,latticed,raild));
else    
    file = importdata(sprintf('../SimLatticeDistance%dLattice%dp5/RailDistance%.1fLattice%dtest%d4000fps.csv',latticed,raild-rem(raild,1),raild,latticed,test));
    numdata = importdata(sprintf('../SimLatticeDistance%dLattice%dp5/ExperimentalSplinedBistableMagneticChainDisplacementsLattice%d.00Rail%dp5.txt',latticed,raild-rem(raild,1),latticed,raild-rem(raild,1)));
end

timestep = 0.25; % 0.25 ms
time = timestep*(0:length(numdata(1,:))-1);
t0 = 0.085;
el = 8;
[b,a] = butter(5,0.01);
expdata = filter(b,a,file.data);

vnum = ArrivalTimeNumerical(numdata,el,latticed,4.5);
[vexp,stdexpv,w,stdw] = WidthVelocity(raild,latticed,test,[1,2,3],3.4);

t = 1000*(expdata(:,1)/4000 + t0);
figure
h = plot(t,0.1*expdata(:,7),'b',time,-numdata(el,:),'b--',t,0.1*expdata(:,17),'m',time,-numdata(el+1,:),'m--',t,0.1*expdata(:,27),'r',time,-numdata(el+2,:),'r--',t,0.1*expdata(:,37),'black',time,-numdata(el+3,:),'black--');
set(h,'LineWidth',1.5)
set(gca,'FontSize',14)
axis([0 900 -5.9 0.5])
xlabel('Time (ms)')
ylabel('Displacement (cm)')
legend('Element 1 (experiment)','Element 1 (numerics)','Element 2 (experiment)','Element 2 (numerics)','Element 3 (experiment)', 'Element 3 (numerics)','Element 4 (experiment)','Element 4 (numerics)')

% plot(time(1:length(time)-1),diff(numdata(el,:))/0.25,time(1:length(time)-1),diff(numdata(el+1,:))/0.25,time(1:length(time)-1),diff(numdata(el+2,:))/0.25,time(1:length(time)-1),diff(numdata(el+3,:))/0.25)
% hold on
% plot(t(1:length(t)-1),diff(-0.1*expdata(:,7))/0.25,t(1:length(t)-1),diff(-0.1*expdata(:,17))/0.25,t(1:length(t)-1),diff(-0.1*expdata(:,27))/0.25,t(1:length(t)-1),diff(-0.1*expdata(:,37))/0.25)
