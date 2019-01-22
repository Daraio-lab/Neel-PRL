% clear all
close all
clc

raild = 22.5;
latticed = 7;
test = 3;

if (rem(raild,2)==0)    
    file = importdata(sprintf('../SimLatticeDistance%dLattice%d/RailDistance%dLattice%dtest%d4000fps.csv',latticed,raild,raild,latticed,test));
    numdata = importdata(sprintf('../SimLatticeDistance%dLattice%d/ExperimentalSplinedBistableMagneticChainDisplacementsLattice%d.00Rail%d.txt',latticed,raild,latticed,raild));
else    
    file = importdata(sprintf('../SimLatticeDistance%dLattice%dp5/RailDistance%.1fLattice%dtest%d4000fps.csv',latticed,raild-rem(raild,1),raild,latticed,test));
    numdata = importdata(sprintf('../SimLatticeDistance%dLattice%dp5/ExperimentalSplinedBistableMagneticChainDisplacementsLattice%d.00Rail%dp5.txt',latticed,raild-rem(raild,1),latticed,raild-rem(raild,1)));
end

timestep = 1/4000; % 0.25 ms
t0 = -0.25;
time = timestep*(0:length(numdata(1,:))-2)+t0;
el = 8;
[b,a] = butter(1,0.1);
expdata = filter(b,a,file.data);

vdata = -(1/timestep)*(10^-2)*[diff(expdata(:,7)),diff(expdata(:,17)),diff(expdata(:,27)),diff(expdata(:,37))];
t = timestep*(0:length(vdata(:,1))-1);

plot(time,(1/timestep)*(10^-1)*diff(numdata(9,:)),'b',t,vdata(:,1),'--b',time,(1/timestep)*(10^-1)*diff(numdata(10,:)),'r',t,vdata(:,2),'--r',...
    time,(1/timestep)*(10^-1)*diff(numdata(11,:)),'k',t,vdata(:,3),'--k',time,(1/timestep)*(10^-1)*diff(numdata(12,:)),'m',t,vdata(:,4),'m--')
axis([0,1,-5,20])

max_index_exp = vdata(:,3)==max(vdata(:,3));
max_index_num = find(diff(numdata(11,:))==max(diff(numdata(11,:))));

x = 1:4;
xexp = 1:4;
particlev = [(1/timestep)*(10^-1)*diff(numdata(9,:));(1/timestep)*(10^-1)*diff(numdata(10,:));(1/timestep)*(10^-1)*diff(numdata(11,:));(1/timestep)*(10^-1)*diff(numdata(12,:))];

[xout1,yout1] = intersections(x,particlev(:,max_index_num(1)),x,max(particlev(:,max_index_num(1)))*ones(length(x),1)/2,1);
[xout2,yout2] = intersections(x,vdata(max_index_exp,:),x,max(vdata(max_index_exp,:))*ones(length(x),1)/2,1);

figure
plot(xexp,particlev(:,max_index_num(1)),'b',xexp,vdata(max_index_exp,:),'--b',xout1,yout1,'vk',xout2,yout2,'^k')

[xout2(2)-xout2(1) xout1(2)-xout1(1)]

[(xout2(2)-xout2(1))/max(vdata(max_index_exp,:)) (xout1(2)-xout1(1))/max(particlev(:,max_index_num(1)))]
