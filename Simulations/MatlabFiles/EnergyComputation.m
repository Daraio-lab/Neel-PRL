%% Open displacements
clc
close all
clear all

latticed = 8;
raild = 21.5;
test = 3;
t0 = 0.078;

if (rem(raild,2)==0)    
    file = importdata(sprintf('../SimLatticeDistance%dLattice%d/RailDistance%dLattice%dtest%d4000fps.csv',latticed,raild,raild,latticed,test));
    data = importdata(sprintf('../SimulationsWithoutExperiments/ExperimentalSplinedBistableMagneticChainDisplacementsLattice%d.0Rail%d.txt',latticed,raild));
else    
    file = importdata(sprintf('../SimLatticeDistance%dLattice%dp5/RailDistance%.1fLattice%dtest%d4000fps.csv',latticed,raild-rem(raild,1),raild,latticed,test));
    data = importdata(sprintf('../SimulationsWithoutExperiments/ExperimentalSplinedBistableMagneticChainDisplacementsLattice%d.0Rail%dp5.txt',latticed,raild-rem(raild,1)));
end

displacements1 = data(1:length(data(:,1)),:);
x = 1:length(displacements1(:,1));
timestep = 0.25; % 0.25 ms
time = timestep*(0:length(displacements1(1,:))-1);
step = 50;

el = 8;

figure
expdata = file.data;
t = 1000*(expdata(:,1)/4000 + t0);

for i=1:60
    expdata1_marked(i) = expdata(50*i,7);
    numdata1_marked(i) = -displacements1(el,50*i);
    expdata2_marked(i) = expdata(50*i,17);
    numdata2_marked(i) = -displacements1(el+1,50*i);
    expdata3_marked(i) = expdata(50*i,27);
    numdata3_marked(i) = -displacements1(el+2,50*i);
    expdata4_marked(i) = expdata(50*i,37);
    numdata4_marked(i) = -displacements1(el+3,50*i);
    t_marked(i) = t(50*i);
    time_marked(i) = time(50*i);
end

h = plot(t_marked,0.1*expdata1_marked,'bd',t_marked,0.1*expdata2_marked,'ms',t_marked,0.1*expdata3_marked,'rp',t_marked,0.1*expdata4_marked,'kh',time_marked,numdata1_marked,'bd',time_marked,numdata2_marked,'ms',time_marked,numdata3_marked,'rp',time_marked,numdata4_marked,'kh',...
t,0.1*expdata(:,7),'b',time,-displacements1(el,:),'b--',t,0.1*expdata(:,17),'m',time,-displacements1(el+1,:),'m--',t,0.1*expdata(:,27),'r',time,-displacements1(el+2,:),'r--',t,0.1*expdata(:,37),'k',time,-displacements1(el+3,:),'k--');
set(h,'LineWidth',2)
set(gca,'FontSize',22)
axis([0 600 -4.5 0.5])
xlabel('Time (ms)')
ylabel('Displacement (cm)')
leg = legend('Element 8', 'Element 9', 'Element 10', 'Element 11');
set(leg,'LineWidth',1)
set(leg,'FontSize',14)

style = '--';
color = zeros(4,3);
for i = 1:4
    color(i,:) = [0 0 0]+0.2*(i-1);
end

figure
contourf(displacements1)

velocity = zeros(length(displacements1(:,1)),length(displacements1(1,:))-1);
energy = zeros(1,length(displacements1(1,:))-1);

for i = 1:length(displacements1(1,:))-1
    velocity(:,i) = (displacements1(:,i+1)-displacements1(:,i))*10^-2/(timestep*10^-3);
    energy(i) = 0.5*sum(velocity(:,i).^2)*latticed*10^-2;
end

mass=21;
[en,momDensity] = EnergyVelocity(displacements1,mass*10^-3/(latticed*10^-2),timestep*10^-3);

timeForEnergy = timestep*10^-3*(0:length(energy)-1);
figure
plot(timeForEnergy,energy*10^3/momDensity,'b','LineWidth',2)
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('Energy per momentum density (mJs/kg)')

figure
plot(x,-displacements1(:,500),'bo-',x,-displacements1(:,1000),'ro-',...
x,-displacements1(:,1500),'ko-','Linewidth',2)
set(gca,'FontSize',18)
axis([1,20,-7,0.1])
xlabel('Element number')
ylabel('Displacements (cm)')
legend('t = 125 ms','t = 250 ms','t = 375 ms')