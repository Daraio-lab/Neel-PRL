function [vmax,stdvmax] = ExperimentalParticleVelocity(raild,latticed)

test = 1:3;
velmax = zeros(3,1);

for i=1:3
if raild==22 
       file1 = importdata(sprintf('../Data/RailDistance%dLattice%dtest%d4000fps.csv',raild,latticed,test(i)));
       numdata = importdata(sprintf('../SimulationsWithoutExperiments/ExperimentalSplinedBistableMagneticChainDisplacementsLattice%.1fRail%.0f.txt',latticed,raild));
else
       file1 = importdata(sprintf('../Data/RailDistance%.1fLattice%dtest%d4000fps.csv',raild,latticed,test(i)));
       numdata = importdata(sprintf('../SimulationsWithoutExperiments/ExperimentalSplinedBistableMagneticChainDisplacementsLattice%.1fRail%.0fp5.txt',latticed,raild-1));
end

pvelocity = zeros(length(numdata(:,1)),length(numdata(1,:))-1);    
for count = 1:length(numdata(1,:))-1
    pvelocity(:,count) = 4000*(numdata(:,count+1)-numdata(:,count))/10;    
end

data1 = file1.data;
expdata1 = [data1(:,7) data1(:,17) data1(:,27) data1(:,37)];

p = 1;
q = 0.01;

% [b,a] = butter(p,q);  
% butterexpdata1 = filter(b,a,expdata1(:,1));
% 
% [b,a] = butter(p,q);  
% butterexpdata2 = filter(b,a,expdata1(:,2));

[b,a] = butter(p,q);  
butterexpdata3 = filter(b,a,expdata1(:,3));

% [b,a] = butter(p,q);  
% butterexpdata4 = filter(b,a,expdata1(:,4));

figure(i)
hold on
plot(pvelocity(10,:),'g')
plot(-4000*diff(butterexpdata3)*10^-2,'r')
legend('numerical','butterworth')
% hold off

velmax(i) = max(-4000*diff(butterexpdata3)*10^-2);

end

vmax = mean(velmax);
stdvmax = std(velmax);
end