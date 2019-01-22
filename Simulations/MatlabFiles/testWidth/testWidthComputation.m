close all
clc

numdata = importdata(sprintf('../../SimulationsWithoutExperiments/ExperimentalSplinedBistableMagneticChainDisplacementsLattice%.1fRail22p5.txt',8));

hold on
plot(-numdata(8,1:3200),'b')
plot(-numdata(9,1:3200),'r')
plot(-numdata(10,1:3200),'k')
plot(-numdata(11,1:3200),'g')