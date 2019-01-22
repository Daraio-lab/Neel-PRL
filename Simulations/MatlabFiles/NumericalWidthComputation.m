function [cm_width, particle_width] = NumericalWidthComputation()

cm_width = zeros(3,21);
particle_width = zeros(3,21);

latticed = 6.0+(0:20)*0.1;

for k = 1:3
for j = 1:21

    if k==1
        numdata = importdata(sprintf('../SimulationsWithoutExperiments/ExperimentalSplinedBistableMagneticChainDisplacementsLattice%.1fRail21p5.txt',latticed(j)));
    elseif k==2
        numdata = importdata(sprintf('../SimulationsWithoutExperiments/ExperimentalSplinedBistableMagneticChainDisplacementsLattice%.1fRail22.txt',latticed(j)));
    elseif k==3
        numdata = importdata(sprintf('../SimulationsWithoutExperiments/ExperimentalSplinedBistableMagneticChainDisplacementsLattice%.1fRail22p5.txt',latticed(j)));
    end
    
strain = zeros(length(numdata(:,1))-1,length(numdata(1,:)));    
for count = 1:length(numdata(:,1))-1
    strain(count,:) = -numdata(count+1,:)+numdata(count,:);    
end

% figure
% hold on
% plot(strain(8,:))
% plot(strain(9,:))
% plot(strain(10,:))

x = 1:5;
T = find(strain(9,:)==max(strain(9,:)));
T = T(1);
strainInSpace = [strain(7,T),strain(8,T),strain(9,T),strain(10,T),strain(11,T)];

[xout,~] = intersections(x,strainInSpace,x,strain(9,T)*ones(length(x),1)/2,1);
% figure
% hold on
% plot(strainInSpace)
% plot(xout,yout,'ro-')

particle_width(k,j) = xout(2)-xout(1);
cm_width(k,j) = latticed(j)*(xout(2)-xout(1));

end
end
end