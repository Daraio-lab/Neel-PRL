function [cm_width, particle_width, vmaxnum] = WidthComputation()

cm_width = zeros(3,21);
particle_width = zeros(3,21);
vmaxnum = zeros(3,21);

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
    
    timestep = 1/4000;
    x = 1:length(numdata(:,1));

    pvelocity = zeros(length(numdata(:,1)),length(numdata(1,:)-1));
    for i = 1:length(numdata(1,:))-1
        pvelocity(:,i) = (numdata(:,i+1)-numdata(:,i))*0.1/timestep;
    end

    element = 10;

    max_index = find(pvelocity(element,:)==max(pvelocity(element,:)));
    vmaxnum(k,j) = max(pvelocity(element,:));

%     figure
%     hold on
%     plot(x,pvelocity(:,max_index(1)),'-*')
%     plot(1:length(numdata(:,1)),max(pvelocity(element,:))*ones(length(x),1)/2,'k')
%     axis([0,26,-5,15])
%     
    [xout,yout] = intersections(x,pvelocity(:,max_index(1)),x,max(pvelocity(element,:))*ones(length(x),1)/2,1);
    plot(xout,yout,'or')

    cm_width(k,j) = latticed(j)*(xout(2)-xout(1));
    particle_width(k,j) = (xout(2)-xout(1));
end
end

% figure
% l = plot(latticed,cm_width(1,:),'b',latticed,cm_width(2,:),'--r',latticed,cm_width(3,:),'-.k');
% set(l,'LineWidth',2)
% set(gca,'FontSize',26)
% axis([6 8 7 15])
% xlabel('Lattice distance (cm)')
% ylabel('FWHM (cm)')
% leg1 = legend('rail distance (21.5 cm)','rail distance (22 cm)','rail distance (22.5 cm)');
% set(leg1,'LineWidth',1)
% set(leg1,'FontSize',18)
% 
% figure
% h = plot(latticed,particle_width(1,:),'b',latticed,particle_width(2,:),'--r',latticed,particle_width(3,:),'-.k');
% set(h,'LineWidth',2)
% set(gca,'FontSize',26)
% axis([6 8 1 2.2])
% xlabel('Lattice distance (cm)')
% ylabel('FWHM (no. of particles)')
% leg2 = legend('rail distance (21.5 cm)','rail distance (22 cm)','rail distance (22.5 cm)');
% set(leg2,'LineWidth',1)
% set(leg2,'FontSize',18)

end