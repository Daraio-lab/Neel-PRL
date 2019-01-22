function [E, rhov] = EnergyVelocity(displacements,mass,timestep)
x = 0:length(displacements(:,1))-1;
time = timestep*(0:length(displacements(1,:))-1);

% figure
% plot(time, displacements(80,:))

En = zeros(length(x),1);

pvelocity = zeros(length(x),length(time)-1);

for i=1:length(displacements(1,:))-1
    pvelocity(:,i) = (displacements(:,i+1)-displacements(:,i))/timestep;
end

index = zeros(1,2499);

for k=1:length(index)
    A = find(pvelocity(:,k)==max(pvelocity(:,k)));
    index(k) = A(1);
end
% figure
% plot(index)

t = timestep*(500:1500);

vel = polyfit(t,index(500:1500),1);

for j=1:length(time)-1
    En(j) = 0.5*mass*trapz(pvelocity(:,j).^2);
end
% figure
% plot(time(1:length(time)-1),En,'linewidth',2)
% set(gca,'fontsize', 24);
% xlabel('Time')
% ylabel('Energy')
% 
% figure
% contourf(time,x,displacements)
% set(gca,'fontsize', 24);
% xlabel('Time')
% ylabel('Nodal position')

E = mean(En(500:1500));

rhov = mass*vel(1);

end