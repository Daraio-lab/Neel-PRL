function [E, rhov1, maxpvel] = EnergyVelocity(displacements,timestep,t1,t2)
x = 0:length(displacements(:,1))-1;
time = timestep*(0:length(displacements(1,:))-1);

En = zeros(length(x),1);

pvelocity = zeros(length(x),length(time)-1);

for i=1:length(displacements(1,:))-1
    pvelocity(:,i) = (displacements(:,i+1)-displacements(:,i))/timestep;
end

maxpvel = max(max(abs(pvelocity)));

index1 = zeros(1,length(displacements(1,:)));

for k=1:length(index1)
    [~, idx1] = min(abs(displacements(:,k)-1.9));
    index1(k) = idx1;
end
% figure
% plot(index1,'r')

t = timestep*(t1:t2);

vel1 = polyfit(t,index1(t1:t2),1);

for j=1:length(time)-1
    En(j) = 0.5*trapz(pvelocity(:,j).^2);
end
% figure
% plot(En)

E = mean(En(t1:t2));

rhov1 = vel1(1);
end