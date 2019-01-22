function [E, Eth] = findEnergy(displacements)
x = 0:length(displacements(:,1))-1;
timestep = 0.1; % 0.25 ms
time = timestep*(0:length(displacements(1,:))-1);
gamma = 1;

psif = -3.09842;
Eth = ones(length(x),1)*-psif/gamma;

figure
plot(time, displacements(80,:))

E = zeros(length(x),1);

pvelocity = zeros(length(x),length(time)-1);

for i=1:length(displacements(1,:))-1
    pvelocity(:,i) = (displacements(:,i+1)-displacements(:,i))/timestep;
end

for j=1:length(x)
    E(j) = timestep*trapz(pvelocity(j,:).^2);
end

end