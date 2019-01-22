%% Finding second equilibrium point
clear all
close all
clc

eqpoint21p5 = 0;
eqpoint22 = 0;
eqpoint22p5 = 0;

for i =3:4
for k = 1:3
    file = importdata(sprintf('../SimLatticeDistance%dLattice21p5/RailDistance21.5Lattice%dtest%d4000fps.csv',2*i,2*i,k));
    data = file.data;
for j = 0:3
    eqpoint21p5 = eqpoint21p5 + mean(data(length(data(:,1))-1000:length(data(:,1)),10*j+7));
end

end
end

eqpoint21p5 = eqpoint21p5/24;

for i =3:4
for k = 1:3
    file = importdata(sprintf('../SimLatticeDistance%dLattice22/RailDistance22Lattice%dtest%d4000fps.csv',2*i,2*i,k));
    data = file.data;
for j = 0:3
    eqpoint22 = eqpoint22 + mean(data(length(data(:,1))-1000:length(data(:,1)),10*j+7));
end

end
end

eqpoint22 = eqpoint22/24;

for i =3:4
for k = 1:3
    file = importdata(sprintf('../SimLatticeDistance%dLattice22p5/RailDistance22.5Lattice%dtest%d4000fps.csv',2*i,2*i,k));
    data = file.data;
for j = 0:3
    eqpoint22p5 = eqpoint22p5 + mean(data(length(data(:,1))-1000:length(data(:,1)),10*j+7));
end

end
end

eqpoint22p5 = eqpoint22p5/24;