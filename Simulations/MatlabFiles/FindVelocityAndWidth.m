%% Computing the values of width and velocity
function [vmat stdvmat wmat stdwmat] = FindVelocityAndWidth(matrix)

Raildistance = [21.5 22 22.5];
Latticedistance = [6 7 8];

wmat = zeros(3,3);
stdwmat = zeros(3,3);

vmat = zeros(3,3);
stdvmat = zeros(3,3);

for i = 1:length(Raildistance)
    for j = 1:length(Latticedistance)
        if (Raildistance(i)==22 && Latticedistance(j)==7) || (Raildistance(i)==21.5 && Latticedistance(j)==7) 
            continue
        end
        [vmat(i,j),stdvmat(i,j),wmat(i,j),stdwmat(i,j)] = WidthVelocity(Raildistance(i),Latticedistance(j),3,matrix((i-1)*3+j,:),3);       
    end
end

figure
hold on
errorbar(Raildistance,10*vmat(:,1),10*stdvmat(:,1),'ob')
errorbar(Raildistance,10*vmat(:,2),10*stdvmat(:,2),'or')
errorbar(Raildistance,10*vmat(:,3),10*stdvmat(:,3),'og')
hold off
axis([21.4,22.6,0,4])
xlabel('Rail distance (cm)')
ylabel('Velocity (m/s)')
legend('lattice distance = 6 cm','lattice distance = 7 cm','lattice distance = 8 cm')

figure
hold on
errorbar(Raildistance,wmat(:,1),stdwmat(:,1),'ob')
errorbar(Raildistance,wmat(:,2),stdwmat(:,2),'or')
errorbar(Raildistance,wmat(:,3),stdwmat(:,3),'og')
hold off
axis([21.4,22.6,0,16])
xlabel('Rail distance (cm)')
ylabel('Width of traveling wave (cm)')
legend('lattice distance = 6 cm','lattice distance = 7 cm','lattice distance = 8 cm')

end