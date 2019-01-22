%% Width and velocity plots for numerical simulations
% clear all
close all
clc

vnumlocal = zeros(3,21);
% vnumoverall = zeros(1,11);
d = zeros(1,21);

for j = 1:3
for i = 0:20
    latticed = 6 + 0.1*i;
    d(i+1) = latticed;
    
    if j==1 
        numdata = importdata(sprintf('../SimulationsWithoutExperiments/ExperimentalSplinedBistableMagneticChainDisplacementsLattice%.1fRail21p5.txt',latticed));
    end
    
    if j==2
        numdata = importdata(sprintf('../SimulationsWithoutExperiments/ExperimentalSplinedBistableMagneticChainDisplacementsLattice%.1fRail22.txt',latticed));
    end
    
    if j==3
        numdata = importdata(sprintf('../SimulationsWithoutExperiments/ExperimentalSplinedBistableMagneticChainDisplacementsLattice%.1fRail22p5.txt',latticed));
    end
    timestep = 0.25;
    time = timestep*(0:length(numdata(1,:))-1);
    el = 8;

    % h = plot(time,-numdata(el,:),'b--',time,-numdata(el+1,:),'m--',time,-numdata(el+2,:),'r--',time,-numdata(el+3,:),'black--');
    % set(h,'LineWidth',1.5)
    % set(gca,'FontSize',14)
    % axis([0 900 -5.9 0.5])
    % xlabel('Time (ms)')
    % ylabel('Displacement (cm)')
    % legend('Element 1 (numerics)','Element 2 (numerics)', 'Element 3 (numerics)','Element 4 (numerics)')

    vnumlocal(j,i+1) = ArrivalTimeNumerical(numdata,el,latticed,3);
    
end
end

% vexp = [0.1845 02787 0.2885 0.1792 0.238 0.1800 0.1203];
% stdvexp = [0.0090 0.0045 0.0100 0.0020 0.0010 0.0081 0.0042];

dexp22p5 = [6 7 8];
dexpother = [6 8];

matrix(1,:) = [1,2,3];
matrix(2,:) = [0,0,0];
matrix(3,:) = [2,3,4];
matrix(4,:) = [1,3,4];
matrix(5,:) = [0,0,0];
matrix(6,:) = [2,3,4];
matrix(7,:) = [1,2,3];
matrix(8,:) = [1,2,3];
matrix(9,:) = [2,3,4];

[vexp stdvexp wexp stdwexp] = FindVelocityAndWidth(matrix);

vexp21p5 = 10*[vexp(1,1) vexp(1,3)];
stdvexp21p5 = 10*[stdvexp(1,1) stdvexp(1,3)];
stdxexp21p5 = 0.1*ones(length(stdvexp21p5),1);

vexp22 = 10*[vexp(2,1) vexp(2,3)];
stdvexp22 = 10*[stdvexp(2,1) stdvexp(2,3)];
stdxexp22 = 0.1*ones(length(stdvexp22),1);

vexp22p5 = 10*vexp(3,:);
stdvexp22p5 = 10*stdvexp(3,:);
stdxexp22p5 = 0.1*ones(length(stdvexp22p5),1);

leg = LEGEND_GENERATOR({'b-','r--','-.k','ob','^r','vk'},{'num (21.5 cm)','num (22 cm)','num (22.5 cm)', 'exp (21.5 cm)', 'exp (22 cm)', 'exp (22.5 cm)'}...
    ,18,8,2);
figure(1)
h = plot(d,10*vnumlocal(1,:),'b',d,10*vnumlocal(2,:),'--r',d,10*vnumlocal(3,:),'-.k');
hold on
p1 = ploterr(dexpother,vexp21p5,stdxexp21p5,stdvexp21p5,'ob');
p2 = ploterr(dexpother,vexp22,stdxexp22,stdvexp22,'^r');
p3 = ploterr(dexp22p5,vexp22p5,stdxexp22p5,stdvexp22p5,'vk');
set(h,'LineWidth',2)
set(p1,'LineWidth',2)
set(p2,'LineWidth',2)
set(p3,'LineWidth',2)
set(h,'LineWidth',2)
set(gca,'Fontsize',26)
set(leg,'Fontsize',16)
xlabel('Lattice distance (cm)')
ylabel('Velocity (m/s)')
% leg = legend('numerics (21.5 cm)','numerics (22 cm)','numerics (22.5 cm)', 'experiment (21.5 cm)', 'experiment (22 cm)', 'experiment (22.5 cm)')
% set(leg,'LineWidth',2)
% set(leg,'FontSize',18)
axis([5.8, 8.2, 0.5, 3])