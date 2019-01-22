%% Wave velocity vs max particle velocity

[~,~,vmaxnum] = WidthComputation();

latticed = [6,6,6,7,8,8,8];
raild = [21.5,22,22.5,22.5,21.5,22,22.5];
vmaxexp = zeros(length(raild),1);
stdvmaxexp = zeros(length(raild),1);

for k = 1:length(raild)
   [vmaxexp(k),stdvmaxexp(k)] = ExperimentalParticleVelocity(raild(k),latticed(k));   
end

WaveVelocityPlot

leg = LEGEND_GENERATOR({'b-','r--','-.k','ob','^r','vk'},{'num (21.5 cm)','num (22 cm)','num (22.5 cm)', 'exp (21.5 cm)', 'exp (22 cm)', 'exp (22.5 cm)'}...
    ,18,8,2);
figure(1)
hold on
h = plot(10*vnumlocal(1,:),vmaxnum(1,:),'-b',10*vnumlocal(2,:),vmaxnum(2,:),'--r',...
    10*vnumlocal(3,:),vmaxnum(3,:),'-.k');
xerr = 0.1*ones(1,length(latticed));
p1 = ploterr(vexp21p5,vmaxexp(1:2),stdvexp21p5,stdvmaxexp(1:2),'ob');
p2 = ploterr(vexp22,vmaxexp(3:4),stdvexp22,stdvmaxexp(3:4),'^r');
p3 = ploterr(vexp22p5,vmaxexp(5:7),stdvexp22p5,stdvmaxexp(5:7),'vk');
set(p1,'LineWidth',2)
set(p2,'LineWidth',2)
set(p3,'LineWidth',2)
set(h,'LineWidth',2)
set(gca,'Fontsize',26)
set(leg,'Fontsize',16)
xlabel('Wave Velocity (m/s)','FontSize',22)
ylabel('Max Particle Velocity (m/s)','FontSize',22)
% leg2 = legend('num (21.5 cm)','num (22 cm)','num (22.5 cm)','exp (21.5 cm)','exp (22 cm)','exp (22.5 cm)');
% set(leg2,'LineWidth',1)
% set(leg2,'FontSize',18)
axis([0.8 3.1 8 28])
