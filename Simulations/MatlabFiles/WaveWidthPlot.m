clc
clear all
close all

[num_cmwidth, num_particlewidth] = NumericalWidthComputation();

latticed_exp = [6,8,6,8,6,7,8];
raild_exp = [21.5,21.5,22,22,22.5,22.5,22.5];
w = zeros(1,length(latticed_exp));
stdw = zeros(1,length(latticed_exp));

for counter = 1:length(w)
   [w(counter),stdw(counter)] = ExperimentalWidthComputation(raild_exp(counter),latticed_exp(counter)); 
end

leg = LEGEND_GENERATOR({'b-','r--','-.k','ob','^r','vk'},{'num (21.5 cm)','num (22 cm)','num (22.5 cm)', 'exp (21.5 cm)', 'exp (22 cm)', 'exp (22.5 cm)'}...
    ,18,8,2);
hold on
latticed_num = 6 + 0.1*(0:20);
h = plot(latticed_num,num_particlewidth(1,:),'-b',latticed_num,num_particlewidth(2,:),'--r',...
    latticed_num,num_particlewidth(3,:),'-.k');

xerr = 0.1*ones(1,length(latticed_exp));
p1 = ploterr(latticed_exp(1:2),w(1:2),xerr(1:2),stdw(1:2),'ob');
p2 = ploterr(latticed_exp(3:4),w(3:4),xerr(3:4),stdw(3:4),'^r');
p3 = ploterr(latticed_exp(5:7),w(5:7),xerr(5:7),stdw(5:7),'vk');
set(p1,'LineWidth',2)
set(p2,'LineWidth',2)
set(p3,'LineWidth',2)
set(h,'LineWidth',2)
set(gca,'Fontsize',26)
set(leg,'Fontsize',16)
xlabel('Lattice distance (cm)','FontSize',26)
ylabel('FWHM (elements)','FontSize',26)
% leg2 = legend('num (21.5 cm)','num (22 cm)','num (22.5 cm)','exp (21.5 cm)','exp (22 cm)','exp (22.5 cm)');
% set(leg2,'LineWidth',1)
% set(leg2,'FontSize',18)
axis([5.8 8.2 0. 2.5])
