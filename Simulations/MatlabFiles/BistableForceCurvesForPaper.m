%% Experimental force plot
clc
clear all
close all

x = [0, 1.75, 5.465 - 0.236, 6.26400];
y = [0, 0.04493, -0.2007, 0];

u1 = x(1)-2:0.001:x(2);
u2 = x(2):0.001:x(3);
u3 = x(3):0.001:x(4)+1.;

F1 = 0. + 0.0513486*u1 - 0.014671*u1.^2;
F2 = -0.204049 + 0.320277*u2 - 0.122133*u2.^2 + 0.0116667*u2.^3;
F3 = 4.92206 - 1.95936*u3 + 0.187356*u3.^2;
u = [u1 u2 u3];
F = [F1 F2 F3];
En = cumtrapz(u,F);

figure
h=ploterr(x,10*y,0.1*ones(length(x),1),0.001*ones(length(x),1),'or');
hold on
q=plot(x,10*y,'or');
[p,a,b]=plotyy(u,10*F,u,0.1*En);%,u2,10*F2,'k',u3,10*F3,'k');
set(h,'LineWidth',1.5)
set(p,'FontSize',18)
a.LineWidth=2;
b.LineWidth=2;
a.Color = 'k';
b.Color = 'w';
a.LineStyle = '--';
set(q,'LineWidth',1.5)
axis([-2.5,7,-2.5,1]);
xlabel('Displacement (cm)')
ylabel(p(1),'Force (N)')
ylabel(p(2),'Energy (mJ)')
legend('experiments','theory')
set(p,{'ycolor'},{'k';'k'})
set(p(1),'Position', [0.13 0.11 0.775-.08 0.815]);
set(p(2),'Position', [0.13 0.11 0.775-.08 0.815]);