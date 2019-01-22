%% Width and Velocity Plot

function [meanv, stdv, meanwidth, stdwidth] = WidthVelocity(Raildistance,Latticedistance,numberofexperiments,membernumber,disp)

if nargin < 4
    membernumber = [1,2,3,4];
end

r = Raildistance; 
d = Latticedistance;
m = membernumber;
v = zeros(1,numberofexperiments);
w = zeros(1,numberofexperiments);

for i=1:numberofexperiments

if (rem(r,2)==0)    
file = importdata(sprintf('../SimLatticeDistance%dLattice%d/RailDistance%dLattice%dtest%d4000fps.csv',d,r,r,d,i));

else    
file = importdata(sprintf('../SimLatticeDistance%dLattice%dp5/RailDistance%.1fLattice%dtest%d4000fps.csv',d,r-rem(r,1),r,d,i));
end

[b,a] = butter(6,0.01);
expdata = filter(b,a,file.data);
t = 1000*(expdata(:,1)/4000);
% figure
% h = plot(t,0.1*expdata(:,7),'b',t,0.1*expdata(:,17),'g',t,0.1*expdata(:,27),'r',t,0.1*expdata(:,37),'k');
% set(h,'LineWidth',1.5)
% set(gca,'FontSize',14)
% axis([0 900 -7.0 0.5])
% xlabel('Time (ms)')
% ylabel('Displacement (cm)')
% legend('Element 1 (experiment)','Element 2 (experiment)','Element 3 (experiment)','Element 4 (experiment)')

index1 = zeros(1,length(m));
timewidth = zeros(1,length(m));

for k=1:length(index1)
    [~, idx1] = min(abs(0.1*expdata(:,7+10*(m(k)-1))+disp));
    meand = mean(expdata(length(expdata(:,1))-100:length(expdata(:,1)),7+10*(m(k)-1)));
    [~, idx2] = min(abs(expdata(:,7+10*(m(k)-1))-0.8*meand));
    [~, idx3] = min(abs(expdata(:,7+10*(m(k)-1))-0.2*meand));
      
    index1(k) = t(idx1);
    timewidth(k) = t(idx2)-t(idx3);
end

x = membernumber*d;

p = polyfit(index1,x,1);

% u = (10:0.001:40)*d;
% figure
% plot(u,polyval(p,u),'b',index1,x,'*r')

v(i) = p(1);
w(i) = p(1)*mean(timewidth);

strain = zeros(length(m)-1,length(expdata(:,1)));

for j = 1:length(m)-1
    strain(j,:)=0.1*expdata(:,7+10*(m(j+1)-1))-0.1*expdata(:,7+10*(m(j)-1));
end

% figure
% plot(t,strain(2,:),'*')

end

meanv = mean(v);
stdv = std(v);

meanwidth = mean(w);
stdwidth = std(w);
end