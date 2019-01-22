function [cmwidth, stdcmwidth, particlewidth, stdparticlewidth, vmaxexp, stdvmaxexp] = WidthComputationExperimental()

cmwidth = zeros(7,1);
particlewidth = zeros(7,1);
stdcmwidth = zeros(7,1);
vmaxexp = zeros(7,1);
stdvmaxexp = zeros(7,1);
stdparticlewidth = zeros(7,1);
F = 4000;
timestep = 1/F;
latticed = [6,8,6,8,6,7,8];
raild =  [21.5,21.5,22,22,22.5,22.5,22.5];

for i = 1:7
   
   if raild(i)==22 
       file1 = importdata(sprintf('../Data/RailDistance%dLattice%dtest%d4000fps.csv',raild(i),latticed(i),1));
       file2 = importdata(sprintf('../Data/RailDistance%dLattice%dtest%d4000fps.csv',raild(i),latticed(i),2));
       file3 = importdata(sprintf('../Data/RailDistance%dLattice%dtest%d4000fps.csv',raild(i),latticed(i),3));
   else
       file1 = importdata(sprintf('../Data/RailDistance%.1fLattice%dtest%d4000fps.csv',raild(i),latticed(i),1));
       file2 = importdata(sprintf('../Data/RailDistance%.1fLattice%dtest%d4000fps.csv',raild(i),latticed(i),2));
       file3 = importdata(sprintf('../Data/RailDistance%.1fLattice%dtest%d4000fps.csv',raild(i),latticed(i),3));
   end
   
   exp1 = file1.data;
   exp2 = file2.data;
   exp3 = file3.data;
   % [b,a] = butter(1,0.1);  
%    expdata1 = zeros(length(exp1(:,1)),4);
%    expdata2 = zeros(length(exp2(:,1)),4);
%    expdata3 = zeros(length(exp3(:,1)),4);
   
   for j = 1:4
       [expdata1(:,j), ~, ~] = SavGolay(exp1(:,(j-1)*10+7),F,8,41);
       [expdata2(:,j), ~, ~] = SavGolay(exp2(:,(j-1)*10+7),F,8,41);
       [expdata3(:,j), ~, ~] = SavGolay(exp3(:,(j-1)*10+7),F,8,41);
   end
   
   vdata1 = -(1/timestep)*(10^-2)*[diff(expdata1(:,1)),diff(expdata1(:,2)),diff(expdata1(:,3)),diff(expdata1(:,4))];
   vdata2 = -(1/timestep)*(10^-2)*[diff(expdata2(:,1)),diff(expdata2(:,2)),diff(expdata2(:,3)),diff(expdata2(:,4))];
   vdata3 = -(1/timestep)*(10^-2)*[diff(expdata3(:,1)),diff(expdata3(:,2)),diff(expdata3(:,3)),diff(expdata3(:,4))];
   
%    plot(vdata1(:,1)-vdata1(:,2))
   
%    figure
%    hold on
%    plot(vdata1(:,1),'b');
%    plot(vdata1(:,2),'r');
%    plot(vdata1(:,3),'g');
%    plot(vdata1(:,4),'k');
%    hold off
   x = 1:4;
   
   el = 3;
   
   max_index1 = vdata1(:,el)==max(vdata1(:,el));
   max_index2 = vdata2(:,el)==max(vdata2(:,el));
   max_index3 = vdata3(:,el)==max(vdata3(:,el));
   
   vmax = zeros(3,1);
   vmax(1) = max(vdata1(:,el));
   vmax(2) = max(vdata2(:,el));
   vmax(3) = max(vdata3(:,el));
   
   [xout1,yout1] = intersections(x,vdata1(max_index1,:),x,max(vdata1(:,el))*ones(length(x),1)/2,1);
   [xout2,yout2] = intersections(x,vdata2(max_index2,:),x,max(vdata2(:,el))*ones(length(x),1)/2,1);
   [xout3,yout3] = intersections(x,vdata3(max_index3,:),x,max(vdata3(:,el))*ones(length(x),1)/2,1);
   
%    figure
%    plot(x,vdata1(max_index1,:),xout1,yout1,'or')
%  
   width = zeros(3,1);
   width(1) = xout1(2)-xout1(1);
   width(2) = xout2(2)-xout2(1);
   width(3) = xout3(2)-xout3(1);
   
   cmwidth(i) = latticed(i)*mean(width);
   particlewidth(i) = mean(width);
   
   stdcmwidth(i) = latticed(i)*std(width);
   stdparticlewidth(i) = std(width);
   
   vmaxexp(i) = mean(vmax);
   stdvmaxexp(i) = std(vmax);
   
   clear('expdata1');
   clear('expdata2');
   clear('expdata3');

end
end