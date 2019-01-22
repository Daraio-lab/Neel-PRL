%% Visualization program

clc
clear all
close all

file = importdata('RailDistance22.5Lattice7test14000fps.csv');
data = file.data;

framesanalyzedpertotal=1;
frameRate = 4000*framesanalyzedpertotal;
time = data(:,1)/frameRate;
 
y = length(data(:,1));
x = 3000;

l1 = mean(-data(x:y,7));
l2 = mean(-data(x:y,17));
l3 = mean(-data(x:y,27));
l4 = mean(-data(x:y,37));

l = [l1 l2 l3 l4]

v1 = mean(-data(x:y,6));
v2 = mean(-data(x:y,16));
v3 = mean(-data(x:y,26));
v4 = mean(-data(x:y,36));

w1 = mean(-data(x:y,5));
w2 = mean(-data(x:y,15));
w3 = mean(-data(x:y,25));
w4 = mean(-data(x:y,35));

mean(v1/l1+v2/l2+v3/l3+v4/l4);
mean(w1/l1+w2/l2+w3/l3+w4/l4);

% figure
% plot(time,data(:,7),time,data(:,6),time,data(:,5))
% title('Displacement vs Time (Rail Distance 22.5 Lattice 6)')
% xlabel('Time (s)')
% ylabel('Displacement (mm)')
% legend('Z (out-of-plane) displacement', 'X displacement', 'Y displacement')




figure
plot(time,data(:,7)/l1,time,data(:,17)/l2,time,data(:,27)/l3,time,data(:,37)/l4)
title('Displacement vs Time (Rail Distance 22.5 Lattice 6)')
xlabel('Time (s)')
ylabel('Displacement (mm)')
legend('element 1', 'element 2', 'element 3', 'element 4')

figure
fig = plot(time,data(:,7),time,data(:,17),time,data(:,27),time,data(:,37));
title('Displacement vs Time (Rail Distance 22.5 Lattice 6)')
xlabel('Time (s)')
ylabel('Displacement (mm) ')
legend('element 1', 'element 2', 'element 3', 'element 4')

% 
nstrain1 = data(:,7)/l1-data(:,17)/l2;
nstrain2 = data(:,17)/l2-data(:,27)/l3;
nstrain3 = data(:,27)/l3-data(:,37)/l4;

figure
plot(time, -nstrain1, time, -nstrain2, time, -nstrain3)
title('Normalized Strain vs Time (Rail Distance 22.5 Lattice 6)')
xlabel('Time (s)')
ylabel('Strain (mm)')
legend('element 1', 'element 2', 'element 3')

strain1 = data(:,7)-data(:,17);
strain2 = data(:,17)-data(:,27);
strain3 = data(:,27)-data(:,37);

figure
plot(time, -strain1, time, -strain2, time, -strain3)
title('Strain vs Time (Rail Distance 22.5 Lattice 6)')
xlabel('Time (s)')
ylabel('Strain (mm)')
legend('element 1', 'element 2', 'element 3')

% 
% velocity1 = diff(data(:,7));
% velocity2 = diff(data(:,17));
% velocity3 = diff(data(:,27));
% velocity4 = diff(data(:,37));
% time2 = data(1:length(data(:,1))-1,1)/frameRate;
% % 
% % figure
% % plot(time2,velocity1)
% 
% l1 = mean(-data(x:y,7));
% l2 = mean(-data(x:y,17));
% l3 = mean(-data(x:y,27));
% l4 = mean(-data(x:y,37));
% 
% 
% disp1 = data(:,7)/l1;
% disp2 = data(:,17)/l2;
% disp3 = data(:,27)/l3;
% disp4 = data(:,37)/l4;
% 
% fn=0.02;
% n=6;
%  
% [bb,a] = butter(n, fn);
%  
% inc_pulse1 = -diff(filtfilt(bb,a,disp1));
% inc_pulse2 = -diff(filtfilt(bb,a,disp2));
% inc_pulse3 = -diff(filtfilt(bb,a,disp3));
% inc_pulse4 = -diff(filtfilt(bb,a,disp4));
% 
% figure
% plot(time, disp1, time, filtfilt(bb,a,disp1))
% figure
% plot(time2, inc_pulse1, time2, inc_pulse2, time2, inc_pulse3, time2, inc_pulse4)
% 
% find(inc_pulse4==max(inc_pulse4))-find(inc_pulse3==max(inc_pulse3))
% find(inc_pulse3==max(inc_pulse3))-find(inc_pulse2==max(inc_pulse2)) 
% find(inc_pulse2==max(inc_pulse2))-find(inc_pulse1==max(inc_pulse1))