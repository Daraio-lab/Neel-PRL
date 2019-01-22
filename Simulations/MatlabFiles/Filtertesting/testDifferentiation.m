% differentiating the noisy signal

close all
clear all
clc

raild = 22.5;
latticed = 7;
test = 1;

if raild==22 
       file1 = importdata(sprintf('../../Data/RailDistance%dLattice%dtest%d4000fps.csv',raild,latticed,test));
       numdata = importdata(sprintf('../../SimulationsWithoutExperiments/ExperimentalSplinedBistableMagneticChainDisplacementsLattice%.1fRail%.0f.txt',latticed,raild));
else
       file1 = importdata(sprintf('../../Data/RailDistance%.1fLattice%dtest%d4000fps.csv',raild,latticed,test));
       numdata = importdata(sprintf('../../SimulationsWithoutExperiments/ExperimentalSplinedBistableMagneticChainDisplacementsLattice%.1fRail%.0fp5.txt',latticed,raild-1));
end

pvelocity = zeros(length(numdata(:,1)),length(numdata(1,:))-1);    
for count = 1:length(numdata(1,:))-1
    pvelocity(:,count) = 4000*(numdata(:,count+1)-numdata(:,count))/10;    
end

data1 = file1.data;
expdata1 = [data1(:,7) data1(:,17) data1(:,27) data1(:,37)];

time = 0:length(expdata1(:,1))-1;

Fs = 4000; % 4000 Hz

p = 1;
q = 0.01;

% figure(1)
% pwelch(expdata1(:,1),[],[],[],Fs)

[b,a] = butter(p,q);  
butterexpdata1 = filter(b,a,expdata1(:,1));
[savgolaywithoutdiffexpdata1, savgolaywithoutdiffvexpdata1, ~] = SavGolay(expdata1(:,1),Fs,8,81);
savgolayexpdata1 = sgolayfilt(expdata1(:,1),6,41);

[b,a] = butter(p,q);  
butterexpdata2 = filter(b,a,expdata1(:,2));
[savgolaywithoutdiffexpdata2, savgolaywithoutdiffvexpdata2, ~] = SavGolay(expdata1(:,2),Fs,8,81);
savgolayexpdata2 = sgolayfilt(expdata1(:,2),6,41);

[b,a] = butter(p,q);  
butterexpdata3 = filter(b,a,expdata1(:,3));
[savgolaywithoutdiffexpdata3, savgolaywithoutdiffvexpdata3, ~] = SavGolay(expdata1(:,3),Fs,8,41);
savgolayexpdata3 = sgolayfilt(expdata1(:,3),6,41);

[b,a] = butter(p,q);  
butterexpdata4 = filter(b,a,expdata1(:,4));
[savgolaywithoutdiffexpdata4, savgolaywithoutdiffvexpdata4, ~] = SavGolay(expdata1(:,4),Fs,8,81);
% [p,q] = butter(1,0.02);  
% savgolaybutterdoublefilt = butter(p,q,savgolaywithoutdiffvexpdata4);
savgolayexpdata4 = sgolayfilt(expdata1(:,4),6,41);

% figure(2)
% plot(time,expdata1(:,1),time,butterexpdata1,time,savgolayexpdata1)
% hold on
% plot(savgolaywithoutdiffexpdata1,'k')
% legend('original','butterworth','savitsy-golay with diff','savitsy-golay without diff')
% 

figure(3)
hold on
plot(pvelocity(10,:),'g')
% plot(-savgolaybutterdoublefilt*10^-2,'b')% plot(-4000*diff(expdata1(:,1))*10^-2,'b')
plot(-4000*diff(butterexpdata3)*10^-2,'r')
% plot(-4000*diff(savgolayexpdata1)*10^-2,'m')
% plot(-savgolaywithoutdiffvexpdata1*10^-2,'k')
legend('numerical','savgolay double filtered','butterworth','savitsy-golay with diff','savitsy-golay without diff')
hold off

stexpdata1 = (expdata1(:,2)-expdata1(:,1))*10^-1;
butterst1 = (butterexpdata2-butterexpdata1)*10^-1;
savgolayst1 = (savgolayexpdata2-savgolayexpdata1)*10^-1;
stexpdata2 = (expdata1(:,3)-expdata1(:,2))*10^-1;
butterst2 = (butterexpdata3-butterexpdata2)*10^-1;
savgolayst2 = (savgolayexpdata3-savgolayexpdata2)*10^-1;
stexpdata3 = (expdata1(:,4)-expdata1(:,3))*10^-1;
butterst3 = (butterexpdata4-butterexpdata3)*10^-1;
savgolayst3 = (savgolayexpdata4-savgolayexpdata3)*10^-1;

% figure(4)
% hold on
% plot((expdata1(:,2)-expdata1(:,1))*10^-1,'b')
% plot((butterexpdata2-butterexpdata1)*10^-1,'r')
% plot((savgolayexpdata2-savgolayexpdata1)*10^-1,'k')
% plot((expdata1(:,3)-expdata1(:,2))*10^-1,'b--')
% plot((butterexpdata3-butterexpdata2)*10^-1,'r--')
% plot((savgolayexpdata3-savgolayexpdata2)*10^-1,'k--')
% plot((expdata1(:,4)-expdata1(:,3))*10^-1,'b-.')
% plot((butterexpdata4-butterexpdata3)*10^-1,'r-.')
% plot((savgolayexpdata4-savgolayexpdata3)*10^-1,'k-.')
% legend('original','butterworth','savitsy-golay with diff','savitsy-golay without diff')
% hold off

Traw = find(stexpdata2==max(stexpdata2));
Tbutter = find(butterst2==max(butterst2));
Tsavgolay = find(savgolayst2==max(savgolayst2));

strainraw = [0,stexpdata1(Traw),stexpdata2(Traw),stexpdata3(Traw),0];
strainbutter = [0,butterst1(Tbutter),butterst2(Tbutter),butterst3(Tbutter)];
strainsavgolay = [0,savgolayst1(Tsavgolay),savgolayst2(Tsavgolay),savgolayst3(Tsavgolay),0];
dispraw = [expdata1(Traw,1),expdata1(Traw,2),expdata1(Traw,3),expdata1(Traw,4)];
dispbutter = [butterexpdata1(Tbutter),butterexpdata2(Tbutter),butterexpdata3(Tbutter),butterexpdata4(Tbutter)];
dispsavgolay = [savgolayexpdata1(Tsavgolay),savgolayexpdata2(Tsavgolay),savgolayexpdata3(Tsavgolay),savgolayexpdata4(Tsavgolay)];

x = 1:5;

[xraw,yraw] = intersections(x,strainraw,x,strainbutter(3)*ones(length(x),1)/2,1);

% figure(5)
% hold on
% plot(strainraw,'b')
% plot(strainbutter,'r')
% plot(strainsavgolay,'k')
% plot(xraw,yraw,'g-')
% 
% figure(6)
% hold on
% plot(dispraw,'b')
% plot(dispbutter,'r')
% plot(dispsavgolay,'k')

xout = xraw(2)-xraw(1);
% Nf = 50;
% Fpass = 1000;
% Fstop = 1200;
% 
% filt = designfilt('differentiatorfir','FilterOrder',Nf, ...
%     'PassbandFrequency',Fpass,'StopbandFrequency',Fstop, ...
%     'SampleRate',Fs);
% 
% fvtool(filt,'MagnitudeDisplay','zero-phase','Fs',Fs)