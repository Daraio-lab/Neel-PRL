% differentiating the noisy signal

function [w, stdw] = ExperimentalWidthComputation(raild,latticed) 
% raild = 21.5;
% latticed = 6;

test = [1,2,3];
width = zeros(length(test),1);

for counter = 1:length(test)

if raild==22 
       file1 = importdata(sprintf('../Data/RailDistance%dLattice%dtest%d4000fps.csv',raild,latticed,test(counter)));
else
       file1 = importdata(sprintf('../Data/RailDistance%.1fLattice%dtest%d4000fps.csv',raild,latticed,test(counter)));
end

data1 = file1.data;
expdata1 = [data1(:,7) data1(:,17) data1(:,27) data1(:,37)];

time = 0:length(expdata1(:,1))-1;

Fs = 4000; % 4000 Hz

% figure(1)
% pwelch(expdata1(:,1),[],[],[],Fs)

[b,a] = butter(1,0.02);  
butterexpdata1 = filter(b,a,expdata1(:,1));
[savgolaywithoutdiffexpdata1, savgolaywithoutdiffvexpdata1, p] = SavGolay(expdata1(:,1),4000,8,81);
savgolayexpdata1 = sgolayfilt(expdata1(:,1),6,41);

[b,a] = butter(1,0.02);  
butterexpdata2 = filter(b,a,expdata1(:,2));
[savgolaywithoutdiffexpdata2, savgolaywithoutdiffvexpdata2, p] = SavGolay(expdata1(:,2),4000,8,81);
savgolayexpdata2 = sgolayfilt(expdata1(:,2),6,41);


[b,a] = butter(1,0.02);  
butterexpdata3 = filter(b,a,expdata1(:,3));
[savgolaywithoutdiffexpdata3, savgolaywithoutdiffvexpdata3, p] = SavGolay(expdata1(:,3),4000,8,41);
savgolayexpdata3 = sgolayfilt(expdata1(:,3),6,41);


[b,a] = butter(1,0.02);  
butterexpdata4 = filter(b,a,expdata1(:,4));
[savgolaywithoutdiffexpdata4, savgolaywithoutdiffvexpdata4, p] = SavGolay(expdata1(:,4),4000,8,41);
savgolayexpdata4 = sgolayfilt(expdata1(:,4),6,41);


% figure(2)
% plot(time,expdata1(:,1),time,butterexpdata1,time,savgolayexpdata1)
% hold on
% plot(savgolaywithoutdiffexpdata1,'k')
% % legend('original',
% legend('butterworth','savitsy-golay with diff','savitsy-golay without diff')
% 
% figure(3)
% hold on
% % plot(-4000*diff(expdata1(:,1))*10^-2,'b')
% plot(-4000*diff(butterexpdata1)*10^-2,'r')
% plot(-4000*diff(savgolayexpdata1)*10^-2,'m')
% plot(-savgolaywithoutdiffvexpdata1*10^-2,'k')
% legend('original','butterworth','savitsy-golay with diff','savitsy-golay without diff')
% hold off

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
% legend('original','butterworth','savitsy-golay with diff')%,'savitsy-golay without diff')
% hold off

Traw = find(stexpdata2==max(stexpdata2));
Tbutter = find(butterst2==max(butterst2));
Tsavgolay = find(savgolayst2==max(savgolayst2));

strainraw = [0,stexpdata1(Traw),stexpdata2(Traw),stexpdata3(Traw),0];
strainbutter = [0,butterst1(Tbutter),butterst2(Tbutter),butterst3(Tbutter),0];
strainsavgolay = [0,savgolayst1(Tsavgolay),savgolayst2(Tsavgolay),savgolayst3(Tsavgolay),0];
dispraw = [expdata1(Traw,1),expdata1(Traw,2),expdata1(Traw,3),expdata1(Traw,4)];
dispbutter = [butterexpdata1(Tbutter),butterexpdata2(Tbutter),butterexpdata3(Tbutter),butterexpdata4(Tbutter)];
dispsavgolay = [savgolayexpdata1(Tsavgolay),savgolayexpdata2(Tsavgolay),savgolayexpdata3(Tsavgolay),savgolayexpdata4(Tsavgolay)];

x = 1:5;

[xraw,~] = intersections(x,strainsavgolay,x,strainsavgolay(3)*ones(length(x),1)/2,1);

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

width(counter) = xraw(2)-xraw(1);

end

w = mean(width);
stdw = std(width);

end
% Nf = 50;
% Fpass = 1000;
% Fstop = 1200;
% 
% filt = designfilt('differentiatorfir','FilterOrder',Nf, ...
%     'PassbandFrequency',Fpass,'StopbandFrequency',Fstop, ...
%     'SampleRate',Fs);
% 
% fvtool(filt,'MagnitudeDisplay','zero-phase','Fs',Fs)