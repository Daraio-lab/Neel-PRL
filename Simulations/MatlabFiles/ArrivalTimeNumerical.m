%% Arrival time measurement of numerical simulations

function vlocal = ArrivalTimeNumerical(numdata,el,latticed,disp)

[idx1, idx1] = min(abs(numdata(el,:)-disp));
[idx2, idx2] = min(abs(numdata(el+1,:)-disp));% max(abs(diff(numdata(el+1,:))));
[idx3, idx3] = min(abs(numdata(el+2,:)-disp));% max(abs(diff(numdata(el+2,:))));
[idx4, idx4] = min(abs(numdata(el+3,:)-disp));% max(abs(diff(numdata(el+3,:))));

% figure
% hold on
% plot(abs(diff(numdata(el,:))),'b')
% plot(abs(diff(numdata(el+1,:))),'g')
% plot(abs(diff(numdata(el+2,:))),'r')
% plot(abs(diff(numdata(el+3,:))),'k')

times = 0.25*[idx1, idx2, idx3, idx4];
position = latticed*[1,2,3,4];

% figure
% plot(times,position,'*')

fit = polyfit(times,position,1);

vlocal = fit(1);
% [E,v1,v2] = EnergyVelocity(numdata,latticed,0.25,500,900);
% voverall = v1;
end