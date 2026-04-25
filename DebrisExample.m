% program DebrisExample1.m

% this example does the following:

% 1) creates a steady state debris-free glacier
% 2) adds debris and runs to a steady state at the same climate
% 3) changes the ELA so that the glacier retreats to a new steady state
% 4) reruns the retreat of (3) but with some cryokarst

% After these computations, the following is plotted:

% a) profiles of glaciers and their debris thicknesses
% b) transient length versus time during the retreat w/ and w/o cryokarst

% To run this example, the following files are required:
% DEBISO.m, config1.m, config2.m, config3.m, config4.m


[stor1] = DEBISO(500,0.01,0,1,'config1');
  
files.H0 = stor1.H(:,end);
files.D = stor1.D(:,end);
 
[stor2] = DEBISO(3000,0.01,1,1,'config2',files);

files.H0 = stor2.H(:,end);
files.D = stor2.D(:,end);
 
[stor3] = DEBISO(1000,0.01,1,1,'config3',files);

[stor4] = DEBISO(1000,0.01,1,1,'config4',files);

figure(1)
clf
subplot(2,1,1)
plot(stor1.x,stor1.H(:,end)+stor1.b,'--m')
hold on
plot(stor1.x,stor2.H(:,end)+stor1.b,'b')
plot(stor1.x,stor3.H(:,end)+stor1.b,'g')
plot(stor1.x,stor1.b,'k')
ylim([1500 3500])
ylabel('Elevation (m)')
legend('c=0.00% ELA = 3000m', 'c=0.25% ELA = 3000m', 'c=0.25% ELA = 3100m')
title('Steady state profiles')

subplot(2,1,2)
plot(stor2.x,stor2.D(:,end),'b')
hold on
plot(stor2.x,stor3.D(:,end),'g')
ylabel('Debris thickness (m)')
xlabel('Distance (m)')

xt = zeros(1,length(stor4.t));
xt_cryo = zeros(1,length(stor4.t));
for j = 1:length(stor4.t)
    ind = find(stor3.H(:,j)<=0,1);
    xt(j) = stor3.x(ind);
    ind = find(stor4.H(:,j)<=0,1);
    xt_cryo(j) = stor4.x(ind);
end

figure(2)
clf
plot(stor4.t,xt,'b')
hold on
plot(stor4.t,xt_cryo,'r')
xlabel('Time (yr)')
ylabel('Extent (m)')
legend('$\lambda_m = 0\%$','$\lambda_m=20\%$','Interpreter','Latex')
title('Retreat with and without cryokarst')

