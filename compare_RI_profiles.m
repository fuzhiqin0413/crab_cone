meanConeL = 470;

n0 = 1.52;

radius = 0:100;

RIConeEqn = n0*sech(pi*radius/2/radius);

RIOrignal = 1.52-0.000004914*radius.^2;

RIEndLongitudional = 1.5+(0.02)-0.00000612*radius.^2;

figure; hold on

plot(radius, RIConeEqn);
plot(radius, RIOrignal);
plot(radius, RIEndLongitudional);
ylim([1.35 1.55])