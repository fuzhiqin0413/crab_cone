meanConeL = 470;
expConeL = 600;

n0 = 1.52;

radius = 0:100;

RIConeEqn_Mean = n0*sech(pi*radius/2/meanConeL);

RIConeEqn_Exp = n0*sech(pi*radius/2/expConeL);

RIOrignal = 1.52-0.000004914*radius.^2;

RIEndLongitudional = 1.5+(0.02)-0.00000612*radius.^2;

figure; hold on

h1 = plot(radius, RIConeEqn_Mean);
h2 = plot(radius, RIConeEqn_Exp);
h3 = plot(radius, RIOrignal);
h4 = plot(radius, RIEndLongitudional);
ylim([1.35 1.55])

legend([h1 h2 h3 h4], 'Lens Cylinder Eqn - 470', 'Lens Cylinder Eqn - 600',...
    'Original radial graident', 'Longitudional at tip + radial gradient')