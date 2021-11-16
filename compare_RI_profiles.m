meanConeL = 470;
expConeL = 600;

n0 = 1.52;

radiusFixed = 0:100;

RIConeEqn_Mean = n0*sech(pi*radiusFixed/2/meanConeL);

RIConeEqn_Exp = n0*sech(pi*radiusFixed/2/expConeL);



% Note - might want to do a range of max radiuses so this pans out

radius = 0:100;

RIOrignal = 1.52-0.000004914*(radius/max(radius)*80).^2;

% RIEndLongitudional = 1.5+(0.02)-0.00000612*(radius/max(radius)*80).^2;

RIEndLongitudional = 1.5+(1*0.01+0.01)-(0.5*1+0.8)*0.00000612*(radius/max(radius)*80).^2;

RIStartLongitudional = 1.5+(0*0.01+0.01)-(0.5*0+0.8)*0.00000612*(radius/max(radius)*80).^2;

figure; hold on

h1 = plot(radiusFixed, RIConeEqn_Mean);
h2 = plot(radiusFixed, RIConeEqn_Exp);
h3 = plot(radius, RIOrignal);
h4 = plot(radius, RIEndLongitudional);
h5 = plot(radius, RIStartLongitudional);
ylim([1.35 1.55])

legend([h1 h2 h3 h4 h5], 'Lens Cylinder Eqn - 470', 'Lens Cylinder Eqn - 600',...
    'Original radial graident', 'Longitudional at tip', 'Longitudional at base')