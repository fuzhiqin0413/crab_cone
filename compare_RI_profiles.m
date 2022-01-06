meanConeL = 470;
expConeL = 600;

n0 = 1.52;

radiusFixed = 0:100;

RIConeEqn_Mean = n0*sech(pi*radiusFixed/2/meanConeL);

RIConeEqn_Exp = n0*sech(pi*radiusFixed/2/expConeL);

radiusFiber = 0:0.005:0.1;

RIFiberEqn_matchFn = sqrt(n0^2*(1-(3.25)^2*(radiusFiber).^2));

RIFiberEqn_matchTip = sqrt(n0^2*(1-(2.2)^2*(radiusFiber).^2));

% Note - might want to do a range of max radiuses so this pans out

radius = 0:80;

RIOrignal = 1.52-0.000004914*(radius/max(radius)*80).^2;

% RIEndLongitudional = 1.5+(0.02)-0.00000612*(radius/max(radius)*80).^2;

RIEndLongitudional = 1.5+(1*0.01+0.01)-(0.5*1+0.8)*0.00000612*(radius/max(radius)*80).^2;

RIEndTest = 1.5+(1*0.01+0.01)-(0.5*1+0.8)*0.00000612*80^2*(radius/max(radius)).^2;

RIStartLongitudional = 1.5+(0*0.01+0.01)-(0.5*0+0.8)*0.00000612*(radius/max(radius)*80).^2;

figure; hold on

h1 = plot(radiusFixed, RIConeEqn_Mean);
h2 = plot(radiusFixed, RIConeEqn_Exp);
h3 = plot(radius, RIOrignal);
h4 = plot(radius, RIEndLongitudional);
h5 = plot(radius, RIStartLongitudional);

h6 = plot(radiusFiber*1000, RIFiberEqn_matchFn);
h7 = plot(radiusFiber*1000, RIFiberEqn_matchTip);

plot = plot(radius, RIEndTest, 'kx');

ylim([1.35 1.55])

legend([h1 h2 h3 h4 h5 h6 h7], 'Lens Cylinder Eqn - 470', 'Lens Cylinder Eqn - 600',...
    'Original radial graident', 'Longitudional at tip', 'Longitudional at base', ...
    'Fiber Eqn', 'Fiber Eqn To Match')