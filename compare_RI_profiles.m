clear; close all; clc

meanConeL = 472;

maxR = 103; %at 425 um
maxRHeight = 425;

endR = 34; %at 15
endRHeight = 15;

figure; set(gcf,'position',[92         418        1522         537])
cols = lines(7);

% At max R 
subplot(1,2,1); hold on

radius = 0:maxR;

radialTopRI = 1.5+(1*0.01+0.01)-(0.5*1+0.8)*0.00000612*(radius/maxR*80).^2;
radialBaseRI = 1.5+(0*0.01+0.01)-(0.5*0+0.8)*0.00000612*(radius/maxR*80).^2;
radialBothRI = 1.5+(((meanConeL-maxRHeight)/meanConeL)*0.01+0.01)-(0.5*((meanConeL-maxRHeight)/meanConeL)+0.8)*0.00000612*(radius/maxR*80).^2;

cylinderRI = 1.52*sech(pi*radius/2/meanConeL);
cylinderRIp36 = 1.52*sech(pi*radius/2/(meanConeL+36));
cylinderRIm36 = 1.52*sech(pi*radius/2/(meanConeL-36));

cylinderNormMax = 1.52*sech(pi*(radius/maxR*maxR)/2/meanConeL);
cylinderNorm80 = 1.52*sech(pi*(radius/maxR*80)/2/meanConeL);

h1 = plot(radius, radialTopRI, 'linewidth',2,'color',cols(1,:));
h2 = plot(radius, radialBaseRI, 'linewidth',2,'color',cols(2,:));
h3 = line([0 110], [1 1]*(1.5+(((meanConeL-maxRHeight)/meanConeL)*0.01+0.01)), 'linewidth',2,'color',cols(3,:));
h4 = plot(radius, radialBothRI, 'linewidth',2,'color',cols(4,:));

h5 = plot(radius, cylinderRI, 'linewidth',2,'color',cols(5,:));
h6 = plot(radius, cylinderRIp36, 'linewidth',2,'color',cols(6,:));
h7 = plot(radius, cylinderRIm36, 'linewidth',2,'color',cols(7,:));

h8 = plot(radius, cylinderNormMax, ':', 'linewidth',2,'color',cols(1,:));
h9 = plot(radius, cylinderNorm80,  ':', 'linewidth',2,'color',cols(2,:));

h10 = line([0 110], [1.52 1.52], 'linewidth',2,'color',cols(3,:),'linestyle', ':');

ylim([1.40 1.55])
xlim([0 maxR])

legend([h1 h2 h3 h4 h5 h6 h7 h8 h9 h10], {'Radial-top', 'Radial-base', 'Linear', 'Radial+linear', ...
    'Theory-equal', 'Theory-plus', 'Theory-minus', 'Normalised-max', 'Normalised-80', 'Uniform'},'Location','southwest');

xlabel('Radius (um)'); ylabel('RI  at 425 um')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20);


% At tip
subplot(1,2,2); hold on
radius = 0:endR;

radialTopRI = 1.5+(1*0.01+0.01)-(0.5*1+0.8)*0.00000612*(radius/endR*80).^2;
radialBaseRI = 1.5+(0*0.01+0.01)-(0.5*0+0.8)*0.00000612*(radius/endR*80).^2;
radialBothRI = 1.5+(((meanConeL-endRHeight)/meanConeL)*0.01+0.01)-(0.5*((meanConeL-endRHeight)/meanConeL)+0.8)*0.00000612*(radius/endR*80).^2;

cylinderRI = 1.52*sech(pi*radius/2/meanConeL);
cylinderRIp36 = 1.52*sech(pi*radius/2/(meanConeL+36));
cylinderRIm36 = 1.52*sech(pi*radius/2/(meanConeL-36));

cylinderNormMax = 1.52*sech(pi*(radius/endR*maxR)/2/meanConeL);
cylinderNorm80 = 1.52*sech(pi*(radius/endR*80)/2/meanConeL);

h1 = plot(radius, radialTopRI, 'linewidth',2,'color',cols(1,:));
h2 = plot(radius, radialBaseRI, 'linewidth',2,'color',cols(2,:));
h3 = line([0 110], [1 1]*(1.5+(((meanConeL-endRHeight)/meanConeL)*0.01+0.01)), 'linewidth',2,'color',cols(3,:));
h4 = plot(radius, radialBothRI, 'linewidth',2,'color',cols(4,:));

h5 = plot(radius, cylinderRI, 'linewidth',2,'color',cols(5,:));
h6 = plot(radius, cylinderRIp36, 'linewidth',2,'color',cols(6,:));
h7 = plot(radius, cylinderRIm36, 'linewidth',2,'color',cols(7,:));

h8 = plot(radius, cylinderNormMax, ':', 'linewidth',2,'color',cols(1,:));
h9 = plot(radius, cylinderNorm80,  ':', 'linewidth',2,'color',cols(2,:));

h10 = line([0 110], [1.52 1.52], 'linewidth',2,'color',cols(3,:),'linestyle', ':');

ylim([1.40 1.55])
xlim([0 endR])

legend([h1 h2 h3 h4 h5 h6 h7 h8 h9 h10], {'Radial-top', 'Radial-base', 'Linear', 'Radial+linear', ...
    'Theory-equal', 'Theory-plus', 'Theory-minus', 'Normalised-max', 'Normalised-80', 'Uniform'},'Location','southwest');

xlabel('Radius (um)'); ylabel('RI at 15 um')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20);

cd('/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Figures')
exportgraphics(gcf,'RIProfiles.pdf','ContentType','vector') 