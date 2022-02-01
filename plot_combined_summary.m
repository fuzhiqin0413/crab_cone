%%% Need to set figure size to get nicely shaped axes...

%%% Check cylinder cylinder for COLC, X Focus and Y Focus height (and radius?)

clear; close all; clc

saveFigures = 1;
saveFolder = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Figures_new/';

% Load reference acceptance function
acceptanceFunctionNight = readmatrix('/Users/gavintaylor/Documents/Matlab/Git_versioned_April_27/Crab_cone_git/Data/Acceptance Function Night.csv');
acceptanceFunctionDay = readmatrix('/Users/gavintaylor/Documents/Matlab/Git_versioned_April_27/Crab_cone_git/Data/Acceptance function Day.csv');

nPlots = 4;
cols = lines(7);
outlineCol = [0.5 0.5 0.5];

% Load data from varying RI series
cd('/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Data_new/Varying RI Profile/')
data_GRIN_both = load('Cone_1000_nm_Cone_0_SD_GRIN_both_SIMDATA.mat');
data_GRIN_radialBase = load('Cone_1000_nm_Cone_0_SD_GRIN_radialBase_SIMDATA.mat');
data_GRIN_radialTop = load('Cone_1000_nm_Cone_0_SD_GRIN_radialTop_SIMDATA.mat');
data_GRIN_linear = load('Cone_1000_nm_Cone_0_SD_GRIN_linear_SIMDATA.mat');

data_GRIN_radialTop_tipCorr = load('Cone_1000_nm_Cone_0_SD_GRIN_radialTop_TipCorrection_SIMDATA.mat');

data_GRIN_cylinder = load('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_SIMDATA.mat');
data_GRIN_cylinder_norm_max = load('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_normalised_max_SIMDATA.mat');
data_GRIN_cylinder_norm_80 = load('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_normalised_80_SIMDATA.mat');

data_GRIN_cylinder_m36 = load('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_m36_SIMDATA.mat');
data_GRIN_cylinder_p36 = load('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_p36_SIMDATA.mat');
data_GRIN_cylinder_shape = load('Cylinder_1000_nm_Cone_0_SD_GRIN_cylinder_SIMDATA.mat');

% Load data from internal structures (w/ varying RI)
cd('/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Data_new/Extra internal structures/Radial RI')
data_GRIN_radialTop_CinC = load('Cone_CinC_1000_nm_Cone_0_SD_GRIN_radialTop_SIMDATA.mat');
data_GRIN_radialTop_EC = load('Cone_EC_1000_nm_Cone_0_SD_GRIN_radialTop_SIMDATA.mat');
data_GRIN_radialTop_CinC_EC = load('Cone_CinC_EC_1000_nm_Cone_0_SD_GRIN_radialTop_SIMDATA.mat');

% Load data from varying SD (unfirom RI) series
cd('/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Data_new/Varying Cone SD')
data_uniform_m4SD = load('Cone_1000_nm_Cone_-4_SD_Uniform_1.52_SIMDATA.mat');
data_uniform_m2SD = load('Cone_1000_nm_Cone_-2_SD_Uniform_1.52_SIMDATA.mat');
data_uniform_0SD = load('Cone_1000_nm_Cone_0_SD_Uniform_1.52_SIMDATA.mat');
data_uniform_p2SD = load('Cone_1000_nm_Cone_2_SD_Uniform_1.52_SIMDATA.mat');
data_uniform_p4SD = load('Cone_1000_nm_Cone_4_SD_Uniform_1.52_SIMDATA.mat');

% Load data from internal structures (w/ uniform RI)
cd('/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Data_new/Extra internal structures/Uniform RI')
data_GRIN_uniform_CinC = load('Cone_CinC_1000_nm_Cone_0_SD_Uniform_1.52_SIMDATA.mat');
data_GRIN_uniform_EC = load('Cone_EC_1000_nm_Cone_0_SD_Uniform_1.52_SIMDATA.mat');
data_GRIN_uniform_CinC_EC = load('Cone_CinC_EC_1000_nm_Cone_0_SD_Uniform_1.52_SIMDATA.mat');

% Load data from the three varying intercone series.
cd('/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Data_new/Variable Intercone/Radial RI')
data_GRIN_radialTop_IC_142 = load('Cone_1000_nm_Cone_0_SD_GRIN_radialTop_IC_142_SIMDATA.mat');
data_GRIN_radialTop_IC_144 = load('Cone_1000_nm_Cone_0_SD_GRIN_radialTop_IC_144_SIMDATA.mat');
data_GRIN_radialTop_IC_146 = load('Cone_1000_nm_Cone_0_SD_GRIN_radialTop_IC_146_SIMDATA.mat');
data_GRIN_radialTop_IC_148 = load('Cone_1000_nm_Cone_0_SD_GRIN_radialTop_IC_148_SIMDATA.mat');

cd('/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Data_new/Variable Intercone/Uniform RI')
data_uniform_0SD_IC_142 = load('Cone_1000_nm_Cone_0_SD_Uniform_1.52_IC_142_SIMDATA.mat');
data_uniform_0SD_IC_144 = load('Cone_1000_nm_Cone_0_SD_Uniform_1.52_IC_144_SIMDATA.mat');
data_uniform_0SD_IC_146 = load('Cone_1000_nm_Cone_0_SD_Uniform_1.52_IC_146_SIMDATA.mat');
data_uniform_0SD_IC_148 = load('Cone_1000_nm_Cone_0_SD_Uniform_1.52_IC_148_SIMDATA.mat');

cd('/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Data_new/Variable Intercone/Theoretical RI')
data_GRIN_cylinder_IC_142 = load('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_IC_142_SIMDATA.mat');
data_GRIN_cylinder_IC_144 = load('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_IC_144_SIMDATA.mat');
data_GRIN_cylinder_IC_146 = load('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_IC_146_SIMDATA.mat');
data_GRIN_cylinder_IC_148 = load('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_IC_148_SIMDATA.mat');

% Load data from varying pigment series
cd('/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Data_new/Variable Pigment/Radial RI')
data_GRIN_radialTop_PG_140 = load('Cone_1000_nm_Cone_0_SD_GRIN_radialTop_PG_140_SIMDATA.mat');
data_GRIN_radialTop_PG_145 = load('Cone_1000_nm_Cone_0_SD_GRIN_radialTop_PG_145_SIMDATA.mat');
data_GRIN_radialTop_PG_150 = load('Cone_1000_nm_Cone_0_SD_GRIN_radialTop_PG_150_SIMDATA.mat');

cd('/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Data_new/Variable Pigment/Uniform RI')
data_uniform_0SD_PG_140 = load('Cone_1000_nm_Cone_0_SD_Uniform_1.52_PG_140_SIMDATA.mat');
data_uniform_0SD_PG_145 = load('Cone_1000_nm_Cone_0_SD_Uniform_1.52_PG_145_SIMDATA.mat');
data_uniform_0SD_PG_150 = load('Cone_1000_nm_Cone_0_SD_Uniform_1.52_PG_150_SIMDATA.mat');

cd('/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Data_new/Variable Pigment/Theoretical RI')
data_GRIN_cylinder_PG_140 = load('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_PG_140_SIMDATA.mat');
data_GRIN_cylinder_PG_145 = load('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_PG_145_SIMDATA.mat');
data_GRIN_cylinder_PG_150 = load('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_PG_150_SIMDATA.mat');

%%
% process acceptance angle from paper
% Fix zero at zero
acceptanceFunctionNight(abs(acceptanceFunctionNight(:,1)) < 0.1,1) = 0;

negInds = find(acceptanceFunctionNight(:,1) <= 0);
posInds = find(acceptanceFunctionNight(:,1) >= 0);

angleSteps = 0:10;
negInterp = interp1(-acceptanceFunctionNight(negInds,1), acceptanceFunctionNight(negInds,2), angleSteps, 'linear', 'extrap');
posInterp = interp1(acceptanceFunctionNight(posInds,1), acceptanceFunctionNight(posInds,2), angleSteps, 'linear', 'extrap');

interpAcceptanceNight = (negInterp+posInterp)/2;

lowInds = find(interpAcceptanceNight > 0.5);
if ~isempty(lowInds)
    % Take last
    lowInds = lowInds(end);
    
    slope = (interpAcceptanceNight(lowInds + 1) - interpAcceptanceNight(lowInds))/...
        (angleSteps(lowInds + 1) - angleSteps(lowInds));
    acceptanceAngleNight = angleSteps(lowInds) + (0.5 - interpAcceptanceNight(lowInds))/slope;
end

% do for day as well
acceptanceFunctionDay(abs(acceptanceFunctionDay(:,1)) < 0.1,1) = 0;

negInds = find(acceptanceFunctionDay(:,1) <= 0);
posInds = find(acceptanceFunctionDay(:,1) >= 0);

angleSteps = 0:10;
negInterp = interp1(-acceptanceFunctionDay(negInds,1), acceptanceFunctionDay(negInds,2), angleSteps, 'linear', 'extrap');
posInterp = interp1(acceptanceFunctionDay(posInds,1), acceptanceFunctionDay(posInds,2), angleSteps, 'linear', 'extrap');

interpAcceptanceDay = (negInterp+posInterp)/2;

lowInds = find(interpAcceptanceDay > 0.5);
if ~isempty(lowInds)
    % Take last
    lowInds = lowInds(end);
    
    slope = (interpAcceptanceDay(lowInds + 1) - interpAcceptanceDay(lowInds))/...
        (angleSteps(lowInds + 1) - angleSteps(lowInds));
    acceptanceAngleDay = angleSteps(lowInds) + (0.5 - interpAcceptanceDay(lowInds))/slope;
end

% figure; hold on
% plot(angleSteps, interpAcceptanceNight, 'k:', 'linewidth',2)
% plot(acceptanceAngleNight, 0.5, 'kd','markersize',10, 'linewidth',2)
% 
% plot(angleSteps, interpAcceptanceDay, 'k:', 'linewidth',2)
% plot(acceptanceAngleDay, 0.5, 'kd','markersize',10, 'linewidth',2)

%% Create legends reference
legendsF = figure; set(gcf, 'position', [-1704          32        1006         675])

subplot(4,nPlots,1); hold on
h1= plot([0 1], [0 1], 'linewidth',2, 'color', cols(1,:));

h2= plot([0 1], [0 1], 'linewidth',2, 'color', cols(4,:));

h3= plot([0 1], [0 1], 'linewidth',2, 'color', cols(5,:));

h4= plot([0 1], [0 1], 'linewidth',2, 'color', cols(6,:));

legend([h1 h2 h3 h4], {'Radial-top','Radial-base','Linear','Radial+linear'},'Location','south')
set(gca,'LineWidth', 1, 'FontSize', 20);
axis off; ylim([-100 100]); xlim([-100 100]);

subplot(4,nPlots,4); hold on
h1= plot([0 1], [0 1], 'linewidth',2, 'color', cols(1,:));

h2= plot([0 1], [0 1], 'linewidth',2, 'color', cols(4,:));

legend([h1 h2], {'Radial-top','W/ Tip Corr.'},'Location','south')
set(gca,'LineWidth', 1, 'FontSize', 20);
axis off; ylim([-100 100]); xlim([-100 100]);

subplot(4,nPlots,3); hold on
h1 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(2,:));

h2 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(4,:));

h3 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(5,:));

h4 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(6,:));

legend([h1 h2 h3 h4], {'Theory-equal', 'Theory-plus', 'Theory-minus', 'Cylinder'},'Location','south');
set(gca,'LineWidth', 1, 'FontSize', 20);
axis off; ylim([-100 100]); xlim([-100 100]);

subplot(4,nPlots,2); hold on
h1 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(1,:));

h2 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(4,:));

h3 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(5,:));

h4 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(6,:));

legend([h1 h2 h3 h4], {'Radial-top', 'CinC', 'EC', 'CinC & EC'},'Location','south')
set(gca,'LineWidth', 1, 'FontSize', 20);
axis off; ylim([-100 100]); xlim([-100 100]);

subplot(4,nPlots,8); hold on
h1 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(1,:));

h2 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(2,:));

h3 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(3,:));

h4 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(4,:));

h5 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(5,:));

legend([h1 h2 h3 h4 h5], {'Radial-top','Theory-equal', 'SD 0', 'Normalised-max', 'Normalised-80'},'Location','south');
set(gca,'LineWidth', 1, 'FontSize', 20);
axis off; ylim([-100 100]); xlim([-100 100]);

% Uniform RI
subplot(4,nPlots,5); hold on
h1 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(4,:));

h2 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(5,:));

h3 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(3,:));

h4 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(6,:));

h5 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(7,:));

legend([h1 h2 h3 h4 h5], {'SD -4', 'SD -2', 'SD 0', 'SD +2', 'SD +4'},'Location','south')
set(gca,'LineWidth', 1, 'FontSize', 20);
axis off; ylim([-100 100]); xlim([-100 100]);

subplot(4,nPlots,6); hold on
h1 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(3,:));

h2 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(4,:));

h3 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(5,:));

h4 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(6,:));

legend([h1 h2 h3 h4], {'SD 0', 'CinC', 'EC', 'CinC & EC'},'Location','south')
set(gca,'LineWidth', 1, 'FontSize', 20);
axis off; ylim([-100 100]); xlim([-100 100]);

% Varying intercone RI
subplot(4,nPlots,9); hold on
h1 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(1,:));

h2 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(4,:));

h3 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(5,:));

h4 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(6,:));

h5 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(7,:));

legend([h1 h2 h3 h4 h5], {'Radial-top', 'IC-1.42', 'IC-1.44', 'IC-1.46', 'IC-1.48'},'Location','south')
set(gca,'LineWidth', 1, 'FontSize', 20);
axis off; ylim([-100 100]); xlim([-100 100]);

subplot(4,nPlots,10); hold on
h1 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(3,:));

h2 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(4,:));

h3 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(5,:));

h4 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(6,:));

h5 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(7,:));

legend([h1 h2 h3 h4 h5], {'SD 0', 'IC-1.42', 'IC-1.44', 'IC-1.46', 'IC-1.48'},'Location','south')
set(gca, 'LineWidth', 1, 'FontSize', 20);
axis off; ylim([-100 100]); xlim([-100 100]);

subplot(4,nPlots,11); hold on
h1 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(2,:));

h2 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(4,:));

h3 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(5,:));

h4 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(6,:));

h5 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(7,:));

legend([h1 h2 h3 h4 h5], {'Theory-equal', 'IC-1.42', 'IC-1.44', 'IC-1.46', 'IC-1.48'},'Location','south');
set(gca,'LineWidth', 1, 'FontSize', 20);
axis off; ylim([-100 100]); xlim([-100 100]);

% Varying pigment RI
subplot(4,nPlots,13); hold on
h1 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(1,:));

h2 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(4,:));

h3 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(5,:));

h4 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(6,:));

legend([h1 h2 h3 h4], {'Radial-top', 'PG-1.40', 'PG-1.45', 'PG-1.50'},'Location','south')
set(gca, 'LineWidth', 1, 'FontSize', 20);
axis off; ylim([-100 100]); xlim([-100 100]);

subplot(4,nPlots,14); hold on
h1 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(3,:));

h2 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(4,:));

h3 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(5,:));

h4 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(6,:));

legend([h1 h2 h3 h4], {'SD 0', 'PG-1.40', 'PG-1.45', 'PG-1.50'},'Location','south')
ylim([0 1]); xlim([0 20]); ylabel('% TIR'); xlabel('Angle (deg)')
set(gca,'LineWidth', 1, 'FontSize', 20);
axis off; ylim([-100 100]); xlim([-100 100]);

subplot(4,nPlots,15); hold on
h1 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(2,:));

h2 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(4,:));

h3 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(5,:));

h4 = plot([0 1], [0 1], 'linewidth',2, 'color', cols(6,:));

lh4 = legend([h1 h2 h3 h4], {'Theory-equal', 'PG-1.40', 'PG-1.45', 'PG-1.50'},'Location','south');
set(gca,'LineWidth', 1, 'FontSize', 20);
axis off; ylim([-100 100]); xlim([-100 100]);
%% Plot acceptance angle for night
acceptanceFNight = figure; set(gcf, 'position', [-1919        -149        1920        1104])

subplot(4,nPlots,1); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.acceptancePercentageNight, 'linewidth',2, 'color', cols(1,:));
plot(data_GRIN_radialTop.acceptanceAngleNight, 0.5, 'o','markersize',5, 'linewidth',2, 'color', cols(1,:))

h2= plot(data_GRIN_radialBase.incidenceAngle, data_GRIN_radialBase.acceptancePercentageNight, 'linewidth',2, 'color', cols(4,:));
plot(data_GRIN_radialBase.acceptanceAngleNight, 0.5, 'o','markersize',5, 'linewidth',2, 'color', cols(4,:))

h3= plot(data_GRIN_linear.incidenceAngle, data_GRIN_linear.acceptancePercentageNight, 'linewidth',2, 'color', cols(5,:));
plot(data_GRIN_linear.acceptanceAngleNight, 0.5, 'o','markersize',5, 'linewidth',2, 'color', cols(5,:))

h4= plot(data_GRIN_both.incidenceAngle, data_GRIN_both.acceptancePercentageNight, 'linewidth',2, 'color', cols(6,:));
plot(data_GRIN_both.acceptanceAngleNight, 0.5, 'o','markersize',5, 'linewidth',2, 'color', cols(6,:))

plot(angleSteps, interpAcceptanceNight, 'k:', 'linewidth',2)
% plot(acceptanceAngleNight, 0.5, 'kd','markersize',10, 'linewidth',2)

ylim([0 1.25]); xlim([0 20])
ylabel('% entering'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

subplot(4,nPlots,4); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.acceptancePercentageNight, 'linewidth',2, 'color', cols(1,:));
plot(data_GRIN_radialTop.acceptanceAngleNight, 0.5, 'o','markersize',5, 'linewidth',2, 'color', cols(1,:))

h2= plot(data_GRIN_radialTop_tipCorr.incidenceAngle, data_GRIN_radialTop_tipCorr.acceptancePercentageNight, 'linewidth',2, 'color', cols(4,:));
plot(data_GRIN_radialTop_tipCorr.acceptanceAngleNight, 0.5, 'o','markersize',5, 'linewidth',2, 'color', cols(4,:))

plot(angleSteps, interpAcceptanceNight, 'k:', 'linewidth',2)
% plot(acceptanceAngleNight, 0.5, 'kd','markersize',10, 'linewidth',2)

ylim([0 1.25]); xlim([0 20])
ylabel('% entering'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

subplot(4,nPlots,3); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.acceptancePercentageNight, 'linewidth',2, 'color', cols(2,:));
plot(data_GRIN_cylinder.acceptanceAngleNight, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(2,:), 'linewidth',2, 'color', cols(2,:))

h2 = plot(data_GRIN_cylinder_p36.incidenceAngle, data_GRIN_cylinder_p36.acceptancePercentageNight, 'linewidth',2, 'color', cols(4,:));
plot(data_GRIN_cylinder_p36.acceptanceAngleNight, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(4,:),'linewidth',2, 'color', cols(4,:))

h3 = plot(data_GRIN_cylinder_m36.incidenceAngle, data_GRIN_cylinder_m36.acceptancePercentageNight, 'linewidth',2, 'color', cols(5,:));
plot(data_GRIN_cylinder_m36.acceptanceAngleNight, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(5,:), 'linewidth',2, 'color', cols(5,:))

h4 = plot(data_GRIN_cylinder_shape.incidenceAngle, data_GRIN_cylinder_shape.acceptancePercentageNight, 'linewidth',2, 'color', cols(6,:));
plot(data_GRIN_cylinder_shape.acceptanceAngleNight, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(6,:), 'linewidth',2, 'color', cols(6,:))

plot(angleSteps, interpAcceptanceNight, 'k:', 'linewidth',2)
% plot(acceptanceAngleNight, 0.5, 'kd','markersize',10, 'linewidth',2)

ylim([0 1.25]); xlim([0 20])
ylabel('% entering'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

subplot(4,nPlots,2); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.acceptancePercentageNight, 'linewidth',2, 'color', cols(1,:));
plot(data_GRIN_radialTop.acceptanceAngleNight, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(1,:), 'linewidth',2, 'color', cols(1,:))

h2 = plot(data_GRIN_radialTop_CinC.incidenceAngle, data_GRIN_radialTop_CinC.acceptancePercentageNight, 'linewidth',2, 'color', cols(4,:));
plot(data_GRIN_radialTop_CinC.acceptanceAngleNight, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(4,:), 'linewidth',2, 'color', cols(4,:))

h3 = plot(data_GRIN_radialTop_EC.incidenceAngle, data_GRIN_radialTop_EC.acceptancePercentageNight, 'linewidth',2, 'color', cols(5,:));
plot(data_GRIN_radialTop_EC.acceptanceAngleNight, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(5,:), 'linewidth',2, 'color', cols(5,:))

h4 = plot(data_GRIN_radialTop_CinC_EC.incidenceAngle, data_GRIN_radialTop_CinC_EC.acceptancePercentageNight, 'linewidth',2, 'color', cols(6,:));
plot(data_GRIN_radialTop_CinC_EC.acceptanceAngleNight, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(6,:), 'linewidth',2, 'color', cols(6,:))

plot(angleSteps, interpAcceptanceNight, 'k:', 'linewidth',2)
% plot(acceptanceAngleNight, 0.5, 'kd','markersize',10, 'linewidth',2)

ylim([0 1.25]); xlim([0 20])
ylabel('% entering'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

subplot(4,nPlots,8); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.acceptancePercentageNight, 'linewidth',2, 'color', cols(1,:));
plot(data_GRIN_radialTop.acceptanceAngleNight, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(1,:), 'linewidth',2, 'color', cols(1,:))

h2 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.acceptancePercentageNight, 'linewidth',2, 'color', cols(2,:));
plot(data_GRIN_cylinder.acceptanceAngleNight, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(2,:),'linewidth',2, 'color', cols(2,:))

h3 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.acceptancePercentageNight, 'linewidth',2, 'color', cols(3,:));
plot(data_uniform_0SD.acceptanceAngleNight, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(3,:),'linewidth',2, 'color', cols(3,:))

h4 = plot(data_GRIN_cylinder_norm_max.incidenceAngle, data_GRIN_cylinder_norm_max.acceptancePercentageNight, 'linewidth',2, 'color', cols(4,:));
plot(data_GRIN_cylinder_norm_max.acceptanceAngleNight, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(4,:), 'linewidth',2, 'color', cols(4,:))

h5 = plot(data_GRIN_cylinder_norm_80.incidenceAngle, data_GRIN_cylinder_norm_80.acceptancePercentageNight, 'linewidth',2, 'color', cols(5,:));
plot(data_GRIN_cylinder_norm_80.acceptanceAngleNight, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(5,:), 'linewidth',2, 'color', cols(5,:))

plot(angleSteps, interpAcceptanceNight, 'k:', 'linewidth',2)
% plot(acceptanceAngleNight, 0.5, 'kd','markersize',10, 'linewidth',2)

ylim([0 1.25]); xlim([0 20])
ylabel('% entering'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

%uniform RI
subplot(4,nPlots,5); hold on
h1 = plot(data_uniform_m4SD.incidenceAngle, data_uniform_m4SD.acceptancePercentageNight, 'linewidth',2, 'color', cols(4,:));
plot(data_uniform_m4SD.acceptanceAngleNight, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(4,:), 'linewidth',2, 'color', cols(4,:))

h2 = plot(data_uniform_m2SD.incidenceAngle, data_uniform_m2SD.acceptancePercentageNight, 'linewidth',2, 'color', cols(5,:));
plot(data_uniform_m2SD.acceptanceAngleNight, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(5,:), 'linewidth',2, 'color', cols(5,:))

h3 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.acceptancePercentageNight, 'linewidth',2, 'color', cols(3,:));
plot(data_uniform_0SD.acceptanceAngleNight, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(3,:), 'linewidth',2, 'color', cols(3,:))

h4 = plot(data_uniform_p2SD.incidenceAngle, data_uniform_p2SD.acceptancePercentageNight, 'linewidth',2, 'color', cols(6,:));
plot(data_uniform_p2SD.acceptanceAngleNight, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(6,:), 'linewidth',2, 'color', cols(6,:))

h5 = plot(data_uniform_p4SD.incidenceAngle, data_uniform_p4SD.acceptancePercentageNight, 'linewidth',2, 'color', cols(7,:));
plot(data_uniform_p4SD.acceptanceAngleNight, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(7,:), 'linewidth',2, 'color', cols(7,:))

plot(angleSteps, interpAcceptanceNight, 'k:', 'linewidth',2)
% plot(acceptanceAngleNight, 0.5, 'kd','markersize',10, 'linewidth',2)

ylim([0 1.25]); xlim([0 20])
ylabel('% entering'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

subplot(4,nPlots,6); hold on
h1 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.acceptancePercentageNight, 'linewidth',2, 'color', cols(3,:));
plot(data_uniform_0SD.acceptanceAngleNight, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(3,:), 'linewidth',2, 'color', cols(3,:))

h2 = plot(data_GRIN_uniform_CinC.incidenceAngle, data_GRIN_uniform_CinC.acceptancePercentageNight, 'linewidth',2, 'color', cols(4,:));
plot(data_GRIN_uniform_CinC.acceptanceAngleNight, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(4,:), 'linewidth',2, 'color', cols(4,:))

h3 = plot(data_GRIN_uniform_EC.incidenceAngle, data_GRIN_uniform_EC.acceptancePercentageNight, 'linewidth',2, 'color', cols(5,:));
plot(data_GRIN_uniform_EC.acceptanceAngleNight, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(5,:), 'linewidth',2, 'color', cols(5,:))

h4 = plot(data_GRIN_uniform_CinC_EC.incidenceAngle, data_GRIN_uniform_CinC_EC.acceptancePercentageNight, 'linewidth',2, 'color', cols(6,:));
plot(data_GRIN_uniform_CinC_EC.acceptanceAngleNight, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(6,:), 'linewidth',2, 'color', cols(6,:))

plot(angleSteps, interpAcceptanceNight, 'k:', 'linewidth',2)
% plot(acceptanceAngleNight, 0.5, 'kd','markersize',10, 'linewidth',2)

ylim([0 1.25]); xlim([0 20])
ylabel('% entering'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

% Varying intercone RI
subplot(4,nPlots,9); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.acceptancePercentageNight, 'linewidth',2, 'color', cols(1,:));
plot(data_GRIN_radialTop.acceptanceAngleNight, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(1,:), 'linewidth',2, 'color', cols(1,:))

h2 = plot(data_GRIN_radialTop_IC_142.incidenceAngle, data_GRIN_radialTop_IC_142.acceptancePercentageNight, 'linewidth',2, 'color', cols(4,:));
plot(data_GRIN_radialTop_IC_142.acceptanceAngleNight, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(4,:), 'linewidth',2, 'color', cols(4,:))

h3 = plot(data_GRIN_radialTop_IC_144.incidenceAngle, data_GRIN_radialTop_IC_144.acceptancePercentageNight, 'linewidth',2, 'color', cols(5,:));
plot(data_GRIN_radialTop_IC_144.acceptanceAngleNight, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(5,:), 'linewidth',2, 'color', cols(5,:))

h4 = plot(data_GRIN_radialTop_IC_146.incidenceAngle, data_GRIN_radialTop_IC_146.acceptancePercentageNight, 'linewidth',2, 'color', cols(6,:));
plot(data_GRIN_radialTop_IC_146.acceptanceAngleNight, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(6,:), 'linewidth',2, 'color', cols(6,:))

h5 = plot(data_GRIN_radialTop_IC_148.incidenceAngle, data_GRIN_radialTop_IC_148.acceptancePercentageNight, 'linewidth',2, 'color', cols(7,:));
plot(data_GRIN_radialTop_IC_148.acceptanceAngleNight, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(7,:), 'linewidth',2, 'color', cols(7,:))

plot(angleSteps, interpAcceptanceNight, 'k:', 'linewidth',2)
% plot(acceptanceAngleNight, 0.5, 'kd','markersize',10, 'linewidth',2)

ylim([0 1.25]); xlim([0 20])
ylabel('% entering'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

subplot(4,nPlots,10); hold on
h1 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.acceptancePercentageNight, 'linewidth',2, 'color', cols(3,:));
plot(data_uniform_0SD.acceptanceAngleNight, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(3,:), 'linewidth',2, 'color', cols(3,:))

h2 = plot(data_uniform_0SD_IC_142.incidenceAngle, data_uniform_0SD_IC_142.acceptancePercentageNight, 'linewidth',2, 'color', cols(4,:));
plot(data_uniform_0SD_IC_142.acceptanceAngleNight, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(4,:), 'linewidth',2, 'color', cols(4,:))

h3 = plot(data_uniform_0SD_IC_144.incidenceAngle, data_uniform_0SD_IC_144.acceptancePercentageNight, 'linewidth',2, 'color', cols(5,:));
plot(data_uniform_0SD_IC_144.acceptanceAngleNight, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(5,:), 'linewidth',2, 'color', cols(5,:))

h4 = plot(data_uniform_0SD_IC_146.incidenceAngle, data_uniform_0SD_IC_146.acceptancePercentageNight, 'linewidth',2, 'color', cols(6,:));
plot(data_uniform_0SD_IC_146.acceptanceAngleNight, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(6,:), 'linewidth',2, 'color', cols(6,:))

h5 = plot(data_uniform_0SD_IC_148.incidenceAngle, data_uniform_0SD_IC_148.acceptancePercentageNight, 'linewidth',2, 'color', cols(7,:));
plot(data_uniform_0SD_IC_148.acceptanceAngleNight, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(7,:), 'linewidth',2, 'color', cols(7,:))

plot(angleSteps, interpAcceptanceNight, 'k:', 'linewidth',2)
% plot(acceptanceAngleNight, 0.5, 'kd','markersize',10, 'linewidth',2)

ylim([0 1.25]); xlim([0 20])
ylabel('% entering'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

subplot(4,nPlots,11); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.acceptancePercentageNight, 'linewidth',2, 'color', cols(2,:));
plot(data_GRIN_cylinder.acceptanceAngleNight, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(2,:), 'linewidth',2, 'color', cols(2,:))

h2 = plot(data_GRIN_cylinder_IC_142.incidenceAngle, data_GRIN_cylinder_IC_142.acceptancePercentageNight, 'linewidth',2, 'color', cols(4,:));
plot(data_GRIN_cylinder_IC_142.acceptanceAngleNight, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(4,:), 'linewidth',2, 'color', cols(4,:))

h3 = plot(data_GRIN_cylinder_IC_144.incidenceAngle, data_GRIN_cylinder_IC_144.acceptancePercentageNight, 'linewidth',2, 'color', cols(5,:));
plot(data_GRIN_cylinder_IC_144.acceptanceAngleNight, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(5,:), 'linewidth',2, 'color', cols(5,:))

h4 = plot(data_GRIN_cylinder_IC_146.incidenceAngle, data_GRIN_cylinder_IC_146.acceptancePercentageNight, 'linewidth',2, 'color', cols(6,:));
plot(data_GRIN_cylinder_IC_146.acceptanceAngleNight, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(6,:), 'linewidth',2, 'color', cols(6,:))

h5 = plot(data_GRIN_cylinder_IC_148.incidenceAngle, data_GRIN_cylinder_IC_148.acceptancePercentageNight, 'linewidth',2, 'color', cols(7,:));
plot(data_GRIN_cylinder_IC_148.acceptanceAngleNight, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(7,:), 'linewidth',2, 'color', cols(7,:))

plot(angleSteps, interpAcceptanceNight, 'k:', 'linewidth',2)
% plot(acceptanceAngleNight, 0.5, 'kd','markersize',10, 'linewidth',2)

ylim([0 1.25]); xlim([0 20])
ylabel('% entering'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

% Varying pigment RI
subplot(4,nPlots,13); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.acceptancePercentageNight, 'linewidth',2, 'color', cols(1,:));
plot(data_GRIN_radialTop.acceptanceAngleNight, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(1,:),'linewidth',2, 'color', cols(1,:))

h2 = plot(data_GRIN_radialTop_PG_140.incidenceAngle, data_GRIN_radialTop_PG_140.acceptancePercentageNight, 'linewidth',2, 'color', cols(4,:));
plot(data_GRIN_radialTop_PG_140.acceptanceAngleNight, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(4,:),'linewidth',2, 'color', cols(4,:))

h3 = plot(data_GRIN_radialTop_PG_145.incidenceAngle, data_GRIN_radialTop_PG_145.acceptancePercentageNight, 'linewidth',2, 'color', cols(5,:));
plot(data_GRIN_radialTop_PG_145.acceptanceAngleNight, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(5,:), 'linewidth',2, 'color', cols(5,:))

h4 = plot(data_GRIN_radialTop_PG_150.incidenceAngle, data_GRIN_radialTop_PG_150.acceptancePercentageNight, 'linewidth',2, 'color', cols(6,:));
plot(data_GRIN_radialTop_PG_150.acceptanceAngleNight, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(6,:),'linewidth',2, 'color', cols(6,:))

plot(angleSteps, interpAcceptanceNight, 'k:', 'linewidth',2)
% plot(acceptanceAngleNight, 0.5, 'kd','markersize',10, 'linewidth',2)

ylim([0 1.25]); xlim([0 20])
ylabel('% entering'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

subplot(4,nPlots,14); hold on
h1 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.acceptancePercentageNight, 'linewidth',2, 'color', cols(3,:));
plot(data_uniform_0SD.acceptanceAngleNight, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(3,:),'linewidth',2, 'color', cols(3,:))

h2 = plot(data_uniform_0SD_PG_140.incidenceAngle, data_uniform_0SD_PG_140.acceptancePercentageNight, 'linewidth',2, 'color', cols(4,:));
plot(data_uniform_0SD_PG_140.acceptanceAngleNight, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(4,:),'linewidth',2, 'color', cols(4,:))

h3 = plot(data_uniform_0SD_PG_145.incidenceAngle, data_uniform_0SD_PG_145.acceptancePercentageNight, 'linewidth',2, 'color', cols(5,:));
plot(data_uniform_0SD_PG_145.acceptanceAngleNight, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(5,:),'linewidth',2, 'color', cols(5,:))

h4 = plot(data_uniform_0SD_PG_150.incidenceAngle, data_uniform_0SD_PG_150.acceptancePercentageNight, 'linewidth',2, 'color', cols(6,:));
plot(data_uniform_0SD_PG_150.acceptanceAngleNight, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(6,:),'linewidth',2, 'color', cols(6,:))

plot(angleSteps, interpAcceptanceNight, 'k:', 'linewidth',2)
% plot(acceptanceAngleNight, 0.5, 'kd','markersize',10, 'linewidth',2)

ylim([0 1.25]); xlim([0 20])
ylabel('% entering'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

subplot(4,nPlots,15); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.acceptancePercentageNight, 'linewidth',2, 'color', cols(2,:));
plot(data_GRIN_cylinder.acceptanceAngleNight, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(2,:),'linewidth',2, 'color', cols(2,:))

h2 = plot(data_GRIN_cylinder_PG_140.incidenceAngle, data_GRIN_cylinder_PG_140.acceptancePercentageNight, 'linewidth',2, 'color', cols(4,:));
plot(data_GRIN_cylinder_PG_140.acceptanceAngleNight, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(4,:),'linewidth',2, 'color', cols(4,:))

h3 = plot(data_GRIN_cylinder_PG_145.incidenceAngle, data_GRIN_cylinder_PG_145.acceptancePercentageNight, 'linewidth',2, 'color', cols(5,:));
plot(data_GRIN_cylinder_PG_145.acceptanceAngleNight, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(5,:),'linewidth',2, 'color', cols(5,:))

h4 = plot(data_GRIN_cylinder_PG_150.incidenceAngle, data_GRIN_cylinder_PG_150.acceptancePercentageNight, 'linewidth',2, 'color', cols(6,:));
plot(data_GRIN_cylinder_PG_150.acceptanceAngleNight, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(6,:),'linewidth',2, 'color', cols(6,:))

plot(angleSteps, interpAcceptanceNight, 'k:', 'linewidth',2)
% plot(acceptanceAngleNight, 0.5, 'kd','markersize',10, 'linewidth',2)

ylim([0 1.25]); xlim([0 20])
ylabel('% entering'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

%% Plot acceptance angle for day
acceptanceFDay = figure; set(gcf, 'position', [-1919        -149        1920        1104])

subplot(4,nPlots,1); hold on
if data_GRIN_radialTop.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h1= plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.acceptancePercentageDay, 'linewidth',2, 'color', cols(1,:), 'linestyle', style);
plot(data_GRIN_radialTop.acceptanceAngleDay, 0.5, 'o','markersize',5, 'linewidth',2, 'color', cols(1,:))

if data_GRIN_radialBase.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h2= plot(data_GRIN_radialBase.incidenceAngle, data_GRIN_radialBase.acceptancePercentageDay, 'linewidth',2, 'color', cols(4,:), 'linestyle', style);
plot(data_GRIN_radialBase.acceptanceAngleDay, 0.5, 'o','markersize',5, 'linewidth',2, 'color', cols(4,:))

if data_GRIN_linear.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h3= plot(data_GRIN_linear.incidenceAngle, data_GRIN_linear.acceptancePercentageDay, 'linewidth',2, 'color', cols(5,:), 'linestyle', style);
plot(data_GRIN_linear.acceptanceAngleDay, 0.5, 'o','markersize',5, 'linewidth',2, 'color', cols(5,:))

if data_GRIN_both.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h4= plot(data_GRIN_both.incidenceAngle, data_GRIN_both.acceptancePercentageDay, 'linewidth',2, 'color', cols(6,:), 'linestyle', style);
plot(data_GRIN_both.acceptanceAngleDay, 0.5, 'o','markersize',5, 'linewidth',2, 'color', cols(6,:))

plot(angleSteps, interpAcceptanceDay, 'k:', 'linewidth',2)
% plot(acceptanceAngleDay, 0.5, 'kd','markersize',10, 'linewidth',2)

ylim([0 1.25]); xlim([0 20])
ylabel('% entering'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

subplot(4,nPlots,4); hold on
if data_GRIN_radialTop.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h1= plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.acceptancePercentageDay, 'linewidth',2, 'color', cols(1,:), 'linestyle', style);
plot(data_GRIN_radialTop.acceptanceAngleDay, 0.5, 'o','markersize',5, 'linewidth',2, 'color', cols(1,:))

if data_GRIN_radialTop_tipCorr.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h2= plot(data_GRIN_radialTop_tipCorr.incidenceAngle, data_GRIN_radialTop_tipCorr.acceptancePercentageDay, 'linewidth',2, 'color', cols(4,:), 'linestyle', style);
plot(data_GRIN_radialTop_tipCorr.acceptanceAngleDay, 0.5, 'o','markersize',5, 'linewidth',2, 'color', cols(4,:))

plot(angleSteps, interpAcceptanceDay, 'k:', 'linewidth',2)
% plot(acceptanceAngleDay, 0.5, 'kd','markersize',10, 'linewidth',2)

ylim([0 1.25]); xlim([0 20])
ylabel('% entering'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

subplot(4,nPlots,3); hold on
if data_GRIN_cylinder.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h1 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.acceptancePercentageDay, 'linewidth',2, 'color', cols(2,:), 'linestyle', style);
plot(data_GRIN_cylinder.acceptanceAngleDay, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(2,:), 'linewidth',2, 'color', cols(2,:))

if data_GRIN_cylinder_p36.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h2 = plot(data_GRIN_cylinder_p36.incidenceAngle, data_GRIN_cylinder_p36.acceptancePercentageDay, 'linewidth',2, 'color', cols(4,:), 'linestyle', style);
plot(data_GRIN_cylinder_p36.acceptanceAngleDay, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(4,:),'linewidth',2, 'color', cols(4,:))

if data_GRIN_cylinder_m36.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h3 = plot(data_GRIN_cylinder_m36.incidenceAngle, data_GRIN_cylinder_m36.acceptancePercentageDay, 'linewidth',2, 'color', cols(5,:), 'linestyle', style);
plot(data_GRIN_cylinder_m36.acceptanceAngleDay, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(5,:), 'linewidth',2, 'color', cols(5,:))

if data_GRIN_cylinder_shape.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h4 = plot(data_GRIN_cylinder_shape.incidenceAngle, data_GRIN_cylinder_shape.acceptancePercentageDay, 'linewidth',2, 'color', cols(6,:), 'linestyle', style);
plot(data_GRIN_cylinder_shape.acceptanceAngleDay, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(6,:), 'linewidth',2, 'color', cols(6,:))

plot(angleSteps, interpAcceptanceDay, 'k:', 'linewidth',2)
% plot(acceptanceAngleDay, 0.5, 'kd','markersize',10, 'linewidth',2)

ylim([0 1.25]); xlim([0 20])
ylabel('% entering'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

subplot(4,nPlots,2); hold on
if data_GRIN_radialTop.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.acceptancePercentageDay, 'linewidth',2, 'color', cols(1,:), 'linestyle', style);
plot(data_GRIN_radialTop.acceptanceAngleDay, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(1,:), 'linewidth',2, 'color', cols(1,:))

if data_GRIN_radialTop_CinC.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h2 = plot(data_GRIN_radialTop_CinC.incidenceAngle, data_GRIN_radialTop_CinC.acceptancePercentageDay, 'linewidth',2, 'color', cols(4,:), 'linestyle', style);
plot(data_GRIN_radialTop_CinC.acceptanceAngleDay, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(4,:), 'linewidth',2, 'color', cols(4,:))

if data_GRIN_radialTop_EC.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h3 = plot(data_GRIN_radialTop_EC.incidenceAngle, data_GRIN_radialTop_EC.acceptancePercentageDay, 'linewidth',2, 'color', cols(5,:), 'linestyle', style);
plot(data_GRIN_radialTop_EC.acceptanceAngleDay, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(5,:), 'linewidth',2, 'color', cols(5,:))

if data_GRIN_radialTop_CinC_EC.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h4 = plot(data_GRIN_radialTop_CinC_EC.incidenceAngle, data_GRIN_radialTop_CinC_EC.acceptancePercentageDay, 'linewidth',2, 'color', cols(6,:), 'linestyle', style);
plot(data_GRIN_radialTop_CinC_EC.acceptanceAngleDay, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(6,:), 'linewidth',2, 'color', cols(6,:))

plot(angleSteps, interpAcceptanceDay, 'k:', 'linewidth',2)
% plot(acceptanceAngleDay, 0.5, 'kd','markersize',10, 'linewidth',2)

ylim([0 1.25]); xlim([0 20])
ylabel('% entering'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

subplot(4,nPlots,8); hold on
if data_GRIN_radialTop.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h1= plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.acceptancePercentageDay, 'linewidth',2, 'color', cols(1,:), 'linestyle', style);
plot(data_GRIN_radialTop.acceptanceAngleDay, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(1,:), 'linewidth',2, 'color', cols(1,:))

if data_GRIN_cylinder.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h2 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.acceptancePercentageDay, 'linewidth',2, 'color', cols(2,:), 'linestyle', style);
plot(data_GRIN_cylinder.acceptanceAngleDay, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(2,:),'linewidth',2, 'color', cols(2,:))

if data_uniform_0SD.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h3 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.acceptancePercentageDay, 'linewidth',2, 'color', cols(3,:), 'linestyle', style);
plot(data_uniform_0SD.acceptanceAngleDay, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(3,:),'linewidth',2, 'color', cols(3,:))

if data_GRIN_cylinder_norm_max.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h4 = plot(data_GRIN_cylinder_norm_max.incidenceAngle, data_GRIN_cylinder_norm_max.acceptancePercentageDay, 'linewidth',2, 'color', cols(4,:), 'linestyle', style);
plot(data_GRIN_cylinder_norm_max.acceptanceAngleDay, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(4,:), 'linewidth',2, 'color', cols(4,:))

if data_GRIN_cylinder_norm_80.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h5 = plot(data_GRIN_cylinder_norm_80.incidenceAngle, data_GRIN_cylinder_norm_80.acceptancePercentageDay, 'linewidth',2, 'color', cols(5,:), 'linestyle', style);
plot(data_GRIN_cylinder_norm_80.acceptanceAngleDay, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(5,:), 'linewidth',2, 'color', cols(5,:))

plot(angleSteps, interpAcceptanceDay, 'k:', 'linewidth',2)
% plot(acceptanceAngleDay, 0.5, 'kd','markersize',10, 'linewidth',2)

ylim([0 1.25]); xlim([0 20])
ylabel('% entering'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

%uniform RI
subplot(4,nPlots,5); hold on
if data_uniform_m4SD.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h1 = plot(data_uniform_m4SD.incidenceAngle, data_uniform_m4SD.acceptancePercentageDay, 'linewidth',2, 'color', cols(4,:), 'linestyle', style);
plot(data_uniform_m4SD.acceptanceAngleDay, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(4,:), 'linewidth',2, 'color', cols(4,:))

if data_uniform_m2SD.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h2 = plot(data_uniform_m2SD.incidenceAngle, data_uniform_m2SD.acceptancePercentageDay, 'linewidth',2, 'color', cols(5,:), 'linestyle', style);
plot(data_uniform_m2SD.acceptanceAngleDay, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(5,:), 'linewidth',2, 'color', cols(5,:))

if data_uniform_0SD.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h3 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.acceptancePercentageDay, 'linewidth',2, 'color', cols(3,:), 'linestyle', style);
plot(data_uniform_0SD.acceptanceAngleDay, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(3,:), 'linewidth',2, 'color', cols(3,:))

if data_uniform_p2SD.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h4 = plot(data_uniform_p2SD.incidenceAngle, data_uniform_p2SD.acceptancePercentageDay, 'linewidth',2, 'color', cols(6,:), 'linestyle', style);
plot(data_uniform_p2SD.acceptanceAngleDay, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(6,:), 'linewidth',2, 'color', cols(6,:))

if data_uniform_p4SD.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h5 = plot(data_uniform_p4SD.incidenceAngle, data_uniform_p4SD.acceptancePercentageDay, 'linewidth',2, 'color', cols(7,:), 'linestyle', style);
plot(data_uniform_p4SD.acceptanceAngleDay, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(7,:), 'linewidth',2, 'color', cols(7,:))

plot(angleSteps, interpAcceptanceDay, 'k:', 'linewidth',2)
% plot(acceptanceAngleDay, 0.5, 'kd','markersize',10, 'linewidth',2)

ylim([0 1.25]); xlim([0 20])
ylabel('% entering'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

subplot(4,nPlots,6); hold on
if data_uniform_0SD.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h1 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.acceptancePercentageDay, 'linewidth',2, 'color', cols(3,:), 'linestyle', style);
plot(data_uniform_0SD.acceptanceAngleDay, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(3,:), 'linewidth',2, 'color', cols(3,:))

if data_GRIN_uniform_CinC.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h2 = plot(data_GRIN_uniform_CinC.incidenceAngle, data_GRIN_uniform_CinC.acceptancePercentageDay, 'linewidth',2, 'color', cols(4,:), 'linestyle', style);
plot(data_GRIN_uniform_CinC.acceptanceAngleDay, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(4,:), 'linewidth',2, 'color', cols(4,:))

if data_GRIN_uniform_EC.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h3 = plot(data_GRIN_uniform_EC.incidenceAngle, data_GRIN_uniform_EC.acceptancePercentageDay, 'linewidth',2, 'color', cols(5,:), 'linestyle', style);
plot(data_GRIN_uniform_EC.acceptanceAngleDay, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(5,:), 'linewidth',2, 'color', cols(5,:))

if data_GRIN_uniform_CinC_EC.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h4 = plot(data_GRIN_uniform_CinC_EC.incidenceAngle, data_GRIN_uniform_CinC_EC.acceptancePercentageDay, 'linewidth',2, 'color', cols(6,:), 'linestyle', style);
plot(data_GRIN_uniform_CinC_EC.acceptanceAngleDay, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(6,:), 'linewidth',2, 'color', cols(6,:))

plot(angleSteps, interpAcceptanceDay, 'k:', 'linewidth',2)
% plot(acceptanceAngleDay, 0.5, 'kd','markersize',10, 'linewidth',2)

ylim([0 1.25]); xlim([0 20])
ylabel('% entering'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

% Varying intercone RI
subplot(4,nPlots,9); hold on
if data_GRIN_radialTop.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.acceptancePercentageDay, 'linewidth',2, 'color', cols(1,:), 'linestyle', style);
plot(data_GRIN_radialTop.acceptanceAngleDay, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(1,:), 'linewidth',2, 'color', cols(1,:))

if data_GRIN_radialTop_IC_142.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h2 = plot(data_GRIN_radialTop_IC_142.incidenceAngle, data_GRIN_radialTop_IC_142.acceptancePercentageDay, 'linewidth',2, 'color', cols(4,:), 'linestyle', style);
plot(data_GRIN_radialTop_IC_142.acceptanceAngleDay, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(4,:), 'linewidth',2, 'color', cols(4,:))

if data_GRIN_radialTop_IC_144.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h3 = plot(data_GRIN_radialTop_IC_144.incidenceAngle, data_GRIN_radialTop_IC_144.acceptancePercentageDay, 'linewidth',2, 'color', cols(5,:), 'linestyle', style);
plot(data_GRIN_radialTop_IC_144.acceptanceAngleDay, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(5,:), 'linewidth',2, 'color', cols(5,:))

if data_GRIN_radialTop_IC_146.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h4 = plot(data_GRIN_radialTop_IC_146.incidenceAngle, data_GRIN_radialTop_IC_146.acceptancePercentageDay, 'linewidth',2, 'color', cols(6,:), 'linestyle', style);
plot(data_GRIN_radialTop_IC_146.acceptanceAngleDay, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(6,:), 'linewidth',2, 'color', cols(6,:))

if data_GRIN_radialTop_IC_148.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h5 = plot(data_GRIN_radialTop_IC_148.incidenceAngle, data_GRIN_radialTop_IC_148.acceptancePercentageDay, 'linewidth',2, 'color', cols(7,:), 'linestyle', style);
plot(data_GRIN_radialTop_IC_148.acceptanceAngleDay, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(7,:), 'linewidth',2, 'color', cols(7,:))

plot(angleSteps, interpAcceptanceDay, 'k:', 'linewidth',2)
% plot(acceptanceAngleDay, 0.5, 'kd','markersize',10, 'linewidth',2)

ylim([0 1.25]); xlim([0 20])
ylabel('% entering'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

subplot(4,nPlots,10); hold on
if data_uniform_0SD.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h1 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.acceptancePercentageDay, 'linewidth',2, 'color', cols(3,:), 'linestyle', style);
plot(data_uniform_0SD.acceptanceAngleDay, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(3,:), 'linewidth',2, 'color', cols(3,:))

if data_uniform_0SD_IC_142.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h2 = plot(data_uniform_0SD_IC_142.incidenceAngle, data_uniform_0SD_IC_142.acceptancePercentageDay, 'linewidth',2, 'color', cols(4,:), 'linestyle', style);
plot(data_uniform_0SD_IC_142.acceptanceAngleDay, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(4,:), 'linewidth',2, 'color', cols(4,:))

if data_uniform_0SD_IC_144.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h3 = plot(data_uniform_0SD_IC_144.incidenceAngle, data_uniform_0SD_IC_144.acceptancePercentageDay, 'linewidth',2, 'color', cols(5,:), 'linestyle', style);
plot(data_uniform_0SD_IC_144.acceptanceAngleDay, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(5,:), 'linewidth',2, 'color', cols(5,:))

if data_uniform_0SD_IC_146.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h4 = plot(data_uniform_0SD_IC_146.incidenceAngle, data_uniform_0SD_IC_146.acceptancePercentageDay, 'linewidth',2, 'color', cols(6,:), 'linestyle', style);
plot(data_uniform_0SD_IC_146.acceptanceAngleDay, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(6,:), 'linewidth',2, 'color', cols(6,:))

if data_uniform_0SD_IC_148.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h5 = plot(data_uniform_0SD_IC_148.incidenceAngle, data_uniform_0SD_IC_148.acceptancePercentageDay, 'linewidth',2, 'color', cols(7,:), 'linestyle', style);
plot(data_uniform_0SD_IC_148.acceptanceAngleDay, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(7,:), 'linewidth',2, 'color', cols(7,:))

plot(angleSteps, interpAcceptanceDay, 'k:', 'linewidth',2)
% plot(acceptanceAngleDay, 0.5, 'kd','markersize',10, 'linewidth',2)

ylim([0 1.25]); xlim([0 20])
ylabel('% entering'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

subplot(4,nPlots,11); hold on
if data_GRIN_cylinder.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h1 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.acceptancePercentageDay, 'linewidth',2, 'color', cols(2,:), 'linestyle', style);
plot(data_GRIN_cylinder.acceptanceAngleDay, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(2,:), 'linewidth',2, 'color', cols(2,:))

if data_GRIN_cylinder_IC_142.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h2 = plot(data_GRIN_cylinder_IC_142.incidenceAngle, data_GRIN_cylinder_IC_142.acceptancePercentageDay, 'linewidth',2, 'color', cols(4,:), 'linestyle', style);
plot(data_GRIN_cylinder_IC_142.acceptanceAngleDay, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(4,:), 'linewidth',2, 'color', cols(4,:))

if data_GRIN_cylinder_IC_144.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h3 = plot(data_GRIN_cylinder_IC_144.incidenceAngle, data_GRIN_cylinder_IC_144.acceptancePercentageDay, 'linewidth',2, 'color', cols(5,:), 'linestyle', style);
plot(data_GRIN_cylinder_IC_144.acceptanceAngleDay, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(5,:), 'linewidth',2, 'color', cols(5,:))

if data_GRIN_cylinder_IC_146.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h4 = plot(data_GRIN_cylinder_IC_146.incidenceAngle, data_GRIN_cylinder_IC_146.acceptancePercentageDay, 'linewidth',2, 'color', cols(6,:), 'linestyle', style);
plot(data_GRIN_cylinder_IC_146.acceptanceAngleDay, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(6,:), 'linewidth',2, 'color', cols(6,:))

if data_GRIN_cylinder_IC_148.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h5 = plot(data_GRIN_cylinder_IC_148.incidenceAngle, data_GRIN_cylinder_IC_148.acceptancePercentageDay, 'linewidth',2, 'color', cols(7,:), 'linestyle', style);
plot(data_GRIN_cylinder_IC_148.acceptanceAngleDay, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(7,:), 'linewidth',2, 'color', cols(7,:))

plot(angleSteps, interpAcceptanceDay, 'k:', 'linewidth',2)
% plot(acceptanceAngleDay, 0.5, 'kd','markersize',10, 'linewidth',2)

ylim([0 1.25]); xlim([0 20])
ylabel('% entering'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

% Varying pigment RI
subplot(4,nPlots,13); hold on
if data_GRIN_radialTop.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.acceptancePercentageDay, 'linewidth',2, 'color', cols(1,:), 'linestyle', style);
plot(data_GRIN_radialTop.acceptanceAngleDay, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(1,:),'linewidth',2, 'color', cols(1,:))

if data_GRIN_radialTop_PG_140.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h2 = plot(data_GRIN_radialTop_PG_140.incidenceAngle, data_GRIN_radialTop_PG_140.acceptancePercentageDay, 'linewidth',2, 'color', cols(4,:), 'linestyle', style);
plot(data_GRIN_radialTop_PG_140.acceptanceAngleDay, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(4,:),'linewidth',2, 'color', cols(4,:))

if data_GRIN_radialTop_PG_145.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h3 = plot(data_GRIN_radialTop_PG_145.incidenceAngle, data_GRIN_radialTop_PG_145.acceptancePercentageDay, 'linewidth',2, 'color', cols(5,:), 'linestyle', style);
plot(data_GRIN_radialTop_PG_145.acceptanceAngleDay, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(5,:), 'linewidth',2, 'color', cols(5,:))

if data_GRIN_radialTop_PG_150.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h4 = plot(data_GRIN_radialTop_PG_150.incidenceAngle, data_GRIN_radialTop_PG_150.acceptancePercentageDay, 'linewidth',2, 'color', cols(6,:), 'linestyle', style);
plot(data_GRIN_radialTop_PG_150.acceptanceAngleDay, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(6,:),'linewidth',2, 'color', cols(6,:))

plot(angleSteps, interpAcceptanceDay, 'k:', 'linewidth',2)
% plot(acceptanceAngleDay, 0.5, 'kd','markersize',10, 'linewidth',2)

ylim([0 1.25]); xlim([0 20])
ylabel('% entering'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

subplot(4,nPlots,14); hold on
if data_uniform_0SD.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h1 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.acceptancePercentageDay, 'linewidth',2, 'color', cols(3,:), 'linestyle', style);
plot(data_uniform_0SD.acceptanceAngleDay, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(3,:),'linewidth',2, 'color', cols(3,:))

if data_uniform_0SD_PG_140.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h2 = plot(data_uniform_0SD_PG_140.incidenceAngle, data_uniform_0SD_PG_140.acceptancePercentageDay, 'linewidth',2, 'color', cols(4,:), 'linestyle', style);
plot(data_uniform_0SD_PG_140.acceptanceAngleDay, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(4,:),'linewidth',2, 'color', cols(4,:))

if data_uniform_0SD_PG_145.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h3 = plot(data_uniform_0SD_PG_145.incidenceAngle, data_uniform_0SD_PG_145.acceptancePercentageDay, 'linewidth',2, 'color', cols(5,:), 'linestyle', style);
plot(data_uniform_0SD_PG_145.acceptanceAngleDay, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(5,:),'linewidth',2, 'color', cols(5,:))

if data_uniform_0SD_PG_150.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h4 = plot(data_uniform_0SD_PG_150.incidenceAngle, data_uniform_0SD_PG_150.acceptancePercentageDay, 'linewidth',2, 'color', cols(6,:), 'linestyle', style);
plot(data_uniform_0SD_PG_150.acceptanceAngleDay, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(6,:),'linewidth',2, 'color', cols(6,:))

plot(angleSteps, interpAcceptanceDay, 'k:', 'linewidth',2)
% plot(acceptanceAngleDay, 0.5, 'kd','markersize',10, 'linewidth',2)

ylim([0 1.25]); xlim([0 20])
ylabel('% entering'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

subplot(4,nPlots,15); hold on
if data_GRIN_cylinder.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h1 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.acceptancePercentageDay, 'linewidth',2, 'color', cols(2,:), 'linestyle', style);
plot(data_GRIN_cylinder.acceptanceAngleDay, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(2,:),'linewidth',2, 'color', cols(2,:))

if data_GRIN_cylinder_PG_140.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h2 = plot(data_GRIN_cylinder_PG_140.incidenceAngle, data_GRIN_cylinder_PG_140.acceptancePercentageDay, 'linewidth',2, 'color', cols(4,:), 'linestyle', style);
plot(data_GRIN_cylinder_PG_140.acceptanceAngleDay, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(4,:),'linewidth',2, 'color', cols(4,:))

if data_GRIN_cylinder_PG_145.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h3 = plot(data_GRIN_cylinder_PG_145.incidenceAngle, data_GRIN_cylinder_PG_145.acceptancePercentageDay, 'linewidth',2, 'color', cols(5,:), 'linestyle', style);
plot(data_GRIN_cylinder_PG_145.acceptanceAngleDay, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(5,:),'linewidth',2, 'color', cols(5,:))

if data_GRIN_cylinder_PG_150.acceptanceNumDay(1) > 1; style = '-'; else; style = ':'; end
h4 = plot(data_GRIN_cylinder_PG_150.incidenceAngle, data_GRIN_cylinder_PG_150.acceptancePercentageDay, 'linewidth',2, 'color', cols(6,:), 'linestyle', style);
plot(data_GRIN_cylinder_PG_150.acceptanceAngleDay, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(6,:),'linewidth',2, 'color', cols(6,:))

plot(angleSteps, interpAcceptanceDay, 'k:', 'linewidth',2)
% plot(acceptanceAngleDay, 0.5, 'kd','markersize',10, 'linewidth',2)

ylim([0 1.25]); xlim([0 20])
ylabel('% entering'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

%% Plot acceptance angle for COLC
acceptanceFCOLC = figure; set(gcf, 'position', [-1919        -149        1920        1104])

subplot(4,nPlots,1); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.acceptancePercentageColc, 'linewidth',2, 'color', cols(1,:));
plot(data_GRIN_radialTop.acceptanceAngleColc, 0.5, 'o','markersize',5, 'linewidth',2, 'color', cols(1,:))

h2= plot(data_GRIN_radialBase.incidenceAngle, data_GRIN_radialBase.acceptancePercentageColc, 'linewidth',2, 'color', cols(4,:));
plot(data_GRIN_radialBase.acceptanceAngleColc, 0.5, 'o','markersize',5, 'linewidth',2, 'color', cols(4,:))

h3= plot(data_GRIN_linear.incidenceAngle, data_GRIN_linear.acceptancePercentageColc, 'linewidth',2, 'color', cols(5,:));
plot(data_GRIN_linear.acceptanceAngleColc, 0.5, 'o','markersize',5, 'linewidth',2, 'color', cols(5,:))

h4= plot(data_GRIN_both.incidenceAngle, data_GRIN_both.acceptancePercentageColc, 'linewidth',2, 'color', cols(6,:));
plot(data_GRIN_both.acceptanceAngleColc, 0.5, 'o','markersize',5, 'linewidth',2, 'color', cols(6,:))

ylim([0 1.25]); xlim([0 20])
ylabel('% exiting'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

subplot(4,nPlots,4); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.acceptancePercentageColc, 'linewidth',2, 'color', cols(1,:));
plot(data_GRIN_radialTop.acceptanceAngleColc, 0.5, 'o','markersize',5, 'linewidth',2, 'color', cols(1,:))

h2= plot(data_GRIN_radialTop_tipCorr.incidenceAngle, data_GRIN_radialTop_tipCorr.acceptancePercentageColc, 'linewidth',2, 'color', cols(4,:));
plot(data_GRIN_radialTop_tipCorr.acceptanceAngleColc, 0.5, 'o','markersize',5, 'linewidth',2, 'color', cols(4,:))

ylim([0 1.25]); xlim([0 20])
ylabel('% exiting'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

subplot(4,nPlots,3); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.acceptancePercentageColc, 'linewidth',2, 'color', cols(2,:));
plot(data_GRIN_cylinder.acceptanceAngleColc, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(2,:), 'linewidth',2, 'color', cols(2,:))

h2 = plot(data_GRIN_cylinder_p36.incidenceAngle, data_GRIN_cylinder_p36.acceptancePercentageColc, 'linewidth',2, 'color', cols(4,:));
plot(data_GRIN_cylinder_p36.acceptanceAngleColc, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(4,:),'linewidth',2, 'color', cols(4,:))

h3 = plot(data_GRIN_cylinder_m36.incidenceAngle, data_GRIN_cylinder_m36.acceptancePercentageColc, 'linewidth',2, 'color', cols(5,:));
plot(data_GRIN_cylinder_m36.acceptanceAngleColc, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(5,:), 'linewidth',2, 'color', cols(5,:))

h4 = plot(data_GRIN_cylinder_shape.incidenceAngle, data_GRIN_cylinder_shape.acceptancePercentageColc, 'linewidth',2, 'color', cols(6,:));
plot(data_GRIN_cylinder_shape.acceptanceAngleColc, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(6,:), 'linewidth',2, 'color', cols(6,:))

ylim([0 1.25]); xlim([0 20])
ylabel('% exiting'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

subplot(4,nPlots,2); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.acceptancePercentageColc, 'linewidth',2, 'color', cols(1,:));
plot(data_GRIN_radialTop.acceptanceAngleColc, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(1,:), 'linewidth',2, 'color', cols(1,:))

h2 = plot(data_GRIN_radialTop_CinC.incidenceAngle, data_GRIN_radialTop_CinC.acceptancePercentageColc, 'linewidth',2, 'color', cols(4,:));
plot(data_GRIN_radialTop_CinC.acceptanceAngleColc, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(4,:), 'linewidth',2, 'color', cols(4,:))

h3 = plot(data_GRIN_radialTop_EC.incidenceAngle, data_GRIN_radialTop_EC.acceptancePercentageColc, 'linewidth',2, 'color', cols(5,:));
plot(data_GRIN_radialTop_EC.acceptanceAngleColc, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(5,:), 'linewidth',2, 'color', cols(5,:))

h4 = plot(data_GRIN_radialTop_CinC_EC.incidenceAngle, data_GRIN_radialTop_CinC_EC.acceptancePercentageColc, 'linewidth',2, 'color', cols(6,:));
plot(data_GRIN_radialTop_CinC_EC.acceptanceAngleColc, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(6,:), 'linewidth',2, 'color', cols(6,:))

ylim([0 1.25]); xlim([0 20])
ylabel('% exiting'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

subplot(4,nPlots,8); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.acceptancePercentageColc, 'linewidth',2, 'color', cols(1,:));
plot(data_GRIN_radialTop.acceptanceAngleColc, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(1,:), 'linewidth',2, 'color', cols(1,:))

h2 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.acceptancePercentageColc, 'linewidth',2, 'color', cols(2,:));
plot(data_GRIN_cylinder.acceptanceAngleColc, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(2,:),'linewidth',2, 'color', cols(2,:))

h3 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.acceptancePercentageColc, 'linewidth',2, 'color', cols(3,:));
plot(data_uniform_0SD.acceptanceAngleColc, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(3,:),'linewidth',2, 'color', cols(3,:))

h4 = plot(data_GRIN_cylinder_norm_max.incidenceAngle, data_GRIN_cylinder_norm_max.acceptancePercentageColc, 'linewidth',2, 'color', cols(4,:));
plot(data_GRIN_cylinder_norm_max.acceptanceAngleColc, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(4,:), 'linewidth',2, 'color', cols(4,:))

h5 = plot(data_GRIN_cylinder_norm_80.incidenceAngle, data_GRIN_cylinder_norm_80.acceptancePercentageColc, 'linewidth',2, 'color', cols(5,:));
plot(data_GRIN_cylinder_norm_80.acceptanceAngleColc, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(5,:), 'linewidth',2, 'color', cols(5,:))

ylim([0 1.25]); xlim([0 20])
ylabel('% exiting'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

%uniform RI
subplot(4,nPlots,5); hold on
h1 = plot(data_uniform_m4SD.incidenceAngle, data_uniform_m4SD.acceptancePercentageColc, 'linewidth',2, 'color', cols(4,:));
plot(data_uniform_m4SD.acceptanceAngleColc, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(4,:), 'linewidth',2, 'color', cols(4,:))

h2 = plot(data_uniform_m2SD.incidenceAngle, data_uniform_m2SD.acceptancePercentageColc, 'linewidth',2, 'color', cols(5,:));
plot(data_uniform_m2SD.acceptanceAngleColc, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(5,:), 'linewidth',2, 'color', cols(5,:))

h3 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.acceptancePercentageColc, 'linewidth',2, 'color', cols(3,:));
plot(data_uniform_0SD.acceptanceAngleColc, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(3,:), 'linewidth',2, 'color', cols(3,:))

h4 = plot(data_uniform_p2SD.incidenceAngle, data_uniform_p2SD.acceptancePercentageColc, 'linewidth',2, 'color', cols(6,:));
plot(data_uniform_p2SD.acceptanceAngleColc, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(6,:), 'linewidth',2, 'color', cols(6,:))

h5 = plot(data_uniform_p4SD.incidenceAngle, data_uniform_p4SD.acceptancePercentageColc, 'linewidth',2, 'color', cols(7,:));
plot(data_uniform_p4SD.acceptanceAngleColc, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(7,:), 'linewidth',2, 'color', cols(7,:))

ylim([0 1.25]); xlim([0 20])
ylabel('% exiting'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

subplot(4,nPlots,6); hold on
h1 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.acceptancePercentageColc, 'linewidth',2, 'color', cols(3,:));
plot(data_uniform_0SD.acceptanceAngleColc, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(3,:), 'linewidth',2, 'color', cols(3,:))

h2 = plot(data_GRIN_uniform_CinC.incidenceAngle, data_GRIN_uniform_CinC.acceptancePercentageColc, 'linewidth',2, 'color', cols(4,:));
plot(data_GRIN_uniform_CinC.acceptanceAngleColc, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(4,:), 'linewidth',2, 'color', cols(4,:))

h3 = plot(data_GRIN_uniform_EC.incidenceAngle, data_GRIN_uniform_EC.acceptancePercentageColc, 'linewidth',2, 'color', cols(5,:));
plot(data_GRIN_uniform_EC.acceptanceAngleColc, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(5,:), 'linewidth',2, 'color', cols(5,:))

h4 = plot(data_GRIN_uniform_CinC_EC.incidenceAngle, data_GRIN_uniform_CinC_EC.acceptancePercentageColc, 'linewidth',2, 'color', cols(6,:));
plot(data_GRIN_uniform_CinC_EC.acceptanceAngleColc, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(6,:), 'linewidth',2, 'color', cols(6,:))

ylim([0 1.25]); xlim([0 20])
ylabel('% exiting'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

% Varying intercone RI
subplot(4,nPlots,9); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.acceptancePercentageColc, 'linewidth',2, 'color', cols(1,:));
plot(data_GRIN_radialTop.acceptanceAngleColc, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(1,:), 'linewidth',2, 'color', cols(1,:))

h2 = plot(data_GRIN_radialTop_IC_142.incidenceAngle, data_GRIN_radialTop_IC_142.acceptancePercentageColc, 'linewidth',2, 'color', cols(4,:));
plot(data_GRIN_radialTop_IC_142.acceptanceAngleColc, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(4,:), 'linewidth',2, 'color', cols(4,:))

h3 = plot(data_GRIN_radialTop_IC_144.incidenceAngle, data_GRIN_radialTop_IC_144.acceptancePercentageColc, 'linewidth',2, 'color', cols(5,:));
plot(data_GRIN_radialTop_IC_144.acceptanceAngleColc, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(5,:), 'linewidth',2, 'color', cols(5,:))

h4 = plot(data_GRIN_radialTop_IC_146.incidenceAngle, data_GRIN_radialTop_IC_146.acceptancePercentageColc, 'linewidth',2, 'color', cols(6,:));
plot(data_GRIN_radialTop_IC_146.acceptanceAngleColc, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(6,:), 'linewidth',2, 'color', cols(6,:))

h5 = plot(data_GRIN_radialTop_IC_148.incidenceAngle, data_GRIN_radialTop_IC_148.acceptancePercentageColc, 'linewidth',2, 'color', cols(7,:));
plot(data_GRIN_radialTop_IC_148.acceptanceAngleColc, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(7,:), 'linewidth',2, 'color', cols(7,:))

ylim([0 1.25]); xlim([0 20])
ylabel('% exiting'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

subplot(4,nPlots,10); hold on
h1 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.acceptancePercentageColc, 'linewidth',2, 'color', cols(3,:));
plot(data_uniform_0SD.acceptanceAngleColc, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(3,:), 'linewidth',2, 'color', cols(3,:))

h2 = plot(data_uniform_0SD_IC_142.incidenceAngle, data_uniform_0SD_IC_142.acceptancePercentageColc, 'linewidth',2, 'color', cols(4,:));
plot(data_uniform_0SD_IC_142.acceptanceAngleColc, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(4,:), 'linewidth',2, 'color', cols(4,:))

h3 = plot(data_uniform_0SD_IC_144.incidenceAngle, data_uniform_0SD_IC_144.acceptancePercentageColc, 'linewidth',2, 'color', cols(5,:));
plot(data_uniform_0SD_IC_144.acceptanceAngleColc, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(5,:), 'linewidth',2, 'color', cols(5,:))

h4 = plot(data_uniform_0SD_IC_146.incidenceAngle, data_uniform_0SD_IC_146.acceptancePercentageColc, 'linewidth',2, 'color', cols(6,:));
plot(data_uniform_0SD_IC_146.acceptanceAngleColc, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(6,:), 'linewidth',2, 'color', cols(6,:))

h5 = plot(data_uniform_0SD_IC_148.incidenceAngle, data_uniform_0SD_IC_148.acceptancePercentageColc, 'linewidth',2, 'color', cols(7,:));
plot(data_uniform_0SD_IC_148.acceptanceAngleColc, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(7,:), 'linewidth',2, 'color', cols(7,:))

ylim([0 1.25]); xlim([0 20])
ylabel('% exiting'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

subplot(4,nPlots,11); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.acceptancePercentageColc, 'linewidth',2, 'color', cols(2,:));
plot(data_GRIN_cylinder.acceptanceAngleColc, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(2,:), 'linewidth',2, 'color', cols(2,:))

h2 = plot(data_GRIN_cylinder_IC_142.incidenceAngle, data_GRIN_cylinder_IC_142.acceptancePercentageColc, 'linewidth',2, 'color', cols(4,:));
plot(data_GRIN_cylinder_IC_142.acceptanceAngleColc, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(4,:), 'linewidth',2, 'color', cols(4,:))

h3 = plot(data_GRIN_cylinder_IC_144.incidenceAngle, data_GRIN_cylinder_IC_144.acceptancePercentageColc, 'linewidth',2, 'color', cols(5,:));
plot(data_GRIN_cylinder_IC_144.acceptanceAngleColc, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(5,:), 'linewidth',2, 'color', cols(5,:))

h4 = plot(data_GRIN_cylinder_IC_146.incidenceAngle, data_GRIN_cylinder_IC_146.acceptancePercentageColc, 'linewidth',2, 'color', cols(6,:));
plot(data_GRIN_cylinder_IC_146.acceptanceAngleColc, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(6,:), 'linewidth',2, 'color', cols(6,:))

h5 = plot(data_GRIN_cylinder_IC_148.incidenceAngle, data_GRIN_cylinder_IC_148.acceptancePercentageColc, 'linewidth',2, 'color', cols(7,:));
plot(data_GRIN_cylinder_IC_148.acceptanceAngleColc, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(7,:), 'linewidth',2, 'color', cols(7,:))

ylim([0 1.25]); xlim([0 20])
ylabel('% exiting'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

% Varying pigment RI
subplot(4,nPlots,13); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.acceptancePercentageColc, 'linewidth',2, 'color', cols(1,:));
plot(data_GRIN_radialTop.acceptanceAngleColc, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(1,:),'linewidth',2, 'color', cols(1,:))

h2 = plot(data_GRIN_radialTop_PG_140.incidenceAngle, data_GRIN_radialTop_PG_140.acceptancePercentageColc, 'linewidth',2, 'color', cols(4,:));
plot(data_GRIN_radialTop_PG_140.acceptanceAngleColc, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(4,:),'linewidth',2, 'color', cols(4,:))

h3 = plot(data_GRIN_radialTop_PG_145.incidenceAngle, data_GRIN_radialTop_PG_145.acceptancePercentageColc, 'linewidth',2, 'color', cols(5,:));
plot(data_GRIN_radialTop_PG_145.acceptanceAngleColc, 0.5, 'o','markersize',5,'MarkerFaceColor', cols(5,:), 'linewidth',2, 'color', cols(5,:))

h4 = plot(data_GRIN_radialTop_PG_150.incidenceAngle, data_GRIN_radialTop_PG_150.acceptancePercentageColc, 'linewidth',2, 'color', cols(6,:));
plot(data_GRIN_radialTop_PG_150.acceptanceAngleColc, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(6,:),'linewidth',2, 'color', cols(6,:))

ylim([0 1.25]); xlim([0 20])
ylabel('% exiting'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

subplot(4,nPlots,14); hold on
h1 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.acceptancePercentageColc, 'linewidth',2, 'color', cols(3,:));
plot(data_uniform_0SD.acceptanceAngleColc, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(3,:),'linewidth',2, 'color', cols(3,:))

h2 = plot(data_uniform_0SD_PG_140.incidenceAngle, data_uniform_0SD_PG_140.acceptancePercentageColc, 'linewidth',2, 'color', cols(4,:));
plot(data_uniform_0SD_PG_140.acceptanceAngleColc, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(4,:),'linewidth',2, 'color', cols(4,:))

h3 = plot(data_uniform_0SD_PG_145.incidenceAngle, data_uniform_0SD_PG_145.acceptancePercentageColc, 'linewidth',2, 'color', cols(5,:));
plot(data_uniform_0SD_PG_145.acceptanceAngleColc, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(5,:),'linewidth',2, 'color', cols(5,:))

h4 = plot(data_uniform_0SD_PG_150.incidenceAngle, data_uniform_0SD_PG_150.acceptancePercentageColc, 'linewidth',2, 'color', cols(6,:));
plot(data_uniform_0SD_PG_150.acceptanceAngleColc, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(6,:),'linewidth',2, 'color', cols(6,:))

ylim([0 1.25]); xlim([0 20])
ylabel('% exiting'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

subplot(4,nPlots,15); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.acceptancePercentageColc, 'linewidth',2, 'color', cols(2,:));
plot(data_GRIN_cylinder.acceptanceAngleColc, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(2,:),'linewidth',2, 'color', cols(2,:))

h2 = plot(data_GRIN_cylinder_PG_140.incidenceAngle, data_GRIN_cylinder_PG_140.acceptancePercentageColc, 'linewidth',2, 'color', cols(4,:));
plot(data_GRIN_cylinder_PG_140.acceptanceAngleColc, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(4,:),'linewidth',2, 'color', cols(4,:))

h3 = plot(data_GRIN_cylinder_PG_145.incidenceAngle, data_GRIN_cylinder_PG_145.acceptancePercentageColc, 'linewidth',2, 'color', cols(5,:));
plot(data_GRIN_cylinder_PG_145.acceptanceAngleColc, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(5,:),'linewidth',2, 'color', cols(5,:))

h4 = plot(data_GRIN_cylinder_PG_150.incidenceAngle, data_GRIN_cylinder_PG_150.acceptancePercentageColc, 'linewidth',2, 'color', cols(6,:));
plot(data_GRIN_cylinder_PG_150.acceptanceAngleColc, 0.5, 'o','markersize',5, 'MarkerFaceColor', cols(6,:),'linewidth',2, 'color', cols(6,:))

ylim([0 1.25]); xlim([0 20])
ylabel('% exiting'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1 1.25],'XTick',[0 5 10 15 20]);

%% Plot TIR
TIRF = figure; set(gcf, 'position', [-1919        -149        1920        1104])

subplot(4,nPlots,1); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.TIRPercengtage, 'linewidth',2, 'color', cols(1,:));

h2= plot(data_GRIN_radialBase.incidenceAngle, data_GRIN_radialBase.TIRPercengtage, 'linewidth',2, 'color', cols(4,:));

h3= plot(data_GRIN_linear.incidenceAngle, data_GRIN_linear.TIRPercengtage, 'linewidth',2, 'color', cols(5,:));

h4= plot(data_GRIN_both.incidenceAngle, data_GRIN_both.TIRPercengtage, 'linewidth',2, 'color', cols(6,:));

ylim([0 1]); xlim([0 20]); ylabel('% TIR'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);

subplot(4,nPlots,4); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.TIRPercengtage, 'linewidth',2, 'color', cols(1,:));

h2= plot(data_GRIN_radialTop_tipCorr.incidenceAngle, data_GRIN_radialTop_tipCorr.TIRPercengtage, 'linewidth',2, 'color', cols(4,:));

ylim([0 1]); xlim([0 20]); ylabel('% TIR'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);

subplot(4,nPlots,3); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.TIRPercengtage, 'linewidth',2, 'color', cols(2,:));

h2 = plot(data_GRIN_cylinder_p36.incidenceAngle, data_GRIN_cylinder_p36.TIRPercengtage, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_cylinder_m36.incidenceAngle, data_GRIN_cylinder_m36.TIRPercengtage, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_cylinder_shape.incidenceAngle, data_GRIN_cylinder_shape.TIRPercengtage, 'linewidth',2, 'color', cols(6,:));

ylim([0 1]); xlim([0 20]); ylabel('% TIR'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);

subplot(4,nPlots,2); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.TIRPercengtage, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_radialTop_CinC.incidenceAngle, data_GRIN_radialTop_CinC.TIRPercengtage, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_radialTop_EC.incidenceAngle, data_GRIN_radialTop_EC.TIRPercengtage, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_radialTop_CinC_EC.incidenceAngle, data_GRIN_radialTop_CinC_EC.TIRPercengtage, 'linewidth',2, 'color', cols(6,:));

ylim([0 1]); xlim([0 20]); ylabel('% TIR'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);

subplot(4,nPlots,8); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.TIRPercengtage, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.TIRPercengtage, 'linewidth',2, 'color', cols(2,:));

h3 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.TIRPercengtage, 'linewidth',2, 'color', cols(3,:));

h4 = plot(data_GRIN_cylinder_norm_max.incidenceAngle, data_GRIN_cylinder_norm_max.TIRPercengtage, 'linewidth',2, 'color', cols(4,:));

h5 = plot(data_GRIN_cylinder_norm_80.incidenceAngle, data_GRIN_cylinder_norm_80.TIRPercengtage, 'linewidth',2, 'color', cols(5,:));

ylim([0 1]); xlim([0 20]); ylabel('% TIR'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);

% Uniform RI
subplot(4,nPlots,5); hold on
h1 = plot(data_uniform_m4SD.incidenceAngle, data_uniform_m4SD.TIRPercengtage, 'linewidth',2, 'color', cols(4,:));

h2 = plot(data_uniform_m2SD.incidenceAngle, data_uniform_m2SD.TIRPercengtage, 'linewidth',2, 'color', cols(5,:));

h3 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.TIRPercengtage, 'linewidth',2, 'color', cols(3,:));

h4 = plot(data_uniform_p2SD.incidenceAngle, data_uniform_p2SD.TIRPercengtage, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_uniform_p4SD.incidenceAngle, data_uniform_p4SD.TIRPercengtage, 'linewidth',2, 'color', cols(7,:));

ylim([0 1]); xlim([0 20]); ylabel('% TIR'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);

subplot(4,nPlots,6); hold on
h1 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.TIRPercengtage, 'linewidth',2, 'color', cols(3,:));

h2 = plot(data_GRIN_uniform_CinC.incidenceAngle, data_GRIN_uniform_CinC.TIRPercengtage, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_uniform_EC.incidenceAngle, data_GRIN_uniform_EC.TIRPercengtage, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_uniform_CinC_EC.incidenceAngle, data_GRIN_uniform_CinC_EC.TIRPercengtage, 'linewidth',2, 'color', cols(6,:));

ylim([0 1]); xlim([0 20]); ylabel('% TIR'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);

% Varying intercone RI
subplot(4,nPlots,9); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.TIRPercengtage, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_radialTop_IC_142.incidenceAngle, data_GRIN_radialTop_IC_142.TIRPercengtage, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_radialTop_IC_144.incidenceAngle, data_GRIN_radialTop_IC_144.TIRPercengtage, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_radialTop_IC_146.incidenceAngle, data_GRIN_radialTop_IC_146.TIRPercengtage, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_GRIN_radialTop_IC_148.incidenceAngle, data_GRIN_radialTop_IC_148.TIRPercengtage, 'linewidth',2, 'color', cols(7,:));

ylim([0 1]); xlim([0 20]); ylabel('% TIR'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);

subplot(4,nPlots,10); hold on
h1 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.TIRPercengtage, 'linewidth',2, 'color', cols(3,:));

h2 = plot(data_uniform_0SD_IC_142.incidenceAngle, data_uniform_0SD_IC_142.TIRPercengtage, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_uniform_0SD_IC_144.incidenceAngle, data_uniform_0SD_IC_144.TIRPercengtage, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_uniform_0SD_IC_146.incidenceAngle, data_uniform_0SD_IC_146.TIRPercengtage, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_uniform_0SD_IC_148.incidenceAngle, data_uniform_0SD_IC_148.TIRPercengtage, 'linewidth',2, 'color', cols(7,:));

ylim([0 1]); xlim([0 20]); ylabel('% TIR'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);

subplot(4,nPlots,11); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.TIRPercengtage, 'linewidth',2, 'color', cols(2,:));

h2 = plot(data_GRIN_cylinder_IC_142.incidenceAngle, data_GRIN_cylinder_IC_142.TIRPercengtage, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_cylinder_IC_144.incidenceAngle, data_GRIN_cylinder_IC_144.TIRPercengtage, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_cylinder_IC_146.incidenceAngle, data_GRIN_cylinder_IC_146.TIRPercengtage, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_GRIN_cylinder_IC_148.incidenceAngle, data_GRIN_cylinder_IC_148.TIRPercengtage, 'linewidth',2, 'color', cols(7,:));

ylim([0 1]); xlim([0 20]); ylabel('% TIR'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);

% Varying pigment RI
subplot(4,nPlots,13); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.TIRPercengtage, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_radialTop_PG_140.incidenceAngle, data_GRIN_radialTop_PG_140.TIRPercengtage, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_radialTop_PG_145.incidenceAngle, data_GRIN_radialTop_PG_145.TIRPercengtage, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_radialTop_PG_150.incidenceAngle, data_GRIN_radialTop_PG_150.TIRPercengtage, 'linewidth',2, 'color', cols(6,:));

ylim([0 1]); xlim([0 20]); ylabel('% TIR'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);

subplot(4,nPlots,14); hold on
h1 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.TIRPercengtage, 'linewidth',2, 'color', cols(3,:));

h2 = plot(data_uniform_0SD_PG_140.incidenceAngle, data_uniform_0SD_PG_140.TIRPercengtage, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_uniform_0SD_PG_145.incidenceAngle, data_uniform_0SD_PG_145.TIRPercengtage, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_uniform_0SD_PG_150.incidenceAngle, data_uniform_0SD_PG_150.TIRPercengtage, 'linewidth',2, 'color', cols(6,:));

ylim([0 1]); xlim([0 20]); ylabel('% TIR'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);

subplot(4,nPlots,15); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.TIRPercengtage, 'linewidth',2, 'color', cols(2,:));

h2 = plot(data_GRIN_cylinder_PG_140.incidenceAngle, data_GRIN_cylinder_PG_140.TIRPercengtage, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_cylinder_PG_145.incidenceAngle, data_GRIN_cylinder_PG_145.TIRPercengtage, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_cylinder_PG_150.incidenceAngle, data_GRIN_cylinder_PG_150.TIRPercengtage, 'linewidth',2, 'color', cols(6,:));

ylim([0 1]); xlim([0 20]); ylabel('% TIR'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);
%% Plot COLC height
COLCHeightF = figure; set(gcf, 'position', [-1919        -149        1920        1104])

subplot(4,nPlots,1); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, (data_GRIN_radialTop.colcHeight - data_GRIN_radialTop.coneTipZ)*1000, 'linewidth',2, 'color', cols(1,:));

h2= plot(data_GRIN_radialBase.incidenceAngle, (data_GRIN_radialBase.colcHeight- data_GRIN_radialBase.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h3= plot(data_GRIN_linear.incidenceAngle, (data_GRIN_linear.colcHeight- data_GRIN_linear.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

h4= plot(data_GRIN_both.incidenceAngle, (data_GRIN_both.colcHeight- data_GRIN_both.coneTipZ)*1000, 'linewidth',2, 'color', cols(6,:));

% line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
ylim([0 50]); ylabel('Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

subplot(4,nPlots,4); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, (data_GRIN_radialTop.colcHeight- data_GRIN_radialTop.coneTipZ)*1000, 'linewidth',2, 'color', cols(1,:));

h2= plot(data_GRIN_radialTop_tipCorr.incidenceAngle, (data_GRIN_radialTop_tipCorr.colcHeight- data_GRIN_radialTop_tipCorr.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

% line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
ylim([0 50]); ylabel('Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

subplot(4,nPlots,3); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, (data_GRIN_cylinder.colcHeight- data_GRIN_cylinder.coneTipZ)*1000, 'linewidth',2, 'color', cols(2,:));

h2 = plot(data_GRIN_cylinder_p36.incidenceAngle, (data_GRIN_cylinder_p36.colcHeight- data_GRIN_cylinder_p36.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_cylinder_m36.incidenceAngle, (data_GRIN_cylinder_m36.colcHeight- data_GRIN_cylinder_m36.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_cylinder_shape.incidenceAngle, (data_GRIN_cylinder_shape.colcHeight- data_GRIN_cylinder_shape.coneTipZ)*1000, 'linewidth',2, 'color', cols(6,:));

% line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
ylim([0 50]); ylabel('Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

subplot(4,nPlots,2); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, (data_GRIN_radialTop.colcHeight- data_GRIN_radialTop.coneTipZ)*1000, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_radialTop_CinC.incidenceAngle, (data_GRIN_radialTop_CinC.colcHeight- data_GRIN_radialTop_CinC.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_radialTop_EC.incidenceAngle, (data_GRIN_radialTop_EC.colcHeight- data_GRIN_radialTop_EC.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_radialTop_CinC_EC.incidenceAngle, (data_GRIN_radialTop_CinC_EC.colcHeight- data_GRIN_radialTop_CinC_EC.coneTipZ)*1000, 'linewidth',2, 'color', cols(6,:));

% line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
ylim([0 50]); ylabel('Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

subplot(4,nPlots,8); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, (data_GRIN_radialTop.colcHeight- data_GRIN_radialTop.coneTipZ)*1000, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_cylinder.incidenceAngle, (data_GRIN_cylinder.colcHeight- data_GRIN_cylinder.coneTipZ)*1000, 'linewidth',2, 'color', cols(2,:));

h3 = plot(data_uniform_0SD.incidenceAngle, (data_uniform_0SD.colcHeight- data_uniform_0SD.coneTipZ)*1000, 'linewidth',2, 'color', cols(3,:));

h4 = plot(data_GRIN_cylinder_norm_max.incidenceAngle, (data_GRIN_cylinder_norm_max.colcHeight- data_GRIN_cylinder_norm_max.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h5 = plot(data_GRIN_cylinder_norm_80.incidenceAngle, (data_GRIN_cylinder_norm_80.colcHeight- data_GRIN_cylinder_norm_80.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

% line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
ylim([0 50]); ylabel('Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

% Uniform RI
subplot(4,nPlots,5); hold on
h1 = plot(data_uniform_m4SD.incidenceAngle, (data_uniform_m4SD.colcHeight- data_uniform_m4SD.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h2 = plot(data_uniform_m2SD.incidenceAngle, (data_uniform_m2SD.colcHeight- data_uniform_m2SD.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

h3 = plot(data_uniform_0SD.incidenceAngle, (data_uniform_0SD.colcHeight- data_uniform_0SD.coneTipZ)*1000, 'linewidth',2, 'color', cols(3,:));

h4 = plot(data_uniform_p2SD.incidenceAngle, (data_uniform_p2SD.colcHeight- data_uniform_p2SD.coneTipZ)*1000, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_uniform_p4SD.incidenceAngle, (data_uniform_p4SD.colcHeight- data_uniform_p4SD.coneTipZ)*1000, 'linewidth',2, 'color', cols(7,:));

% line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
ylim([0 50]); ylabel('Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

subplot(4,nPlots,6); hold on
h1 = plot(data_uniform_0SD.incidenceAngle, (data_uniform_0SD.colcHeight- data_uniform_0SD.coneTipZ)*1000, 'linewidth',2, 'color', cols(3,:));

h2 = plot(data_GRIN_uniform_CinC.incidenceAngle, (data_GRIN_uniform_CinC.colcHeight- data_GRIN_uniform_CinC.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_uniform_EC.incidenceAngle, (data_GRIN_uniform_EC.colcHeight- data_GRIN_uniform_EC.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_uniform_CinC_EC.incidenceAngle, (data_GRIN_uniform_CinC_EC.colcHeight- data_GRIN_uniform_CinC_EC.coneTipZ)*1000, 'linewidth',2, 'color', cols(6,:));

% line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
ylim([0 50]); ylabel('Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

% Varying intercone RI
subplot(4,nPlots,9); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, (data_GRIN_radialTop.colcHeight- data_GRIN_radialTop.coneTipZ)*1000, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_radialTop_IC_142.incidenceAngle, (data_GRIN_radialTop_IC_142.colcHeight- data_GRIN_radialTop_IC_142.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_radialTop_IC_144.incidenceAngle, (data_GRIN_radialTop_IC_144.colcHeight- data_GRIN_radialTop_IC_144.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_radialTop_IC_146.incidenceAngle, (data_GRIN_radialTop_IC_146.colcHeight- data_GRIN_radialTop_IC_146.coneTipZ)*1000, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_GRIN_radialTop_IC_148.incidenceAngle, (data_GRIN_radialTop_IC_148.colcHeight- data_GRIN_radialTop_IC_148.coneTipZ)*1000, 'linewidth',2, 'color', cols(7,:));

% line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
ylim([0 50]); ylabel('Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

subplot(4,nPlots,10); hold on
h1 = plot(data_uniform_0SD.incidenceAngle, (data_uniform_0SD.colcHeight- data_uniform_0SD.coneTipZ)*1000, 'linewidth',2, 'color', cols(3,:));

h2 = plot(data_uniform_0SD_IC_142.incidenceAngle, (data_uniform_0SD_IC_142.colcHeight- data_uniform_0SD_IC_142.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_uniform_0SD_IC_144.incidenceAngle, (data_uniform_0SD_IC_144.colcHeight- data_uniform_0SD_IC_144.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_uniform_0SD_IC_146.incidenceAngle, (data_uniform_0SD_IC_146.colcHeight- data_uniform_0SD_IC_146.coneTipZ)*1000, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_uniform_0SD_IC_148.incidenceAngle, (data_uniform_0SD_IC_148.colcHeight- data_uniform_0SD_IC_148.coneTipZ)*1000, 'linewidth',2, 'color', cols(7,:));

% line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
ylim([0 50]); ylabel('Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

subplot(4,nPlots,11); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, (data_GRIN_cylinder.colcHeight- data_GRIN_cylinder.coneTipZ)*1000, 'linewidth',2, 'color', cols(2,:));

h2 = plot(data_GRIN_cylinder_IC_142.incidenceAngle, (data_GRIN_cylinder_IC_142.colcHeight- data_GRIN_cylinder_IC_142.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_cylinder_IC_144.incidenceAngle, (data_GRIN_cylinder_IC_144.colcHeight- data_GRIN_cylinder_IC_144.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_cylinder_IC_146.incidenceAngle, (data_GRIN_cylinder_IC_146.colcHeight- data_GRIN_cylinder_IC_146.coneTipZ)*1000, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_GRIN_cylinder_IC_148.incidenceAngle, (data_GRIN_cylinder_IC_148.colcHeight- data_GRIN_cylinder_IC_148.coneTipZ)*1000, 'linewidth',2, 'color', cols(7,:));

% line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
ylim([0 50]); ylabel('Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

% Varying pigment RI
subplot(4,nPlots,13); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, (data_GRIN_radialTop.colcHeight- data_GRIN_radialTop.coneTipZ)*1000, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_radialTop_PG_140.incidenceAngle, (data_GRIN_radialTop_PG_140.colcHeight- data_GRIN_radialTop_PG_140.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_radialTop_PG_145.incidenceAngle, (data_GRIN_radialTop_PG_145.colcHeight- data_GRIN_radialTop_PG_145.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_radialTop_PG_150.incidenceAngle, (data_GRIN_radialTop_PG_150.colcHeight- data_GRIN_radialTop_PG_150.coneTipZ)*1000, 'linewidth',2, 'color', cols(6,:));

% line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
ylim([0 50]); ylabel('Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

subplot(4,nPlots,14); hold on
h1 = plot(data_uniform_0SD.incidenceAngle, (data_uniform_0SD.colcHeight- data_uniform_0SD.coneTipZ)*1000, 'linewidth',2, 'color', cols(3,:));

h2 = plot(data_uniform_0SD_PG_140.incidenceAngle, (data_uniform_0SD_PG_140.colcHeight- data_uniform_0SD_PG_140.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_uniform_0SD_PG_145.incidenceAngle, (data_uniform_0SD_PG_145.colcHeight- data_uniform_0SD_PG_145.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_uniform_0SD_PG_150.incidenceAngle, (data_uniform_0SD_PG_150.colcHeight- data_uniform_0SD_PG_150.coneTipZ)*1000, 'linewidth',2, 'color', cols(6,:));

% line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
ylim([0 50]); ylabel('Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

subplot(4,nPlots,15); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, (data_GRIN_cylinder.colcHeight- data_GRIN_cylinder.coneTipZ)*1000, 'linewidth',2, 'color', cols(2,:));

h2 = plot(data_GRIN_cylinder_PG_140.incidenceAngle, (data_GRIN_cylinder_PG_140.colcHeight- data_GRIN_cylinder_PG_140.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_cylinder_PG_145.incidenceAngle, (data_GRIN_cylinder_PG_145.colcHeight- data_GRIN_cylinder_PG_145.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_cylinder_PG_150.incidenceAngle, (data_GRIN_cylinder_PG_150.colcHeight- data_GRIN_cylinder_PG_150.coneTipZ)*1000, 'linewidth',2, 'color', cols(6,:));

% line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
ylim([0 50]); ylabel('Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

%% Plot COLC radius
COLCRadiusF = figure; set(gcf, 'position', [-1919        -149        1920        1104])

subplot(4,nPlots,1); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.colcRadius*1000, 'linewidth',2, 'color', cols(1,:));

h2= plot(data_GRIN_radialBase.incidenceAngle, data_GRIN_radialBase.colcRadius*1000, 'linewidth',2, 'color', cols(4,:));

h3= plot(data_GRIN_linear.incidenceAngle, data_GRIN_linear.colcRadius*1000, 'linewidth',2, 'color', cols(5,:));

h4= plot(data_GRIN_both.incidenceAngle, data_GRIN_both.colcRadius*1000, 'linewidth',2, 'color', cols(6,:));

ylim([0 50]); ylabel('Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

subplot(4,nPlots,4); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.colcRadius*1000, 'linewidth',2, 'color', cols(1,:));

h2= plot(data_GRIN_radialTop_tipCorr.incidenceAngle, data_GRIN_radialTop_tipCorr.colcRadius*1000, 'linewidth',2, 'color', cols(4,:));

ylim([0 50]); ylabel('Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

subplot(4,nPlots,3); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.colcRadius*1000, 'linewidth',2, 'color', cols(2,:));

h2 = plot(data_GRIN_cylinder_p36.incidenceAngle, data_GRIN_cylinder_p36.colcRadius*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_cylinder_m36.incidenceAngle, data_GRIN_cylinder_m36.colcRadius*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_cylinder_shape.incidenceAngle, data_GRIN_cylinder_shape.colcRadius*1000, 'linewidth',2, 'color', cols(6,:));

ylim([0 50]); ylabel('Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

subplot(4,nPlots,2); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.colcRadius*1000, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_radialTop_CinC.incidenceAngle, data_GRIN_radialTop_CinC.colcRadius*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_radialTop_EC.incidenceAngle, data_GRIN_radialTop_EC.colcRadius*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_radialTop_CinC_EC.incidenceAngle, data_GRIN_radialTop_CinC_EC.colcRadius*1000, 'linewidth',2, 'color', cols(6,:));

ylim([0 50]); ylabel('Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

subplot(4,nPlots,8); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.colcRadius*1000, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.colcRadius*1000, 'linewidth',2, 'color', cols(2,:));

h3 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.colcRadius*1000, 'linewidth',2, 'color', cols(3,:));

h4 = plot(data_GRIN_cylinder_norm_max.incidenceAngle, data_GRIN_cylinder_norm_max.colcRadius*1000, 'linewidth',2, 'color', cols(4,:));

h5 = plot(data_GRIN_cylinder_norm_80.incidenceAngle, data_GRIN_cylinder_norm_80.colcRadius*1000, 'linewidth',2, 'color', cols(5,:));

ylim([0 50]); ylabel('Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

% Uniform
subplot(4,nPlots,5); hold on
h1 = plot(data_uniform_m4SD.incidenceAngle, data_uniform_m4SD.colcRadius*1000, 'linewidth',2, 'color', cols(4,:));

h2 = plot(data_uniform_m2SD.incidenceAngle, data_uniform_m2SD.colcRadius*1000, 'linewidth',2, 'color', cols(5,:));

h3 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.colcRadius*1000, 'linewidth',2, 'color', cols(3,:));

h4 = plot(data_uniform_p2SD.incidenceAngle, data_uniform_p2SD.colcRadius*1000, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_uniform_p4SD.incidenceAngle, data_uniform_p4SD.colcRadius*1000, 'linewidth',2, 'color', cols(7,:));

ylim([0 50]); ylabel('Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

subplot(4,nPlots,6); hold on
h1 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.colcRadius*1000, 'linewidth',2, 'color', cols(3,:));

h2 = plot(data_GRIN_uniform_CinC.incidenceAngle, data_GRIN_uniform_CinC.colcRadius*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_uniform_EC.incidenceAngle, data_GRIN_uniform_EC.colcRadius*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_uniform_CinC_EC.incidenceAngle, data_GRIN_uniform_CinC_EC.colcRadius*1000, 'linewidth',2, 'color', cols(6,:));

ylim([0 50]); ylabel('Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

% Varying intercone RI
subplot(4,nPlots,9); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.colcRadius*1000, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_radialTop_IC_142.incidenceAngle, data_GRIN_radialTop_IC_142.colcRadius*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_radialTop_IC_144.incidenceAngle, data_GRIN_radialTop_IC_144.colcRadius*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_radialTop_IC_146.incidenceAngle, data_GRIN_radialTop_IC_146.colcRadius*1000, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_GRIN_radialTop_IC_148.incidenceAngle, data_GRIN_radialTop_IC_148.colcRadius*1000, 'linewidth',2, 'color', cols(7,:));

ylim([0 50]); ylabel('Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

subplot(4,nPlots,10); hold on
h1 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.colcRadius*1000, 'linewidth',2, 'color', cols(3,:));

h2 = plot(data_uniform_0SD_IC_142.incidenceAngle, data_uniform_0SD_IC_142.colcRadius*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_uniform_0SD_IC_144.incidenceAngle, data_uniform_0SD_IC_144.colcRadius*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_uniform_0SD_IC_146.incidenceAngle, data_uniform_0SD_IC_146.colcRadius*1000, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_uniform_0SD_IC_148.incidenceAngle, data_uniform_0SD_IC_148.colcRadius*1000, 'linewidth',2, 'color', cols(7,:));

ylim([0 50]); ylabel('Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

subplot(4,nPlots,11); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.colcRadius*1000, 'linewidth',2, 'color', cols(2,:));

h2 = plot(data_GRIN_cylinder_IC_142.incidenceAngle, data_GRIN_cylinder_IC_142.colcRadius*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_cylinder_IC_144.incidenceAngle, data_GRIN_cylinder_IC_144.colcRadius*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_cylinder_IC_146.incidenceAngle, data_GRIN_cylinder_IC_146.colcRadius*1000, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_GRIN_cylinder_IC_148.incidenceAngle, data_GRIN_cylinder_IC_148.colcRadius*1000, 'linewidth',2, 'color', cols(7,:));

ylim([0 50]); ylabel('Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

% Varying pigment RI
subplot(4,nPlots,13); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.colcRadius*1000, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_radialTop_PG_140.incidenceAngle, data_GRIN_radialTop_PG_140.colcRadius*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_radialTop_PG_145.incidenceAngle, data_GRIN_radialTop_PG_145.colcRadius*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_radialTop_PG_150.incidenceAngle, data_GRIN_radialTop_PG_150.colcRadius*1000, 'linewidth',2, 'color', cols(6,:));

ylim([0 50]); ylabel('Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

subplot(4,nPlots,14); hold on
h1 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.colcRadius*1000, 'linewidth',2, 'color', cols(3,:));

h2 = plot(data_uniform_0SD_PG_140.incidenceAngle, data_uniform_0SD_PG_140.colcRadius*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_uniform_0SD_PG_145.incidenceAngle, data_uniform_0SD_PG_145.colcRadius*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_uniform_0SD_PG_150.incidenceAngle, data_uniform_0SD_PG_150.colcRadius*1000, 'linewidth',2, 'color', cols(6,:));

ylim([0 50]); ylabel('Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

subplot(4,nPlots,15); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.colcRadius*1000, 'linewidth',2, 'color', cols(2,:));

h2 = plot(data_GRIN_cylinder_PG_140.incidenceAngle, data_GRIN_cylinder_PG_140.colcRadius*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_cylinder_PG_145.incidenceAngle, data_GRIN_cylinder_PG_145.colcRadius*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_cylinder_PG_150.incidenceAngle, data_GRIN_cylinder_PG_150.colcRadius*1000, 'linewidth',2, 'color', cols(6,:));

ylim([0 50]); ylabel('Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

%% Plot COLC divergence
COLCDivergenceF = figure; set(gcf, 'position', [-1919        -149        1920        1104])

subplot(4,nPlots,1); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.beamDivergenceMax, 'linewidth',2, 'color', cols(1,:));

h2= plot(data_GRIN_radialBase.incidenceAngle, data_GRIN_radialBase.beamDivergenceMax, 'linewidth',2, 'color', cols(4,:));

h3= plot(data_GRIN_linear.incidenceAngle, data_GRIN_linear.beamDivergenceMax, 'linewidth',2, 'color', cols(5,:));

h4= plot(data_GRIN_both.incidenceAngle, data_GRIN_both.beamDivergenceMax, 'linewidth',2, 'color', cols(6,:));

ylim([0 90]); ylabel('Divergence (deg)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:30:90,'XTick',[0 5 10 15 20]);

subplot(4,nPlots,4); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.beamDivergenceMax, 'linewidth',2, 'color', cols(1,:));

h2= plot(data_GRIN_radialTop_tipCorr.incidenceAngle, data_GRIN_radialTop_tipCorr.beamDivergenceMax, 'linewidth',2, 'color', cols(4,:));

ylim([0 90]); ylabel('Divergence (deg)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:30:90,'XTick',[0 5 10 15 20]);

subplot(4,nPlots,3); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.beamDivergenceMax, 'linewidth',2, 'color', cols(2,:));

h2 = plot(data_GRIN_cylinder_p36.incidenceAngle, data_GRIN_cylinder_p36.beamDivergenceMax, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_cylinder_m36.incidenceAngle, data_GRIN_cylinder_m36.beamDivergenceMax, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_cylinder_shape.incidenceAngle, data_GRIN_cylinder_shape.beamDivergenceMax, 'linewidth',2, 'color', cols(6,:));

ylim([0 90]); ylabel('Divergence (deg)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:30:90,'XTick',[0 5 10 15 20]);

subplot(4,nPlots,2); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.beamDivergenceMax, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_radialTop_CinC.incidenceAngle, data_GRIN_radialTop_CinC.beamDivergenceMax, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_radialTop_EC.incidenceAngle, data_GRIN_radialTop_EC.beamDivergenceMax, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_radialTop_CinC_EC.incidenceAngle, data_GRIN_radialTop_CinC_EC.beamDivergenceMax, 'linewidth',2, 'color', cols(6,:));

ylim([0 90]); ylabel('Divergence (deg)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:30:90,'XTick',[0 5 10 15 20]);

subplot(4,nPlots,8); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.beamDivergenceMax, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.beamDivergenceMax, 'linewidth',2, 'color', cols(2,:));

h3 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.beamDivergenceMax, 'linewidth',2, 'color', cols(3,:));

h4 = plot(data_GRIN_cylinder_norm_max.incidenceAngle, data_GRIN_cylinder_norm_max.beamDivergenceMax, 'linewidth',2, 'color', cols(4,:));

h5 = plot(data_GRIN_cylinder_norm_80.incidenceAngle, data_GRIN_cylinder_norm_80.beamDivergenceMax, 'linewidth',2, 'color', cols(5,:));

ylim([0 90]); ylabel('Divergence (deg)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:30:90,'XTick',[0 5 10 15 20]);

% Uniform
subplot(4,nPlots,5); hold on
h1 = plot(data_uniform_m4SD.incidenceAngle, data_uniform_m4SD.beamDivergenceMax, 'linewidth',2, 'color', cols(4,:));

h2 = plot(data_uniform_m2SD.incidenceAngle, data_uniform_m2SD.beamDivergenceMax, 'linewidth',2, 'color', cols(5,:));

h3 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.beamDivergenceMax, 'linewidth',2, 'color', cols(3,:));

h4 = plot(data_uniform_p2SD.incidenceAngle, data_uniform_p2SD.beamDivergenceMax, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_uniform_p4SD.incidenceAngle, data_uniform_p4SD.beamDivergenceMax, 'linewidth',2, 'color', cols(7,:));

ylim([0 90]); ylabel('Divergence (deg)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:30:90,'XTick',[0 5 10 15 20]);

subplot(4,nPlots,6); hold on
h1 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.beamDivergenceMax, 'linewidth',2, 'color', cols(3,:));

h2 = plot(data_GRIN_uniform_CinC.incidenceAngle, data_GRIN_uniform_CinC.beamDivergenceMax, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_uniform_EC.incidenceAngle, data_GRIN_uniform_EC.beamDivergenceMax, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_uniform_CinC_EC.incidenceAngle, data_GRIN_uniform_CinC_EC.beamDivergenceMax, 'linewidth',2, 'color', cols(6,:));

ylim([0 90]); ylabel('Divergence (deg)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:30:90,'XTick',[0 5 10 15 20]);

% Varying intercone RI
subplot(4,nPlots,9); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.beamDivergenceMax, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_radialTop_IC_142.incidenceAngle, data_GRIN_radialTop_IC_142.beamDivergenceMax, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_radialTop_IC_144.incidenceAngle, data_GRIN_radialTop_IC_144.beamDivergenceMax, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_radialTop_IC_146.incidenceAngle, data_GRIN_radialTop_IC_146.beamDivergenceMax, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_GRIN_radialTop_IC_148.incidenceAngle, data_GRIN_radialTop_IC_148.beamDivergenceMax, 'linewidth',2, 'color', cols(7,:));

ylim([0 90]); ylabel('Divergence (deg)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:30:90,'XTick',[0 5 10 15 20]);

subplot(4,nPlots,10); hold on
h1 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.beamDivergenceMax, 'linewidth',2, 'color', cols(3,:));

h2 = plot(data_uniform_0SD_IC_142.incidenceAngle, data_uniform_0SD_IC_142.beamDivergenceMax, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_uniform_0SD_IC_144.incidenceAngle, data_uniform_0SD_IC_144.beamDivergenceMax, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_uniform_0SD_IC_146.incidenceAngle, data_uniform_0SD_IC_146.beamDivergenceMax, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_uniform_0SD_IC_148.incidenceAngle, data_uniform_0SD_IC_148.beamDivergenceMax, 'linewidth',2, 'color', cols(7,:));

ylim([0 90]); ylabel('Divergence (deg)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:30:90,'XTick',[0 5 10 15 20]);

subplot(4,nPlots,11); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.beamDivergenceMax, 'linewidth',2, 'color', cols(2,:));

h2 = plot(data_GRIN_cylinder_IC_142.incidenceAngle, data_GRIN_cylinder_IC_142.beamDivergenceMax, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_cylinder_IC_144.incidenceAngle, data_GRIN_cylinder_IC_144.beamDivergenceMax, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_cylinder_IC_146.incidenceAngle, data_GRIN_cylinder_IC_146.beamDivergenceMax, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_GRIN_cylinder_IC_148.incidenceAngle, data_GRIN_cylinder_IC_148.beamDivergenceMax, 'linewidth',2, 'color', cols(7,:));

ylim([0 90]); ylabel('Divergence (deg)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:30:90,'XTick',[0 5 10 15 20]);

% Varying pigment RI
subplot(4,nPlots,13); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.beamDivergenceMax, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_radialTop_PG_140.incidenceAngle, data_GRIN_radialTop_PG_140.beamDivergenceMax, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_radialTop_PG_145.incidenceAngle, data_GRIN_radialTop_PG_145.beamDivergenceMax, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_radialTop_PG_150.incidenceAngle, data_GRIN_radialTop_PG_150.beamDivergenceMax, 'linewidth',2, 'color', cols(6,:));

ylim([0 90]); ylabel('Divergence (deg)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:30:90,'XTick',[0 5 10 15 20]);

subplot(4,nPlots,14); hold on
h1 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.beamDivergenceMax, 'linewidth',2, 'color', cols(3,:));

h2 = plot(data_uniform_0SD_PG_140.incidenceAngle, data_uniform_0SD_PG_140.beamDivergenceMax, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_uniform_0SD_PG_145.incidenceAngle, data_uniform_0SD_PG_145.beamDivergenceMax, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_uniform_0SD_PG_150.incidenceAngle, data_uniform_0SD_PG_150.beamDivergenceMax, 'linewidth',2, 'color', cols(6,:));

ylim([0 90]); ylabel('Divergence (deg)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:30:90,'XTick',[0 5 10 15 20]);

subplot(4,nPlots,15); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.beamDivergenceMax, 'linewidth',2, 'color', cols(2,:));

h2 = plot(data_GRIN_cylinder_PG_140.incidenceAngle, data_GRIN_cylinder_PG_140.beamDivergenceMax, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_cylinder_PG_145.incidenceAngle, data_GRIN_cylinder_PG_145.beamDivergenceMax, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_cylinder_PG_150.incidenceAngle, data_GRIN_cylinder_PG_150.beamDivergenceMax, 'linewidth',2, 'color', cols(6,:));

ylim([0 90]); ylabel('Divergence (deg)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:30:90,'XTick',[0 5 10 15 20]);

%% 
if saveFigures
    currentDirectory = pwd; 
    
    cd(saveFolder)
    exportgraphics(legendsF,'Legends_Summary.pdf','ContentType','vector') 
    
    exportgraphics(acceptanceFNight,'AcceptanceFunction_Night_Summary.pdf','ContentType','vector') 
    exportgraphics(acceptanceFDay,'AcceptanceFunction_Day_Summary.pdf','ContentType','vector') 
    exportgraphics(acceptanceFCOLC,'AcceptanceFunction_COLC_Summary.pdf','ContentType','vector') 
    
    exportgraphics(TIRF,'COLCTIRRatio_Summary.pdf','ContentType','vector') 
    
    exportgraphics(COLCHeightF,'COLCHeight_Summary.pdf','ContentType','vector') 
    exportgraphics(COLCRadiusF,'COLCRadius_Summary.pdf','ContentType','vector') 
    exportgraphics(COLCDivergenceF,'COLCDivergence_Summary.pdf','ContentType','vector') 
    
    cd(currentDirectory)
end