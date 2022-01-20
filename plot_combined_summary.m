%%% Need to set figure size to get nicely shaped axes...

%%% Check cylinder cylinder for COLC, X Focus and Y Focus height (and radius?)

clear; close all; clc

saveFigures = 0;
saveFolder = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Figures/';

% Load reference acceptance function
acceptanceAngle = readmatrix('/Users/gavintaylor/Documents/Matlab/Git_versioned_April_27/Crab_cone_git/Data/Acceptacne Function.csv');

nPlots = 4;
cols = lines(7);
outlineCol = [0.5 0.5 0.5];

% Load data from varying RI series
cd('/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Data/Varying RI Profile/')
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
cd('/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Data/Extra internal structures/Radial RI')
data_GRIN_radialTop_CinC = load('Cone_CinC_1000_nm_Cone_0_SD_GRIN_radialTop_SIMDATA.mat');
data_GRIN_radialTop_EC = load('Cone_EC_1000_nm_Cone_0_SD_GRIN_radialTop_SIMDATA.mat');
data_GRIN_radialTop_CinC_EC = load('Cone_CinC_EC_1000_nm_Cone_0_SD_GRIN_radialTop_SIMDATA.mat');

% Load data from varying SD (unfirom RI) series
cd('/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Data/Varying Cone SD')
data_uniform_m4SD = load('Cone_1000_nm_Cone_-4_SD_Uniform_1.52_SIMDATA.mat');
data_uniform_m2SD = load('Cone_1000_nm_Cone_-2_SD_Uniform_1.52_SIMDATA.mat');
data_uniform_0SD = load('Cone_1000_nm_Cone_0_SD_Uniform_1.52_SIMDATA.mat');
data_uniform_p2SD = load('Cone_1000_nm_Cone_2_SD_Uniform_1.52_SIMDATA.mat');
data_uniform_p4SD = load('Cone_1000_nm_Cone_4_SD_Uniform_1.52_SIMDATA.mat');

% Load data from internal structures (w/ uniform RI)
cd('/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Data/Extra internal structures/Uniform RI')
data_GRIN_uniform_CinC = load('Cone_CinC_1000_nm_Cone_0_SD_Uniform_1.52_SIMDATA.mat');
data_GRIN_uniform_EC = load('Cone_EC_1000_nm_Cone_0_SD_Uniform_1.52_SIMDATA.mat');
data_GRIN_uniform_CinC_EC = load('Cone_CinC_EC_1000_nm_Cone_0_SD_Uniform_1.52_SIMDATA.mat');

% Load data from the three varying intercone series.
cd('/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Data/Variable Intercone/Radial RI')
data_GRIN_radialTop_IC_142 = load('Cone_1000_nm_Cone_0_SD_GRIN_radialTop_IC_142_SIMDATA.mat');
data_GRIN_radialTop_IC_144 = load('Cone_1000_nm_Cone_0_SD_GRIN_radialTop_IC_144_SIMDATA.mat');
data_GRIN_radialTop_IC_146 = load('Cone_1000_nm_Cone_0_SD_GRIN_radialTop_IC_146_SIMDATA.mat');
data_GRIN_radialTop_IC_148 = load('Cone_1000_nm_Cone_0_SD_GRIN_radialTop_IC_148_SIMDATA.mat');

cd('/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Data/Variable Intercone/Uniform RI')
data_uniform_0SD_IC_142 = load('Cone_1000_nm_Cone_0_SD_Uniform_1.52_IC_142_SIMDATA.mat');
data_uniform_0SD_IC_144 = load('Cone_1000_nm_Cone_0_SD_Uniform_1.52_IC_144_SIMDATA.mat');
data_uniform_0SD_IC_146 = load('Cone_1000_nm_Cone_0_SD_Uniform_1.52_IC_146_SIMDATA.mat');
data_uniform_0SD_IC_148 = load('Cone_1000_nm_Cone_0_SD_Uniform_1.52_IC_148_SIMDATA.mat');

cd('/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Data/Variable Intercone/Theoretical RI')
data_GRIN_cylinder_IC_142 = load('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_IC_142_SIMDATA.mat');
data_GRIN_cylinder_IC_144 = load('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_IC_144_SIMDATA.mat');
data_GRIN_cylinder_IC_146 = load('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_IC_146_SIMDATA.mat');
data_GRIN_cylinder_IC_148 = load('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_IC_148_SIMDATA.mat');

% process acceptance angle from paper
% Fix zero at zero
acceptanceAngle(abs(acceptanceAngle(:,1)) < 0.1,1) = 0;

negInds = find(acceptanceAngle(:,1) <= 0);
posInds = find(acceptanceAngle(:,1) >= 0);

angleSteps = 0:10;
negInterp = interp1(-acceptanceAngle(negInds,1), acceptanceAngle(negInds,2), angleSteps, 'linear', 'extrap');
posInterp = interp1(acceptanceAngle(posInds,1), acceptanceAngle(posInds,2), angleSteps, 'linear', 'extrap');

interpAcceptance = (negInterp+posInterp)/2;

lowInds = find(interpAcceptance > 0.5);
if ~isempty(lowInds)
    % Take last
    lowInds = lowInds(end);
    
    slope = (interpAcceptance(lowInds + 1) - interpAcceptance(lowInds))/...
        (angleSteps(lowInds + 1) - angleSteps(lowInds));
    acceptanceAngle = angleSteps(lowInds) + (0.5 - interpAcceptance(lowInds))/slope;
end

% figure; hold on
% plot(angleSteps, interpAcceptance, 'k-', 'linewidth',2)
% plot(acceptanceAngle, 0.5, 'kd','markersize',10, 'linewidth',2)


%% Plot acceptance angle
acceptanceF = figure; set(gcf, 'position', [1 85 1680 870])

subplot(3,nPlots,1); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.acceptancePercentage, 'linewidth',2, 'color', cols(1,:));
plot(data_GRIN_radialTop.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(1,:))

h2= plot(data_GRIN_radialBase.incidenceAngle, data_GRIN_radialBase.acceptancePercentage, 'linewidth',2, 'color', cols(4,:));
plot(data_GRIN_radialBase.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(4,:))

h3= plot(data_GRIN_linear.incidenceAngle, data_GRIN_linear.acceptancePercentage, 'linewidth',2, 'color', cols(5,:));
plot(data_GRIN_linear.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(5,:))

h4= plot(data_GRIN_both.incidenceAngle, data_GRIN_both.acceptancePercentage, 'linewidth',2, 'color', cols(6,:));
plot(data_GRIN_both.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(6,:))

legend([h1 h2 h3 h4], {'Radial-top','Radial-base','Linear','Radial+linear'},'Location','southwest')
ylabel('% entering receptor'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);

subplot(3,nPlots,4); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.acceptancePercentage, 'linewidth',2, 'color', cols(1,:));
plot(data_GRIN_radialTop.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(1,:))

h2= plot(data_GRIN_radialTop_tipCorr.incidenceAngle, data_GRIN_radialTop_tipCorr.acceptancePercentage, 'linewidth',2, 'color', cols(4,:));
plot(data_GRIN_radialTop_tipCorr.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(4,:))

legend([h1 h2], {'Radial-top','W/ Tip Corr.'},'Location','southwest')
ylabel('% entering receptor'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);

subplot(3,nPlots,3); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.acceptancePercentage, 'linewidth',2, 'color', cols(2,:));
plot(data_GRIN_cylinder.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(2,:))

h2 = plot(data_GRIN_cylinder_p36.incidenceAngle, data_GRIN_cylinder_p36.acceptancePercentage, 'linewidth',2, 'color', cols(4,:));
plot(data_GRIN_cylinder_p36.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(4,:))

h3 = plot(data_GRIN_cylinder_m36.incidenceAngle, data_GRIN_cylinder_m36.acceptancePercentage, 'linewidth',2, 'color', cols(5,:));
plot(data_GRIN_cylinder_m36.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(5,:))

h4 = plot(data_GRIN_cylinder_shape.incidenceAngle, data_GRIN_cylinder_shape.acceptancePercentage, 'linewidth',2, 'color', cols(6,:));
plot(data_GRIN_cylinder_shape.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(6,:))

lh = legend([h1 h2 h3 h4], {'Theory-equal', 'Theory-plus', 'Theory-minus', 'Cylinder'},'Location','southwest');
ylabel('% entering receptor'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);
set(lh,'Position',[0.5774    0.5470    0.0905    0.1075])

subplot(3,nPlots,2); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.acceptancePercentage, 'linewidth',2, 'color', cols(1,:));
plot(data_GRIN_radialTop.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(1,:))

h2 = plot(data_GRIN_radialTop_CinC.incidenceAngle, data_GRIN_radialTop_CinC.acceptancePercentage, 'linewidth',2, 'color', cols(4,:));
plot(data_GRIN_radialTop_CinC.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(4,:))

h3 = plot(data_GRIN_radialTop_EC.incidenceAngle, data_GRIN_radialTop_EC.acceptancePercentage, 'linewidth',2, 'color', cols(5,:));
plot(data_GRIN_radialTop_EC.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(5,:))

h4 = plot(data_GRIN_radialTop_CinC_EC.incidenceAngle, data_GRIN_radialTop_CinC_EC.acceptancePercentage, 'linewidth',2, 'color', cols(6,:));
plot(data_GRIN_radialTop_CinC_EC.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(6,:))

legend([h1 h2 h3 h4], {'Radial-top', 'CinC', 'EC', 'CinC & EC'},'Location','southwest')
ylabel('% entering receptor'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);

subplot(3,nPlots,8); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.acceptancePercentage, 'linewidth',2, 'color', cols(1,:));
plot(data_GRIN_radialTop.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(1,:))

h2 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.acceptancePercentage, 'linewidth',2, 'color', cols(2,:));
plot(data_GRIN_cylinder.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(2,:))

h3 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.acceptancePercentage, 'linewidth',2, 'color', cols(3,:));
plot(data_uniform_0SD.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(3,:))

h4 = plot(data_GRIN_cylinder_norm_max.incidenceAngle, data_GRIN_cylinder_norm_max.acceptancePercentage, 'linewidth',2, 'color', cols(4,:));
plot(data_GRIN_cylinder_norm_max.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(4,:))

h5 = plot(data_GRIN_cylinder_norm_80.incidenceAngle, data_GRIN_cylinder_norm_80.acceptancePercentage, 'linewidth',2, 'color', cols(5,:));
plot(data_GRIN_cylinder_norm_80.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(5,:))

lh2 = legend([h1 h2 h3 h4 h5], {'Radial-top','Theory-equal', 'SD 0', 'Normalised-max', 'Normalised-80'},'Location','southwest');
ylabel('% entering receptor'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);
set(lh2,'Position',[0.7756    0.2459    0.1042    0.1075])

%uniform RI
subplot(3,nPlots,5); hold on
h1 = plot(data_uniform_m4SD.incidenceAngle, data_uniform_m4SD.acceptancePercentage, 'linewidth',2, 'color', cols(4,:));
plot(data_uniform_m4SD.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(4,:))

h2 = plot(data_uniform_m2SD.incidenceAngle, data_uniform_m2SD.acceptancePercentage, 'linewidth',2, 'color', cols(5,:));
plot(data_uniform_m2SD.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(5,:))

h3 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.acceptancePercentage, 'linewidth',2, 'color', cols(3,:));
plot(data_uniform_0SD.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(3,:))

h4 = plot(data_uniform_p2SD.incidenceAngle, data_uniform_p2SD.acceptancePercentage, 'linewidth',2, 'color', cols(6,:));
plot(data_uniform_p2SD.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(6,:))

h5 = plot(data_uniform_p4SD.incidenceAngle, data_uniform_p4SD.acceptancePercentage, 'linewidth',2, 'color', cols(7,:));
plot(data_uniform_p4SD.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(7,:))

legend([h1 h2 h3 h4 h5], {'SD -4', 'SD -2', 'SD 0', 'SD +2', 'SD +4'},'Location','southwest')
ylabel('% entering receptor'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);

subplot(3,nPlots,6); hold on
h1 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.acceptancePercentage, 'linewidth',2, 'color', cols(3,:));
plot(data_uniform_0SD.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(3,:))

h2 = plot(data_GRIN_uniform_CinC.incidenceAngle, data_GRIN_uniform_CinC.acceptancePercentage, 'linewidth',2, 'color', cols(4,:));
plot(data_GRIN_uniform_CinC.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(4,:))

h3 = plot(data_GRIN_uniform_EC.incidenceAngle, data_GRIN_uniform_EC.acceptancePercentage, 'linewidth',2, 'color', cols(5,:));
plot(data_GRIN_uniform_EC.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(5,:))

h4 = plot(data_GRIN_uniform_CinC_EC.incidenceAngle, data_GRIN_uniform_CinC_EC.acceptancePercentage, 'linewidth',2, 'color', cols(6,:));
plot(data_GRIN_uniform_CinC_EC.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(6,:))

legend([h1 h2 h3 h4], {'SD 0', 'CinC', 'EC', 'CinC & EC'},'Location','southwest')
ylabel('% entering receptor'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);

% Varying intercone RI
subplot(3,nPlots,9); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.acceptancePercentage, 'linewidth',2, 'color', cols(1,:));
plot(data_GRIN_radialTop.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(1,:))

h2 = plot(data_GRIN_radialTop_IC_142.incidenceAngle, data_GRIN_radialTop_IC_142.acceptancePercentage, 'linewidth',2, 'color', cols(4,:));
plot(data_GRIN_radialTop_IC_142.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(4,:))

h3 = plot(data_GRIN_radialTop_IC_144.incidenceAngle, data_GRIN_radialTop_IC_144.acceptancePercentage, 'linewidth',2, 'color', cols(5,:));
plot(data_GRIN_radialTop_IC_144.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(5,:))

h4 = plot(data_GRIN_radialTop_IC_146.incidenceAngle, data_GRIN_radialTop_IC_146.acceptancePercentage, 'linewidth',2, 'color', cols(6,:));
plot(data_GRIN_radialTop_IC_146.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(6,:))

h5 = plot(data_GRIN_radialTop_IC_148.incidenceAngle, data_GRIN_radialTop_IC_148.acceptancePercentage, 'linewidth',2, 'color', cols(7,:));
plot(data_GRIN_radialTop_IC_148.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(7,:))

legend([h1 h2 h3 h4 h5], {'Radial-top', 'IC-1.42', 'IC-1.44', 'IC-1.46', 'IC-1.48'},'Location','southwest')
ylabel('% entering receptor'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);

subplot(3,nPlots,10); hold on
h1 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.acceptancePercentage, 'linewidth',2, 'color', cols(3,:));
plot(data_uniform_0SD.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(3,:))

h2 = plot(data_uniform_0SD_IC_142.incidenceAngle, data_uniform_0SD_IC_142.acceptancePercentage, 'linewidth',2, 'color', cols(4,:));
plot(data_uniform_0SD_IC_142.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(4,:))

h3 = plot(data_uniform_0SD_IC_144.incidenceAngle, data_uniform_0SD_IC_144.acceptancePercentage, 'linewidth',2, 'color', cols(5,:));
plot(data_uniform_0SD_IC_144.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(5,:))

h4 = plot(data_uniform_0SD_IC_146.incidenceAngle, data_uniform_0SD_IC_146.acceptancePercentage, 'linewidth',2, 'color', cols(6,:));
plot(data_uniform_0SD_IC_146.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(6,:))

h5 = plot(data_uniform_0SD_IC_148.incidenceAngle, data_uniform_0SD_IC_148.acceptancePercentage, 'linewidth',2, 'color', cols(7,:));
plot(data_uniform_0SD_IC_148.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(7,:))

legend([h1 h2 h3 h4 h5], {'SD 0', 'IC-1.42', 'IC-1.44', 'IC-1.46', 'IC-1.48'},'Location','southwest')
ylabel('% entering receptor'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);

subplot(3,nPlots,11); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.acceptancePercentage, 'linewidth',2, 'color', cols(2,:));
plot(data_GRIN_cylinder.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(2,:))

h2 = plot(data_GRIN_cylinder_IC_142.incidenceAngle, data_GRIN_cylinder_IC_142.acceptancePercentage, 'linewidth',2, 'color', cols(4,:));
plot(data_GRIN_cylinder_IC_142.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(4,:))

h3 = plot(data_GRIN_cylinder_IC_144.incidenceAngle, data_GRIN_cylinder_IC_144.acceptancePercentage, 'linewidth',2, 'color', cols(5,:));
plot(data_GRIN_cylinder_IC_144.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(5,:))

h4 = plot(data_GRIN_cylinder_IC_146.incidenceAngle, data_GRIN_cylinder_IC_146.acceptancePercentage, 'linewidth',2, 'color', cols(6,:));
plot(data_GRIN_cylinder_IC_146.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(6,:))

h5 = plot(data_GRIN_cylinder_IC_148.incidenceAngle, data_GRIN_cylinder_IC_148.acceptancePercentage, 'linewidth',2, 'color', cols(7,:));
plot(data_GRIN_cylinder_IC_148.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(7,:))

lh3 =legend([h1 h2 h3 h4 h5], {'Theory-equal', 'IC-1.42', 'IC-1.44', 'IC-1.46', 'IC-1.48'},'Location','northwest');
ylabel('% entering receptor'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);
set(lh3,'Position',[0.6083    0.3001    0.0881    0.1333])
%% Plot TIR
TIRF = figure; set(gcf, 'position', [1 85 1680 870])

subplot(3,nPlots,1); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.TIRPercengtage, 'linewidth',2, 'color', cols(1,:));

h2= plot(data_GRIN_radialBase.incidenceAngle, data_GRIN_radialBase.TIRPercengtage, 'linewidth',2, 'color', cols(4,:));

h3= plot(data_GRIN_linear.incidenceAngle, data_GRIN_linear.TIRPercengtage, 'linewidth',2, 'color', cols(5,:));

h4= plot(data_GRIN_both.incidenceAngle, data_GRIN_both.TIRPercengtage, 'linewidth',2, 'color', cols(6,:));

legend([h1 h2 h3 h4], {'Radial-top','Radial-base','Linear','Radial+linear'},'Location','southeast')
ylim([0 1]); xlim([0 20]); ylabel('% TIR in COLC'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);

subplot(3,nPlots,4); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.TIRPercengtage, 'linewidth',2, 'color', cols(1,:));

h2= plot(data_GRIN_radialTop_tipCorr.incidenceAngle, data_GRIN_radialTop_tipCorr.TIRPercengtage, 'linewidth',2, 'color', cols(4,:));

legend([h1 h2], {'Radial-top','W/ Tip Corr.'},'Location','southeast')
ylim([0 1]); xlim([0 20]); ylabel('% TIR in COLC'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);

subplot(3,nPlots,3); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.TIRPercengtage, 'linewidth',2, 'color', cols(2,:));

h2 = plot(data_GRIN_cylinder_p36.incidenceAngle, data_GRIN_cylinder_p36.TIRPercengtage, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_cylinder_m36.incidenceAngle, data_GRIN_cylinder_m36.TIRPercengtage, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_cylinder_shape.incidenceAngle, data_GRIN_cylinder_shape.TIRPercengtage, 'linewidth',2, 'color', cols(6,:));

lh = legend([h1 h2 h3 h4], {'Theory-equal', 'Theory-plus', 'Theory-minus', 'Cylinder'},'Location','southeast');
ylim([0 1]); xlim([0 20]); ylabel('% TIR in COLC'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);
set(lh,'Position',[0.5774    0.5470    0.0905    0.1075])

subplot(3,nPlots,2); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.TIRPercengtage, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_radialTop_CinC.incidenceAngle, data_GRIN_radialTop_CinC.TIRPercengtage, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_radialTop_EC.incidenceAngle, data_GRIN_radialTop_EC.TIRPercengtage, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_radialTop_CinC_EC.incidenceAngle, data_GRIN_radialTop_CinC_EC.TIRPercengtage, 'linewidth',2, 'color', cols(6,:));

legend([h1 h2 h3 h4], {'Radial-top', 'CinC', 'EC', 'CinC & EC'},'Location','southeast')
ylim([0 1]); xlim([0 20]); ylabel('% TIR in COLC'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);

subplot(3,nPlots,8); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.TIRPercengtage, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.TIRPercengtage, 'linewidth',2, 'color', cols(2,:));

h3 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.TIRPercengtage, 'linewidth',2, 'color', cols(3,:));

h4 = plot(data_GRIN_cylinder_norm_max.incidenceAngle, data_GRIN_cylinder_norm_max.TIRPercengtage, 'linewidth',2, 'color', cols(4,:));

h5 = plot(data_GRIN_cylinder_norm_80.incidenceAngle, data_GRIN_cylinder_norm_80.TIRPercengtage, 'linewidth',2, 'color', cols(5,:));

lh2 = legend([h1 h2 h3 h4 h5], {'Radial-top','Theory-equal', 'SD 0', 'Normalised-max', 'Normalised-80'},'Location','southeast');
ylim([0 1]); xlim([0 20]); ylabel('% TIR in COLC'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);
set(lh2,'Position',[0.7756    0.2459    0.1042    0.1075])

% Uniform RI
subplot(3,nPlots,5); hold on
h1 = plot(data_uniform_m4SD.incidenceAngle, data_uniform_m4SD.TIRPercengtage, 'linewidth',2, 'color', cols(4,:));

h2 = plot(data_uniform_m2SD.incidenceAngle, data_uniform_m2SD.TIRPercengtage, 'linewidth',2, 'color', cols(5,:));

h3 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.TIRPercengtage, 'linewidth',2, 'color', cols(3,:));

h4 = plot(data_uniform_p2SD.incidenceAngle, data_uniform_p2SD.TIRPercengtage, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_uniform_p4SD.incidenceAngle, data_uniform_p4SD.TIRPercengtage, 'linewidth',2, 'color', cols(7,:));

legend([h1 h2 h3 h4 h5], {'SD -4', 'SD -2', 'SD 0', 'SD +2', 'SD +4'},'Location','southeast')
ylim([0 1]); xlim([0 20]); ylabel('% TIR in COLC'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);

subplot(3,nPlots,6); hold on
h1 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.TIRPercengtage, 'linewidth',2, 'color', cols(3,:));

h2 = plot(data_GRIN_uniform_CinC.incidenceAngle, data_GRIN_uniform_CinC.TIRPercengtage, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_uniform_EC.incidenceAngle, data_GRIN_uniform_EC.TIRPercengtage, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_uniform_CinC_EC.incidenceAngle, data_GRIN_uniform_CinC_EC.TIRPercengtage, 'linewidth',2, 'color', cols(6,:));

legend([h1 h2 h3 h4], {'SD 0', 'CinC', 'EC', 'CinC & EC'},'Location','southeast')
ylim([0 1]); xlim([0 20]); ylabel('% TIR in COLC'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);

% Varying intercone RI
subplot(3,nPlots,9); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.TIRPercengtage, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_radialTop_IC_142.incidenceAngle, data_GRIN_radialTop_IC_142.TIRPercengtage, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_radialTop_IC_144.incidenceAngle, data_GRIN_radialTop_IC_144.TIRPercengtage, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_radialTop_IC_146.incidenceAngle, data_GRIN_radialTop_IC_146.TIRPercengtage, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_GRIN_radialTop_IC_148.incidenceAngle, data_GRIN_radialTop_IC_148.TIRPercengtage, 'linewidth',2, 'color', cols(7,:));

legend([h1 h2 h3 h4 h5], {'Radial-top', 'IC-1.42', 'IC-1.44', 'IC-1.46', 'IC-1.48'},'Location','southeast')
ylim([0 1]); xlim([0 20]); ylabel('% TIR in COLC'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);

subplot(3,nPlots,10); hold on
h1 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.TIRPercengtage, 'linewidth',2, 'color', cols(3,:));

h2 = plot(data_uniform_0SD_IC_142.incidenceAngle, data_uniform_0SD_IC_142.TIRPercengtage, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_uniform_0SD_IC_144.incidenceAngle, data_uniform_0SD_IC_144.TIRPercengtage, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_uniform_0SD_IC_146.incidenceAngle, data_uniform_0SD_IC_146.TIRPercengtage, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_uniform_0SD_IC_148.incidenceAngle, data_uniform_0SD_IC_148.TIRPercengtage, 'linewidth',2, 'color', cols(7,:));

legend([h1 h2 h3 h4 h5], {'SD 0', 'IC-1.42', 'IC-1.44', 'IC-1.46', 'IC-1.48'},'Location','southeast')
ylim([0 1]); xlim([0 20]); ylabel('% TIR in COLC'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);

subplot(3,nPlots,11); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.TIRPercengtage, 'linewidth',2, 'color', cols(2,:));

h2 = plot(data_GRIN_cylinder_IC_142.incidenceAngle, data_GRIN_cylinder_IC_142.TIRPercengtage, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_cylinder_IC_144.incidenceAngle, data_GRIN_cylinder_IC_144.TIRPercengtage, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_cylinder_IC_146.incidenceAngle, data_GRIN_cylinder_IC_146.TIRPercengtage, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_GRIN_cylinder_IC_148.incidenceAngle, data_GRIN_cylinder_IC_148.TIRPercengtage, 'linewidth',2, 'color', cols(7,:));

lh3 = legend([h1 h2 h3 h4 h5], {'Theory-equal', 'IC-1.42', 'IC-1.44', 'IC-1.46', 'IC-1.48'},'Location','southeast');
ylim([0 1]); xlim([0 20]); ylabel('% TIR in COLC'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);
set(lh3,'Position',[0.6083    0.3001    0.0881    0.1333])
%% Plot COLC height
COLCHeightF = figure; set(gcf, 'position', [1 85 1680 870])

subplot(3,nPlots,1); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, (data_GRIN_radialTop.colcHeight - data_GRIN_radialTop.coneTipZ)*1000, 'linewidth',2, 'color', cols(1,:));

h2= plot(data_GRIN_radialBase.incidenceAngle, (data_GRIN_radialBase.colcHeight- data_GRIN_radialBase.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h3= plot(data_GRIN_linear.incidenceAngle, (data_GRIN_linear.colcHeight- data_GRIN_linear.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

h4= plot(data_GRIN_both.incidenceAngle, (data_GRIN_both.colcHeight- data_GRIN_both.coneTipZ)*1000, 'linewidth',2, 'color', cols(6,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
legend([h1 h2 h3 h4], {'Radial-top','Radial-base','Linear','Radial+linear'},'Location','southwest')
ylim([-150 50]); ylabel('COLC Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -150:50:50,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,4); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, (data_GRIN_radialTop.colcHeight- data_GRIN_radialTop.coneTipZ)*1000, 'linewidth',2, 'color', cols(1,:));

h2= plot(data_GRIN_radialTop_tipCorr.incidenceAngle, (data_GRIN_radialTop_tipCorr.colcHeight- data_GRIN_radialTop_tipCorr.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
legend([h1 h2], {'Radial-top','W/ Tip Corr.'},'Location','southwest')
ylim([-150 50]); ylabel('COLC Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -150:50:50,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,3); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, (data_GRIN_cylinder.colcHeight- data_GRIN_cylinder.coneTipZ)*1000, 'linewidth',2, 'color', cols(2,:));

h2 = plot(data_GRIN_cylinder_p36.incidenceAngle, (data_GRIN_cylinder_p36.colcHeight- data_GRIN_cylinder_p36.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_cylinder_m36.incidenceAngle, (data_GRIN_cylinder_m36.colcHeight- data_GRIN_cylinder_m36.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_cylinder_shape.incidenceAngle, (data_GRIN_cylinder_shape.colcHeight- data_GRIN_cylinder_shape.coneTipZ)*1000, 'linewidth',2, 'color', cols(6,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
lh = legend([h1 h2 h3 h4], {'Theory-equal', 'Theory-plus', 'Theory-minus', 'Cylinder'},'Location','southwest');
ylim([-150 50]); ylabel('COLC Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -150:50:50,'XTick',[0 5 10 15 20]);
set(lh,'Position',[0.5774    0.5470    0.0905    0.1075])

subplot(3,nPlots,2); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, (data_GRIN_radialTop.colcHeight- data_GRIN_radialTop.coneTipZ)*1000, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_radialTop_CinC.incidenceAngle, (data_GRIN_radialTop_CinC.colcHeight- data_GRIN_radialTop_CinC.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_radialTop_EC.incidenceAngle, (data_GRIN_radialTop_EC.colcHeight- data_GRIN_radialTop_EC.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_radialTop_CinC_EC.incidenceAngle, (data_GRIN_radialTop_CinC_EC.colcHeight- data_GRIN_radialTop_CinC_EC.coneTipZ)*1000, 'linewidth',2, 'color', cols(6,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
legend([h1 h2 h3 h4], {'Radial-top', 'CinC', 'EC', 'CinC & EC'},'Location','southwest')
ylim([-150 50]); ylabel('COLC Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -150:50:50,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,8); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, (data_GRIN_radialTop.colcHeight- data_GRIN_radialTop.coneTipZ)*1000, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_cylinder.incidenceAngle, (data_GRIN_cylinder.colcHeight- data_GRIN_cylinder.coneTipZ)*1000, 'linewidth',2, 'color', cols(2,:));

h3 = plot(data_uniform_0SD.incidenceAngle, (data_uniform_0SD.colcHeight- data_uniform_0SD.coneTipZ)*1000, 'linewidth',2, 'color', cols(3,:));

h4 = plot(data_GRIN_cylinder_norm_max.incidenceAngle, (data_GRIN_cylinder_norm_max.colcHeight- data_GRIN_cylinder_norm_max.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h5 = plot(data_GRIN_cylinder_norm_80.incidenceAngle, (data_GRIN_cylinder_norm_80.colcHeight- data_GRIN_cylinder_norm_80.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
lh2 = legend([h1 h2 h3 h4 h5], {'Radial-top','Theory-equal', 'SD 0', 'Normalised-max', 'Normalised-80'},'Location','southwest');
ylim([-150 50]); ylabel('COLC Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -150:50:50,'XTick',[0 5 10 15 20]);
set(lh2,'Position',[0.7756    0.2459    0.1042    0.1075])

% Uniform RI
subplot(3,nPlots,5); hold on
h1 = plot(data_uniform_m4SD.incidenceAngle, (data_uniform_m4SD.colcHeight- data_uniform_m4SD.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h2 = plot(data_uniform_m2SD.incidenceAngle, (data_uniform_m2SD.colcHeight- data_uniform_m2SD.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

h3 = plot(data_uniform_0SD.incidenceAngle, (data_uniform_0SD.colcHeight- data_uniform_0SD.coneTipZ)*1000, 'linewidth',2, 'color', cols(3,:));

h4 = plot(data_uniform_p2SD.incidenceAngle, (data_uniform_p2SD.colcHeight- data_uniform_p2SD.coneTipZ)*1000, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_uniform_p4SD.incidenceAngle, (data_uniform_p4SD.colcHeight- data_uniform_p4SD.coneTipZ)*1000, 'linewidth',2, 'color', cols(7,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
legend([h1 h2 h3 h4 h5], {'SD -4', 'SD -2', 'SD 0', 'SD +2', 'SD +4'},'Location','southwest')
ylim([-150 50]); ylabel('COLC Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -150:50:50,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,6); hold on
h1 = plot(data_uniform_0SD.incidenceAngle, (data_uniform_0SD.colcHeight- data_uniform_0SD.coneTipZ)*1000, 'linewidth',2, 'color', cols(3,:));

h2 = plot(data_GRIN_uniform_CinC.incidenceAngle, (data_GRIN_uniform_CinC.colcHeight- data_GRIN_uniform_CinC.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_uniform_EC.incidenceAngle, (data_GRIN_uniform_EC.colcHeight- data_GRIN_uniform_EC.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_uniform_CinC_EC.incidenceAngle, (data_GRIN_uniform_CinC_EC.colcHeight- data_GRIN_uniform_CinC_EC.coneTipZ)*1000, 'linewidth',2, 'color', cols(6,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
legend([h1 h2 h3 h4], {'SD 0', 'CinC', 'EC', 'CinC & EC'},'Location','southwest')
ylim([-150 50]); ylabel('COLC Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -150:50:50,'XTick',[0 5 10 15 20]);

% Varying intercone RI
subplot(3,nPlots,9); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, (data_GRIN_radialTop.colcHeight- data_GRIN_radialTop.coneTipZ)*1000, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_radialTop_IC_142.incidenceAngle, (data_GRIN_radialTop_IC_142.colcHeight- data_GRIN_radialTop_IC_142.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_radialTop_IC_144.incidenceAngle, (data_GRIN_radialTop_IC_144.colcHeight- data_GRIN_radialTop_IC_144.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_radialTop_IC_146.incidenceAngle, (data_GRIN_radialTop_IC_146.colcHeight- data_GRIN_radialTop_IC_146.coneTipZ)*1000, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_GRIN_radialTop_IC_148.incidenceAngle, (data_GRIN_radialTop_IC_148.colcHeight- data_GRIN_radialTop_IC_148.coneTipZ)*1000, 'linewidth',2, 'color', cols(7,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
legend([h1 h2 h3 h4 h5], {'Radial-top', 'IC-1.42', 'IC-1.44', 'IC-1.46', 'IC-1.48'},'Location','southwest')
ylim([-150 50]); ylabel('COLC Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -150:50:50,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,10); hold on
h1 = plot(data_uniform_0SD.incidenceAngle, (data_uniform_0SD.colcHeight- data_uniform_0SD.coneTipZ)*1000, 'linewidth',2, 'color', cols(3,:));

h2 = plot(data_uniform_0SD_IC_142.incidenceAngle, (data_uniform_0SD_IC_142.colcHeight- data_uniform_0SD_IC_142.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_uniform_0SD_IC_144.incidenceAngle, (data_uniform_0SD_IC_144.colcHeight- data_uniform_0SD_IC_144.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_uniform_0SD_IC_146.incidenceAngle, (data_uniform_0SD_IC_146.colcHeight- data_uniform_0SD_IC_146.coneTipZ)*1000, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_uniform_0SD_IC_148.incidenceAngle, (data_uniform_0SD_IC_148.colcHeight- data_uniform_0SD_IC_148.coneTipZ)*1000, 'linewidth',2, 'color', cols(7,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
legend([h1 h2 h3 h4 h5], {'SD 0', 'IC-1.42', 'IC-1.44', 'IC-1.46', 'IC-1.48'},'Location','southwest')
ylim([-150 50]); ylabel('COLC Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -150:50:50,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,11); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, (data_GRIN_cylinder.colcHeight- data_GRIN_cylinder.coneTipZ)*1000, 'linewidth',2, 'color', cols(2,:));

h2 = plot(data_GRIN_cylinder_IC_142.incidenceAngle, (data_GRIN_cylinder_IC_142.colcHeight- data_GRIN_cylinder_IC_142.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_cylinder_IC_144.incidenceAngle, (data_GRIN_cylinder_IC_144.colcHeight- data_GRIN_cylinder_IC_144.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_cylinder_IC_146.incidenceAngle, (data_GRIN_cylinder_IC_146.colcHeight- data_GRIN_cylinder_IC_146.coneTipZ)*1000, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_GRIN_cylinder_IC_148.incidenceAngle, (data_GRIN_cylinder_IC_148.colcHeight- data_GRIN_cylinder_IC_148.coneTipZ)*1000, 'linewidth',2, 'color', cols(7,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
lh3 =legend([h1 h2 h3 h4 h5], {'Theory-equal', 'IC-1.42', 'IC-1.44', 'IC-1.46', 'IC-1.48'},'Location','southwest');
ylim([-150 50]); ylabel('COLC Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -150:50:50,'XTick',[0 5 10 15 20]);
% set(lh3,'Position',[0.6083    0.3001    0.0881    0.1333])
%% Plot COLC radius
COLCRadiusF = figure; set(gcf, 'position', [1 85 1680 870])

subplot(3,nPlots,1); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.colcRadius*1000, 'linewidth',2, 'color', cols(1,:));

h2= plot(data_GRIN_radialBase.incidenceAngle, data_GRIN_radialBase.colcRadius*1000, 'linewidth',2, 'color', cols(4,:));

h3= plot(data_GRIN_linear.incidenceAngle, data_GRIN_linear.colcRadius*1000, 'linewidth',2, 'color', cols(5,:));

h4= plot(data_GRIN_both.incidenceAngle, data_GRIN_both.colcRadius*1000, 'linewidth',2, 'color', cols(6,:));

ylim([0 50]); ylabel('COLC Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,4); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.colcRadius*1000, 'linewidth',2, 'color', cols(1,:));

h2= plot(data_GRIN_radialTop_tipCorr.incidenceAngle, data_GRIN_radialTop_tipCorr.colcRadius*1000, 'linewidth',2, 'color', cols(4,:));

ylim([0 50]); ylabel('COLC Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,3); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.colcRadius*1000, 'linewidth',2, 'color', cols(2,:));

h2 = plot(data_GRIN_cylinder_p36.incidenceAngle, data_GRIN_cylinder_p36.colcRadius*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_cylinder_m36.incidenceAngle, data_GRIN_cylinder_m36.colcRadius*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_cylinder_shape.incidenceAngle, data_GRIN_cylinder_shape.colcRadius*1000, 'linewidth',2, 'color', cols(6,:));

ylim([0 50]); ylabel('COLC Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,2); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.colcRadius*1000, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_radialTop_CinC.incidenceAngle, data_GRIN_radialTop_CinC.colcRadius*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_radialTop_EC.incidenceAngle, data_GRIN_radialTop_EC.colcRadius*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_radialTop_CinC_EC.incidenceAngle, data_GRIN_radialTop_CinC_EC.colcRadius*1000, 'linewidth',2, 'color', cols(6,:));

ylim([0 50]); ylabel('COLC Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,8); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.colcRadius*1000, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.colcRadius*1000, 'linewidth',2, 'color', cols(2,:));

h3 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.colcRadius*1000, 'linewidth',2, 'color', cols(3,:));

h4 = plot(data_GRIN_cylinder_norm_max.incidenceAngle, data_GRIN_cylinder_norm_max.colcRadius*1000, 'linewidth',2, 'color', cols(4,:));

h5 = plot(data_GRIN_cylinder_norm_80.incidenceAngle, data_GRIN_cylinder_norm_80.colcRadius*1000, 'linewidth',2, 'color', cols(5,:));

ylim([0 50]); ylabel('COLC Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

% Uniform
subplot(3,nPlots,5); hold on
h1 = plot(data_uniform_m4SD.incidenceAngle, data_uniform_m4SD.colcRadius*1000, 'linewidth',2, 'color', cols(4,:));

h2 = plot(data_uniform_m2SD.incidenceAngle, data_uniform_m2SD.colcRadius*1000, 'linewidth',2, 'color', cols(5,:));

h3 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.colcRadius*1000, 'linewidth',2, 'color', cols(3,:));

h4 = plot(data_uniform_p2SD.incidenceAngle, data_uniform_p2SD.colcRadius*1000, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_uniform_p4SD.incidenceAngle, data_uniform_p4SD.colcRadius*1000, 'linewidth',2, 'color', cols(7,:));

ylim([0 50]); ylabel('COLC Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,6); hold on
h1 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.colcRadius*1000, 'linewidth',2, 'color', cols(3,:));

h2 = plot(data_GRIN_uniform_CinC.incidenceAngle, data_GRIN_uniform_CinC.colcRadius*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_uniform_EC.incidenceAngle, data_GRIN_uniform_EC.colcRadius*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_uniform_CinC_EC.incidenceAngle, data_GRIN_uniform_CinC_EC.colcRadius*1000, 'linewidth',2, 'color', cols(6,:));

ylim([0 50]); ylabel('COLC Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

% Varying intercone RI
subplot(3,nPlots,9); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.colcRadius*1000, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_radialTop_IC_142.incidenceAngle, data_GRIN_radialTop_IC_142.colcRadius*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_radialTop_IC_144.incidenceAngle, data_GRIN_radialTop_IC_144.colcRadius*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_radialTop_IC_146.incidenceAngle, data_GRIN_radialTop_IC_146.colcRadius*1000, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_GRIN_radialTop_IC_148.incidenceAngle, data_GRIN_radialTop_IC_148.colcRadius*1000, 'linewidth',2, 'color', cols(7,:));

ylim([0 50]); ylabel('COLC Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,10); hold on
h1 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.colcRadius*1000, 'linewidth',2, 'color', cols(3,:));

h2 = plot(data_uniform_0SD_IC_142.incidenceAngle, data_uniform_0SD_IC_142.colcRadius*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_uniform_0SD_IC_144.incidenceAngle, data_uniform_0SD_IC_144.colcRadius*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_uniform_0SD_IC_146.incidenceAngle, data_uniform_0SD_IC_146.colcRadius*1000, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_uniform_0SD_IC_148.incidenceAngle, data_uniform_0SD_IC_148.colcRadius*1000, 'linewidth',2, 'color', cols(7,:));

ylim([0 50]); ylabel('COLC Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,11); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.colcRadius*1000, 'linewidth',2, 'color', cols(2,:));

h2 = plot(data_GRIN_cylinder_IC_142.incidenceAngle, data_GRIN_cylinder_IC_142.colcRadius*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_cylinder_IC_144.incidenceAngle, data_GRIN_cylinder_IC_144.colcRadius*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_cylinder_IC_146.incidenceAngle, data_GRIN_cylinder_IC_146.colcRadius*1000, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_GRIN_cylinder_IC_148.incidenceAngle, data_GRIN_cylinder_IC_148.colcRadius*1000, 'linewidth',2, 'color', cols(7,:));

ylim([0 50]); ylabel('COLC Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);
%% Plot X Focus Height
XFocusHeightF = figure; set(gcf, 'position', [1 85 1680 870])

subplot(3,nPlots,1); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, (data_GRIN_radialTop.focusXHeight- data_GRIN_radialTop.coneTipZ)*1000, 'linewidth',2, 'color', cols(1,:));

h2= plot(data_GRIN_radialBase.incidenceAngle, (data_GRIN_radialBase.focusXHeight- data_GRIN_radialBase.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h3= plot(data_GRIN_linear.incidenceAngle, (data_GRIN_linear.focusXHeight- data_GRIN_linear.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

h4= plot(data_GRIN_both.incidenceAngle, (data_GRIN_both.focusXHeight- data_GRIN_both.coneTipZ)*1000, 'linewidth',2, 'color', cols(6,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
legend([h1 h2 h3 h4], {'Radial-top','Radial-base','Linear','Radial+linear'},'Location','southwest')
ylim([-250 50]); ylabel('X-Focus Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -250:50:50,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,4); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, (data_GRIN_radialTop.focusXHeight- data_GRIN_radialTop.coneTipZ)*1000, 'linewidth',2, 'color', cols(1,:));

h2= plot(data_GRIN_radialTop_tipCorr.incidenceAngle, (data_GRIN_radialTop_tipCorr.focusXHeight- data_GRIN_radialTop_tipCorr.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
legend([h1 h2], {'Radial-top','W/ Tip Corr.'},'Location','southwest')
ylim([-250 50]); ylabel('X-Focus Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -250:50:50,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,3); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, (data_GRIN_cylinder.focusXHeight- data_GRIN_cylinder.coneTipZ)*1000, 'linewidth',2, 'color', cols(2,:));

h2 = plot(data_GRIN_cylinder_p36.incidenceAngle, (data_GRIN_cylinder_p36.focusXHeight- data_GRIN_cylinder_p36.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_cylinder_m36.incidenceAngle, (data_GRIN_cylinder_m36.focusXHeight- data_GRIN_cylinder_m36.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_cylinder_shape.incidenceAngle, (data_GRIN_cylinder_shape.focusXHeight- data_GRIN_cylinder_shape.coneTipZ)*1000, 'linewidth',2, 'color', cols(6,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
lh = legend([h1 h2 h3 h4], {'Theory-equal', 'Theory-plus', 'Theory-minus', 'Cylinder'},'Location','southwest');
ylim([-250 50]); ylabel('X-Focus Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -250:50:50,'XTick',[0 5 10 15 20]);
% set(lh,'Position',[0.5774    0.5470    0.0905    0.1075])

subplot(3,nPlots,2); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, (data_GRIN_radialTop.focusXHeight- data_GRIN_radialTop.coneTipZ)*1000, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_radialTop_CinC.incidenceAngle, (data_GRIN_radialTop_CinC.focusXHeight- data_GRIN_radialTop_CinC.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_radialTop_EC.incidenceAngle, (data_GRIN_radialTop_EC.focusXHeight- data_GRIN_radialTop_EC.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_radialTop_CinC_EC.incidenceAngle, (data_GRIN_radialTop_CinC_EC.focusXHeight- data_GRIN_radialTop_CinC_EC.coneTipZ)*1000, 'linewidth',2, 'color', cols(6,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
legend([h1 h2 h3 h4], {'Radial-top', 'CinC', 'EC', 'CinC & EC'},'Location','southwest')
ylim([-250 50]); ylabel('X-Focus Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -250:50:50,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,8); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, (data_GRIN_radialTop.focusXHeight- data_GRIN_radialTop.coneTipZ)*1000, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_cylinder.incidenceAngle, (data_GRIN_cylinder.focusXHeight- data_GRIN_cylinder.coneTipZ)*1000, 'linewidth',2, 'color', cols(2,:));

h3 = plot(data_uniform_0SD.incidenceAngle, (data_uniform_0SD.focusXHeight- data_uniform_0SD.coneTipZ)*1000, 'linewidth',2, 'color', cols(3,:));

h4 = plot(data_GRIN_cylinder_norm_max.incidenceAngle, (data_GRIN_cylinder_norm_max.focusXHeight- data_GRIN_cylinder_norm_max.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h5 = plot(data_GRIN_cylinder_norm_80.incidenceAngle, (data_GRIN_cylinder_norm_80.focusXHeight- data_GRIN_cylinder_norm_80.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
lh2 = legend([h1 h2 h3 h4 h5], {'Radial-top','Theory-equal', 'SD 0', 'Normalised-max', 'Normalised-80'},'Location','southwest');
ylim([-250 50]); ylabel('X-Focus Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -250:50:50,'XTick',[0 5 10 15 20]);
% set(lh2,'Position',[0.7756    0.2459    0.1042    0.1075])

% Uniform
subplot(3,nPlots,5); hold on
h1 = plot(data_uniform_m4SD.incidenceAngle, (data_uniform_m4SD.focusXHeight- data_uniform_m4SD.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h2 = plot(data_uniform_m2SD.incidenceAngle, (data_uniform_m2SD.focusXHeight- data_uniform_m2SD.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

h3 = plot(data_uniform_0SD.incidenceAngle, (data_uniform_0SD.focusXHeight- data_uniform_0SD.coneTipZ)*1000, 'linewidth',2, 'color', cols(3,:));

h4 = plot(data_uniform_p2SD.incidenceAngle, (data_uniform_p2SD.focusXHeight- data_uniform_p2SD.coneTipZ)*1000, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_uniform_p4SD.incidenceAngle, (data_uniform_p4SD.focusXHeight- data_uniform_p4SD.coneTipZ)*1000, 'linewidth',2, 'color', cols(7,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
legend([h1 h2 h3 h4 h5], {'SD -4', 'SD -2', 'SD 0', 'SD +2', 'SD +4'},'Location','southwest')
ylim([-250 50]); ylabel('X-Focus Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -250:50:50,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,6); hold on
h1 = plot(data_uniform_0SD.incidenceAngle, (data_uniform_0SD.focusXHeight- data_uniform_0SD.coneTipZ)*1000, 'linewidth',2, 'color', cols(3,:));

h2 = plot(data_GRIN_uniform_CinC.incidenceAngle, (data_GRIN_uniform_CinC.focusXHeight- data_GRIN_uniform_CinC.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_uniform_EC.incidenceAngle, (data_GRIN_uniform_EC.focusXHeight- data_GRIN_uniform_EC.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_uniform_CinC_EC.incidenceAngle, (data_GRIN_uniform_CinC_EC.focusXHeight- data_GRIN_uniform_CinC_EC.coneTipZ)*1000, 'linewidth',2, 'color', cols(6,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
legend([h1 h2 h3 h4], {'SD 0', 'CinC', 'EC', 'CinC & EC'},'Location','southwest')
ylim([-250 50]); ylabel('X-Focus Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -250:50:50,'XTick',[0 5 10 15 20]);

% Varying intercone RI
subplot(3,nPlots,9); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, (data_GRIN_radialTop.focusXHeight- data_GRIN_radialTop.coneTipZ)*1000, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_radialTop_IC_142.incidenceAngle, (data_GRIN_radialTop_IC_142.focusXHeight- data_GRIN_radialTop_IC_142.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_radialTop_IC_144.incidenceAngle, (data_GRIN_radialTop_IC_144.focusXHeight- data_GRIN_radialTop_IC_144.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_radialTop_IC_146.incidenceAngle, (data_GRIN_radialTop_IC_146.focusXHeight- data_GRIN_radialTop_IC_146.coneTipZ)*1000, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_GRIN_radialTop_IC_148.incidenceAngle, (data_GRIN_radialTop_IC_148.focusXHeight- data_GRIN_radialTop_IC_148.coneTipZ)*1000, 'linewidth',2, 'color', cols(7,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
legend([h1 h2 h3 h4 h5], {'Radial-top', 'IC-1.42', 'IC-1.44', 'IC-1.46', 'IC-1.48'},'Location','southwest')
ylim([-250 50]); ylabel('X-Focus Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -250:50:50,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,10); hold on
h1 = plot(data_uniform_0SD.incidenceAngle, (data_uniform_0SD.focusXHeight- data_uniform_0SD.coneTipZ)*1000, 'linewidth',2, 'color', cols(3,:));

h2 = plot(data_uniform_0SD_IC_142.incidenceAngle, (data_uniform_0SD_IC_142.focusXHeight- data_uniform_0SD_IC_142.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_uniform_0SD_IC_144.incidenceAngle, (data_uniform_0SD_IC_144.focusXHeight- data_uniform_0SD_IC_144.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_uniform_0SD_IC_146.incidenceAngle, (data_uniform_0SD_IC_146.focusXHeight- data_uniform_0SD_IC_146.coneTipZ)*1000, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_uniform_0SD_IC_148.incidenceAngle, (data_uniform_0SD_IC_148.focusXHeight- data_uniform_0SD_IC_148.coneTipZ)*1000, 'linewidth',2, 'color', cols(7,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
legend([h1 h2 h3 h4 h5], {'SD 0', 'IC-1.42', 'IC-1.44', 'IC-1.46', 'IC-1.48'},'Location','southwest')
ylim([-250 50]); ylabel('X-Focus Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -250:50:50,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,11); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, (data_GRIN_cylinder.focusXHeight- data_GRIN_cylinder.coneTipZ)*1000, 'linewidth',2, 'color', cols(2,:));

h2 = plot(data_GRIN_cylinder_IC_142.incidenceAngle, (data_GRIN_cylinder_IC_142.focusXHeight- data_GRIN_cylinder_IC_142.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_cylinder_IC_144.incidenceAngle, (data_GRIN_cylinder_IC_144.focusXHeight- data_GRIN_cylinder_IC_144.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_cylinder_IC_146.incidenceAngle, (data_GRIN_cylinder_IC_146.focusXHeight- data_GRIN_cylinder_IC_146.coneTipZ)*1000, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_GRIN_cylinder_IC_148.incidenceAngle, (data_GRIN_cylinder_IC_148.focusXHeight- data_GRIN_cylinder_IC_148.coneTipZ)*1000, 'linewidth',2, 'color', cols(7,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
lh3 =legend([h1 h2 h3 h4 h5], {'Theory-equal', 'IC-1.42', 'IC-1.44', 'IC-1.46', 'IC-1.48'},'Location','southwest');
ylim([-250 50]); ylabel('X-Focus Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -250:50:50,'XTick',[0 5 10 15 20]);
% set(lh3,'Position',[0.6083    0.3001    0.0881    0.1333])
%% Plot X Focus Radius
XFocusRadiusF = figure; set(gcf, 'position', [1 85 1680 870])

subplot(3,nPlots,1); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.focusXRadius*1000, 'linewidth',2, 'color', cols(1,:));

h2= plot(data_GRIN_radialBase.incidenceAngle, data_GRIN_radialBase.focusXRadius*1000, 'linewidth',2, 'color', cols(4,:));

h3= plot(data_GRIN_linear.incidenceAngle, data_GRIN_linear.focusXRadius*1000, 'linewidth',2, 'color', cols(5,:));

h4= plot(data_GRIN_both.incidenceAngle, data_GRIN_both.focusXRadius*1000, 'linewidth',2, 'color', cols(6,:));

ylim([0 100]); ylabel('X-Focus Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:25:100,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,4); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.focusXRadius*1000, 'linewidth',2, 'color', cols(1,:));

h2= plot(data_GRIN_radialTop_tipCorr.incidenceAngle, data_GRIN_radialTop_tipCorr.focusXRadius*1000, 'linewidth',2, 'color', cols(4,:));

ylim([0 100]); ylabel('X-Focus Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:25:100,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,3); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.focusXRadius*1000, 'linewidth',2, 'color', cols(2,:));

h2 = plot(data_GRIN_cylinder_p36.incidenceAngle, data_GRIN_cylinder_p36.focusXRadius*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_cylinder_m36.incidenceAngle, data_GRIN_cylinder_m36.focusXRadius*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_cylinder_shape.incidenceAngle, data_GRIN_cylinder_shape.focusXRadius*1000, 'linewidth',2, 'color', cols(6,:));

ylim([0 100]); ylabel('X-Focus Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:25:100,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,2); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.focusXRadius*1000, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_radialTop_CinC.incidenceAngle, data_GRIN_radialTop_CinC.focusXRadius*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_radialTop_EC.incidenceAngle, data_GRIN_radialTop_EC.focusXRadius*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_radialTop_CinC_EC.incidenceAngle, data_GRIN_radialTop_CinC_EC.focusXRadius*1000, 'linewidth',2, 'color', cols(6,:));

ylim([0 100]); ylabel('X-Focus Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:25:100,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,8); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.focusXRadius*1000, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.focusXRadius*1000, 'linewidth',2, 'color', cols(2,:));

h3 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.focusXRadius*1000, 'linewidth',2, 'color', cols(3,:));

h4 = plot(data_GRIN_cylinder_norm_max.incidenceAngle, data_GRIN_cylinder_norm_max.focusXRadius*1000, 'linewidth',2, 'color', cols(4,:));

h5 = plot(data_GRIN_cylinder_norm_80.incidenceAngle, data_GRIN_cylinder_norm_80.focusXRadius*1000, 'linewidth',2, 'color', cols(5,:));

ylim([0 100]); ylabel('X-Focus Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:25:100,'XTick',[0 5 10 15 20]);

% Uniform
subplot(3,nPlots,5); hold on
h1 = plot(data_uniform_m4SD.incidenceAngle, data_uniform_m4SD.focusXRadius*1000, 'linewidth',2, 'color', cols(4,:));

h2 = plot(data_uniform_m2SD.incidenceAngle, data_uniform_m2SD.focusXRadius*1000, 'linewidth',2, 'color', cols(5,:));

h3 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.focusXRadius*1000, 'linewidth',2, 'color', cols(3,:));

h4 = plot(data_uniform_p2SD.incidenceAngle, data_uniform_p2SD.focusXRadius*1000, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_uniform_p4SD.incidenceAngle, data_uniform_p4SD.focusXRadius*1000, 'linewidth',2, 'color', cols(7,:));

ylim([0 100]); ylabel('X-Focus Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:25:100,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,6); hold on
h1 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.focusXRadius*1000, 'linewidth',2, 'color', cols(3,:));

h2 = plot(data_GRIN_uniform_CinC.incidenceAngle, data_GRIN_uniform_CinC.focusXRadius*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_uniform_EC.incidenceAngle, data_GRIN_uniform_EC.focusXRadius*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_uniform_CinC_EC.incidenceAngle, data_GRIN_uniform_CinC_EC.focusXRadius*1000, 'linewidth',2, 'color', cols(6,:));

ylim([0 100]); ylabel('X-Focus Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:25:100,'XTick',[0 5 10 15 20]);

% Varying intercone RI
subplot(3,nPlots,9); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.focusXRadius*1000, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_radialTop_IC_142.incidenceAngle, data_GRIN_radialTop_IC_142.focusXRadius*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_radialTop_IC_144.incidenceAngle, data_GRIN_radialTop_IC_144.focusXRadius*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_radialTop_IC_146.incidenceAngle, data_GRIN_radialTop_IC_146.focusXRadius*1000, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_GRIN_radialTop_IC_148.incidenceAngle, data_GRIN_radialTop_IC_148.focusXRadius*1000, 'linewidth',2, 'color', cols(7,:));

ylim([0 100]); ylabel('X-Focus Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:25:100,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,10); hold on
h1 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.focusXRadius*1000, 'linewidth',2, 'color', cols(3,:));

h2 = plot(data_uniform_0SD_IC_142.incidenceAngle, data_uniform_0SD_IC_142.focusXRadius*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_uniform_0SD_IC_144.incidenceAngle, data_uniform_0SD_IC_144.focusXRadius*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_uniform_0SD_IC_146.incidenceAngle, data_uniform_0SD_IC_146.focusXRadius*1000, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_uniform_0SD_IC_148.incidenceAngle, data_uniform_0SD_IC_148.focusXRadius*1000, 'linewidth',2, 'color', cols(7,:));

ylim([0 100]); ylabel('X-Focus Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:25:100,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,11); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.focusXRadius*1000, 'linewidth',2, 'color', cols(2,:));

h2 = plot(data_GRIN_cylinder_IC_142.incidenceAngle, data_GRIN_cylinder_IC_142.focusXRadius*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_cylinder_IC_144.incidenceAngle, data_GRIN_cylinder_IC_144.focusXRadius*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_cylinder_IC_146.incidenceAngle, data_GRIN_cylinder_IC_146.focusXRadius*1000, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_GRIN_cylinder_IC_148.incidenceAngle, data_GRIN_cylinder_IC_148.focusXRadius*1000, 'linewidth',2, 'color', cols(7,:));

ylim([0 100]); ylabel('X-Focus Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:25:100,'XTick',[0 5 10 15 20]);
%% Plot Y Focus Height
YFocusHeightF = figure; set(gcf, 'position', [1 85 1680 870])

subplot(3,nPlots,1); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, (data_GRIN_radialTop.focusYHeight- data_GRIN_radialTop.coneTipZ)*1000, 'linewidth',2, 'color', cols(1,:));

h2= plot(data_GRIN_radialBase.incidenceAngle, (data_GRIN_radialBase.focusYHeight- data_GRIN_radialBase.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h3= plot(data_GRIN_linear.incidenceAngle, (data_GRIN_linear.focusYHeight- data_GRIN_linear.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

h4= plot(data_GRIN_both.incidenceAngle, (data_GRIN_both.focusYHeight- data_GRIN_both.coneTipZ)*1000, 'linewidth',2, 'color', cols(6,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
legend([h1 h2 h3 h4], {'Radial-top','Radial-base','Linear','Radial+linear'},'Location','southwest')
ylim([-250 50]); ylabel('Y-Focus Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -250:50:50,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,4); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, (data_GRIN_radialTop.focusYHeight- data_GRIN_radialTop.coneTipZ)*1000, 'linewidth',2, 'color', cols(1,:));

h2= plot(data_GRIN_radialTop_tipCorr.incidenceAngle, (data_GRIN_radialTop_tipCorr.focusYHeight- data_GRIN_radialTop_tipCorr.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
legend([h1 h2], {'Radial-top','W/ Tip Corr.'},'Location','southwest')
ylim([-250 50]); ylabel('Y-Focus Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -250:50:50,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,3); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, (data_GRIN_cylinder.focusYHeight- data_GRIN_cylinder.coneTipZ)*1000, 'linewidth',2, 'color', cols(2,:));

h2 = plot(data_GRIN_cylinder_p36.incidenceAngle, (data_GRIN_cylinder_p36.focusYHeight- data_GRIN_cylinder_p36.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_cylinder_m36.incidenceAngle, (data_GRIN_cylinder_m36.focusYHeight- data_GRIN_cylinder_m36.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_cylinder_shape.incidenceAngle, (data_GRIN_cylinder_shape.focusYHeight- data_GRIN_cylinder_shape.coneTipZ)*1000, 'linewidth',2, 'color', cols(6,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
lh = legend([h1 h2 h3 h4], {'Theory-equal', 'Theory-plus', 'Theory-minus', 'Cylinder'},'Location','southwest');
ylim([-250 50]); ylabel('Y-Focus Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -250:50:50,'XTick',[0 5 10 15 20]);
% set(lh,'Position',[0.5774    0.5470    0.0905    0.1075])

subplot(3,nPlots,2); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, (data_GRIN_radialTop.focusYHeight- data_GRIN_radialTop.coneTipZ)*1000, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_radialTop_CinC.incidenceAngle, (data_GRIN_radialTop_CinC.focusYHeight- data_GRIN_radialTop_CinC.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_radialTop_EC.incidenceAngle, (data_GRIN_radialTop_EC.focusYHeight- data_GRIN_radialTop_EC.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_radialTop_CinC_EC.incidenceAngle, (data_GRIN_radialTop_CinC_EC.focusYHeight- data_GRIN_radialTop_CinC_EC.coneTipZ)*1000, 'linewidth',2, 'color', cols(6,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
legend([h1 h2 h3 h4], {'Radial-top', 'CinC', 'EC', 'CinC & EC'},'Location','southwest')
ylim([-250 50]); ylabel('Y-Focus Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -250:50:50,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,8); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, (data_GRIN_radialTop.focusYHeight- data_GRIN_radialTop.coneTipZ)*1000, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_cylinder.incidenceAngle, (data_GRIN_cylinder.focusYHeight- data_GRIN_cylinder.coneTipZ)*1000, 'linewidth',2, 'color', cols(2,:));

h3 = plot(data_uniform_0SD.incidenceAngle, (data_uniform_0SD.focusYHeight- data_uniform_0SD.coneTipZ)*1000, 'linewidth',2, 'color', cols(3,:));

h4 = plot(data_GRIN_cylinder_norm_max.incidenceAngle, (data_GRIN_cylinder_norm_max.focusYHeight- data_GRIN_cylinder_norm_max.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h5 = plot(data_GRIN_cylinder_norm_80.incidenceAngle, (data_GRIN_cylinder_norm_80.focusYHeight- data_GRIN_cylinder_norm_80.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
lh2 = legend([h1 h2 h3 h4 h5], {'Radial-top','Theory-equal', 'SD 0', 'Normalised-max', 'Normalised-80'},'Location','southwest');
ylim([-250 50]); ylabel('Y-Focus Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -250:50:50,'XTick',[0 5 10 15 20]);
% set(lh2,'Position',[0.7756    0.2459    0.1042    0.1075])

% uniform
subplot(3,nPlots,5); hold on
h1 = plot(data_uniform_m4SD.incidenceAngle, (data_uniform_m4SD.focusYHeight- data_uniform_m4SD.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h2 = plot(data_uniform_m2SD.incidenceAngle, (data_uniform_m2SD.focusYHeight- data_uniform_m2SD.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

h3 = plot(data_uniform_0SD.incidenceAngle, (data_uniform_0SD.focusYHeight- data_uniform_0SD.coneTipZ)*1000, 'linewidth',2, 'color', cols(3,:));

h4 = plot(data_uniform_p2SD.incidenceAngle, (data_uniform_p2SD.focusYHeight- data_uniform_p2SD.coneTipZ)*1000, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_uniform_p4SD.incidenceAngle, (data_uniform_p4SD.focusYHeight- data_uniform_p4SD.coneTipZ)*1000, 'linewidth',2, 'color', cols(7,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
legend([h1 h2 h3 h4 h5], {'SD -4', 'SD -2', 'SD 0', 'SD +2', 'SD +4'},'Location','southwest')
ylim([-250 50]); ylabel('Y-Focus Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -250:50:50,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,6); hold on
h1 = plot(data_uniform_0SD.incidenceAngle, (data_uniform_0SD.focusYHeight- data_uniform_0SD.coneTipZ)*1000, 'linewidth',2, 'color', cols(3,:));

h2 = plot(data_GRIN_uniform_CinC.incidenceAngle, (data_GRIN_uniform_CinC.focusYHeight- data_GRIN_uniform_CinC.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_uniform_EC.incidenceAngle, (data_GRIN_uniform_EC.focusYHeight- data_GRIN_uniform_EC.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_uniform_CinC_EC.incidenceAngle, (data_GRIN_uniform_CinC_EC.focusYHeight- data_GRIN_uniform_CinC_EC.coneTipZ)*1000, 'linewidth',2, 'color', cols(6,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
legend([h1 h2 h3 h4], {'SD 0', 'CinC', 'EC', 'CinC & EC'},'Location','southwest')
ylim([-250 50]); ylabel('Y-Focus Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -250:50:50,'XTick',[0 5 10 15 20]);

% Varying intercone RI
subplot(3,nPlots,9); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, (data_GRIN_radialTop.focusYHeight- data_GRIN_radialTop.coneTipZ)*1000, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_radialTop_IC_142.incidenceAngle, (data_GRIN_radialTop_IC_142.focusYHeight- data_GRIN_radialTop_IC_142.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_radialTop_IC_144.incidenceAngle, (data_GRIN_radialTop_IC_144.focusYHeight- data_GRIN_radialTop_IC_144.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_radialTop_IC_146.incidenceAngle, (data_GRIN_radialTop_IC_146.focusYHeight- data_GRIN_radialTop_IC_146.coneTipZ)*1000, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_GRIN_radialTop_IC_148.incidenceAngle, (data_GRIN_radialTop_IC_148.focusYHeight- data_GRIN_radialTop_IC_148.coneTipZ)*1000, 'linewidth',2, 'color', cols(7,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
legend([h1 h2 h3 h4 h5], {'Radial-top', 'IC-1.42', 'IC-1.44', 'IC-1.46', 'IC-1.48'},'Location','southwest')
ylim([-250 50]); ylabel('Y-Focus Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -250:50:50,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,10); hold on
h1 = plot(data_uniform_0SD.incidenceAngle, (data_uniform_0SD.focusYHeight- data_uniform_0SD.coneTipZ)*1000, 'linewidth',2, 'color', cols(3,:));

h2 = plot(data_uniform_0SD_IC_142.incidenceAngle, (data_uniform_0SD_IC_142.focusYHeight- data_uniform_0SD_IC_142.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_uniform_0SD_IC_144.incidenceAngle, (data_uniform_0SD_IC_144.focusYHeight- data_uniform_0SD_IC_144.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_uniform_0SD_IC_146.incidenceAngle, (data_uniform_0SD_IC_146.focusYHeight- data_uniform_0SD_IC_146.coneTipZ)*1000, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_uniform_0SD_IC_148.incidenceAngle, (data_uniform_0SD_IC_148.focusYHeight- data_uniform_0SD_IC_148.coneTipZ)*1000, 'linewidth',2, 'color', cols(7,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
legend([h1 h2 h3 h4 h5], {'SD 0', 'IC-1.42', 'IC-1.44', 'IC-1.46', 'IC-1.48'},'Location','southwest')
ylim([-250 50]); ylabel('Y-Focus Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -250:50:50,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,11); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, (data_GRIN_cylinder.focusYHeight- data_GRIN_cylinder.coneTipZ)*1000, 'linewidth',2, 'color', cols(2,:));

h2 = plot(data_GRIN_cylinder_IC_142.incidenceAngle, (data_GRIN_cylinder_IC_142.focusYHeight- data_GRIN_cylinder_IC_142.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_cylinder_IC_144.incidenceAngle, (data_GRIN_cylinder_IC_144.focusYHeight- data_GRIN_cylinder_IC_144.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_cylinder_IC_146.incidenceAngle, (data_GRIN_cylinder_IC_146.focusYHeight- data_GRIN_cylinder_IC_146.coneTipZ)*1000, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_GRIN_cylinder_IC_148.incidenceAngle, (data_GRIN_cylinder_IC_148.focusYHeight- data_GRIN_cylinder_IC_148.coneTipZ)*1000, 'linewidth',2, 'color', cols(7,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
lh3 =legend([h1 h2 h3 h4 h5], {'Theory-equal', 'IC-1.42', 'IC-1.44', 'IC-1.46', 'IC-1.48'},'Location','southwest');
ylim([-250 50]); ylabel('Y-Focus Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -250:50:50,'XTick',[0 5 10 15 20]);
% set(lh3,'Position',[0.6083    0.3001    0.0881    0.1333])
%% Plot Y Focus Radius
YFocusRadiusF = figure; set(gcf, 'position', [1 85 1680 870])

subplot(3,nPlots,1); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.focusYRadius*1000, 'linewidth',2, 'color', cols(1,:));

h2= plot(data_GRIN_radialBase.incidenceAngle, data_GRIN_radialBase.focusYRadius*1000, 'linewidth',2, 'color', cols(4,:));

h3= plot(data_GRIN_linear.incidenceAngle, data_GRIN_linear.focusYRadius*1000, 'linewidth',2, 'color', cols(5,:));

h4= plot(data_GRIN_both.incidenceAngle, data_GRIN_both.focusYRadius*1000, 'linewidth',2, 'color', cols(6,:));

ylim([0 100]); ylabel('Y-Focus Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:25:100,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,4); hold on

h1= plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.focusYRadius*1000, 'linewidth',2, 'color', cols(1,:));

h2= plot(data_GRIN_radialTop_tipCorr.incidenceAngle, data_GRIN_radialTop_tipCorr.focusYRadius*1000, 'linewidth',2, 'color', cols(4,:));

ylim([0 100]); ylabel('Y-Focus Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:25:100,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,3); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.focusYRadius*1000, 'linewidth',2, 'color', cols(2,:));

h2 = plot(data_GRIN_cylinder_p36.incidenceAngle, data_GRIN_cylinder_p36.focusYRadius*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_cylinder_m36.incidenceAngle, data_GRIN_cylinder_m36.focusYRadius*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_cylinder_shape.incidenceAngle, data_GRIN_cylinder_shape.focusYRadius*1000, 'linewidth',2, 'color', cols(6,:));

ylim([0 100]); ylabel('Y-Focus Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:25:100,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,2); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.focusYRadius*1000, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_radialTop_CinC.incidenceAngle, data_GRIN_radialTop_CinC.focusYRadius*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_radialTop_EC.incidenceAngle, data_GRIN_radialTop_EC.focusYRadius*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_radialTop_CinC_EC.incidenceAngle, data_GRIN_radialTop_CinC_EC.focusYRadius*1000, 'linewidth',2, 'color', cols(6,:));

ylim([0 100]); ylabel('Y-Focus Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:25:100,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,8); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.focusYRadius*1000, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.focusYRadius*1000, 'linewidth',2, 'color', cols(2,:));

h3 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.focusYRadius*1000, 'linewidth',2, 'color', cols(3,:));

h4 = plot(data_GRIN_cylinder_norm_max.incidenceAngle, data_GRIN_cylinder_norm_max.focusYRadius*1000, 'linewidth',2, 'color', cols(4,:));

h5 = plot(data_GRIN_cylinder_norm_80.incidenceAngle, data_GRIN_cylinder_norm_80.focusYRadius*1000, 'linewidth',2, 'color', cols(5,:));

ylim([0 100]); ylabel('Y-Focus Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:25:100,'XTick',[0 5 10 15 20]);

% uniform
subplot(3,nPlots,5); hold on
h1 = plot(data_uniform_m4SD.incidenceAngle, data_uniform_m4SD.focusYRadius*1000, 'linewidth',2, 'color', cols(4,:));

h2 = plot(data_uniform_m2SD.incidenceAngle, data_uniform_m2SD.focusYRadius*1000, 'linewidth',2, 'color', cols(5,:));

h3 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.focusYRadius*1000, 'linewidth',2, 'color', cols(3,:));

h4 = plot(data_uniform_p2SD.incidenceAngle, data_uniform_p2SD.focusYRadius*1000, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_uniform_p4SD.incidenceAngle, data_uniform_p4SD.focusYRadius*1000, 'linewidth',2, 'color', cols(7,:));

ylim([0 100]); ylabel('Y-Focus Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:25:100,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,6); hold on
h1 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.focusYRadius*1000, 'linewidth',2, 'color', cols(3,:));

h2 = plot(data_GRIN_uniform_CinC.incidenceAngle, data_GRIN_uniform_CinC.focusYRadius*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_uniform_EC.incidenceAngle, data_GRIN_uniform_EC.focusYRadius*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_uniform_CinC_EC.incidenceAngle, data_GRIN_uniform_CinC_EC.focusYRadius*1000, 'linewidth',2, 'color', cols(6,:));

ylim([0 100]); ylabel('Y-Focus Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:25:100,'XTick',[0 5 10 15 20]);

% Varying intercone RI
subplot(3,nPlots,9); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.focusYRadius*1000, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_radialTop_IC_142.incidenceAngle, data_GRIN_radialTop_IC_142.focusYRadius*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_radialTop_IC_144.incidenceAngle, data_GRIN_radialTop_IC_144.focusYRadius*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_radialTop_IC_146.incidenceAngle, data_GRIN_radialTop_IC_146.focusYRadius*1000, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_GRIN_radialTop_IC_148.incidenceAngle, data_GRIN_radialTop_IC_148.focusYRadius*1000, 'linewidth',2, 'color', cols(7,:));

ylim([0 100]); ylabel('Y-Focus Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:25:100,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,10); hold on
h1 = plot(data_uniform_0SD.incidenceAngle, data_uniform_0SD.focusYRadius*1000, 'linewidth',2, 'color', cols(3,:));

h2 = plot(data_uniform_0SD_IC_142.incidenceAngle, data_uniform_0SD_IC_142.focusYRadius*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_uniform_0SD_IC_144.incidenceAngle, data_uniform_0SD_IC_144.focusYRadius*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_uniform_0SD_IC_146.incidenceAngle, data_uniform_0SD_IC_146.focusYRadius*1000, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_uniform_0SD_IC_148.incidenceAngle, data_uniform_0SD_IC_148.focusYRadius*1000, 'linewidth',2, 'color', cols(7,:));

ylim([0 100]); ylabel('Y-Focus Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:25:100,'XTick',[0 5 10 15 20]);

subplot(3,nPlots,11); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.focusYRadius*1000, 'linewidth',2, 'color', cols(2,:));

h2 = plot(data_GRIN_cylinder_IC_142.incidenceAngle, data_GRIN_cylinder_IC_142.focusYRadius*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_cylinder_IC_144.incidenceAngle, data_GRIN_cylinder_IC_144.focusYRadius*1000, 'linewidth',2, 'color', cols(5,:));

h4 = plot(data_GRIN_cylinder_IC_146.incidenceAngle, data_GRIN_cylinder_IC_146.focusYRadius*1000, 'linewidth',2, 'color', cols(6,:));

h5 = plot(data_GRIN_cylinder_IC_148.incidenceAngle, data_GRIN_cylinder_IC_148.focusYRadius*1000, 'linewidth',2, 'color', cols(7,:));

ylim([0 100]); ylabel('Y-Focus Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:25:100,'XTick',[0 5 10 15 20]);
%% 
if saveFigures
    currentDirectory = pwd; 
    cd(saveFolder)
    
    exportgraphics(acceptanceF,'AcceptanceAngleSummary.pdf','ContentType','vector') 
    exportgraphics(TIRF,'TIRRatioSummary.pdf','ContentType','vector') 
    
    exportgraphics(COLCHeightF,'COLCHeightSummary.pdf','ContentType','vector') 
    exportgraphics(COLCRadiusF,'COLCRadiusSummary.pdf','ContentType','vector') 
    
    exportgraphics(XFocusHeightF,'XFocusHeightSummary.pdf','ContentType','vector') 
    exportgraphics(XFocusRadiusF,'XFocusRadiusSummary.pdf','ContentType','vector') 
    
    exportgraphics(YFocusHeightF,'YFocusHeightSummary.pdf','ContentType','vector') 
    exportgraphics(YFocusRadiusF,'YFocusRadiusSummary.pdf','ContentType','vector') 
    
    cd(currentDirectory)
end