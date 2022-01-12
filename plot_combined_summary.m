%%% Need to set figure size to get nicely shaped axes...

%%% Check cylinder cylinder for COLC, X Focus and Y Focus height (and radius?)

clear; close all; clc

% Load reference acceptance function
acceptanceAngle = readmatrix('/Users/gavintaylor/Documents/Matlab/Git_versioned_April_27/Crab_cone_git/Data/Acceptacne Function.csv');

nPlots = 4;
cols = lines(5);
outlineCol = [0.5 0.5 0.5];

% Load data from varying RI series
cd('/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Data/Varying RI Profile/')
data_GRIN_both = load('Cone_1000_nm_Cone_0_SD_GRIN_both_SIMDATA.mat');
data_GRIN_radialBase = load('Cone_1000_nm_Cone_0_SD_GRIN_radialBase_SIMDATA.mat');
data_GRIN_radialTop = load('Cone_1000_nm_Cone_0_SD_GRIN_radialTop_SIMDATA.mat');
data_GRIN_linear = load('Cone_1000_nm_Cone_0_SD_GRIN_linear_SIMDATA.mat');

data_GRIN_radialTop_tipCorr = load('Cone_1000_nm_Cone_0_SD_GRIN_radialTop_TipCorrection_SIMDATA.mat');

data_GRIN_cylinder = load('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_SIMDATA.mat');
data_GRIN_cylinder_m36 = load('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_m36_SIMDATA.mat');
data_GRIN_cylinder_p36 = load('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_p36_SIMDATA.mat');
data_GRIN_cylinder_shape = load('Cylinder_1000_nm_Cone_0_SD_GRIN_cylinder_SIMDATA.mat');

% Load data from internal structures (w/ varying RI)
cd('/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Data/Extra internal structures/Radial RI')
data_GRIN_radialTop_CinC = load('Cone_CinC_1000_nm_Cone_0_SD_GRIN_radialTop_SIMDATA.mat');
data_GRIN_radialTop_EC = load('Cone_EC_1000_nm_Cone_0_SD_GRIN_radialTop_SIMDATA.mat');
data_GRIN_radialTop_CinC_EC = load('Cone_CinC_EC_1000_nm_Cone_0_SD_GRIN_radialTop_SIMDATA.mat');

%%% Need to add varying shape and constant RI w/ other structs



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
acceptanceF = figure;

subplot(1,nPlots,1); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.acceptancePercentage, 'linewidth',2, 'color', cols(1,:));
plot(data_GRIN_radialTop.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(1,:))

h2= plot(data_GRIN_radialBase.incidenceAngle, data_GRIN_radialBase.acceptancePercentage, 'linewidth',2, 'color', cols(2,:));
plot(data_GRIN_radialBase.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(2,:))

h3= plot(data_GRIN_linear.incidenceAngle, data_GRIN_linear.acceptancePercentage, 'linewidth',2, 'color', cols(3,:));
plot(data_GRIN_linear.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(3,:))

h4= plot(data_GRIN_both.incidenceAngle, data_GRIN_both.acceptancePercentage, 'linewidth',2, 'color', cols(4,:));
plot(data_GRIN_both.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(4,:))

legend([h1 h2 h3 h4], {'Radial-top','Radial-base','Linear','Radial+linear'},'Location','southwest')
ylabel('% entering receptor'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);

subplot(1,nPlots,2); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.acceptancePercentage, 'linewidth',2, 'color', cols(1,:));
plot(data_GRIN_radialTop.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(1,:))

h2= plot(data_GRIN_radialTop_tipCorr.incidenceAngle, data_GRIN_radialTop_tipCorr.acceptancePercentage, 'linewidth',2, 'color', cols(2,:));
plot(data_GRIN_radialTop_tipCorr.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(2,:))

legend([h1 h2], {'Radial-top','W/ Tip Corr.'},'Location','southwest')
ylabel('% entering receptor'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);

subplot(1,nPlots,3); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.acceptancePercentage, 'linewidth',2, 'color', cols(2,:));
plot(data_GRIN_cylinder.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(2,:))

h2 = plot(data_GRIN_cylinder_p36.incidenceAngle, data_GRIN_cylinder_p36.acceptancePercentage, 'linewidth',2, 'color', cols(4,:));
plot(data_GRIN_cylinder_p36.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(4,:))

h3 = plot(data_GRIN_cylinder_m36.incidenceAngle, data_GRIN_cylinder_m36.acceptancePercentage, 'linewidth',2, 'color', cols(3,:));
plot(data_GRIN_cylinder_m36.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(3,:))

h4 = plot(data_GRIN_cylinder_shape.incidenceAngle, data_GRIN_cylinder_shape.acceptancePercentage, 'linewidth',2, 'color', cols(5,:));
plot(data_GRIN_cylinder_shape.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(5,:))

legend([h1 h2 h3 h4], {'Theory-equal', 'Theory-plus', 'Theory-minus', 'Cylinder'},'Location','southwest')
ylabel('% entering receptor'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);

subplot(1,nPlots,4); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.acceptancePercentage, 'linewidth',2, 'color', cols(1,:));
plot(data_GRIN_radialTop.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(1,:))

h2 = plot(data_GRIN_radialTop_CinC.incidenceAngle, data_GRIN_radialTop_CinC.acceptancePercentage, 'linewidth',2, 'color', cols(2,:));
plot(data_GRIN_radialTop_CinC.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(2,:))

h3 = plot(data_GRIN_radialTop_EC.incidenceAngle, data_GRIN_radialTop_EC.acceptancePercentage, 'linewidth',2, 'color', cols(3,:));
plot(data_GRIN_radialTop_EC.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(3,:))

h4 = plot(data_GRIN_radialTop_CinC_EC.incidenceAngle, data_GRIN_radialTop_CinC_EC.acceptancePercentage, 'linewidth',2, 'color', cols(4,:));
plot(data_GRIN_radialTop_CinC_EC.acceptanceAngle, 0.5, 'd','markersize',10, 'linewidth',2, 'color', cols(4,:))

legend([h1 h2 h3 h4], {'Radial-top', 'CinC', 'EC', 'CinC & EC'},'Location','southwest')
ylabel('% entering receptor'); 
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);

%% Plot TIR
TIRF = figure;

subplot(1,nPlots,1); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.TIRPercengtage, 'linewidth',2, 'color', cols(1,:));

h2= plot(data_GRIN_radialBase.incidenceAngle, data_GRIN_radialBase.TIRPercengtage, 'linewidth',2, 'color', cols(2,:));

h3= plot(data_GRIN_linear.incidenceAngle, data_GRIN_linear.TIRPercengtage, 'linewidth',2, 'color', cols(3,:));

h4= plot(data_GRIN_both.incidenceAngle, data_GRIN_both.TIRPercengtage, 'linewidth',2, 'color', cols(4,:));

legend([h1 h2 h3 h4], {'Radial-top','Radial-base','Linear','Radial+linear'},'Location','southwest')
ylim([0 1]); xlim([0 20]); ylabel('% TIR in COLC'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);

subplot(1,nPlots,2); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.TIRPercengtage, 'linewidth',2, 'color', cols(1,:));

h2= plot(data_GRIN_radialTop_tipCorr.incidenceAngle, data_GRIN_radialTop_tipCorr.TIRPercengtage, 'linewidth',2, 'color', cols(2,:));

legend([h1 h2], {'Radial-top','W/ Tip Corr.'},'Location','southwest')
ylim([0 1]); xlim([0 20]); ylabel('% TIR in COLC'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);

subplot(1,nPlots,3); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.TIRPercengtage, 'linewidth',2, 'color', cols(2,:));

h2 = plot(data_GRIN_cylinder_p36.incidenceAngle, data_GRIN_cylinder_p36.TIRPercengtage, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_cylinder_m36.incidenceAngle, data_GRIN_cylinder_m36.TIRPercengtage, 'linewidth',2, 'color', cols(3,:));

h4 = plot(data_GRIN_cylinder_shape.incidenceAngle, data_GRIN_cylinder_shape.TIRPercengtage, 'linewidth',2, 'color', cols(5,:));

legend([h1 h2 h3 h4], {'Theory-equal', 'Theory-plus', 'Theory-minus', 'Cylinder'},'Location','southwest')
ylim([0 1]); xlim([0 20]); ylabel('% TIR in COLC'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);


subplot(1,nPlots,4); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.TIRPercengtage, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_radialTop_CinC.incidenceAngle, data_GRIN_radialTop_CinC.TIRPercengtage, 'linewidth',2, 'color', cols(2,:));

h3 = plot(data_GRIN_radialTop_EC.incidenceAngle, data_GRIN_radialTop_EC.TIRPercengtage, 'linewidth',2, 'color', cols(3,:));

h4 = plot(data_GRIN_radialTop_CinC_EC.incidenceAngle, data_GRIN_radialTop_CinC_EC.TIRPercengtage, 'linewidth',2, 'color', cols(4,:));

legend([h1 h2 h3 h4], {'Radial-top', 'CinC', 'EC', 'CinC & EC'},'Location','southwest')
ylim([0 1]); xlim([0 20]); ylabel('% TIR in COLC'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);

%% Plot COLC height
COLCHeightF = figure;

subplot(1,nPlots,1); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, (data_GRIN_radialTop.colcHeight - data_GRIN_radialTop.coneTipZ)*1000, 'linewidth',2, 'color', cols(1,:));

h2= plot(data_GRIN_radialBase.incidenceAngle, (data_GRIN_radialBase.colcHeight- data_GRIN_radialBase.coneTipZ)*1000, 'linewidth',2, 'color', cols(2,:));

h3= plot(data_GRIN_linear.incidenceAngle, (data_GRIN_linear.colcHeight- data_GRIN_linear.coneTipZ)*1000, 'linewidth',2, 'color', cols(3,:));

h4= plot(data_GRIN_both.incidenceAngle, (data_GRIN_both.colcHeight- data_GRIN_both.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
legend([h1 h2 h3 h4], {'Radial-top','Radial-base','Linear','Radial+linear'},'Location','southwest')
ylim([-150 50]); ylabel('COLC Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -150:50:50,'XTick',[0 5 10 15 20]);

subplot(1,nPlots,2); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, (data_GRIN_radialTop.colcHeight- data_GRIN_radialTop.coneTipZ)*1000, 'linewidth',2, 'color', cols(1,:));

h2= plot(data_GRIN_radialTop_tipCorr.incidenceAngle, (data_GRIN_radialTop_tipCorr.colcHeight- data_GRIN_radialTop_tipCorr.coneTipZ)*1000, 'linewidth',2, 'color', cols(2,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
legend([h1 h2], {'Radial-top','W/ Tip Corr.'},'Location','southwest')
ylim([-150 50]); ylabel('COLC Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -150:50:50,'XTick',[0 5 10 15 20]);

subplot(1,nPlots,3); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, (data_GRIN_cylinder.colcHeight- data_GRIN_cylinder.coneTipZ)*1000, 'linewidth',2, 'color', cols(2,:));

h2 = plot(data_GRIN_cylinder_p36.incidenceAngle, (data_GRIN_cylinder_p36.colcHeight- data_GRIN_cylinder_p36.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_cylinder_m36.incidenceAngle, (data_GRIN_cylinder_m36.colcHeight- data_GRIN_cylinder_m36.coneTipZ)*1000, 'linewidth',2, 'color', cols(3,:));

h4 = plot(data_GRIN_cylinder_shape.incidenceAngle, (data_GRIN_cylinder_shape.colcHeight- data_GRIN_cylinder_shape.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
legend([h1 h2 h3 h4], {'Theory-equal', 'Theory-plus', 'Theory-minus', 'Cylinder'},'Location','southwest')
ylim([-150 50]); ylabel('COLC Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -150:50:50,'XTick',[0 5 10 15 20]);

subplot(1,nPlots,4); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, (data_GRIN_radialTop.colcHeight- data_GRIN_radialTop.coneTipZ)*1000, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_radialTop_CinC.incidenceAngle, (data_GRIN_radialTop_CinC.colcHeight- data_GRIN_radialTop_CinC.coneTipZ)*1000, 'linewidth',2, 'color', cols(2,:));

h3 = plot(data_GRIN_radialTop_EC.incidenceAngle, (data_GRIN_radialTop_EC.colcHeight- data_GRIN_radialTop_EC.coneTipZ)*1000, 'linewidth',2, 'color', cols(3,:));

h4 = plot(data_GRIN_radialTop_CinC_EC.incidenceAngle, (data_GRIN_radialTop_CinC_EC.colcHeight- data_GRIN_radialTop_CinC_EC.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
legend([h1 h2 h3 h4], {'Radial-top', 'CinC', 'EC', 'CinC & EC'},'Location','southwest')
ylim([-150 50]); ylabel('COLC Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -150:50:50,'XTick',[0 5 10 15 20]);

%% Plot COLC radius
COLCRadiusF = figure;

subplot(1,nPlots,1); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.colcRadius*1000, 'linewidth',2, 'color', cols(1,:));

h2= plot(data_GRIN_radialBase.incidenceAngle, data_GRIN_radialBase.colcRadius*1000, 'linewidth',2, 'color', cols(2,:));

h3= plot(data_GRIN_linear.incidenceAngle, data_GRIN_linear.colcRadius*1000, 'linewidth',2, 'color', cols(3,:));

h4= plot(data_GRIN_both.incidenceAngle, data_GRIN_both.colcRadius*1000, 'linewidth',2, 'color', cols(4,:));

legend([h1 h2 h3 h4], {'Radial-top','Radial-base','Linear','Radial+linear'},'Location','southwest')
ylim([0 50]); ylabel('COLC Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

subplot(1,nPlots,2); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.colcRadius*1000, 'linewidth',2, 'color', cols(1,:));

h2= plot(data_GRIN_radialTop_tipCorr.incidenceAngle, data_GRIN_radialTop_tipCorr.colcRadius*1000, 'linewidth',2, 'color', cols(2,:));

legend([h1 h2], {'Radial-top','W/ Tip Corr.'},'Location','southwest')
ylim([0 50]); ylabel('COLC Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

subplot(1,nPlots,3); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.colcRadius*1000, 'linewidth',2, 'color', cols(2,:));

h2 = plot(data_GRIN_cylinder_p36.incidenceAngle, data_GRIN_cylinder_p36.colcRadius*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_cylinder_m36.incidenceAngle, data_GRIN_cylinder_m36.colcRadius*1000, 'linewidth',2, 'color', cols(3,:));

h4 = plot(data_GRIN_cylinder_shape.incidenceAngle, data_GRIN_cylinder_shape.colcRadius*1000, 'linewidth',2, 'color', cols(5,:));

legend([h1 h2 h3 h4], {'Theory-equal', 'Theory-plus', 'Theory-minus', 'Cylinder'},'Location','southwest')
ylim([0 50]); ylabel('COLC Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

subplot(1,nPlots,4); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.colcRadius*1000, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_radialTop_CinC.incidenceAngle, data_GRIN_radialTop_CinC.colcRadius*1000, 'linewidth',2, 'color', cols(2,:));

h3 = plot(data_GRIN_radialTop_EC.incidenceAngle, data_GRIN_radialTop_EC.colcRadius*1000, 'linewidth',2, 'color', cols(3,:));

h4 = plot(data_GRIN_radialTop_CinC_EC.incidenceAngle, data_GRIN_radialTop_CinC_EC.colcRadius*1000, 'linewidth',2, 'color', cols(4,:));

legend([h1 h2 h3 h4], {'Radial-top', 'CinC', 'EC', 'CinC & EC'},'Location','southwest')
ylim([0 50]); ylabel('COLC Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:10:50,'XTick',[0 5 10 15 20]);

%% Plot X Focus Height
XFocusHeightF = figure;

subplot(1,nPlots,1); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, (data_GRIN_radialTop.focusXHeight- data_GRIN_radialTop.coneTipZ)*1000, 'linewidth',2, 'color', cols(1,:));

h2= plot(data_GRIN_radialBase.incidenceAngle, (data_GRIN_radialBase.focusXHeight- data_GRIN_radialBase.coneTipZ)*1000, 'linewidth',2, 'color', cols(2,:));

h3= plot(data_GRIN_linear.incidenceAngle, (data_GRIN_linear.focusXHeight- data_GRIN_linear.coneTipZ)*1000, 'linewidth',2, 'color', cols(3,:));

h4= plot(data_GRIN_both.incidenceAngle, (data_GRIN_both.focusXHeight- data_GRIN_both.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
legend([h1 h2 h3 h4], {'Radial-top','Radial-base','Linear','Radial+linear'},'Location','southwest')
ylim([-250 50]); ylabel('X-Focus Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -250:50:50,'XTick',[0 5 10 15 20]);

subplot(1,nPlots,2); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, (data_GRIN_radialTop.focusXHeight- data_GRIN_radialTop.coneTipZ)*1000, 'linewidth',2, 'color', cols(1,:));

h2= plot(data_GRIN_radialTop_tipCorr.incidenceAngle, (data_GRIN_radialTop_tipCorr.focusXHeight- data_GRIN_radialTop_tipCorr.coneTipZ)*1000, 'linewidth',2, 'color', cols(2,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
legend([h1 h2], {'Radial-top','W/ Tip Corr.'},'Location','southwest')
ylim([-250 50]); ylabel('X-Focus Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -250:50:50,'XTick',[0 5 10 15 20]);

subplot(1,nPlots,3); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, (data_GRIN_cylinder.focusXHeight- data_GRIN_cylinder.coneTipZ)*1000, 'linewidth',2, 'color', cols(2,:));

h2 = plot(data_GRIN_cylinder_p36.incidenceAngle, (data_GRIN_cylinder_p36.focusXHeight- data_GRIN_cylinder_p36.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_cylinder_m36.incidenceAngle, (data_GRIN_cylinder_m36.focusXHeight- data_GRIN_cylinder_m36.coneTipZ)*1000, 'linewidth',2, 'color', cols(3,:));

h4 = plot(data_GRIN_cylinder_shape.incidenceAngle, (data_GRIN_cylinder_shape.focusXHeight- data_GRIN_cylinder_shape.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
legend([h1 h2 h3 h4], {'Theory-equal', 'Theory-plus', 'Theory-minus', 'Cylinder'},'Location','southwest')
ylim([-250 50]); ylabel('X-Focus Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -250:50:50,'XTick',[0 5 10 15 20]);

subplot(1,nPlots,4); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, (data_GRIN_radialTop.focusXHeight- data_GRIN_radialTop.coneTipZ)*1000, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_radialTop_CinC.incidenceAngle, (data_GRIN_radialTop_CinC.focusXHeight- data_GRIN_radialTop_CinC.coneTipZ)*1000, 'linewidth',2, 'color', cols(2,:));

h3 = plot(data_GRIN_radialTop_EC.incidenceAngle, (data_GRIN_radialTop_EC.focusXHeight- data_GRIN_radialTop_EC.coneTipZ)*1000, 'linewidth',2, 'color', cols(3,:));

h4 = plot(data_GRIN_radialTop_CinC_EC.incidenceAngle, (data_GRIN_radialTop_CinC_EC.focusXHeight- data_GRIN_radialTop_CinC_EC.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
legend([h1 h2 h3 h4], {'Radial-top', 'CinC', 'EC', 'CinC & EC'},'Location','southwest')
ylim([-250 50]); ylabel('X-Focus Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -250:50:50,'XTick',[0 5 10 15 20]);

%% Plot X Focus Radius
XFocusRadiusF = figure;

subplot(1,nPlots,1); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.focusXRadius*1000, 'linewidth',2, 'color', cols(1,:));

h2= plot(data_GRIN_radialBase.incidenceAngle, data_GRIN_radialBase.focusXRadius*1000, 'linewidth',2, 'color', cols(2,:));

h3= plot(data_GRIN_linear.incidenceAngle, data_GRIN_linear.focusXRadius*1000, 'linewidth',2, 'color', cols(3,:));

h4= plot(data_GRIN_both.incidenceAngle, data_GRIN_both.focusXRadius*1000, 'linewidth',2, 'color', cols(4,:));

legend([h1 h2 h3 h4], {'Radial-top','Radial-base','Linear','Radial+linear'},'Location','southwest')
ylim([0 75]); ylabel('X-Focus Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:25:75,'XTick',[0 5 10 15 20]);

subplot(1,nPlots,2); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.focusXRadius*1000, 'linewidth',2, 'color', cols(1,:));

h2= plot(data_GRIN_radialTop_tipCorr.incidenceAngle, data_GRIN_radialTop_tipCorr.focusXRadius*1000, 'linewidth',2, 'color', cols(2,:));

legend([h1 h2], {'Radial-top','W/ Tip Corr.'},'Location','southwest')
ylim([0 75]); ylabel('X-Focus Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:25:75,'XTick',[0 5 10 15 20]);

subplot(1,nPlots,3); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.focusXRadius*1000, 'linewidth',2, 'color', cols(2,:));

h2 = plot(data_GRIN_cylinder_p36.incidenceAngle, data_GRIN_cylinder_p36.focusXRadius*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_cylinder_m36.incidenceAngle, data_GRIN_cylinder_m36.focusXRadius*1000, 'linewidth',2, 'color', cols(3,:));

h4 = plot(data_GRIN_cylinder_shape.incidenceAngle, data_GRIN_cylinder_shape.focusXRadius*1000, 'linewidth',2, 'color', cols(5,:));

legend([h1 h2 h3 h4], {'Theory-equal', 'Theory-plus', 'Theory-minus', 'Cylinder'},'Location','southwest')
ylim([0 75]); ylabel('X-Focus Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:25:75,'XTick',[0 5 10 15 20]);

subplot(1,nPlots,4); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.focusXRadius*1000, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_radialTop_CinC.incidenceAngle, data_GRIN_radialTop_CinC.focusXRadius*1000, 'linewidth',2, 'color', cols(2,:));

h3 = plot(data_GRIN_radialTop_EC.incidenceAngle, data_GRIN_radialTop_EC.focusXRadius*1000, 'linewidth',2, 'color', cols(3,:));

h4 = plot(data_GRIN_radialTop_CinC_EC.incidenceAngle, data_GRIN_radialTop_CinC_EC.focusXRadius*1000, 'linewidth',2, 'color', cols(4,:));

legend([h1 h2 h3 h4], {'Radial-top', 'CinC', 'EC', 'CinC & EC'},'Location','southwest')
ylim([0 75]); ylabel('X-Focus Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:25:75,'XTick',[0 5 10 15 20]);

%% Plot Y Focus Height
YFocusHeightF = figure;

subplot(1,nPlots,1); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, (data_GRIN_radialTop.focusYHeight- data_GRIN_radialTop.coneTipZ)*1000, 'linewidth',2, 'color', cols(1,:));

h2= plot(data_GRIN_radialBase.incidenceAngle, (data_GRIN_radialBase.focusYHeight- data_GRIN_radialBase.coneTipZ)*1000, 'linewidth',2, 'color', cols(2,:));

h3= plot(data_GRIN_linear.incidenceAngle, (data_GRIN_linear.focusYHeight- data_GRIN_linear.coneTipZ)*1000, 'linewidth',2, 'color', cols(3,:));

h4= plot(data_GRIN_both.incidenceAngle, (data_GRIN_both.focusYHeight- data_GRIN_both.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
legend([h1 h2 h3 h4], {'Radial-top','Radial-base','Linear','Radial+linear'},'Location','southwest')
ylim([-250 50]); ylabel('Y-Focus Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -250:50:50,'XTick',[0 5 10 15 20]);

subplot(1,nPlots,2); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, (data_GRIN_radialTop.focusYHeight- data_GRIN_radialTop.coneTipZ)*1000, 'linewidth',2, 'color', cols(1,:));

h2= plot(data_GRIN_radialTop_tipCorr.incidenceAngle, (data_GRIN_radialTop_tipCorr.focusYHeight- data_GRIN_radialTop_tipCorr.coneTipZ)*1000, 'linewidth',2, 'color', cols(2,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
legend([h1 h2], {'Radial-top','W/ Tip Corr.'},'Location','southwest')
ylim([-250 50]); ylabel('Y-Focus Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -250:50:50,'XTick',[0 5 10 15 20]);

subplot(1,nPlots,3); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, (data_GRIN_cylinder.focusYHeight- data_GRIN_cylinder.coneTipZ)*1000, 'linewidth',2, 'color', cols(2,:));

h2 = plot(data_GRIN_cylinder_p36.incidenceAngle, (data_GRIN_cylinder_p36.focusYHeight- data_GRIN_cylinder_p36.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_cylinder_m36.incidenceAngle, (data_GRIN_cylinder_m36.focusYHeight- data_GRIN_cylinder_m36.coneTipZ)*1000, 'linewidth',2, 'color', cols(3,:));

h4 = plot(data_GRIN_cylinder_shape.incidenceAngle, (data_GRIN_cylinder_shape.focusYHeight- data_GRIN_cylinder_shape.coneTipZ)*1000, 'linewidth',2, 'color', cols(5,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
legend([h1 h2 h3 h4], {'Theory-equal', 'Theory-plus', 'Theory-minus', 'Cylinder'},'Location','southwest')
ylim([-250 50]); ylabel('Y-Focus Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -250:50:50,'XTick',[0 5 10 15 20]);

subplot(1,nPlots,4); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, (data_GRIN_radialTop.focusYHeight- data_GRIN_radialTop.coneTipZ)*1000, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_radialTop_CinC.incidenceAngle, (data_GRIN_radialTop_CinC.focusYHeight- data_GRIN_radialTop_CinC.coneTipZ)*1000, 'linewidth',2, 'color', cols(2,:));

h3 = plot(data_GRIN_radialTop_EC.incidenceAngle, (data_GRIN_radialTop_EC.focusYHeight- data_GRIN_radialTop_EC.coneTipZ)*1000, 'linewidth',2, 'color', cols(3,:));

h4 = plot(data_GRIN_radialTop_CinC_EC.incidenceAngle, (data_GRIN_radialTop_CinC_EC.focusYHeight- data_GRIN_radialTop_CinC_EC.coneTipZ)*1000, 'linewidth',2, 'color', cols(4,:));

line([0 20], [0 0], 'color', outlineCol, 'linewidth',2,'linestyle',':')
legend([h1 h2 h3 h4], {'Radial-top', 'CinC', 'EC', 'CinC & EC'},'Location','southwest')
ylim([-250 50]); ylabel('Y-Focus Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -250:50:50,'XTick',[0 5 10 15 20]);

%% Plot Y Focus Radius
YFocusRadiusF = figure;

subplot(1,nPlots,1); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.focusYRadius*1000, 'linewidth',2, 'color', cols(1,:));

h2= plot(data_GRIN_radialBase.incidenceAngle, data_GRIN_radialBase.focusYRadius*1000, 'linewidth',2, 'color', cols(2,:));

h3= plot(data_GRIN_linear.incidenceAngle, data_GRIN_linear.focusYRadius*1000, 'linewidth',2, 'color', cols(3,:));

h4= plot(data_GRIN_both.incidenceAngle, data_GRIN_both.focusYRadius*1000, 'linewidth',2, 'color', cols(4,:));

legend([h1 h2 h3 h4], {'Radial-top','Radial-base','Linear','Radial+linear'},'Location','southwest')
ylim([0 75]); ylabel('Y-Focus Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:25:75,'XTick',[0 5 10 15 20]);

subplot(1,nPlots,2); hold on
h1= plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.focusYRadius*1000, 'linewidth',2, 'color', cols(1,:));

h2= plot(data_GRIN_radialTop_tipCorr.incidenceAngle, data_GRIN_radialTop_tipCorr.focusYRadius*1000, 'linewidth',2, 'color', cols(2,:));

legend([h1 h2], {'Radial-top','W/ Tip Corr.'},'Location','southwest')
ylim([0 75]); ylabel('Y-Focus Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:25:75,'XTick',[0 5 10 15 20]);

subplot(1,nPlots,3); hold on
h1 = plot(data_GRIN_cylinder.incidenceAngle, data_GRIN_cylinder.focusYRadius*1000, 'linewidth',2, 'color', cols(2,:));

h2 = plot(data_GRIN_cylinder_p36.incidenceAngle, data_GRIN_cylinder_p36.focusYRadius*1000, 'linewidth',2, 'color', cols(4,:));

h3 = plot(data_GRIN_cylinder_m36.incidenceAngle, data_GRIN_cylinder_m36.focusYRadius*1000, 'linewidth',2, 'color', cols(3,:));

h4 = plot(data_GRIN_cylinder_shape.incidenceAngle, data_GRIN_cylinder_shape.focusYRadius*1000, 'linewidth',2, 'color', cols(5,:));

legend([h1 h2 h3 h4], {'Theory-equal', 'Theory-plus', 'Theory-minus', 'Cylinder'},'Location','southwest')
ylim([0 75]); ylabel('Y-Focus Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:25:75,'XTick',[0 5 10 15 20]);

subplot(1,nPlots,4); hold on
h1 = plot(data_GRIN_radialTop.incidenceAngle, data_GRIN_radialTop.focusYRadius*1000, 'linewidth',2, 'color', cols(1,:));

h2 = plot(data_GRIN_radialTop_CinC.incidenceAngle, data_GRIN_radialTop_CinC.focusYRadius*1000, 'linewidth',2, 'color', cols(2,:));

h3 = plot(data_GRIN_radialTop_EC.incidenceAngle, data_GRIN_radialTop_EC.focusYRadius*1000, 'linewidth',2, 'color', cols(3,:));

h4 = plot(data_GRIN_radialTop_CinC_EC.incidenceAngle, data_GRIN_radialTop_CinC_EC.focusYRadius*1000, 'linewidth',2, 'color', cols(4,:));

legend([h1 h2 h3 h4], {'Radial-top', 'CinC', 'EC', 'CinC & EC'},'Location','southwest')
ylim([0 75]); ylabel('Y-Focus Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', 0:25:75,'XTick',[0 5 10 15 20]);