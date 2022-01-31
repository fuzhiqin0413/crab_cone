clear; 

saveFlag = 0;

% Varying RI series
dataDirectory = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Data/Varying RI Profile/';
imageDirectory = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Figures/Varying RI Profile/';
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_both_SIMDATA.mat', dataDirectory, imageDirectory, saveFlag)
%%

loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_radialBase_SIMDATA.mat', dataDirectory, imageDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_radialTop_SIMDATA.mat', dataDirectory, imageDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_linear_SIMDATA.mat', dataDirectory, imageDirectory, saveFlag)

loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_radialTop_TipCorrection_SIMDATA.mat', dataDirectory, imageDirectory, saveFlag)

loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_SIMDATA.mat', dataDirectory, imageDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_normalised_max_SIMDATA.mat', dataDirectory, imageDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_normalised_80_SIMDATA.mat', dataDirectory, imageDirectory, saveFlag)

loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_m36_SIMDATA.mat', dataDirectory, imageDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_p36_SIMDATA.mat', dataDirectory, imageDirectory, saveFlag)
loadReplot('Cylinder_1000_nm_Cone_0_SD_GRIN_cylinder_SIMDATA.mat', dataDirectory, imageDirectory, saveFlag)

% Internal structures (w/ varying RI)
dataDirectory = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Data/Extra internal structures/Radial RI';
imageDirectory = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Figures/Extra internal structures/Radial RI';
loadReplot('Cone_CinC_1000_nm_Cone_0_SD_GRIN_radialTop_SIMDATA.mat', dataDirectory, imageDirectory, saveFlag)
loadReplot('Cone_EC_1000_nm_Cone_0_SD_GRIN_radialTop_SIMDATA.mat', dataDirectory, imageDirectory, saveFlag)
loadReplot('Cone_CinC_EC_1000_nm_Cone_0_SD_GRIN_radialTop_SIMDATA.mat', dataDirectory, imageDirectory, saveFlag)

% Varying SD (unfirom RI) series
dataDirectory = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Data/Varying Cone SD';
imageDirectory = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Figures/Varying Cone SD';
loadReplot('Cone_1000_nm_Cone_-4_SD_Uniform_1.52_SIMDATA.mat', dataDirectory, imageDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_-2_SD_Uniform_1.52_SIMDATA.mat', dataDirectory, imageDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_Uniform_1.52_SIMDATA.mat', dataDirectory, imageDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_2_SD_Uniform_1.52_SIMDATA.mat', dataDirectory, imageDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_4_SD_Uniform_1.52_SIMDATA.mat', dataDirectory, imageDirectory, saveFlag)

% Internal structures (w/ uniform RI)
dataDirectory = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Data/Extra internal structures/Uniform RI';
imageDirectory = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Figures/Extra internal structures/Uniform RI';
loadReplot('Cone_CinC_1000_nm_Cone_0_SD_Uniform_1.52_SIMDATA.mat', dataDirectory, imageDirectory, saveFlag)
loadReplot('Cone_EC_1000_nm_Cone_0_SD_Uniform_1.52_SIMDATA.mat', dataDirectory, imageDirectory, saveFlag)
loadReplot('Cone_CinC_EC_1000_nm_Cone_0_SD_Uniform_1.52_SIMDATA.mat', dataDirectory, imageDirectory, saveFlag)

% Varying intercone
dataDirectory = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Data/Variable Intercone/Radial RI';
imageDirectory = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Figures/Variable Intercone/Radial RI';
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_radialTop_IC_142_SIMDATA.mat', dataDirectory, imageDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_radialTop_IC_144_SIMDATA.mat', dataDirectory, imageDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_radialTop_IC_146_SIMDATA.mat', dataDirectory, imageDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_radialTop_IC_148_SIMDATA.mat', dataDirectory, imageDirectory, saveFlag)

dataDirectory = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Data/Variable Intercone/Uniform RI';
imageDirectory = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Figures/Variable Intercone/Uniform RI';
loadReplot('Cone_1000_nm_Cone_0_SD_Uniform_1.52_IC_142_SIMDATA.mat', dataDirectory, imageDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_Uniform_1.52_IC_144_SIMDATA.mat', dataDirectory, imageDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_Uniform_1.52_IC_146_SIMDATA.mat', dataDirectory, imageDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_Uniform_1.52_IC_148_SIMDATA.mat', dataDirectory, imageDirectory, saveFlag)

dataDirectory = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Data/Variable Intercone/Theoretical RI';
imageDirectory = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Figures/Variable Intercone/Theoretical RI';
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_IC_142_SIMDATA.mat', dataDirectory, imageDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_IC_144_SIMDATA.mat', dataDirectory, imageDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_IC_146_SIMDATA.mat', dataDirectory, imageDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_IC_148_SIMDATA.mat', dataDirectory, imageDirectory, saveFlag)

% Varying pigment
dataDirectory = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Data/Variable Pigment/Radial RI';
imageDirectory = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Figures/Variable Pigment/Radial RI';
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_radialTop_PG_140_SIMDATA.mat', dataDirectory, imageDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_radialTop_PG_145_SIMDATA.mat', dataDirectory, imageDirectory, saveFlag)
% loadReplot(filename, dataDirectory, imageDirectory, saveFlag)

dataDirectory = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Data/Variable Pigment/Uniform RI';
imageDirectory = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Figures/Variable Pigment/Uniform RI';
% loadReplot(filename, dataDirectory, imageDirectory, saveFlag)
% loadReplot(filename, dataDirectory, imageDirectory, saveFlag)
% loadReplot(filename, dataDirectory, imageDirectory, saveFlag)

dataDirectory = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Data/Variable Pigment/Theoretical RI';
imageDirectory = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/Figures/Variable Pigment/Theoretical RI ';
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_PG_140_SIMDATA.mat', dataDirectory, imageDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_PG_145_SIMDATA.mat', dataDirectory, imageDirectory, saveFlag)
% loadReplot(filename, dataDirectory, imageDirectory, saveFlag)


function loadReplot(filename, dataDirectory, imageDirectory, saveFlag)
close all; clc

cd(dataDirectory)
data = load(filename);

% Pull fields from data structure to direct workspace variables
fields = fieldnames(data);

for i = 1:length(fields)
    eval(sprintf('%s = data.%s;', fields{i}, fields{i}))
end

receptorRadius = NaN;
acceptanceUsingReceptor = NaN;

limitCOLCToExterior = 1;

% All from schematic Chamberlain and Barlow 1987
receptorRadiusDay = 30/2;
receptorRadiusNight = 50/2;

receptorDistanceDay = 70;
receptorDistanceNight = 5;

apertureRadiusDay = 5/2;
apertureRadiusNight = 60/2;

dayNarrowStart = 10;
dayNarrowEnd = 55;

%%% Do values need to be stripped out of data strucutre?

plot_cone_data

if saveFlag
   %%% Add specific new/revised values for COLC, X/Y Focus and divergence to data files
    
    cd(imageDirectory)

    %%% Plot over exisiting figures
    % Remove
    %acceptanceAngle, acceptancePercentage, receptorRadius, acceptanceUsingReceptor
    
    % Add - variables above, and:
    % acceptancePercentageNight, acceptanceAngleNight
    % acceptancePercentageDay, acceptanceAngleDay
    % acceptancePercentageColc, acceptanceAngleColc
    % beamDivergenceMax, beamDivergenceAverage
end
end