clear; 

saveFlag = 1;

rootDirectory = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/';

%%
% Varying RI series
seriesDirectory = '/Varying RI Profile/';
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_both_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_radialBase_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_radialTop_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_linear_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)

loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_radialTop_TipCorrection_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)

loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_normalised_max_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_normalised_80_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)

loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_m36_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_p36_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)
loadReplot('Cylinder_1000_nm_Cone_0_SD_GRIN_cylinder_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)

% Internal structures (w/ varying RI)
seriesDirectory = '/Extra internal structures/Radial RI/';
loadReplot('Cone_CinC_1000_nm_Cone_0_SD_GRIN_radialTop_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)
loadReplot('Cone_EC_1000_nm_Cone_0_SD_GRIN_radialTop_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)
loadReplot('Cone_CinC_EC_1000_nm_Cone_0_SD_GRIN_radialTop_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)

% Varying SD (unfirom RI) series
seriesDirectory = '/Varying Cone SD/';
loadReplot('Cone_1000_nm_Cone_-4_SD_Uniform_1.52_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_-2_SD_Uniform_1.52_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_Uniform_1.52_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_2_SD_Uniform_1.52_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_4_SD_Uniform_1.52_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)

% Internal structures (w/ uniform RI)
seriesDirectory = '/Extra internal structures/Uniform RI/';
loadReplot('Cone_CinC_1000_nm_Cone_0_SD_Uniform_1.52_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)
loadReplot('Cone_EC_1000_nm_Cone_0_SD_Uniform_1.52_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)
loadReplot('Cone_CinC_EC_1000_nm_Cone_0_SD_Uniform_1.52_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)

% Varying intercone
seriesDirectory = '/Variable Intercone/Radial RI/';
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_radialTop_IC_142_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_radialTop_IC_144_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_radialTop_IC_146_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_radialTop_IC_148_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)

seriesDirectory = '/Variable Intercone/Uniform RI/';
loadReplot('Cone_1000_nm_Cone_0_SD_Uniform_1.52_IC_142_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_Uniform_1.52_IC_144_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_Uniform_1.52_IC_146_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_Uniform_1.52_IC_148_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)

seriesDirectory = '/Variable Intercone/Theoretical RI/';
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_IC_142_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_IC_144_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_IC_146_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_IC_148_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)

% Varying pigment
seriesDirectory = '/Variable Pigment/Radial RI/';
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_radialTop_PG_140_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_radialTop_PG_145_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_radialTop_PG_150_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)

seriesDirectory = '/Variable Pigment/Uniform RI/';
loadReplot('Cone_1000_nm_Cone_0_SD_Uniform_1.52_PG_140_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_Uniform_1.52_PG_145_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_Uniform_1.52_PG_150_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)

seriesDirectory = '/Variable Pigment/Theoretical RI/';
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_PG_140_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_PG_145_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)
loadReplot('Cone_1000_nm_Cone_0_SD_GRIN_cylinder_PG_150_SIMDATA.mat', rootDirectory, seriesDirectory, saveFlag)


function loadReplot(filename, rootDirectory, seriesDirectory, saveFlag)
    close all; clc

    cd(sprintf('%sData%s',rootDirectory,seriesDirectory))
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

    % check if some variables that were not included on all files need to be added
    if exist('externalPigmentValue','var') == 0 
        externalPigmentValue = NaN;
    end
    
    if exist('firstIntersectCells','var') == 0
        firstIntersectCells = NaN;
    end
    
    plot_cone_data

    if saveFlag

        % save figures
        cd(sprintf('%sFigures_new%s',rootDirectory,seriesDirectory))

        imageFile = sprintf('%sSpotDiagram.pdf', filename(1:end-11));
        exportgraphics(spotF,imageFile,'ContentType','vector') 

        imageFile = sprintf('%sRaysX.pdf', filename(1:end-11));
        exportgraphics(rayFX,imageFile,'ContentType','vector')

        imageFile = sprintf('%sRaysY.pdf', filename(1:end-11));
        exportgraphics(rayFY,imageFile,'ContentType','vector')

        imageFile = sprintf('%sCOLCLines.pdf', filename(1:end-11));
        exportgraphics(tipF,imageFile,'ContentType','vector')

        imageFile = sprintf('%sSummary.pdf', filename(1:end-11));
        exportgraphics(sumF,imageFile,'ContentType','vector')

        % save data    
        cd(sprintf('%sData_new%s',rootDirectory,seriesDirectory))

        save(filename, 'rayPathCells', 'finalIntersectCells', 'finalRayCells', 'finalRayTRefractCells', 'rayReverseCells', 'timePerAngle', 'TIRFlagCells', 'firstIntersectCells', ...
            'dataFile', 'metaFile', 'incidenceAngle', 'receptorAcceptance', 'exposedHeight', 'blockExposedRentry', ...
            'xSpacing', 'plotSpacing', 'interpType', 'tolerance', 'initialDeltaS', 'iterativeFinal', 'epsilon',  'interfaceRefraction', 'blockMultipleExits', 'limitToConeBase', 'clearReverseRays', 'trace3D', ...
            'alphaForIntercone', 'alphaForCone', 'alphaForCinC', 'dilateBorderRadius', 'flipVolume', 'flipSurfaces', ...
            'focusXHeight', 'focusXRadius', 'focusYHeight', 'focusYRadius', 'colcHeight', 'colcRadius', 'rayReverseNum', 'TIRPercengtage', ...
            'nOrigins', 'plotLineCols', 'colorTIR', 'plotOrigins', 'scaleBarsOnRayDiagram', 'plotColorsOnSummary', 'justPlotCenterRays', 'raysOnXZPlane', 'raysOnYZPlane', 'plotRIImageOnRayDiagram', ...
            'coneTipZ', 'corneaZ', 'epicorneaZ', 'coneBaseZ', 'interConeValueUsed', 'externalPigmentValue', 'zSteps', ...
            'coneProfileR', 'coneProfileZ', 'cInCProfileR', 'cInCProfileZ', 'epicorneaProfileR', 'epicorneaProfileZ', 'interconeProfileR', 'interconeProfileZ', ...
            'limitCOLCToExterior', 'receptorRadiusDay', 'receptorRadiusNight', 'receptorDistanceDay', 'receptorDistanceNight', 'apertureRadiusDay', 'apertureRadiusNight', 'dayNarrowStart', 'dayNarrowEnd', ...
            'acceptancePercentageNight', 'acceptanceAngleNight', 'acceptancePercentageDay', 'acceptanceAngleDay', 'acceptancePercentageColc', 'acceptanceAngleColc', 'beamDivergenceMax', 'beamDivergenceAverage', ...
            'acceptanceNumNight', 'acceptanceNumDay', 'acceptanceNumColC', ... 
            'voxelSize', 'volumeSize', 'xStartPointsOrig', 'rayOrigins', 'metaData')
    end
end