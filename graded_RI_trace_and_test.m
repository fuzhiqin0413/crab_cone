% Set up to do both tracing on data and have testing mode on either
% lundaberg lens or graded fibre, as in nishdate 2011 papers
    
% clear; 
clearvars -except dataFolder metaFile interConeValueToUse externalPigmentValue
clc; close all

%% Set parameters

useRealData = 1;

%%% Change back
incidenceAngle = 0:2.5:20; [0 5 10];  % deg, in XZ - plane    

% interConeValueToUse = 1.40;

%%% Change back
xSpacing = 5; % in voxels...
    plotSpacing = 2; % Integer mutiple of xSpaing to plot 

if useRealData

    analysisFolder = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisResults/';
    saveData = 1; %%%
    saveFigures = 1; %%%
    
    %%% For 2 micron data
    % Varying RI
%     dataFolder = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisVolumes/2 micron/Varying RI Profile/';
%     
%     metaFile = 'Cone_1000_nm_Cone_0_SD_GRIN_both.mat';
%     metaFile = 'Cone_1000_nm_Cone_0_SD_GRIN_cylinder.mat';
%     metaFile = 'Cone_1000_nm_Cone_0_SD_GRIN_linear.mat';
%     metaFile = 'Cone_1000_nm_Cone_0_SD_GRIN_radialTop.mat';
%     metaFile = 'Cone_1000_nm_Cone_0_SD_GRIN_radialTop_TipCorrection.mat';

    %%% For 5 micron testing
%     dataFolder = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisVolumes/5 micron/';

%     metaFile = 'Cylinder_1000_nm_Cone_0_SD_GRIN_cylinder.mat';
%     metaFile = 'Cone_1000_nm_Cone_0_SD_GRIN_both.mat';
%     metaFile = 'Cone_1000_nm_Cone_0_SD_Uniform_1.52.mat';
%     metaFile = 'Cone_CinC_EC_1000_nm_Cone_0_SD_GRIN_both.mat';

    dataFile = sprintf('Volume_%s', metaFile);
    
    % Need to set so these don't have holes - note they are in voxels
        % A low value also often causes planes across middle of cone
        % For 5 micron test, 10, 50 and 10 were ok
        % Going to 2 micron doubling seems ok. Could also be influenced by angle step on profiles (2 for all)
    alphaForIntercone = 20;
    alphaForCone = 100;
    alphaForCinC = 20;
    
    dilateBorderRadius = 3; %1 - was ok for 5 micron
    
    % Usually saved pointing down
    flipVolume = 1;    
    flipSurfaces = 1;
    
    createGradedFiber = 0;
    createLunebergLens = 0;
    plotReferenceFiber = 0;
    plotRIImageOnRayDiagram = 0;
    useTestData = 0;
    
    receptorRadius = 30; % micron - night basically 0 out, 5 at 70 out for for day
    receptorAcceptance = 90; % degrees
    
    exposedHeight = 15; % micron - night, maybe around 20 for day light
    blockExposedRentry = 1;
        
    testOrigin = []; [5 6 7 8];
else
    useTestData = 1;
    
    radius = 1; %mm - used for both lens and fiber
    voxelSize = 0.065; % mm

    exteriorRI = 1;

    createLunebergLens = 1;
        
    createGradedFiber = 0;
        %%% Add flag to create gaussian profile then update solution and period
        fiberLength = 10;
        n0 = sqrt(2.5); 
        alpha = 1/sqrt(2.5); 0.7; 
        beta = alpha*n0;
        plotReferenceFiber = 0; % will plot data from Nishidate and only trace 0.5
        
    blockExposedRentry = 0;
end

trace3D = 1; % Otherwise just line across X

% Options, 1, 4S, 4RKN, 5RKN, 45RKN 
% Difference 4 to 5 is quite small
interpType = '4S'; %'4S';
    %This doesn't have big effect on error < 10^-9, ok even at 10^-6 but then collapses
        % May effect peripheral rays more
    tolerance = 10^-9; % Nishidate uses 10^-12, seems a bit too tight
    
    % Seems good to start a bit lower than the general deltaS
        % in fixed step, 10^-3 vs. 10^-4 gives rough order of magnitude error
%     initialDeltaS = 0.5*10^-3;
    % Need for uniform
    initialDeltaS = 0.5*10^-3;

    % This keeps spot size surpsingly small, even on loose tolerance
        % (I guess it ends up stepping further at loose tolerance...)
    iterativeFinal = 1;  
    
%error term, sqrt of machine precision, from Nishidate paper sqrt(1.49e-8)
    % Take geometric mean of single and double precision, double takes ages to evaluate
epsilon = sqrt(geomean([1.19e-7 2.22e-16])); sqrt(1.49e-8); %sqrt(geomean([1.19e-7 2.22e-16])); 

% Should be on, but can switch off to test    
interfaceRefraction = 1;  

% Ray that has exited will error on re-entry
blockMultipleExits = 0;

% Rays won't continue if it hits intercone base - effectively makes base of cone the aperture
limitToConeBase = 1;

% extend on final plots - only added for real data
%%% Now done as extension of zsteps
% extendRayLength = 1; % mm

% Best to do this as I haven't really delt with them yet
clearReverseRays = 1;

plotLineCols = 0;
colorTIR = 1; % Cant use both but not tested...
plotColorsOnSummary = 0;

justPlotCenterRays = 0;
scaleBarsOnRayDiagram = 1;
acceptanceUsingReceptor = 1;

testPlot = 0;
testPlotSurface = 0;

loopLim = 10000;

%% Load or create test data     
if useRealData
    temp = load(sprintf('%s%s',dataFolder, dataFile));
    lensRIVolume = temp.volume;
    
    clear temp
    
    metaData = load(sprintf('%s%s',dataFolder, metaFile));
    voxelSize = metaData.voxSize3D/1000; %um to mm
    
    % Flip to change z direction
    if flipVolume
        lensRIVolume = permute(flipud(permute(lensRIVolume,[3 1 2])), [2 3 1]);
    end
    
    volumeSize = size(lensRIVolume);
    
    % Using negative coding for GRIN
    if metaData.interconeValue > 0
        RIFlagVolume = lensRIVolume < 0;
        interConeValueUsed = metaData.interconeValue;
    else
        RIFlagVolume = lensRIVolume < -1;
        
        lensRIVolume(lensRIVolume == metaData.interconeValue) = interConeValueToUse;
        interConeValueUsed = interConeValueToUse;
    end
    
    if isnumeric(metaData.coneValue)
       fixedConeRI = 1;
    else
       fixedConeRI = 0; 
    end
    
    lensRIVolume(RIFlagVolume) = -lensRIVolume(RIFlagVolume);

    % make meshes for each layer
    meshAngles = (0:2:360)/180*pi;

    %%% Note, border cones won't be defined using this process
        % Add an error if they are present (rotated volume won't work anyway)
    
    % For cornea - flat
    % single Z at end
    corneaZ = (metaData.coneLengthToUse + metaData.outerCorneaLengthToUse + metaData.epicorneaLengthToUse)* ...
        metaData.voxSize/metaData.voxSize3D;

    [xGrid, yGrid, zGrid] = meshgrid(1:volumeSize(1), 1:volumeSize(2), corneaZ);
    xGrid = xGrid(:); yGrid = yGrid(:); zGrid = zGrid(:);
    
    corneaVertices = [xGrid, yGrid, zGrid];
    corneaVertices(:,3) = corneaVertices(:,3) + metaData.tipOffset*metaData.voxSize/metaData.voxSize3D;
    
    % Get 2D triangulation
    tempTriangulation = delaunayTriangulation(corneaVertices(:,1:2));

    corneaSurface.faces = tempTriangulation.ConnectivityList;
    corneaSurface.vertices = corneaVertices;

    corneaZ = corneaZ + metaData.tipOffset*metaData.voxSize/metaData.voxSize3D;
    
    if flipSurfaces
        corneaSurface.vertices(:,3) = -(corneaSurface.vertices(:,3)-volumeSize(3)/2)+volumeSize(3)/2;
        
        corneaZ = -(corneaZ-volumeSize(3)/2)+volumeSize(3)/2;
    end
    corneaBorderVolume = polygon2voxel(corneaSurface, volumeSize, 'none');
    
    % for top of outer cornea - can be curved
    if metaData.displayProfiles.EpicorneaCone

        % has a curve
        epicorneaProfileR = metaData.epicorneaProfileToUse*metaData.voxSize/metaData.voxSize3D;
        epicorneaProfileZ = metaData.corneaXRef*metaData.voxSize/metaData.voxSize3D;
        
        epicorneaProfileZ(isnan(epicorneaProfileR)) = [];
        epicorneaProfileR(isnan(epicorneaProfileR)) = [];

        % Adjust to actual length
        epicorneaProfileZ(1) = 0;
        epicorneaProfileZ = epicorneaProfileZ + (metaData.coneLengthToUse + metaData.outerCorneaLengthToUse)*metaData.voxSize/metaData.voxSize3D;
        
        epicorneaVertices = zeros(length(epicorneaProfileR)*length(meshAngles),3);

        profileMax = max(epicorneaProfileR);
        
        for i = 1:length(epicorneaProfileR)
            for j = 1:length(meshAngles)
                epicorneaVertices((i-1)*length(meshAngles)+j,:) = [epicorneaProfileR(i)*sin(meshAngles(j)), ...
                    epicorneaProfileR(i)*cos(meshAngles(j)), epicorneaProfileZ(i)];
            end
        end
        
        % Get cap to pad front
        [xGrid, yGrid, zGrid] = meshgrid(-floor(profileMax):ceil(profileMax), -floor(profileMax):ceil(profileMax), epicorneaProfileZ(end));

        xGrid = xGrid(:); yGrid = yGrid(:); zGrid = zGrid(:);
        rGrid = sqrt(xGrid.^2 + yGrid.^2);

        frontInds = find(rGrid < epicorneaProfileR(end));

        epicorneaVertices = [epicorneaVertices', [xGrid(frontInds), yGrid(frontInds), zGrid(frontInds)]']';

        epicorneaVertices(:,1) = epicorneaVertices(:,1) + volumeSize(1)/2;
        epicorneaVertices(:,2) = epicorneaVertices(:,2) + volumeSize(2)/2;
        epicorneaVertices(:,3) = epicorneaVertices(:,3) + metaData.tipOffset*metaData.voxSize/metaData.voxSize3D;
        
        % Get ring to pad back
        [xGrid, yGrid, zGrid] = meshgrid(1:volumeSize(1), 1:volumeSize(2), epicorneaProfileZ(1) + metaData.tipOffset*metaData.voxSize/metaData.voxSize3D);
        xGrid = xGrid(:); yGrid = yGrid(:); zGrid = zGrid(:);

        rGrid = sqrt((xGrid-volumeSize(1)/2).^2 + (yGrid-volumeSize(2)/2).^2);
    
        endInds = find(rGrid > epicorneaProfileR(1));
        
        epicorneaVertices = [[xGrid(endInds), yGrid(endInds), zGrid(endInds)]', epicorneaVertices']';
        
        % Make another layer at base to force 3D
        [xGrid, yGrid, zGrid] = meshgrid(1:volumeSize(1), 1:volumeSize(2), epicorneaProfileZ(1)-1 + metaData.tipOffset*metaData.voxSize/metaData.voxSize3D);
        xGrid = xGrid(:); yGrid = yGrid(:); zGrid = zGrid(:);
        
        epicorneaVertices = [epicorneaVertices' [xGrid, yGrid, zGrid]' ]';
        forcedInds = find(epicorneaVertices(:,3) == epicorneaProfileZ(1)-1 + metaData.tipOffset*metaData.voxSize/metaData.voxSize3D);
        
        tempShape = alphaShape(epicorneaVertices, 10);
        [tempFaces, tempVertices] = boundaryFacets(tempShape);
        
        if ~isempty(forcedInds)
            % remove forced inds at base
            % not all of the original forced inds are used in the final trinagulation
            usedForcedInds = find(tempVertices(:,3) == epicorneaVertices(forcedInds(1),3));

            facesToRemove = zeros(length(tempFaces),1, 'logical');
            for i = 1:length(usedForcedInds)
                % Flag for each column
                facesToRemove(tempFaces(:,1) == usedForcedInds(i)) = 1;
                facesToRemove(tempFaces(:,2) == usedForcedInds(i)) = 1;
                facesToRemove(tempFaces(:,3) == usedForcedInds(i)) = 1;
            end

            tempFaces(facesToRemove,:) = [];
            
            % Cant remove from vertices or it will screw up face referencing
        %     tempVertices(usedForcedInds,:) = [];

            epicorneaVertices(forcedInds,:) = [];
        end
        
        epicorneaSurface.faces = tempFaces;
        epicorneaSurface.vertices = tempVertices;
        
        epicorneaProfileZ = epicorneaProfileZ + metaData.tipOffset*metaData.voxSize/metaData.voxSize3D;
        
        if flipSurfaces
            epicorneaProfileZ = -(epicorneaProfileZ-volumeSize(3)/2)+volumeSize(3)/2;
        end
        
    else
       % just flat, can copy cornea with adjusted Z 
       epicorneaVertices = corneaVertices; 
       epicorneaVertices(:,3) = epicorneaVertices(:,3) - metaData.epicorneaLengthToUse* metaData.voxSize/metaData.voxSize3D;
       
       epicorneaSurface = corneaSurface;
       epicorneaSurface.vertices = epicorneaVertices;
       
       epicorneaProfileR = [];
       epicorneaProfileZ = [];
    end
    
    epicorneaZ = (metaData.coneLengthToUse + metaData.outerCorneaLengthToUse)* metaData.voxSize/metaData.voxSize3D;
           
    epicorneaZ = epicorneaZ + metaData.tipOffset*metaData.voxSize/metaData.voxSize3D;
    
    if flipSurfaces
        epicorneaSurface.vertices(:,3) = -(epicorneaSurface.vertices(:,3)-volumeSize(3)/2)+volumeSize(3)/2;
        
        epicorneaZ = -(epicorneaZ-volumeSize(3)/2)+volumeSize(3)/2;
    end
    
    epicorneaBorderVolume = polygon2voxel(epicorneaSurface, volumeSize, 'none');

    
    % for base outer cornea/intercone around cone -  flat
    % Make flat grid as for cornea
    coneBaseZ = (metaData.coneLengthToUse)* metaData.voxSize/metaData.voxSize3D;

    [xGrid, yGrid, zGrid] = meshgrid((1:volumeSize(1))-volumeSize(1)/2, (1:volumeSize(2))-volumeSize(2)/2, coneBaseZ);
    xGrid = xGrid(:); yGrid = yGrid(:); zGrid = zGrid(:);
    
    % Shrink vertices that are inside radius at base of cone
    rGrid = sqrt(xGrid.^2 + yGrid.^2);
    
    if ~isnan(metaData.coneProfileToUse(end)) & ~metaData.createCylinder
        coneBaseRadius = metaData.coneProfileToUse(end)*metaData.voxSize/metaData.voxSize3D;
        
    elseif isnan(metaData.coneProfileToUse(end)) & ~metaData.createCylinder
        error('Need to get last good value in cone profile')
    elseif metaData.createCylinder
        coneBaseRadius = max(metaData.coneProfileToUse)*metaData.voxSize/metaData.voxSize3D;
    end
    
    verticesToRemove = find(rGrid < coneBaseRadius);
    
    interconeBaseVertices = [xGrid, yGrid, zGrid];
    
    % Step in slightly so these don't get incorperated into faces along base of cone
    interconeBaseVertices(verticesToRemove, 1:2) = interconeBaseVertices(verticesToRemove, 1:2)*0.9;
    
    % Add ring of vertices at cone base
    interconeBaseVertices = [interconeBaseVertices' [coneBaseRadius*sin(meshAngles)' coneBaseRadius*cos(meshAngles)' coneBaseZ*ones(length(meshAngles),1)]']';
    
    interconeBaseVertices(:,1) = interconeBaseVertices(:,1) + volumeSize(1)/2;
    interconeBaseVertices(:,2) = interconeBaseVertices(:,2) + volumeSize(2)/2;
    interconeBaseVertices(:,3) = interconeBaseVertices(:,3) + metaData.tipOffset*metaData.voxSize/metaData.voxSize3D;
   
    % Get 2D triangulation
    tempTriangulation = delaunayTriangulation(interconeBaseVertices(:,1:2));
    
    interconeBaseSurface.faces = tempTriangulation.ConnectivityList;
    interconeBaseSurface.vertices = interconeBaseVertices;
    
    facesToRemove = zeros(length(interconeBaseSurface.faces),1,'logical');
    for i = 1:length(verticesToRemove)
        % Flag for each column
        facesToRemove(interconeBaseSurface.faces(:,1) == verticesToRemove(i)) = 1;
        facesToRemove(interconeBaseSurface.faces(:,2) == verticesToRemove(i)) = 1;
        facesToRemove(interconeBaseSurface.faces(:,3) == verticesToRemove(i)) = 1;
    end
    
    interconeBaseSurface.faces(facesToRemove,:) = [];

    coneBaseZ = coneBaseZ + metaData.tipOffset*metaData.voxSize/metaData.voxSize3D;
    
    if flipSurfaces
        interconeBaseSurface.vertices(:,3) = -(interconeBaseSurface.vertices(:,3)-volumeSize(3)/2)+volumeSize(3)/2;
        
        coneBaseZ = -(coneBaseZ-volumeSize(3)/2)+volumeSize(3)/2;
    end
    interconeBaseBorderVolume = polygon2voxel(interconeBaseSurface, volumeSize, 'none');


    
    if metaData.displayProfiles.CinC
        if metaData.createCylinder
            error('Not coded') 
        end
        
        % Make CinC profile at base of cone
        cInCProfileR = metaData.CinCProfileToUse*metaData.voxSize/metaData.voxSize3D;
        cInCProfileZ = metaData.coneXRef*metaData.voxSize/metaData.voxSize3D;

        cInCProfileZ(isnan(cInCProfileR)) = [];
        cInCProfileR(isnan(cInCProfileR)) = [];

        % Adjust to actual length
        cInCProfileZ(end) = metaData.coneLengthToUse*metaData.voxSize/metaData.voxSize3D;

        cInCVertices = zeros(length(cInCProfileR)*length(meshAngles),3);

        for i = 1:length(cInCProfileR)
            for j = 1:length(meshAngles)
                cInCVertices((i-1)*length(meshAngles)+j,:) = [cInCProfileR(i)*sin(meshAngles(j)), ...
                    cInCProfileR(i)*cos(meshAngles(j)), cInCProfileZ(i)];
            end  
        end
        
        % Cap to pad front and suround sides of cinc (within cone base)
        [xGrid, yGrid, zGrid] = meshgrid(-coneBaseRadius:coneBaseRadius, -coneBaseRadius:coneBaseRadius, 1);

        xGrid = xGrid(:); yGrid = yGrid(:); zGrid = zGrid(:);
        rGrid = sqrt(xGrid.^2 + yGrid.^2);
        
        firstInds = find(rGrid < cInCProfileR(1));
        cInCVertices = [[xGrid(firstInds), yGrid(firstInds), (cInCProfileZ(1))*zGrid(firstInds)]', cInCVertices']'; 
        
        sideInds = find(rGrid > cInCProfileR(end) & rGrid < coneBaseRadius);
        cInCVertices = [[xGrid(sideInds), yGrid(sideInds), (cInCProfileZ(end))*zGrid(sideInds)]', cInCVertices']'; 
        
        % Add ring of vertices at cone base
        cInCVertices = [cInCVertices' [coneBaseRadius*sin(meshAngles)' coneBaseRadius*cos(meshAngles)' cInCProfileZ(end)*ones(length(meshAngles),1)]']';

        % Add ring and base stepped down
        baseInds = find(rGrid < coneBaseRadius);
        cInCVertices = [[xGrid(baseInds), yGrid(baseInds), (cInCProfileZ(end)+1)*zGrid(baseInds)]', cInCVertices']'; 
        
        cInCVertices = [cInCVertices' [coneBaseRadius*sin(meshAngles)' coneBaseRadius*cos(meshAngles)' (cInCProfileZ(end)+1)*ones(length(meshAngles),1)]']';
        forcedInds = find(cInCVertices(:,3) == (cInCProfileZ(end)+1));
        
        cInCVertices(:,1) = cInCVertices(:,1) + volumeSize(1)/2;
        cInCVertices(:,2) = cInCVertices(:,2) + volumeSize(2)/2;
        cInCVertices(:,3) = cInCVertices(:,3) + metaData.tipOffset*metaData.voxSize/metaData.voxSize3D;
        
        tempShape = alphaShape(cInCVertices, alphaForCinC);
        [tempFaces, tempVertices] = boundaryFacets(tempShape);
        
        if ~isempty(forcedInds)
            % remove forced inds at base
            % not all of the original forced inds are used in the final trinagulation
            usedForcedInds = find(tempVertices(:,3) == cInCVertices(forcedInds(1),3));

            facesToRemove = zeros(length(tempFaces),1, 'logical');
            for i = 1:length(usedForcedInds)
                % Flag for each column
                facesToRemove(tempFaces(:,1) == usedForcedInds(i)) = 1;
                facesToRemove(tempFaces(:,2) == usedForcedInds(i)) = 1;
                facesToRemove(tempFaces(:,3) == usedForcedInds(i)) = 1;
            end

            tempFaces(facesToRemove,:) = [];

            cInCVertices(forcedInds,:) = [];
        end
        
        cInCSurface.faces = tempFaces;
        cInCSurface.vertices = tempVertices;
        
        cInCProfileZ = cInCProfileZ + metaData.tipOffset*metaData.voxSize/metaData.voxSize3D;
        
        if flipSurfaces
            cInCProfileZ = -(cInCProfileZ-volumeSize(3)/2)+volumeSize(3)/2;
        end
        
    else
        % use a flat cInC - basically as for intercone base
        tempZ = (metaData.coneLengthToUse)* metaData.voxSize/metaData.voxSize3D;
        
        [xGrid, yGrid, zGrid] = meshgrid((1:volumeSize(1))-volumeSize(1)/2, (1:volumeSize(2))-volumeSize(2)/2, tempZ);
        xGrid = xGrid(:); yGrid = yGrid(:); zGrid = zGrid(:);

        % Shrink vertices that are inside radius at base of cone
        rGrid = sqrt(xGrid.^2 + yGrid.^2);

        % Remove some at start so this doesn't spread to wide
        cInCVertices = [xGrid, yGrid, zGrid]; 
        cInCVertices(rGrid > coneBaseRadius*1.25, :) = []; rGrid(rGrid > coneBaseRadius*1.25) = [];
        
        verticesToRemove = find(rGrid > coneBaseRadius);
        
        % Step out slightly so these don't get incorperated into faces along base of cone
        cInCVertices(verticesToRemove, 1:2) = cInCVertices(verticesToRemove, 1:2)*1.1;

        % Add ring of vertices at cone base
        cInCVertices = [cInCVertices' [coneBaseRadius*sin(meshAngles)' coneBaseRadius*cos(meshAngles)' tempZ*ones(length(meshAngles),1)]']';

        cInCVertices(:,1) = cInCVertices(:,1) + volumeSize(1)/2;
        cInCVertices(:,2) = cInCVertices(:,2) + volumeSize(2)/2;
        cInCVertices(:,3) = cInCVertices(:,3) + metaData.tipOffset*metaData.voxSize/metaData.voxSize3D;
   
        % Get 2D triangulation
        tempTriangulation = delaunayTriangulation(cInCVertices(:,1:2));

        cInCSurface.faces = tempTriangulation.ConnectivityList;
        cInCSurface.vertices = cInCVertices;

        facesToRemove = zeros(length(cInCSurface.faces),1,'logical');
        for i = 1:length(verticesToRemove)
            % Flag for each column
            facesToRemove(cInCSurface.faces(:,1) == verticesToRemove(i)) = 1;
            facesToRemove(cInCSurface.faces(:,2) == verticesToRemove(i)) = 1;
            facesToRemove(cInCSurface.faces(:,3) == verticesToRemove(i)) = 1;
        end

        cInCSurface.faces(facesToRemove,:) = [];
        
        cInCProfileR = [];
        cInCProfileZ = [];
    end
    
    if flipSurfaces
        cInCSurface.vertices(:,3) = -(cInCSurface.vertices(:,3)-volumeSize(3)/2)+volumeSize(3)/2;
    end
    
    
    
    %%% Would be good to make a funciton for doing the rotation in these...
    % for cone - is curved note voxSize is actually 2D pixel size
    coneProfileR = metaData.coneProfileToUse*metaData.voxSize/metaData.voxSize3D;
    coneProfileZ = metaData.coneXRef*metaData.voxSize/metaData.voxSize3D;

    coneProfileZ(isnan(coneProfileR)) = [];
    coneProfileR(isnan(coneProfileR)) = [];
    
    % Adjust to actual length
    coneProfileZ(end) = metaData.coneLengthToUse*metaData.voxSize/metaData.voxSize3D;
    
    profileMax = max(coneProfileR);

    coneVertices = zeros(length(coneProfileR)*length(meshAngles),3);

    for i = 1:length(coneProfileR)
        for j = 1:length(meshAngles)
            if ~metaData.createCylinder
                coneVertices((i-1)*length(meshAngles)+j,:) = [coneProfileR(i)*sin(meshAngles(j)), ...
                    coneProfileR(i)*cos(meshAngles(j)), coneProfileZ(i)];
            else
                coneVertices((i-1)*length(meshAngles)+j,:) = [profileMax*sin(meshAngles(j)), ...
                    profileMax*cos(meshAngles(j)), coneProfileZ(i)];
            end
        end
    end
    
    % Get caps to pad front and back, otherwise there are degenerate faces
    [xGrid, yGrid, zGrid] = meshgrid(-floor(profileMax):ceil(profileMax), -floor(profileMax):ceil(profileMax), 1);

    xGrid = xGrid(:); yGrid = yGrid(:); zGrid = zGrid(:);
    rGrid = sqrt(xGrid.^2 + yGrid.^2);

    if ~metaData.createCylinder
        firstInds = find(rGrid < coneProfileR(1));
    else
        firstInds = find(rGrid < profileMax);
    end
    
    coneVertices = [ [xGrid(firstInds), yGrid(firstInds), coneProfileZ(1)*zGrid(firstInds)]', coneVertices']'; 
    
    if ~metaData.createCylinder
        endInds = find(rGrid < coneProfileR(end));
    else
        endInds = find(rGrid < profileMax);
    end
    
    coneVertices = [ [xGrid(endInds), yGrid(endInds), (coneProfileZ(end)+1)*zGrid(endInds)]', coneVertices']'; 
    
    forcedInds = find(coneVertices(:,3) == (coneProfileZ(end)+1));
    
    coneVertices(:,1) = coneVertices(:,1) + volumeSize(1)/2;
    coneVertices(:,2) = coneVertices(:,2) + volumeSize(2)/2;
    coneVertices(:,3) = coneVertices(:,3) + metaData.tipOffset*metaData.voxSize/metaData.voxSize3D;

    % Switch to alpha shape to prevent degenerate triangles
    tempShape = alphaShape(coneVertices, alphaForCone);
    [tempFaces, tempVertices] = boundaryFacets(tempShape);
    
%     tempTriangulation = delaunayTriangulation(coneVertices);
%     [tempFaces, tempVertices] = freeBoundary(tempTriangulation); % this seems to flip vertices in z

    if ~isempty(forcedInds)
        % remove forced inds at base
        % not all of the original forced inds are used in the final trinagulation
        usedForcedInds = find(tempVertices(:,3) == coneVertices(forcedInds(1),3));

        facesToRemove = zeros(length(tempFaces),1, 'logical');
        for i = 1:length(usedForcedInds)
            % Flag for each column
            facesToRemove(tempFaces(:,1) == usedForcedInds(i)) = 1;
            facesToRemove(tempFaces(:,2) == usedForcedInds(i)) = 1;
            facesToRemove(tempFaces(:,3) == usedForcedInds(i)) = 1;
        end

        tempFaces(facesToRemove,:) = [];

        coneVertices(forcedInds,:) = [];
    end

    grinSurface.faces = tempFaces;
    grinSurface.vertices = tempVertices;

    coneProfileZ = coneProfileZ + metaData.tipOffset*metaData.voxSize/metaData.voxSize3D;
    
    % Get cone border volume
    if flipSurfaces
        grinSurface.vertices(:,3) = -(grinSurface.vertices(:,3)-volumeSize(3)/2)+volumeSize(3)/2;
        
        coneProfileZ = -(coneProfileZ-volumeSize(3)/2)+volumeSize(3)/2;
    end
    
    % Combine cone surfaces and really hope they are water tight...
    grinSurface.vertices = vertcat(grinSurface.vertices, cInCSurface.vertices);
    
    % Reverse order of CinC faces so they point in same polarity (at base so will point inward)
    if metaData.displayProfiles.CinC
        cInCSurface.faces = cInCSurface.faces(:,[3 2 1]);
    end

    % update the face list
    faceUpdate = cInCSurface.faces+max(max(grinSurface.faces));
    grinSurface.faces = vertcat(grinSurface.faces,faceUpdate);

    coneBorderVolume = polygon2voxel(grinSurface, volumeSize, 'none');
    

    % for inner intercone - is curved
    if any(metaData.exposedInterconeProfileToUseLeft - metaData.exposedInterconeProfileToUseRight)
        error('Right and left intercone profiles arent equal - rotated volume also wont work')
    end
    
    profileToUse = find(~isnan(metaData.exposedInterconeProfileToUseLeft));
    
    if ~metaData.createCylinder
        interconeProfileR = (metaData.exposedInterconeProfileToUseLeft(profileToUse)+metaData.coneProfileToUse(profileToUse))*metaData.voxSize/metaData.voxSize3D;
    
        % Shrink 1st in profile to cone
        % Tried adding a step before hand, but always caused problems somewhere
        interconeProfileR(1) = metaData.coneProfileToUse(profileToUse(1))*metaData.voxSize/metaData.voxSize3D;
    else
        interconeProfileR = (metaData.exposedInterconeProfileToUseLeft(profileToUse)+max(metaData.coneProfileToUse))*metaData.voxSize/metaData.voxSize3D;
    
        interconeProfileR(1) = max(metaData.coneProfileToUse)*metaData.voxSize/metaData.voxSize3D;
    end
    
    interconeProfileZ = metaData.coneXRef(profileToUse)*metaData.voxSize/metaData.voxSize3D;
    
    interconeVertices = zeros(length(interconeProfileR)*length(meshAngles),3);

    profileMax = max(interconeProfileR);  
    
    for i = 1:length(interconeProfileR)
        for j = 1:length(meshAngles)
            interconeVertices((i-1)*length(meshAngles)+j,:) = [interconeProfileR(i)*sin(meshAngles(j)), ...
                interconeProfileR(i)*cos(meshAngles(j)), interconeProfileZ(i)];
        end
    end
    
    % Get fake cap to pad top
    [xGrid, yGrid, zGrid] = meshgrid(-floor(profileMax):ceil(profileMax), -floor(profileMax):ceil(profileMax), (interconeProfileZ(1)-1));

    xGrid = xGrid(:); yGrid = yGrid(:); zGrid = zGrid(:);
    rGrid = sqrt(xGrid.^2 + yGrid.^2);

    startInds = find(rGrid < interconeProfileR(1));

    interconeVertices = [ [xGrid(startInds), yGrid(startInds), zGrid(startInds)]', interconeVertices']';

    forcedCapInds = find(interconeVertices(:,3) == interconeProfileZ(1)-1);
    % Need a lower multiple for alpha shape
    interconeVertices(forcedCapInds, 1:2) = interconeVertices(forcedCapInds, 1:2)*0.75;
    
    interconeVertices(:,1) = interconeVertices(:,1) + volumeSize(1)/2;
    interconeVertices(:,2) = interconeVertices(:,2) + volumeSize(2)/2;
    interconeVertices(:,3) = interconeVertices(:,3) + metaData.tipOffset*metaData.voxSize/metaData.voxSize3D;
    
    % Now get a ring to force base
    [xGrid, yGrid, zGrid] = meshgrid(1:volumeSize(1), 1:volumeSize(2), interconeProfileZ(end)+metaData.tipOffset*metaData.voxSize/metaData.voxSize3D);
    xGrid = xGrid(:); yGrid = yGrid(:); zGrid = zGrid(:);
    
    rGrid = sqrt((xGrid-volumeSize(1)/2).^2 + (yGrid-volumeSize(2)/2).^2);
    endInds = find(rGrid > interconeProfileR(end));
    
    interconeVertices = [interconeVertices', [xGrid(endInds), yGrid(endInds), zGrid(endInds)]']';
    
    % Make another layer at base to force 3D
    [xGrid, yGrid, zGrid] = meshgrid(1:volumeSize(1), 1:volumeSize(2), interconeProfileZ(end)+1+metaData.tipOffset*metaData.voxSize/metaData.voxSize3D);
    xGrid = xGrid(:); yGrid = yGrid(:); zGrid = zGrid(:);

    interconeVertices = [interconeVertices', [xGrid, yGrid, zGrid]']';
    forcedGridInds = find(interconeVertices(:,3) == interconeProfileZ(end)+1+metaData.tipOffset*metaData.voxSize/metaData.voxSize3D);

%     tempTriangulation = delaunayTriangulation(interconeVertices);
%     [tempFaces, tempVertices] = freeBoundary(tempTriangulation); 
    
    % as exposed cone is concave has to be done with alpha - but leads to some flattening around border
    tempShape = alphaShape(interconeVertices, alphaForIntercone);
    [tempFaces, tempVertices] = boundaryFacets(tempShape);

    % remove forced inds at cap and base
    usedForcedCapInds = find(tempVertices(:,3) == interconeVertices(forcedCapInds(1),3));
    usedForcedBaseInds = find(tempVertices(:,3) == interconeVertices(forcedGridInds(1),3));
    
    facesToRemove = zeros(length(tempFaces),1,'logical');
    for i = 1:length(usedForcedBaseInds)
        % Flag for each column
        facesToRemove(tempFaces(:,1) == usedForcedBaseInds(i)) = 1;
        facesToRemove(tempFaces(:,2) == usedForcedBaseInds(i)) = 1;
        facesToRemove(tempFaces(:,3) == usedForcedBaseInds(i)) = 1;
    end
    
    for i = 1:length(usedForcedCapInds)
        % Flag for each column
        facesToRemove(tempFaces(:,1) == usedForcedCapInds(i)) = 1;
        facesToRemove(tempFaces(:,2) == usedForcedCapInds(i)) = 1;
        facesToRemove(tempFaces(:,3) == usedForcedCapInds(i)) = 1;
    end

    tempFaces(facesToRemove,:) = [];

    interconeVertices([forcedCapInds' forcedGridInds']',:) = [];
    
    interconeSurface.faces = tempFaces;
    interconeSurface.vertices = tempVertices;

    interconeProfileZ = interconeProfileZ + metaData.tipOffset*metaData.voxSize/metaData.voxSize3D;
    
    if flipSurfaces
        interconeSurface.vertices(:,3) = -(interconeSurface.vertices(:,3)-volumeSize(3)/2)+volumeSize(3)/2;
        
        interconeProfileZ = -(interconeProfileZ-volumeSize(3)/2)+volumeSize(3)/2;
    end
    interconeBorderVolume = polygon2voxel(interconeSurface, volumeSize, 'none');

    % plots 
    if testPlotSurface
        figure;
    %     plot3(coneVertices(:,1), coneVertices(:,2), coneVertices(:,3), '.')
        hold on; axis equal
        trisurf(grinSurface.faces, grinSurface.vertices(:,1), grinSurface.vertices(:,2), grinSurface.vertices(:,3), ...
           'FaceColor','g','FaceAlpha',0.1);

    %     plot3(corneaVertices(:,1), corneaVertices(:,2), corneaVertices(:,3), '.')
        trisurf(corneaSurface.faces, corneaSurface.vertices(:,1), corneaSurface.vertices(:,2), corneaSurface.vertices(:,3), ...
           'FaceColor','m','FaceAlpha',0.8);

    %    plot3(epicorneaVertices(:,1), epicorneaVertices(:,2), epicorneaVertices(:,3), '.')
       trisurf(epicorneaSurface.faces, epicorneaSurface.vertices(:,1), epicorneaSurface.vertices(:,2), epicorneaSurface.vertices(:,3), ...
           'FaceColor','c','FaceAlpha',0.8);

    %    plot3(interconeBaseVertices(:,1), interconeBaseVertices(:,2), interconeBaseVertices(:,3), '.')
       trisurf(interconeBaseSurface.faces, interconeBaseSurface.vertices(:,1), interconeBaseSurface.vertices(:,2), interconeBaseSurface.vertices(:,3), ...
           'FaceColor','b','FaceAlpha',0.8);

    %    plot3(cInCVertices(:,1), cInCVertices(:,2), cInCVertices(:,3), '.')
    %    trisurf(cInCSurface.faces, cInCSurface.vertices(:,1), cInCSurface.vertices(:,2), cInCSurface.vertices(:,3), ...
    %        'FaceColor','r','FaceAlpha',0.8);

    %    plot3(interconeVertices(:,1), interconeVertices(:,2), interconeVertices(:,3), '.')
       trisurf(interconeSurface.faces, interconeSurface.vertices(:,1), interconeSurface.vertices(:,2), interconeSurface.vertices(:,3), ...
           'FaceColor','b','FaceAlpha',0.8);
    end
   
   figure;
   subplot(1,2,1)
   tempVolume = coneBorderVolume + corneaBorderVolume + epicorneaBorderVolume + interconeBaseBorderVolume + interconeBorderVolume;
   imshow(permute(tempVolume(round(volumeSize(1)/2),:,:), [2 3 1])');
   
   subplot(1,2,2)
   imshow(permute((lensRIVolume(round(volumeSize(1)/2),:,:)-1.35)/(1.54-1.35), [2 3 1])');

elseif useTestData
    if ~xor(createLunebergLens, createGradedFiber)
        error('Both or neither lens/fiber flags set. Check.')
    end
    
    if createLunebergLens
        volumeSize = ceil(radius*2*1.5/voxelSize*[1 1 1]);
        
        lensRIVolume = createluneberglens(radius/voxelSize, volumeSize);
    elseif createGradedFiber
        volumeSize = ceil([radius*2 radius*2 fiberLength]*1.5/voxelSize);
        
        lensRIVolume = creategradedfiber(radius/voxelSize, fiberLength/voxelSize, volumeSize, ...
            n0, alpha, voxelSize);
    end
    
    RIFlagVolume = ~isnan(lensRIVolume);

    lensRIVolume(isnan(lensRIVolume)) = exteriorRI;
    
    coneBorderVolume = imdilate(~RIFlagVolume, strel('sphere', 1)) & RIFlagVolume;
    
    % Note that meshes are used in test cases, refraction is based on an analytical solution
end

%% Do general set up

if useRealData
    coneTipZ = coneProfileZ(1)*voxelSize;
    
    % Connect GRIN vertexes to an exterior RI
    tempRIVolume = lensRIVolume;

    tempRIVolume(RIFlagVolume) = 0;

    tempRIVolume = imdilate(tempRIVolume, strel('Sphere',1));

    
    vertexIndexes = sub2ind(volumeSize, round(grinSurface.vertices(:,1)), round(grinSurface.vertices(:,2)), round(grinSurface.vertices(:,3)));

    vertexExteriorRI = tempRIVolume(vertexIndexes);

    % If any are missing loop until they are caught
    while any(vertexExteriorRI == 0)
        tempRIVolume2 = imdilate(tempRIVolume, strel('Sphere',1));

        vertexExteriorRI(vertexExteriorRI == 0) = tempRIVolume2(vertexIndexes(vertexExteriorRI == 0));
    end
    
    exposedConeFlag = vertexExteriorRI == metaData.innerValue;
    
    % Change RI of exterior at border of exposed cone
    if ~isempty(externalPigmentValue)
        vertexExteriorRI( grinSurface.vertices(:,3)*voxelSize < coneTipZ - exposedHeight/1000 & vertexExteriorRI == metaData.innerValue) = externalPigmentValue;
    end
    
    % Adjust vertices
    grinSurface.vertices = grinSurface.vertices*voxelSize;

    corneaSurface.vertices = corneaSurface.vertices*voxelSize;
     
    epicorneaSurface.vertices = epicorneaSurface.vertices*voxelSize;
      
    interconeBaseSurface.vertices = interconeBaseSurface.vertices*voxelSize;
       
    interconeSurface.vertices = interconeSurface.vertices*voxelSize;
    
    
    % Get normals for all
    grinNormals = meshFaceNormals(grinSurface.vertices, grinSurface.faces);
    
    corneaNormals = meshFaceNormals(corneaSurface.vertices, corneaSurface.faces);
    
    epicorneaNormals = meshFaceNormals(epicorneaSurface.vertices, epicorneaSurface.faces);
    
    interconeBaseNormals = meshFaceNormals(interconeBaseSurface.vertices, interconeBaseSurface.faces);
    
    interconeNormals = meshFaceNormals(interconeSurface.vertices, interconeSurface.faces);
    
    if testPlotSurface
        %%% Direction of open surfaces does not seem to matter, but cone normals all need to point in
        figure; hold on; axis equal
        trisurf(grinSurface.faces, grinSurface.vertices(:,1), grinSurface.vertices(:,2), grinSurface.vertices(:,3), ...
               'FaceColor','g','FaceAlpha',0.8);
        quiver3(grinSurface.vertices(grinSurface.faces(:,1),1), grinSurface.vertices(grinSurface.faces(:,1),2), grinSurface.vertices(grinSurface.faces(:,1),3), ...
            grinNormals(:,1)*5, grinNormals(:,2)*5, grinNormals(:,3)*5);

        trisurf(corneaSurface.faces, corneaSurface.vertices(:,1), corneaSurface.vertices(:,2), corneaSurface.vertices(:,3), ...
           'FaceColor','m','FaceAlpha',0.8);
       quiver3(corneaSurface.vertices(corneaSurface.faces(:,1),1), corneaSurface.vertices(corneaSurface.faces(:,1),2), corneaSurface.vertices(corneaSurface.faces(:,1),3), ...
            corneaNormals(:,1)*5, corneaNormals(:,2)*5, corneaNormals(:,3)*5);
           
        trisurf(epicorneaSurface.faces, epicorneaSurface.vertices(:,1), epicorneaSurface.vertices(:,2), epicorneaSurface.vertices(:,3), ...
           'FaceColor','c','FaceAlpha',0.8);
       quiver3(epicorneaSurface.vertices(epicorneaSurface.faces(:,1),1), epicorneaSurface.vertices(epicorneaSurface.faces(:,1),2), epicorneaSurface.vertices(epicorneaSurface.faces(:,1),3), ...
            epicorneaNormals(:,1)*5, epicorneaNormals(:,2)*5, epicorneaNormals(:,3)*5);
       
       trisurf(interconeBaseSurface.faces, interconeBaseSurface.vertices(:,1), interconeBaseSurface.vertices(:,2), interconeBaseSurface.vertices(:,3), ...
           'FaceColor','b','FaceAlpha',0.8);
       quiver3(interconeBaseSurface.vertices(interconeBaseSurface.faces(:,1),1), interconeBaseSurface.vertices(interconeBaseSurface.faces(:,1),2), interconeBaseSurface.vertices(interconeBaseSurface.faces(:,1),3), ...
            interconeBaseNormals(:,1)*5, interconeBaseNormals(:,2)*5, interconeBaseNormals(:,3)*5);
       
       trisurf(interconeSurface.faces, interconeSurface.vertices(:,1), interconeSurface.vertices(:,2), interconeSurface.vertices(:,3), ...
           'FaceColor','b','FaceAlpha',0.8);
       quiver3(interconeSurface.vertices(interconeSurface.faces(:,1),1), interconeSurface.vertices(interconeSurface.faces(:,1),2), interconeSurface.vertices(interconeSurface.faces(:,1),3), ...
            interconeNormals(:,1)*5, interconeNormals(:,2)*5, interconeNormals(:,3)*5);
    end
    
    % Correct normal direction - should point out...
    
    % It seems that some vertices can end up outside of the border volume..
    % Could try add border subscripts as well?
    coneBorderVolume = imdilate(coneBorderVolume, strel('sphere',dilateBorderRadius));
    
%     tempV = RIFlagVolume*2 + coneBorderVolume;
%     figure;   imshow(permute(tempV(round(volumeSize(1)/2),:,:)/3, [2 3 1])');
end

% Do plotting - for test volume
if useTestData  
    % Get surface voxels of RI volume
    surfaceInds = find(coneBorderVolume);

    [surfaceX, surfaceY, surfaceZ] = ind2sub(volumeSize, surfaceInds);

    bottomZ = min(surfaceZ);

    topZ = max(surfaceZ);
    
    figure;
    subplot(1,3,1);
    imshow((lensRIVolume(:,:,round(volumeSize(3)/2)))/(max(lensRIVolume(:))))
    
    subplot(1,3,2);
    imshow(permute((lensRIVolume(:,round(volumeSize(2)/2),:))/...
        (max(lensRIVolume(:))), [1 3 2]))

    subplot(1,3,3); hold on
    plot(lensRIVolume(:,round(volumeSize(2)/2),bottomZ),'b');
    
    centreLine = find(RIFlagVolume(:,round(volumeSize(2)/2),bottomZ));
    tempRadius = sqrt((centreLine - volumeSize(1)/2).^2 + ...
        (round(volumeSize(2)/2) - volumeSize(2)/2)^2 );
    
    % Nishidate eqn
    if createGradedFiber
        plot(centreLine, sqrt(2.5-(tempRadius*voxelSize/radius).^2), 'rx');

        plot(centreLine, sqrt(n0^2*(1-alpha^2*(tempRadius*voxelSize).^2)), 'g');
    end
end

% Search is limited to good in points
RIFlagInds = find(RIFlagVolume);

if useRealData
    plotInds = 1:50:size(grinSurface.vertices,1);
    
elseif useTestData
    if createLunebergLens
        plotInds = 1:length(surfaceX);

    elseif createGradedFiber
        
        plotInds = 1:50:length(surfaceX);
    end
end

% Get voxel coordinates
[tempX, tempY, tempZ] = meshgrid(1:volumeSize(1), 1:volumeSize(2), 1:volumeSize(3));

volInds = sub2ind(volumeSize, tempX(:), tempY(:), tempZ(:));

volCoords = ([tempX(:), tempY(:), tempZ(:)])*voxelSize;

zSteps = (1:volumeSize(3)*1.5)*voxelSize;

if useRealData
    
    xStartPoints = volumeSize(1)/2:xSpacing:volumeSize(1)*1.5; 
    xStartPointsOrig = xStartPoints;

    xStartPoints = [volumeSize(1)-fliplr(xStartPoints(2:end)) xStartPoints];

    if trace3D
        [xGrid, yGrid, zGrid] = meshgrid(xStartPoints, xStartPoints, 1);
        
        rayOrigins = [xGrid(:), yGrid(:), zGrid(:)]*voxelSize; 
        
    else
        rayOrigins = [xStartPoints(:), ones(numel(xStartPoints),1)*volumeSize(2)/2,  ...
            ones(numel(xStartPoints),1)]*voxelSize; 
    end
    
    if ~isempty(testOrigin)
        rayOrigins = rayOrigins(testOrigin,:);
    end
else
    %%% Does in 2D, could add trace3D as well (need to update solver)
    
    if ~plotReferenceFiber | createLunebergLens
        xStartPoints = volumeSize(1)/2:xSpacing:volumeSize(1); 
        
        xStartPoints = [(xStartPoints(1:end-1)-volumeSize(1)/2+1) xStartPoints];
    else
        xStartPoints = volumeSize(1)/2+radius/voxelSize*0.5;
    end

    rayOrigins = [xStartPoints(:), ones(numel(xStartPoints),1)*volumeSize(2)/2,  ...
        ones(numel(xStartPoints),1)]*voxelSize; 
end

raysOnXZPlane = rayOrigins(:,2)/voxelSize == volumeSize(2)/2;

raysOnYZPlane = rayOrigins(:,1)/voxelSize == volumeSize(1)/2;

nOrigins = size(rayOrigins,1);

rayCols = lines(nOrigins);

%% Run ray tracer

volumeCenter = volumeSize/2*voxelSize;

% Get border test values

if useRealData
    % Calculating inpolyhedron is very slow, so test on volume first
    intersectionFn = @(x1, x0)(surfaceIntersectFunction(RIFlagVolume, coneBorderVolume, x1, x0, voxelSize, grinSurface));
       
    %%% Correctly calculated as ray distance
    lambdaFn = @(x1, x0)surfaceLambdaFunction(x1, x0, grinSurface);

elseif useTestData
    if createLunebergLens
        % Get sphere intersection
       intersectionFn = @(x1, x0)(sqrt(sum((x1 - volumeCenter).^2)) <= radius);

       % Calculate as radius of ratios.
            %%% As below, want dist to border along ray-axis instead of z dist
       warning('Lambada calculation needs correction')
       
       lambdaFn = @(x1, x0)((radius-sqrt(sum((x0 - volumeCenter).^2)))/...
           (sqrt(sum((x1 - volumeCenter).^2))-sqrt(sum((x0 - volumeCenter).^2))));

    elseif createGradedFiber
       % Just checks X position, not Y
       intersectionFn = @(x1, x0)((x1(3) <= volumeCenter(3) + fiberLength/2 & ...
            x1(3) >= volumeCenter(3) - fiberLength/2  & ...
            x1(1) >= volumeCenter(1) - radius & x1(1) <= volumeCenter(1) + radius));
        
        %%% Lambda is not entirely correct as we want to border along ray-axis instead of z dist. 
            %%% But probably negligible for small distances
        warning('Lambada calculation needs correction')
            
       lambdaFn = @(x1, x0)(((volumeCenter(3) + fiberLength/2) - x0(3))/(x1(3)-x0(3)));
    end
end

rayPathCells = cell(length(incidenceAngle), 1);
finalIntersectCells = cell(length(incidenceAngle), 1); 
finalRayCells = cell(length(incidenceAngle), 1);
finalRayTRefractCells = cell(length(incidenceAngle), 1);
firstIntersectCells = cell(length(incidenceAngle), 1); 
TIRFlagCells = cell(length(incidenceAngle), 1);
rayReverseCells = cell(length(incidenceAngle), 1);

timePerAngle = zeros(length(incidenceAngle),1);

for aAngle = 1:length(incidenceAngle);
    
    tic
    
    startRayT = [sin(incidenceAngle(aAngle)/180*pi), 0, cos(incidenceAngle(aAngle)/180*pi)];    

    rayPathArray = zeros(3, length(zSteps), nOrigins)*NaN;

    rayPathLengthArray = zeros(length(zSteps), nOrigins)*NaN;

    if testPlot
        figure; hold on; axis equal; set(gca,'Clipping','off')
        view(0, 0)
        if useRealData
            plot3(grinSurface.vertices(plotInds,1), grinSurface.vertices(plotInds,2),...
                grinSurface.vertices(plotInds,3), '.')
        elseif useTestData
            plot3(surfaceX(plotInds)*voxelSize, surfaceY(plotInds)*voxelSize,...
                surfaceZ(plotInds)*voxelSize, '.')
        end
    end

    % Turn of warning on singular matrix, otherwise backslash will slow down a lot
    warning('off', 'MATLAB:nearlySingularMatrix')

    minDeltaS = ones(nOrigins, 1);

    firstIntersect = zeros(nOrigins, 3)*NaN;
    finalIntersect = zeros(nOrigins, 3)*NaN;

    finalPathLength = zeros(nOrigins, 1);
    finalRayT = zeros(nOrigins, 3)*NaN;

    finalRayTRefract = zeros(nOrigins, 3)*NaN;
    finalRay = zeros(nOrigins, 3)*NaN;
    
    TIRFlag = zeros(nOrigins, 2);
    rayReversed = zeros(nOrigins, 1);
    
    if createGradedFiber
       % Period is for parabolic profile, gaussian period is different and independent of entry.
       %%% This is approx as specific value varies depending on entry RI
       halfPeriodInFiber = 2*pi*n0/beta*0.5;   

       numPeriods = floor(fiberLength/halfPeriodInFiber);

       periodSpotsArray = zeros(3, numPeriods, nOrigins);
    end

    for iOrigin = 1:nOrigins;
        
        tempTime = toc;
        [aAngle iOrigin round(tempTime/60)]
        
        %%% Could set up to do for series of ray angles
        rayT = startRayT; %3x1 : ray will move from bottom to top

        % In phyisical coordiantes
        %%% Was in voxel coordinates
        rayX = rayOrigins(iOrigin,:);

        inGraded = 0;

        %%% Was ray x in phyiscial coordinates
        voxelX = round(rayOrigins(iOrigin,:)/voxelSize);

        lastNearbyX = [rayX(1), rayX(2), -1];

        go = 1;

        deltaS = initialDeltaS;

        currentStep = 1;

        rayPathArray(:, currentStep, iOrigin) = rayX;

        intersectResult = 0;

        pathLength = 0; % in grin area

        propogateFinalRay = 0;

        numberOfExits = 0;
        
        hitsCornea = 0;

        previousIntersectValue = 0; minIntersect = 0;
        
        currentPeriod = 1; % only used for graded fiber

        loopSteps = 0;
        
        lastTIRStep = 0;
        
        rayDiff = [1 1 1];
        
        while go

            % If outside, extend ray to grin area
            if ~inGraded
                
                x0 = rayX;
                t0 = rayT;

                % entering, need to find intersect to graded volume
                if useRealData
                    previousIntersectValue = minIntersect;
                    
                    % X and T combined on input, T needs to point back
                    lineDef = [rayX rayT];

                    % Check intersect against all surfaces
                    [intersectPointsGrin, intersectDistanceGrin, intersectFacesGrin] = intersectLineMesh3d(lineDef, grinSurface.vertices, grinSurface.faces);

                    [intersectPointsCornea, intersectDistanceCornea, intersectFacesCornea] = intersectLineMesh3d(lineDef, corneaSurface.vertices, corneaSurface.faces);

                    [intersectPointsEpicornea, intersectDistanceEpicornea, intersectFacesEpicornea] = intersectLineMesh3d(lineDef, epicorneaSurface.vertices, epicorneaSurface.faces);

                    [intersectPointsCinC, intersectDistanceCinC, intersectFacesCinC] = intersectLineMesh3d(lineDef, interconeBaseSurface.vertices, interconeBaseSurface.faces);

                    [intersectPointsIntercone, intersectDistanceIntercone, intersectFacesIntercone] = intersectLineMesh3d(lineDef, interconeSurface.vertices, interconeSurface.faces);

                    % put miniminum intersection distances into array - could be better way to do this
                    if ~isempty(intersectDistanceGrin)
                        % Remove points with dist less than ~zero - should modify function to do this...
                        % Or outside of bounds
                        indsOut = find(intersectDistanceGrin <= 1e-6 | ...
                            intersectPointsGrin(:,1) < voxelSize | intersectPointsGrin(:,1) > volumeSize(1)*voxelSize | ...
                            intersectPointsGrin(:,2) < voxelSize | intersectPointsGrin(:,2) > volumeSize(2)*voxelSize | ...
                            intersectPointsGrin(:,3) < voxelSize | intersectPointsGrin(:,3) > volumeSize(3)*voxelSize);

                        intersectPointsGrin(indsOut, :) = [];
                        intersectFacesGrin(indsOut) = [];
                        intersectDistanceGrin(indsOut) = [];

                    end
                    if ~isempty(intersectDistanceGrin); minDistArray = min(intersectDistanceGrin); else; minDistArray = Inf; end

                    if ~isempty(intersectDistanceCornea);  
                        intersectPointsCornea(intersectDistanceCornea <= 1e-6, :) = [];
                        intersectFacesCornea(intersectDistanceCornea <= 1e-6) = [];
                        intersectDistanceCornea(intersectDistanceCornea <= 1e-6) = [];
                    end
                    if ~isempty(intersectDistanceCornea); minDistArray = [minDistArray min(intersectDistanceCornea)]; else; minDistArray = [minDistArray Inf]; end

                    if ~isempty(intersectDistanceEpicornea);  
                        intersectPointsEpicornea(intersectDistanceEpicornea <= 1e-6, :) = [];
                        intersectFacesEpicornea(intersectDistanceEpicornea <= 1e-6) = [];
                        intersectDistanceEpicornea(intersectDistanceEpicornea <= 1e-6) = [];
                    end
                    if ~isempty(intersectDistanceEpicornea); minDistArray = [minDistArray min(intersectDistanceEpicornea)]; else; minDistArray = [minDistArray Inf]; end

                    if ~isempty(intersectDistanceCinC);  
                        intersectPointsCinC(intersectDistanceCinC <= 1e-6, :) = [];
                        intersectFacesCinC(intersectDistanceCinC <= 1e-6) = [];
                        intersectDistanceCinC(intersectDistanceCinC <= 1e-6) = [];
                    end
                    if ~isempty(intersectDistanceCinC);   minDistArray = [minDistArray min(intersectDistanceCinC)]; else; minDistArray = [minDistArray Inf]; end

                    if ~isempty(intersectDistanceIntercone);  
                        intersectPointsIntercone(intersectDistanceIntercone <= 1e-6, :) = [];
                        intersectFacesIntercone(intersectDistanceIntercone <= 1e-6) = [];
                        intersectDistanceIntercone(intersectDistanceIntercone <= 1e-6) = [];
                    end
                    if ~isempty(intersectDistanceIntercone); minDistArray = [minDistArray min(intersectDistanceIntercone)]; else; minDistArray = [minDistArray Inf]; end

                    if any(~isinf(minDistArray))
                        [~, minIntersect] = min(minDistArray);

                        switch minIntersect

                            case 1 % Cone is first

                                inds2Use = find(intersectDistanceGrin == min(intersectDistanceGrin));

                                rayX = mean(intersectPointsGrin(inds2Use,:),1);

                                faceIndices = intersectFacesGrin(inds2Use);

                                intersectResult = 1;

                            case 2 % cornea is first
                                inds2Use = find(intersectDistanceCornea == min(intersectDistanceCornea));

                                rayX = mean(intersectPointsCornea(inds2Use,:),1);

                                faceIndices = intersectFacesCornea(inds2Use);
                                
                                hitsCornea = 1;

                            case 3 % epicornea is first
                                inds2Use = find(intersectDistanceEpicornea == min(intersectDistanceEpicornea));

                                rayX = mean(intersectPointsEpicornea(inds2Use,:),1);

                                faceIndices = intersectFacesEpicornea(inds2Use);    

                            case 4 % CinC is first
                                inds2Use = find(intersectDistanceCinC == min(intersectDistanceCinC));

                                rayX = mean(intersectPointsCinC(inds2Use,:),1);

                                faceIndices = intersectFacesCinC(inds2Use);

                            case 5 % intercone is first  
                                inds2Use = find(intersectDistanceIntercone == min(intersectDistanceIntercone));

                                rayX = mean(intersectPointsIntercone(inds2Use,:),1);

                                faceIndices = intersectFacesIntercone(inds2Use);

                                if limitToConeBase & numberOfExits == 0 
                                    % Block ray as it can't have entered cone yet 
                                    go = 0;
                                end
                        end

                        %%% Seems like there could be a better place for this...
                        propogateFinalRay = 1*go;
                    else
                        go = 0;
                    end

                elseif useTestData
                    origin = rayX - volumeSize/2*voxelSize;
                    minIntersect = -1;

                    if sqrt(origin(1)^2 + origin(2)^2) < radius
                        if createLunebergLens                 
                            % Get sphere intersection
                            centre = 0;
                            originalDistance = centre - origin;

                            t_ca = dot(originalDistance,rayT);
                            d = sqrt(dot(originalDistance,originalDistance)-t_ca^2);
                            t_hc = sqrt(radius^2-d^2);
                            deltaTemp = t_ca - t_hc;

                        elseif createGradedFiber
                            % Get z distance past bottom of gradient section
                            deltaTemp = ((volumeCenter(3) - fiberLength/2) - rayX(3))/rayT(3);
                        end

                        % This is mostly to stop ray jumping back at the end
                        if deltaTemp > 0
                            % Steps ray back to intersection
                            rayX = rayX + rayT*deltaTemp; %3x1 coordinates

                            intersectResult = 1;

                            propogateFinalRay = 1;
                        else
                           go = 0; 
                        end

                    else
                        go = 0;
                    end
                end

                if go
                    voxelX = round(rayX/voxelSize); 
                    
                    % extend ray path up to intersection - just for recording/plotting
                    misingSteps = find(zSteps > x0(3) & zSteps < rayX(3));

                    if ~isempty(misingSteps)
                        tempT = rayT/rayT(3);
                        for i = misingSteps
                            rayPathArray(:, i, iOrigin) = x0 + tempT.*(zSteps(i) - x0(3));
                        end

                        if testPlot
                            plot3(rayPathArray(1, misingSteps, iOrigin), rayPathArray(2, misingSteps, iOrigin), ...
                                rayPathArray(3, misingSteps, iOrigin), '.', 'color', rayCols(iOrigin,:));
                        end

                        currentStep = misingSteps(end);
                    end

                    if useTestData | (useRealData & minIntersect == 1) | (useRealData & fixedConeRI & minIntersect == 1 & previousIntersectValue ~= 1)
                       if useRealData
                            enteringFromExposed = any(exposedConeFlag(grinSurface.faces(faceIndices,:)));
                       else
                          enteringFromExposed = 0; 
                       end
                    end    
                        
                    % Entering into graded region
                    if useTestData | (useRealData & minIntersect == 1) & ~fixedConeRI
                        
                        % Check conditions are ok
                       if ~(blockMultipleExits & numberOfExits > 0) & ~(blockExposedRentry & enteringFromExposed) & (hitsCornea | useTestData)
                            % Get interior RI at entry point
                            [~, rOut] = numerical_dT_dt(rayX, volCoords(RIFlagInds,:), RIFlagInds, lensRIVolume, NaN);

                            if interfaceRefraction
                                % Get surface normal and exterior RI at entry point
                                if useRealData
                                    rIn = vertexExteriorRI(grinSurface.faces(faceIndices,:));
                                    rIn = mean(rIn(:));    

                                    surfaceNormal = mean(grinNormals(faceIndices,:),1);
                                    
                                elseif useTestData
                                    if createLunebergLens
                                        surfaceNormal = rayX - volumeSize/2*voxelSize;

                                        % Should be very close to 1, but allowing this removes notch in error profile
                                            % However, it does increase error spot slightly
                                        %rIn = 1;

                                    elseif createGradedFiber
                                        % Just points backwards
                                        surfaceNormal = [0 0 -1];

                                        halfPeriodInFiber = 2*pi*rOut/beta*0.5; 
                                    end

                                    rIn = exteriorRI;
                                end

                                surfaceNormal = surfaceNormal/norm(surfaceNormal);

                                % Test normal points against ray
                                if acos(dot(rayT, surfaceNormal)/(norm(surfaceNormal)*norm(rayT))) < pi/2
                                    surfaceNormal = -surfaceNormal;
                                end

                                % Calculate entry refraction
                                nRatio = rIn/rOut;
                                cosI = -dot(surfaceNormal, rayT);
                                sinT2 = nRatio^2*(1-cosI^2);
                                cosT = sqrt(1-sinT2);

                                if sinT2 < 1
                                    % Assuming all refracted, none reflected
                                    rayT = nRatio*rayT + (nRatio*cosI-cosT)*surfaceNormal;

                                    inGraded = 1;
                                    
                                    % now multiplied by initial RI
                                    rayT = rayT/norm(rayT)*rOut;
                                    
                                    if all(isnan(firstIntersect))
                                        firstIntersect(iOrigin,:) = rayX;
                                    end

                                    if numberOfExits > 0 
                                        if testPlot
                                           plot3(rayX(1), rayX(2), rayX(3), 'ro', 'markersize', 8)

                                           plot3(grinSurface.vertices(grinSurface.faces(faceIndices,:),1), grinSurface.vertices(grinSurface.faces(faceIndices,:),2), ...
                                               grinSurface.vertices(grinSurface.faces(faceIndices,:),3), 'rx')

                                           trisurf(grinSurface.faces(faceIndices,:), grinSurface.vertices(:,1), grinSurface.vertices(:,2), grinSurface.vertices(:,3), ...
                                                'FaceColor','r');
                                        end

                                       % Changing epsilon or deltaS may help 
                                       warning('ray is reentering') 
                                    end

                                    if testPlot; plot3(rayX(1), rayX(2), rayX(3), 'go'); end
                                else
                                    % TIR
                                    rayT = rayT-2*dot(surfaceNormal,rayT)*surfaceNormal;

                                    rayT = rayT/norm(rayT);
                                    
                                    inGraded = 0;
                                    
                                    intersectResult = 0;
                                    
                                    warning('TIR away from GRIN region')

                                    if testPlot; plot3(rayX(1), rayX(2), rayX(3), 'g*'); end
                                end

                            else
                                inGraded = 1; 

                                rayT = rayT*rOut;
                                
                                if testPlot; plot3(rayX(1), rayX(2), rayX(3), 'gd'); end
                            end

                       else
                           propogateFinalRay = 0;

                           go = 0;
                       end

                    elseif (useRealData & minIntersect > 1) | (useRealData & fixedConeRI)

                        if interfaceRefraction
                            % do refraction between layers

                            %%% currently sets RI assuming ray comes from base
                                % Could test direction to normal and invert if negative?
                            switch minIntersect
                                case 1 % cone
                                    % coming fromt outside
                                    if previousIntersectValue ~= 1
                                        rIn = vertexExteriorRI(grinSurface.faces(faceIndices,:));
                                        rIn = mean(rIn(:)); 
                                        
                                        rOut = metaData.coneValue;
                                        
                                        if ~(~(blockMultipleExits & numberOfExits > 0) & ~(blockExposedRentry & enteringFromExposed) & hitsCornea)
                                            propogateFinalRay = 0;

                                            go = 0;
                                        end
                                    else
                                       % coming from inside 
                                       rIn = metaData.coneValue;
                                       
                                       rOut = vertexExteriorRI(grinSurface.faces(faceIndices,:));
                                       rOut = mean(rOut(:)); 
                                    end
                                    
                                    surfaceNormal = mean(grinNormals(faceIndices,:),1);
                                    
                                case 2 % cornea
                                    rIn = metaData.outerValue;

                                    rOut = metaData.epicorneaValue;

                                    surfaceNormal = mean(corneaNormals(faceIndices,:),1);

                                case 3 % epicornea 
                                    rIn = metaData.epicorneaValue;

                                    rOut = metaData.outerCorneaValue;

                                    surfaceNormal = mean(epicorneaNormals(faceIndices,:),1);

                                case 4 % CinC (into intercone, not cone)
                                    rIn = metaData.outerCorneaValue;

                                    rOut = interConeValueUsed;

                                    surfaceNormal = mean(interconeBaseNormals(faceIndices,:),1);

                                    if limitToConeBase

                                        propogateFinalRay = 0;

                                        go = 0;
                                    end

                                case 5 % intercone 
                                    rIn = interConeValueUsed;

                                    if isempty(externalPigmentValue)
                                        rOut = metaData.innerValue;
                                    else
                                        % may not be neccersary, I exepect all intercone is below this
                                        if rayX(3) < coneTipZ - exposedHeight/1000
                                            rOut = externalPigmentValue; 
                                        else
                                            rOut = metaData.innerValue;
                                        end
                                    end
                                    surfaceNormal = mean(interconeNormals(faceIndices,:),1);
                            end

                            surfaceNormal = surfaceNormal/norm(surfaceNormal);
                            
                            % Test normal points against ray
                            if acos(dot(rayT, surfaceNormal)/(norm(surfaceNormal)*norm(rayT))) < pi/2
                                surfaceNormal = -surfaceNormal;
                            end

                            % Calculate entry refraction
                            nRatio = rIn/rOut;
                            cosI = -dot(surfaceNormal, rayT);
                            sinT2 = nRatio^2*(1-cosI^2);
                            cosT = sqrt(1-sinT2);

                            if sinT2 < 1
                                % Ray is leaving from inside of grin
                                if minIntersect == 1 & previousIntersectValue == 1 
                                   finalIntersect(iOrigin,:) = rayX;

                                   finalRayT(iOrigin,:) = rayT;
                                   
                                   numberOfExits = numberOfExits + 1;
                                    
                                % Ray is entering grin    
                                elseif minIntersect == 1 & previousIntersectValue ~= 1
                                    if all(isnan(firstIntersect(iOrigin,:)))
                                        firstIntersect(iOrigin,:) = rayX;
                                    end

                                    if numberOfExits > 0 
                                        if testPlot
                                           plot3(rayX(1), rayX(2), rayX(3), 'ro', 'markersize', 8)

                                           plot3(grinSurface.vertices(grinSurface.faces(faceIndices,:),1), grinSurface.vertices(grinSurface.faces(faceIndices,:),2), ...
                                               grinSurface.vertices(grinSurface.faces(faceIndices,:),3), 'rx')

                                           trisurf(grinSurface.faces(faceIndices,:), grinSurface.vertices(:,1), grinSurface.vertices(:,2), grinSurface.vertices(:,3), ...
                                                'FaceColor','r');
                                        end

                                       % Changing epsilon or deltaS may help 
                                       warning('ray is reentering') 
                                    end
                                end
                                
                                % Assuming all refracted, none reflected
                                rayT = nRatio*rayT + (nRatio*cosI-cosT)*surfaceNormal;
                                    
                                if minIntersect == 1 & previousIntersectValue == 1 
                                    finalRayTRefract(iOrigin,:) = rayT;
                                end
                                
                                if testPlot; plot3(rayX(1), rayX(2), rayX(3), 'ko'); end
                            else
                                % TIR
                                rayT = rayT-2*dot(surfaceNormal,rayT)*surfaceNormal;

                                if minIntersect == 5
                                    % Check for TIR from intercone 
                                    TIRFlag(iOrigin,2) = TIRFlag(iOrigin,2) + 1;
                                elseif minIntersect == 1 & previousIntersectValue == 1 
                                    TIRFlag(iOrigin,1) = TIRFlag(iOrigin,1) + 1;
                                end
                                
                                if testPlot; plot3(rayX(1), rayX(2), rayX(3), 'k*'); end
                            end

                            rayT = rayT/norm(rayT);
                        else
                            if testPlot; plot3(rayX(1), rayX(2), rayX(3), 'kd'); end
                        end
                    end
                end
            end

            if inGraded
                %check if x has moved to include new voxels
                if (norm(rayX-lastNearbyX) > voxelSize) & inGraded
                    %if so, reselect nearby points
                    lastNearbyX = rayX;

                    % Get array of voxels around current
                    [tempX, tempY, tempZ] = meshgrid((voxelX(1)-4:voxelX(1)+4)', (voxelX(2)-4:voxelX(2)+4)', ...
                        (voxelX(3)-4:voxelX(3)+4)');

                    testCoords= [tempX(:), tempY(:), tempZ(:)];

                    % remove outside of bounds
                    tempInd = find(testCoords(:,1) < 1 | testCoords(:,2) < 1 | testCoords(:,3) < 1 | ...
                            testCoords(:,1) > volumeSize(1) | testCoords(:,2) > volumeSize(2) | testCoords(:,3) > volumeSize(3));

                    testCoords(tempInd, :) = [];

                    testInds = sub2ind(volumeSize, testCoords(:,1), testCoords(:,2), testCoords(:,3));

                    % check for intersect with good volume
                    [testInds, tempInds] = intersect(testInds, RIFlagInds);

                    testCoords = testCoords(tempInds,:)*voxelSize;
                    
                    % only nearest 128 used, but seems to need some extra
                    [~, nearInds] = sort(sqrt((testCoords(:,1)-rayX(1)).^2+(testCoords(:,2)-rayX(2)).^2+...
                        (testCoords(:,3)-rayX(3)).^2));

                    if length(nearInds) > 200
                        nearInds = nearInds(1:200);
                    end

                    testInds = testInds(nearInds);

                    testCoords = testCoords(nearInds,:);

                    % If all voxels equal in patch can flag to just use that
                    if all(lensRIVolume(testInds) == lensRIVolume(testInds(1)))
                        RIToUse = lensRIVolume(testInds(1));
                    else
                        RIToUse = NaN;
                    end
                        
                    %%% If doing direction updating, calculate RI for current direction
                    
%                     [rayX rayT deltaS*10^6 length(testInds)]
                end
                
                x0 = rayX;
                t0 = rayT;

                % Calc is fixed to isotropic RI
                [x, t, deltaS] = ray_interpolation(interpType, 'iso', x0', t0', deltaS, testCoords, ...
                    testInds, lensRIVolume, tolerance, RIToUse);

                rayX = x'; rayT = t';

                pathLength = pathLength + sqrt((x0(1) - rayX(1))^2 + (x0(2) - rayX(2))^2 + (x0(3) - rayX(3))^2);

                if deltaS < minDeltaS(iOrigin)
                    minDeltaS(iOrigin) = deltaS;
                end

                if createGradedFiber
                    % test if focal point is passed and capture if so
                    if rayX(3)-(volumeSize(3)/2*voxelSize - fiberLength/2) > halfPeriodInFiber*(currentPeriod-0.5)

                        % Step back - this is just a linear interpolation up to point
                        deltaTemp = (halfPeriodInFiber*(currentPeriod-0.5) - (rayX(3)-(volumeSize(3)/2*voxelSize - fiberLength/2)))/rayT(3);

                        periodSpotsArray(:,currentPeriod,iOrigin) = rayX + deltaTemp*rayT; 

                        if testPlot;
                            plot3(periodSpotsArray(1,currentPeriod,iOrigin), periodSpotsArray(2,currentPeriod,iOrigin), periodSpotsArray(3,currentPeriod,iOrigin), 'xr');
                        end
                        
                        currentPeriod = currentPeriod + 1;
                    end
                end

                intersectResult = intersectionFn(rayX, x0);

                % Test if leaving a GRIN area
                if ~intersectResult

                    % remove last pathlength step
                    pathLength = pathLength - sqrt((x0(1) - rayX(1))^2 + (x0(2) - rayX(2))^2 + (x0(3) - rayX(3))^2);

                    % Resort voxels
                    % Get array of voxels around current - using old voxelX
                    [tempX, tempY, tempZ] = meshgrid((voxelX(1)-5:voxelX(1)+5)', (voxelX(2)-5:voxelX(2)+5)', ...
                        (voxelX(3)-5:voxelX(3)+5)');

                    testCoords= [tempX(:), tempY(:), tempZ(:)];

                    % remove outside of bounds
                    tempInd = find(testCoords(:,1) < 1 | testCoords(:,2) < 1 | testCoords(:,3) < 1 | ...
                            testCoords(:,1) > volumeSize(1) | testCoords(:,2) > volumeSize(2) | testCoords(:,3) > volumeSize(3));

                    testCoords(tempInd, :) = [];

                    testInds = sub2ind(volumeSize, testCoords(:,1), testCoords(:,2), testCoords(:,3));
                    
                    % check for intersect with good volume
                    [testInds, tempInds] = intersect(testInds, RIFlagInds);

                    testCoords = testCoords(tempInds,:)*voxelSize;

                    % Sort forwards for closest 128
                    [~, nearInds] = sort(sqrt((testCoords(:,1)-x0(1)).^2+(testCoords(:,2)-x0(2)).^2+...
                        (testCoords(:,3)-x0(3)).^2));

                    if length(nearInds) > 400
                        nearInds = nearInds(1:400);
                    end

                    testInds = testInds(nearInds);
                    testCoords = testCoords(nearInds,:);

                    if all(lensRIVolume(testInds) == lensRIVolume(testInds(1)))
                        RIToUse = lensRIVolume(testInds(1));
                    else
                        RIToUse = NaN;
                    end
                    
                    lambda0 = lambdaFn(rayX, x0);

                    if iterativeFinal && ~isnan(lambda0)
                        % Test intersection, leaving could just be a rounding point           
                        deltaS0_final = deltaS;

                        its = 0;

                        % Top test could just be written as lambda0 < 1+1e-9
                        while deltaS0_final*norm(t0) > epsilon && abs(lambda0-1.0) > 1e-9 
                            %%% using 10^-3 as test for zero in tests above and below

                            SF = 1.0; %A from paper - saftey factor

                            if lambda0 < 1e-9
                               rayT = t0; 
                               rayX = x0;
                               break
                            end

                            its = its + 1;

                            %this loop steps back from overshoot
                            while ~intersectionFn(rayX, x0)
                                if SF > 0.11
                                   SF = SF - 0.1; 
                                else
                                   SF = 0.1*SF;
                                end

                                deltaS_final = SF*lambda0*deltaS0_final;

                                % step ray backward with RK5 with fixed step
                                [x, t] = ray_interpolation('5RKN', 'iso', x0', t0', deltaS_final, testCoords, ...
                                    testInds, lensRIVolume, tolerance, RIToUse);

                                rayX = x'; rayT = t';
                            end

                            deltaS_final = (1-SF)*lambda0*deltaS0_final;

                            % This loop steps forward to find new contact position just before overshoot.
                            while intersectionFn(rayX, x0)
                                x0 = rayX; %think this is the case, we are stepping forward
                                t0 = rayT;

                                [x, t] = ray_interpolation('5RKN', 'iso', x0', t0', deltaS_final, testCoords, ...
                                    testInds, lensRIVolume, tolerance, RIToUse);

                                rayX = x'; rayT = t';

                                pathLength = pathLength + sqrt((x0(1) - rayX(1))^2 + (x0(2) - rayX(2))^2 + (x0(3) - rayX(3))^2);
                            end

                            %update parameters for next loop
                            lambda0 = lambdaFn(rayX, x0);

                            if isnan(lambda0)
                               break;
                            end

                            deltaS0_final = deltaS_final;
                        end

                    elseif ~iterativeFinal && ~isnan(lambda0)
                        [x, t] = ray_interpolation(interpType, 'iso', x0', t0', deltaS*lambda0, testCoords_forwards, ...
                            testInds_forwards, lensRIVolume, tolerance, RIToUse);

                       rayX = x';
                       rayT = t';

                       pathLength = pathLength + sqrt((x0(1) - rayX(1))^2 + (x0(2) - rayX(2))^2 + (x0(3) - rayX(3))^2);
                    end

                    % Get surface normal and RI at exit point
                    [~, rIn] = numerical_dT_dt(rayX, volCoords(RIFlagInds,:), RIFlagInds, lensRIVolume, NaN); 

                    % In test case normalization seems to be very similar to dividing by RI at exit
                    rayT = rayT/rIn;

                    % Do refraction at border
                    if interfaceRefraction && ~isnan(lambda0)

                        if useRealData
                            lineDef = [rayX rayT];

                            [intersectPoints, intersectDistance, intersectFaces] = intersectLineMesh3d(lineDef, grinSurface.vertices, grinSurface.faces);

                            inds2Use = find(abs(intersectDistance) == min(abs(intersectDistance)));

                            %%% Considered adjusting rayX to be intersect value but decided against 
                             % changes in X are all dealt with via numerical integration

                            faceIndices = intersectFaces(inds2Use);

                            surfaceNormal = mean(grinNormals(faceIndices,:),1);

                            surfaceNormal = surfaceNormal/norm(surfaceNormal);

                            rOut = vertexExteriorRI(grinSurface.faces(faceIndices,:));
                            rOut = mean(rOut(:));  

                        elseif useTestData
                            if createLunebergLens
                                surfaceNormal = -(rayX - volumeSize/2*voxelSize);

                                surfaceNormal = surfaceNormal/norm(surfaceNormal);

                                % Usually not exactly 1, this allows refraction but actually decreases error on exit trajectory
            %                   rOut = 1;

                            elseif createGradedFiber
                                % Just points backwards
                                surfaceNormal = [0 0 -1];   
                            end

                            rOut = exteriorRI;
                        end

                        if acos(dot(rayT, surfaceNormal)/(norm(surfaceNormal)*norm(rayT))) < pi/2
                            surfaceNormal = -surfaceNormal;
                        end

                        % Calculate exit refraction
                        nRatio = rIn/rOut;
                        cosI = -dot(surfaceNormal, rayT);
                        sinT2 = nRatio^2*(1-cosI^2);
                        cosT = sqrt(1-sinT2);

                        if sinT2 < 1
                           % leaving so store values                     
                           finalIntersect(iOrigin,:) = rayX;

                           finalPathLength(iOrigin,:) = pathLength;

                           finalRayT(iOrigin,:) = rayT*rIn; % was divided previously

                           % Assuming all refracted, none reflected
                           rayT = nRatio*rayT + (nRatio*cosI-cosT)*surfaceNormal;

                           if testPlot; plot3(rayX(1), rayX(2), rayX(3), 'mo'); end 

                           inGraded = 0;

                           numberOfExits = numberOfExits + 1;

                           rayT = rayT/norm(rayT);

                           finalRayTRefract(iOrigin,:) = rayT;
                        else
                            % Step back a bit so that rayX will be just inside mesh instead of just outside
                            if lambda0 > 1e-9
                                if useRealData
                                    SF = 1.0;
                                    % Note sure if SF test is requied, but safer
                                    while ~intersectionFn(rayX, x0) & SF >= 1e-9
                                        if SF > 0.11
                                           SF = SF - 0.1; 
                                        else
                                           SF = 0.1*SF;
                                        end

                                        % note that t0 wasn't scaled back...
                                        [x, t] = ray_interpolation('5RKN', 'iso', x0', t0', lambda0*deltaS_final*SF, testCoords, ...
                                            testInds, lensRIVolume, tolerance, RIToUse);

                                        rayX = x'; rayT = t';
                                    end

                                else
                                    error('Add for treatment test cases')
                                end

                                if ~intersectionFn(rayX, x0)
                                    %%% Untested error
                                    error('TIR ray didnt step back enough - check this')
                                end

                                %%% Could update rIn, rOut and surfaceNormal and refraction calc after step back 
                                    % Unlikely these have changed much? - probably just if triangle intersect has changed...
                                lineDef = [rayX rayT];

                                [~, tempDistance, tempFaces] = intersectLineMesh3d(lineDef, grinSurface.vertices, grinSurface.faces);

                                inds2Use = find(abs(tempDistance) == min(abs(tempDistance)));

                                %%% Considered adjusting rayX to be intersect value but decided against 
                                 % changes in X are all dealt with via numerical integration

                                tempIndices = tempFaces(inds2Use);

                                if ~isempty(setdiff(faceIndices, tempIndices))

                                    %%% Untested error
                                    error('Different triangles used')
                                end
                            else
                               % Just set to start point as in original iterative final
                               rayT = t0; 
                               rayX = x0;
                            end

                            rayT = rayT/rIn;

                            % TIR
                            rayT = rayT-2*dot(surfaceNormal,rayT)*surfaceNormal;
                            rayT = rayT/norm(rayT)*rIn; % rIn is scaled back

                            if testPlot; plot3(rayX(1), rayX(2), rayX(3), 'm*'); end 

                            % Check good triangle intersect is used 
                            [dist, proj] = distancePointMesh(rayX, grinSurface.vertices, grinSurface.faces);

                            if abs(dist - min(abs(intersectDistance))) > 0.05e-3
                                % if this is too far off it could be a problem - check after 50 nm

                                if testPlot
                                    plot3(rayX(:,1), rayX(:,2), rayX(:,3), 'ko')
                                    plot3(proj(:,1), proj(:,2), proj(:,3), 'k*')
                                    trisurf(grinSurface.faces(faceIndices,:), grinSurface.vertices(:,1), grinSurface.vertices(:,2), grinSurface.vertices(:,3), ...
                                        'FaceColor','k','FaceAlpha',0.8);
                                end
                                [acos(dot((proj-rayX), rayT) / (norm(proj-rayX) * norm(rayT)))/pi*180 abs(dist - min(abs(intersectDistance)))*1e6]

                                % May need to find new points after step back (if used)
                                %%% Untested error
                                error('Difference between nearest point and intersect')
                            end

                            % Catch just in case...
                            if ~intersectionFn(rayX, x0)

                                error('Should be inside after TIR')
                            end

                            TIRFlag(iOrigin,1) = TIRFlag(iOrigin,1) + 1;

                            % Does not exit by default
                            inGraded = 1;

                            % check if ray is oscillating without position changing
                            if all(rayX == x0) & lastTIRStep == loopSteps-1 & all((rayT-t0)+rayDiff < 1e-9)
                                % note really sure how this can happen
                                    % occurs when lambda0 equals zero and ray actually doesn't move

                                if testPlot; plot3(rayX(1), rayX(2), rayX(3), 'rd', 'markersize',10); end

                                go = 0;

                                propogateFinalRay = 0;

                                warning('Ray direction was oscillating at TIR')
                            end

                            lastTIRStep = loopSteps;

                            rayDiff = rayT-t0;
                        end
                    elseif ~interfaceRefraction && ~isnan(lambda0)
                       if testPlot; plot3(rayX(1), rayX(2), rayX(3), 'mo'); end

                       inGraded = 0;

                       numberOfExits = numberOfExits + 1;
                    end
                   
                    if isnan(lambda0)
                       go = 0; 

                       propogateFinalRay = 0; 
                       
                       if testPlot; plot3(rayX(1), rayX(2), rayX(3), 'mx', 'markersize', 8); end
                    end
                end

                % At each z step, record path for plotting later 
                if rayX(3) >= zSteps(currentStep+1)
                    if currentStep + 1 < length(zSteps)
                        currentStep = currentStep + 1;

                        rayPathArray(:, currentStep, iOrigin) = rayX;

                        rayPathLengthArray(currentStep, iOrigin) = pathLength;

                    else
                        % Break ray if it has reach end of z steps
                        go = 0;
                    end

                    if testPlot; plot3(rayX(1), rayX(2), rayX(3), '.', 'color', rayCols(iOrigin,:)); end
                end
            end

            if rayT(3) < 0
               % kill ray if it reverses
               %%% Should flag for final plotting
               go = 0; 

               propogateFinalRay = 0;

               if testPlot; plot3(rayX(1), rayX(2), rayX(3), 'rd', 'markersize', 8); end

               warning('Ray has reversed')

               rayReversed(iOrigin) = 1;
            end
            
            % check if sim seems ot have gone too long
            if loopSteps > loopLim
                go = 0;
                
                propogateFinalRay = 0;
                
                % could change to a warning but probably should investigate
                error('Loop steps exceeded limit')
            end
            
            voxelX = round(rayX/voxelSize);

            % Check ray has exitied volume
            if any(voxelX > volumeSize) | any(voxelX < 1)
               go = 0; 
            end
            
            loopSteps = loopSteps + 1;
        end
  
        if numberOfExits == 0
            propogateFinalRay = 0;
        end
        
        % extend ray path to end of z steps
        if propogateFinalRay | useTestData
            finalRay(iOrigin,:) = rayT;

            if currentStep < length(zSteps)
                rayT = rayT/rayT(3);
                
                rayPathArray(1, currentStep+1:end, iOrigin) = rayX(1) + rayT(1).*(zSteps(currentStep+1:end) - rayX(3));
                rayPathArray(2, currentStep+1:end, iOrigin) = rayX(2) + rayT(2).*(zSteps(currentStep+1:end) - rayX(3));
                rayPathArray(3, currentStep+1:end, iOrigin) = rayX(3) + rayT(3).*(zSteps(currentStep+1:end) - rayX(3));

                if testPlot
                    plot3(rayPathArray(1, currentStep+1:length(zSteps), iOrigin), rayPathArray(2, currentStep+1:length(zSteps), iOrigin), ...
                        rayPathArray(3, currentStep+1:length(zSteps), iOrigin), '.', 'color', rayCols(iOrigin,:));
                end
            end
        elseif clearReverseRays

            finalRay(iOrigin,:) = NaN;
            finalRayTRefract(iOrigin,:) = NaN;

            finalIntersect(iOrigin,:) = NaN;
            firstIntersect(iOrigin,:) = NaN;
            rayPathArray(:,:,iOrigin) = NaN;
        end

        if testPlot
            pause(0.01)
        end
    end
    
    rayPathCells{aAngle} = rayPathArray; 
    finalIntersectCells{aAngle} = finalIntersect; 
    finalRayCells{aAngle} = finalRay; 
    finalRayTRefractCells{aAngle} = finalRayTRefract;
    rayReverseCells{aAngle} = rayReversed;
    TIRFlagCells{aAngle} = TIRFlag;
    firstIntersectCells{aAngle} = firstIntersect;
    
    timePerAngle(aAngle) = toc;
    
    warning('on', 'MATLAB:nearlySingularMatrix')

    if strcmp(interpType, '45RKN') 
        tolerance
        minDeltaS(minDeltaS < 1)*10^6
        mean(minDeltaS(minDeltaS < 1))*10^6
    end

end

round(timePerAngle/60)'
%% Plot results
% close all

if useRealData
    
    plot_cone_data
    
    if interConeValueToUse ~= 1.40
        ICAddition = sprintf('_IC_%i', interConeValueToUse*100);
    else
        ICAddition = '';
    end
    
    if ~isempty(externalPigmentValue)
        PGAddition = sprintf('_PG_%i', externalPigmentValue*100);
    else
        PGAddition = '';
    end
        
    % Do some saving
    if saveData
        saveFile = sprintf('%s/Data/%s%s%s_SIMDATA.mat', analysisFolder, metaFile(1:end-4), ICAddition,PGAddition);
        save(saveFile, 'rayPathCells', 'finalIntersectCells', 'finalRayCells', 'finalRayTRefractCells', 'rayReverseCells', 'timePerAngle', 'TIRFlagCells', 'firstIntersectCells', ...
            'dataFile', 'metaFile', 'incidenceAngle', 'receptorRadius', 'receptorAcceptance', 'exposedHeight', 'blockExposedRentry', ...
            'xSpacing', 'plotSpacing', 'interpType', 'tolerance', 'initialDeltaS', 'iterativeFinal', 'epsilon',  'interfaceRefraction', 'blockMultipleExits', 'limitToConeBase', 'clearReverseRays', 'trace3D', ...
            'alphaForIntercone', 'alphaForCone', 'alphaForCinC', 'dilateBorderRadius', 'flipVolume', 'flipSurfaces', ...
            'acceptancePercentage', 'focusXHeight', 'focusXRadius', 'focusYHeight', 'focusYRadius', 'colcHeight', 'colcRadius', 'rayReverseNum', 'TIRPercengtage', 'acceptanceAngle', ...
            'nOrigins', 'plotLineCols', 'colorTIR', 'plotOrigins', 'scaleBarsOnRayDiagram', 'plotColorsOnSummary', 'justPlotCenterRays', 'raysOnXZPlane', 'raysOnYZPlane', 'plotRIImageOnRayDiagram', 'acceptanceUsingReceptor', ...
            'coneTipZ', 'corneaZ', 'epicorneaZ', 'coneBaseZ', 'interConeValueUsed', 'externalPigmentValue', 'zSteps', ...
            'coneProfileR', 'coneProfileZ', 'cInCProfileR', 'cInCProfileZ', 'epicorneaProfileR', 'epicorneaProfileZ', 'interconeProfileR', 'interconeProfileZ', ...
            'voxelSize', 'volumeSize', 'xStartPointsOrig', 'rayOrigins', 'metaData')
    end
    
    if saveFigures
        imageFile = sprintf('%s/Figures/%s%s_SpotDiagram.pdf', analysisFolder, metaFile(1:end-4), ICAddition);
        exportgraphics(spotF,imageFile,'ContentType','vector') 
        
        imageFile = sprintf('%s/Figures/%s%s_RaysX.pdf', analysisFolder, metaFile(1:end-4), ICAddition);
        exportgraphics(rayFX,imageFile,'ContentType','vector')
        
        imageFile = sprintf('%s/Figures/%s%s_RaysY.pdf', analysisFolder, metaFile(1:end-4), ICAddition);
        exportgraphics(rayFY,imageFile,'ContentType','vector')
        
        imageFile = sprintf('%s/Figures/%s%s_COLCLines.pdf', analysisFolder, metaFile(1:end-4), ICAddition);
        exportgraphics(tipF,imageFile,'ContentType','vector')
        
        imageFile = sprintf('%s/Figures/%s%s_Summary.pdf', analysisFolder, metaFile(1:end-4), ICAddition);
        exportgraphics(sumF,imageFile,'ContentType','vector')
    end
    
elseif useTestData
    
    if createLunebergLens
        % Plot ray path in object
        figure; subplot(2,2,1); hold on; axis equal;
        view(0, 0)
        for iOrigin = 1:nOrigins
            rayPath = permute(rayPathArray(:, :, iOrigin), [2 1]);

            plot3(rayPath(:,1), rayPath(:,2), rayPath(:,3), 'color', rayCols(iOrigin,:));
        end
        % create and plot circle
        circleX = sin(0:(2*pi/100):2*pi)+volumeSize(1)/2*voxelSize;
        circleY = cos(0:(2*pi/100):2*pi)+volumeSize(2)/2*voxelSize;
        
        plot3(circleX, zeros(length(circleY),1), circleY, 'b')
        
        zlim([0 3])
        title('Ray path')
        
        % Plot analytic ray path - for on axis rays along central y axis
        subplot(2,2,2); hold on; axis equal;
        view(0, 0)
        trueRayPathArray = zeros(3, volumeSize(3), nOrigins)*NaN;
        
        lensCentre = volumeSize/2*voxelSize;
        
        startZPerOrigin = zeros(nOrigins,1);
        
        for iOrigin = 1:nOrigins
            tempRayPath = zeros(3, volumeSize(3))*NaN;
            
            if rayOrigins(iOrigin,2) ~= lensCentre(2) | ... 
                    any(startRayT - [0, 0, 1])
                %%% Can generalize to X component and off-axis rays from Eq 2 in Babayigit ... Turduev 2019
                
                error('Analytic solution only set up for Y on midline and parallel rays')
            end
                       
            % Check it's within fiber radius
            if abs(rayOrigins(iOrigin,1)-lensCentre(1)) < radius

                %Find first pixel in lens for ray
                startZ = find(RIFlagVolume(round(rayOrigins(iOrigin,1)/voxelSize),...
                    round(rayOrigins(iOrigin,2)/voxelSize), :));

                if ~isempty(startZ)
                    startZ = startZ(1);
                    startZPerOrigin(iOrigin) = startZ;

                    % Just a line up to enterance
                    tempRayPath(1,1:startZ) = rayOrigins(iOrigin,1)-lensCentre(1); % offset to sphere centre
                    tempRayPath(2,:) = rayOrigins(iOrigin,2) - lensCentre(2);
    %                 tempRayPath(3,:) = zSteps - lensCentre(3);
                    % Use actual z values
                    tempRayPath(3,:) = rayPathArray(3, :, iOrigin) - lensCentre(3);

                    % get non stepped z values for start and end
                    tempRayPath(3,startZ) = firstIntersect(iOrigin,3)- lensCentre(3);
                    tempRayPath(3,topZ) = finalIntersect(iOrigin,3)- lensCentre(3);

                    %%% Would be good to have solution as function of optical path length

                    % From eqn 3 in Turduev cloaking paper (eqn 2 for non parallel rays)
                    tempRayPath(1, startZ:topZ) = tempRayPath(1,startZ)*(tempRayPath(3,startZ)*tempRayPath(3,startZ:topZ) + ...
                        radius*sqrt(radius^2 + tempRayPath(3,startZ)^2 - tempRayPath(3,startZ:topZ).^2))/...
                        (tempRayPath(3,startZ)^2 + radius^2);

                    finalR = radius;
                    finalX = tempRayPath(1,startZ)*(tempRayPath(3,startZ)*finalR + ...
                        radius*sqrt(radius^2 + tempRayPath(3,startZ)^2 - finalR.^2))/...
                        (tempRayPath(3,startZ)^2 + radius^2);

                    % From eqn 4 in Turdev, note simplified given theta is 0
                        % Note error in differentiation, final part of sqrt in denominator should be +2*x0^2
                    finalP = 2*tempRayPath(1,startZ)*tempRayPath(3,startZ)/(2*tempRayPath(3,startZ)^2+2*radius^2) + ...
                        2*sqrt(2)*radius*finalR*tempRayPath(1,startZ)*-1/ ...
                        ((2*tempRayPath(3,startZ)^2+2*radius^2)*sqrt(2*radius^2-2*finalR^2 + 2*tempRayPath(3,startZ)^2));
                    finalT = [finalP, 0, 1];

                    % Propogate ray
                    tempRayPath(1, topZ+1:end) = (tempRayPath(3, topZ+1:end)-finalR)*finalT(1) + finalX;

                    % Add sphere centre back
                    tempRayPath(1,:) = tempRayPath(1,:) + lensCentre(1);
                    tempRayPath(2,:) = tempRayPath(2,:) + lensCentre(2);
                    tempRayPath(3,:) = tempRayPath(3,:) + lensCentre(3);   
                else
                   startZPerOrigin(iOrigin) = NaN; 
                end
            else
                % outside fiber, just propogate along z
                
                tempRayPath(1,:) = rayOrigins(iOrigin,1); 
                tempRayPath(2,:) = rayOrigins(iOrigin,2);
                tempRayPath(3,:) = zSteps;
            end
            
            plot3(tempRayPath(1,:), tempRayPath(2,:), tempRayPath(3,:), 'color', rayCols(iOrigin,:));
                
            trueRayPathArray(:,:,iOrigin) = tempRayPath;
        end
        plot3(circleX, zeros(length(circleY),1), circleY, 'b')
        zlim([0 3])
        
        title('Plot true ray path')

        % Plot X error along Z
        subplot(2,2,3); hold on;
        for iOrigin = 1:nOrigins
            if minDeltaS(iOrigin) < 1
                rayPath = permute(rayPathArray(:, :, iOrigin), [2 1]);
                
                tempRayPath = permute(trueRayPathArray(:, :, iOrigin), [2 1]);

                if ~isnan(startZPerOrigin(iOrigin))
                    % Put actual start and end points into recorded path
                    rayPath(startZPerOrigin(iOrigin),:) = firstIntersect(iOrigin,:);
                    rayPath(topZ,:) = finalIntersect(iOrigin,:);

                    plot((tempRayPath(:,1) - rayPath(:,1))*10^6, zSteps, 'color', rayCols(iOrigin,:))

                    plot((tempRayPath(startZPerOrigin(iOrigin),1) - rayPath(startZPerOrigin(iOrigin),1))*10^6, zSteps(startZPerOrigin(iOrigin)), 'rx')
                    plot((tempRayPath(topZ,1) - rayPath(topZ,1))*10^6, zSteps(topZ), 'rx')
                end
            end
        end
        ylim([0 3])
%         xlim([-50 50])

        xlabel('Path error in nm')
        title('Ray error')
        
        % plot spot diagram
        subplot(2,2,4); hold on; axis equal;
        centerRef = volumeSize/2*voxelSize;
        for iOrigin = 1:nOrigins
            if finalIntersect(iOrigin,3) ~= 0 & ~isnan(startZPerOrigin(iOrigin))
                plot((finalIntersect(iOrigin,1)-centerRef(1))*10^6, (finalIntersect(iOrigin,3)-radius-centerRef(3))*10^6, 'o',  'color', rayCols(iOrigin,:));
            end
        end

%         xlim([-50 50]); 
%         ylim([-50 50])
        xlabel('Spot error in nm')
        title('Spot diagram')
    else
        if plotReferenceFiber
            cd('/Users/gavintaylor/Documents/Matlab/Git_versioned_April_27/Crab_cone_git/Data')
            
            % These are for x start 0.5, given sqrt(2.5-(r/R)^2) RI profile 
            
            % For iterative solution
            nishidatError = readmatrix('Iterative.csv');
            
            % vDash is supposed to be refracted ray, but I think v is.
            nishidatPath = readmatrix('v.csv');
        end
        
        % Plot ray path in object
        figure; subplot(2,2,1); hold on; axis equal;
        view(0, 0)
        
        for iOrigin = 1:nOrigins
            rayPath = permute(rayPathArray(:, :, iOrigin), [2 1]);

            plot3(rayPath(:,1), rayPath(:,2), rayPath(:,3), 'color', rayCols(iOrigin,:));
        end
        
        if plotReferenceFiber
            
            % Align bottom inds in X
            zPoints = nishidatPath(:,1)+volumeSize(3)/2*voxelSize-fiberLength/2;
            xPoints = nishidatPath(:,2)+volumeSize(1)/2*voxelSize;
            
            zInds = find(zPoints < volumeSize(3)/2*voxelSize-fiberLength/2);
            
            xPoints = xPoints - mean(xPoints(zInds)) + rayOrigins(1);
            
            nishidatPath(:,1) = xPoints;
            nishidatPath(:,2) = zPoints;
            
            plot3(nishidatPath(:,1), ones(size(nishidatPath,1),1)*volumeSize(2)/2*voxelSize,...
                nishidatPath(:,2), 'r')
        end
        
        title('Plot ray path')
        
        % plot lines at top and bottom
        line([-1 1]*radius+volumeSize(1)/2*voxelSize, [-1 1]*radius+volumeSize(1)/2*voxelSize, ...
            [-1 -1]*fiberLength/2 + volumeSize(3)/2*voxelSize, 'color', 'b')
        
        line([-1 1]*radius+volumeSize(1)/2*voxelSize, [-1 1]*radius+volumeSize(1)/2*voxelSize, ...
            [1 1]*fiberLength/2 + volumeSize(3)/2*voxelSize, 'color', 'b')
        
        % Plot analytic ray path - for on axis rays along central y axis
        subplot(2,2,2); hold on; axis equal;
        view(0, 0)
        
        trueRayPathArray = zeros(3, volumeSize(3), nOrigins)*NaN;
        
        fiberCentre = volumeSize(1:3)/2*voxelSize;
        
        warning('off', 'MATLAB:nearlySingularMatrix')
        
        % Test compared to paper       
        for iOrigin = 1:nOrigins
            tempRayPath = zeros(3, volumeSize(3))*NaN;
            
            if rayOrigins(iOrigin,2) ~= fiberCentre(2) | ... 
                    any(startRayT/norm(startRayT) - [0, 0, 1])
                %%% Can gernalize to X component and off-axis/helical rays from Section 5.3 in Merchland 1978
                
                error('Analytic solution only set up for Y on midline and parallel rays')
            end
            
            % Check it's within fiber radius
            if abs(rayOrigins(iOrigin,1)-fiberCentre(1)) < radius
            
                % Just a line up to enterance
                tempRayPath(1,1:bottomZ) = rayOrigins(iOrigin,1)-fiberCentre(1); % offset to fiber centre
                tempRayPath(2,:) = rayOrigins(iOrigin,2)-fiberCentre(2);
                % Use actual z value after step.
                tempRayPath(3,:) = rayPathArray(3, :, iOrigin) - (fiberCentre(3) - fiberLength/2);

                % get non stepped z values for start and end
                tempRayPath(3,bottomZ) = firstIntersect(iOrigin,3) - (fiberCentre(3) - fiberLength/2);
                tempRayPath(3,topZ) = finalIntersect(iOrigin,3) - (fiberCentre(3) - fiberLength/2);
                
                % Ray path in fiber - Nishidate 2011, eqn 31 
                if plotReferenceFiber
                    nishidateSolution = tempRayPath(1, :);
                    nIn = sqrt(2.5-0.5^2);
                    nishidateSolution(bottomZ:topZ) = 0.5*cos(tempRayPath(3,bottomZ:topZ)/nIn);
                    
                    yc = 0.5*cos(fiberLength/nIn);
                    qc = -0.5*sin(fiberLength/nIn);
                    nEnd = sqrt(2.5-yc^2);
                    
                    % This solutions seems correct, but is more like v than v' in figure 8
                    nishidateSolution(topZ+1:end) = nEnd*qc/sqrt(nIn^2+(1-nEnd^2)*qc^2)*(tempRayPath(3, topZ+1:end)-fiberLength) + yc;

                end
                
                % Gives correct result From Merchland 5.41 and other refs 
                    %%% Was previously not dividing by RI/n0
                firstCoord = [rayOrigins(iOrigin,1), fiberCentre(2), fiberCentre(3) - fiberLength/2];
                [~, rStart] = numerical_dT_dt(firstCoord, volCoords(RIFlagInds,:), RIFlagInds, lensRIVolume, NaN);  

                tempRayPath(1, bottomZ:topZ) = tempRayPath(1,bottomZ)*cos(tempRayPath(3,bottomZ:topZ)*beta/rStart);
                
                % Analytic expresion for ray after exit
                % Get final Coord and RI
                finalX = tempRayPath(1,bottomZ)*cos(fiberLength*beta/rStart);
                finalCoord = [finalX + fiberCentre(1), fiberCentre(2), fiberCentre(3) + fiberLength/2];
                
                [~, rEnd] = numerical_dT_dt(finalCoord, volCoords(RIFlagInds,:), RIFlagInds, lensRIVolume, NaN);  

                % Get final vector and normalize - 5.43 in Merchland
                finalP = -beta*tempRayPath(1,bottomZ)*sin(fiberLength*beta/rStart);
                
                % Note rStart is divded by rEnd before refraction... But matches what I did in ray-tracing?
                finalT = [finalP, 0, rStart]/rEnd;
                % Division by rEnd seems to effectively normalize
%                 finalT = finalT/norm(finalT);
                
                surfaceNormal = [0 0 -1];
                rOut = exteriorRI;
                
                % Refract
                nRatio = rEnd/rOut;
                cosI = -dot(surfaceNormal, finalT);
                sinT2 = nRatio^2*(1-cosI^2);
                cosT = sqrt(1-sinT2);

                if sinT2 < 1
                    % Assuming all refracted, none reflected
                    refractT = nRatio*finalT + (nRatio*cosI-cosT)*surfaceNormal;
                else
                    % TIR
                    refractT = finalT-2*dot(surfaceNormal,finalT)*surfaceNormal;
                end
                
                refractT = refractT/norm(refractT);
                
                % Propogate ray
                % Normalize to ray z as steps are along z
                refractT = refractT/refractT(3);
                tempRayPath(1, topZ+1:end) = (tempRayPath(3, topZ+1:end)-fiberLength)*refractT(1) + finalX;
                
%                % Add fiber centre back for X
                 tempRayPath(1,:) = tempRayPath(1,:) + fiberCentre(1);
                 tempRayPath(2,:) = tempRayPath(2,:) + fiberCentre(2);
                 tempRayPath(3,:) = tempRayPath(3,:) + (fiberCentre(3) - fiberLength/2);
                 
            else
                % outside fiber, just propogate along z
                
                tempRayPath(1,:) = rayOrigins(iOrigin,1); 
                tempRayPath(2,:) = rayOrigins(iOrigin,2);
                tempRayPath(3,:) = zSteps;
            end
            
            plot3(tempRayPath(1,:), tempRayPath(2,:), tempRayPath(3,:), 'color', rayCols(iOrigin,:)); 
            
            trueRayPathArray(:,:,iOrigin) = tempRayPath;
        end   
        
        if plotReferenceFiber
            plot3(nishidateSolution+fiberCentre(1), tempRayPath(2,:), tempRayPath(3,:), 'r--')
        end
        
        
        warning('on', 'MATLAB:nearlySingularMatrix')
        
        % plot lines at top and bottom
        line([-1 1]*radius+volumeSize(1)/2*voxelSize, [-1 1]*radius+volumeSize(1)/2*voxelSize, ...
            [-1 -1]*fiberLength/2 + volumeSize(3)/2*voxelSize, 'color', 'b')
        
        line([-1 1]*radius+volumeSize(1)/2*voxelSize, [-1 1]*radius+volumeSize(1)/2*voxelSize, ...
            [1 1]*fiberLength/2 + volumeSize(3)/2*voxelSize, 'color', 'b')
        
        title('Plot true ray path')
        
        % Plot X error along Z
        subplot(2,2,3); hold on; %axis equal;

        for iOrigin = 1:nOrigins
            if minDeltaS(iOrigin) < 1
                rayPath = permute(rayPathArray(:, :, iOrigin), [2 1]);

%                 rayPath(bottomZ,:) = firstIntersect(iOrigin,:);
%                 rayPath(topZ,:) = finalIntersect(iOrigin,:);
                
                tempRayPath = permute(trueRayPathArray(:, :, iOrigin), [2 1]);

                plot((tempRayPath(:,1) - rayPath(:,1))*10^6, zSteps, 'color', rayCols(iOrigin,:))
                
                plot((tempRayPath(bottomZ,1) - rayPath(bottomZ,1))*10^6, zSteps(bottomZ), 'rx')
                plot((tempRayPath(topZ,1) - rayPath(topZ,1))*10^6, zSteps(topZ), 'rx')
            end
        end
       
       if plotReferenceFiber
          plot(nishidatError(:,2), nishidatError(:,1)+volumeSize(3)/2*voxelSize-fiberLength/2, 'r') 
          
          % Recalcualte error for traced vs calculated
          % First interp to same zsteps
          interpNishidatPath = interp1(nishidatPath(:,2), nishidatPath(:,1), tempRayPath(:,3), 'linear', 0);
          
          interpNishidatPath(interpNishidatPath == 0) = NaN;
          
          interpNishidatPath = interpNishidatPath - volumeSize(1)/2*voxelSize;
          
%           plot((interpNishidatPath - nishidateSolution')*10^6, zSteps, 'r--')

       end
       
        xlim([-1 1]*100)
        title('Ray error')
        
        subplot(2,2,4); hold on;
        
        colsPeriod = winter(numPeriods);
        
        for iOrigin = 1:nOrigins
            for jPeriod = 1:numPeriods
                if minDeltaS(iOrigin) < 1
                    plot((periodSpotsArray(1, jPeriod, iOrigin)-fiberCentre(1))*10^6, periodSpotsArray(3, jPeriod, iOrigin),...
                        'o', 'color', rayCols(iOrigin,:));
                end
            end
        end
        title('Spot diagram')
    end
end

%% 
function  intersect = surfaceIntersectFunction(volumeFull, volumeBorder, x, x0, scale, surface) 
    % Return 1 if inside
    
   voxelX = round(x/scale);
    
   %%% All mesh functions (distance, intersect and inpoly) are slow
        % inpoly is slowest but accepts mesh normals - provide to speed up
   
   if ~( any(voxelX > size(volumeFull)) | any(voxelX < 1) )
        % Check if on border voxel and fine test required
        if volumeBorder(voxelX(1), voxelX(2), voxelX(3))

            % Check points aren't coincident
            if norm(x-x0) > 0
                
                % Check angle to surface
                %%% This takes nearly as long as line call 
                    %%% could possibly speed up if vector to surface is pre-calculated by voxel - accuracy?
                [~, proj] = distancePointMesh(x, surface.vertices, surface.faces);
                
                angle = acos(dot((proj-x), (x-x0)) / (norm(proj-x) * norm(x-x0)))/pi*180;
                
%              [abs(angle - 90) inpolyhedron(surface, x, 'tol', 0, 'flipNormals',true)]
                
                %%% If this too low, will probably get lambda error with all inds on one side
                if  abs(angle - 90) > 3      
                % If large use line intersects - faster
                    
                    lineDef = [x (x-x0)/norm((x-x0))];  

                    % tolerance of 1e-15 sometimes misses true co-planar intersections
                        % 1e-21 catches these but could cause problems else where

                    [tempPoints, intersectDistance] = intersectLineMesh3d(lineDef, surface.vertices, surface.faces);  

                    % for debug            
        %             plot3(tempPoints(:,1), tempPoints(:,2), tempPoints(:,3), 'rx')

                    if length(intersectDistance) > 1
                        zeroInds = find(intersectDistance == 0);
                        
                        % If any inds are zero this probably wont work
                        if isempty(zeroInds)
                            
                            backInds = find(intersectDistance < 0);
                            frontInds = find(intersectDistance > 0);

                            if isempty(backInds) | isempty(frontInds)
                                % all inds on one side, so should be not in volume
                                %  but this is not always reliable for coplanar intersections
                                intersect = 0;

                            elseif ~isempty(backInds) & ~isempty(frontInds)
                               % must be in volume as intersects on both sides
                               intersect = 1;
                            end
                        else
                            intersect = inpolyhedron(surface, x, 'tol', 0, 'flipNormals',true);
                        end

                    elseif length(intersectDistance) == 1
                        % can't be resolved unambigiously by distances
                            % This is really slow! Note matlab normals seem to be opposite to convention
                            % Tolerance cant be used - lambda will fall short by roughly tolerance amount
                        intersect = inpolyhedron(surface, x, 'tol', 0, 'flipNormals',true);

                    elseif isempty(intersectDistance)
                       % Not intersect so can't be in GRIN volume (?)

                       % Test this just in case
                       intersect = inpolyhedron(surface, x, 'tol', 0, 'flipNormals',true);
                       if intersect ~= 0
                          error('No line intersect but inside surface?') 
                       end
                    end
                    
                else
                   % If small angle problems can occur with line intersects
                    intersect = inpolyhedron(surface, x, 'tol', 0, 'flipNormals',true);
                end
            else
                % Points are equal (can't make line) so just test if inside polyhedron
                intersect = inpolyhedron(surface, x, 'tol', 0, 'flipNormals',true);
%               warning('Query points are equal in intersect')
            end
            
        elseif volumeFull(voxelX(1), voxelX(2), voxelX(3))
             % just inside    
            intersect = 1; 
        else
            % just outside
            intersect = 0;
        end
   else
      intersect = 0; 
   end
end

function lambda = surfaceLambdaFunction(x1, x0, surface)

    if norm(x1-x0) > 0
        lineDef = [x0 (x1-x0)/norm(x1-x0)];  

        [intersectPoints, intersectDistance] = intersectLineMesh3d(lineDef, surface.vertices, surface.faces);
        % plot3(tempPoints(:,1), tempPoints(:,2), tempPoints(:,3), 'rx')

        % Have some threshold for zero...
            % Note that lambda0 threshold is 1e-9. 
            % In iterative final delta can go to 1e-6 to 1e-12 - threshold here can make main loop threshold redundant
            %%% Used 1e-9 for both, cylinder cone, (linear and on?), 1e-15 for cylinder cylinder
        zeroInd = find(abs(intersectDistance) < 1e-9);
        
        % Sort out intersections
        if length(intersectDistance) > 1
            
            if isempty(zeroInd)
                % If no zeros there should be ind/s in front and ind/s behind
                backInds = find(intersectDistance < 0);
                frontInds = find(intersectDistance > 0);

%                 [~, closestBackInd] = min(abs(intersectDistance(backInds)));
                [~, closestFrontInd] = min(intersectDistance(frontInds));

                if isempty(backInds) | isempty(frontInds)
                    % This should be called while x0 is inside mesh and x1 is outside
                    %  but this is not always reliable for coplanar intersections, so check with angle calc
                    
                    [closestDist, proj] = distancePointMesh(x0, surface.vertices, surface.faces);
                
                    angle = acos(dot((proj-x0), (x1-x0)) / (norm(proj-x0) * norm(x1-x0)))/pi*180;
                    if  abs(angle - 90) < 3
                        % don't really know actualy distance ...
                        
                        if min(abs(intersectDistance)) < norm(x1-x0)
                            % Use intersect distance if less than normal
                            nearestDist = min(abs(intersectDistance));
                            
                            warning('Distance to intersect unknown, using absolute in lambda0 calc')
                        elseif closestDist < norm(x1-x0)
                            % Use closest distance if less than normal
                            nearestDist = closestDist;
                            
                            warning('Distance to intersect unknown, using closest in lambda0 calc')
                        else
                            %%% Could just use 0 or 0.5?
                            error('Unsure which distance to use in lambda0 calc')
                        end
                        
                    else
                        % error has been a sign of problems with intersection calculation, sometimes several steps before
                        [inpolyhedron(surface, x1, 'tol', 0, 'flipNormals',true) ...
                            inpolyhedron(surface, x1, 'tol', 0, 'flipNormals',true) norm(x1-x0) min(abs(intersectDistance))]
                        
                        error('Check usage - all inds either in front or behind and angle is large')
                    end
                else
                    nearestDist = intersectDistance(frontInds(closestFrontInd));
                end
            else
                % if there is a zero distance lambda will be zero
                nearestDist = 0;
            end

        elseif length(intersectDistance) == 1 
            if isempty(zeroInd)
                if intersectDistance >= 0
                    nearestDist = intersectDistance(1);
                else
                    if x0(3) < 0.193
                        %%% A few missing faces at very base of cone
                            % Ray exit should be calculated properly, but won't be able to get lambda.
                            %%% Remove if faces fixed...
                        nearestDist = NaN;
                        warning('Ray went through missing face at cone base')
                    else
                        error('Only one intersect and behind x0') 
                    end
                end
            else
                nearestDist = 0;
            end

        elseif isempty(intersectDistance)
            error('No intersect')
        end

        % Same as distance value in this context
%         lambda = norm(intersectPoints(nearestInd,:)-x0)/norm(x1-x0);
        lambda = nearestDist/norm(x1-x0);
        
        % Seem to need a threshold here as well...
        if lambda > 1 && lambda < 1 + 1e-9
            lambda = 1;
        elseif lambda < 0 && lambda > -1e-9
            lambda = 0;
        end
        
    else
       % points are equal, set to zero to break
       lambda = 0;
%        warning('Query points are equal in lambda')
    end
    
    if lambda > 1
        [lambda inpolyhedron(surface, x0, 'tol', 0, 'flipNormals',true) ...
            inpolyhedron(surface, x1, 'tol', 0, 'flipNormals',true) norm(x1-x0) min(abs(intersectDistance))]
        
        error('Large lambda step, intersect was triggered too early, check threshold in intersection function')
    elseif lambda < 0
       error('Check this') 
    end
end
