% Set up to do both tracing on data and have testing mode on either
% lundaberg lens or graded fibre, as in nishdate 2011 papers
    
clear; 
clc; close all

%% Set parameters

useRealData = 1;

if useRealData
    dataFolder = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/AnalysisVolumes/';

    % New radial
%     dataFile = 'Volume_Cone_1000_nm_Cone_0_SD_GRIN_radial.mat';
%     metaFile = 'Cone_1000_nm_Cone_0_SD_GRIN_radial.mat';
    
    % Cylinder
%     dataFile = 'Volume_Cylinder_1000_nm_Cone_0_SD_GRIN_cylinder.mat';
%     metaFile = 'Cylinder_1000_nm_Cone_0_SD_GRIN_cylinder.mat';
    
    % Orignal Radial
    dataFile = 'Volume_ConeOriginal_1000_nm_Cone_0_SD_GRIN_radial.mat';
    metaFile = 'ConeOriginal_1000_nm_Cone_0_SD_GRIN_radial.mat';

    % Usually saved pointing down
    flipVolume = 1;    
    flipSurfaces = 1;
    
    createGradedFiber = 0;
    createLunebergLens = 0;
    plotReferenceFiber = 0;
    useTestData = 0;
    
    testOrigin = []; [5 6 7 8];
else
    useTestData = 1;
    
    radius = 1; %mm - used for both lens and fiber
    voxelSize = 0.065; % mm

    exteriorRI = 1;

    createLunebergLens = 0;
        
    createGradedFiber = 1;
        %%% Add flag to create gaussian profile then update solution and period
        fiberLength = 10;
        n0 = sqrt(2.5); 
        alpha = 1/sqrt(2.5); 0.7; 
        beta = alpha*n0;
        plotReferenceFiber = 0; % will plot data from Nishidate and only trace 0.5
end

xSpacing = 10; % in voxels...

% Options, 1, 4S, 4RKN, 5RKN, 45RKN 
% Difference 4 to 5 is quite small
interpType = '4S'; %'4S';
    %This doesn't have big effect on error < 10^-9, ok even at 10^-6 but then collapses
        % May effect peripheral rays more
    tolerance = 10^-9; % Nishidate uses 10^-12, seems a bit too tight
    
    % Seems good to start a bit lower than the general deltaS
        % in fixed step, 10^-3 vs. 10^-4 gives rough order of magnitude error
    initialDeltaS = 0.5*10^-3;
    
    % This keeps spot size surpsingly small, even on loose tolerance
        % (I guess it ends up stepping further at loose tolerance...)
    iterativeFinal = 0;  
    
incidenceAngle = 0; % deg, in XZ - plane    
    
%error term, sqrt of machine precision, from Nishidate paper sqrt(1.49e-8)
    % Take geometric mean of single and double precision, double takes ages to evaluate
epsilon = sqrt(geomean([1.19e-7 2.22e-16])); 

% Should be on, but can switch off to test    
interfaceRefraction = 1;  

% Ray that has exited will error on re-entry
blockMultipleExits = 1;

% Rays won't continue if it hits intercone base
limitToConeBase = 1;

% extend on final plots - only added for real data
extendRayLength = 2; % mm
    
testPlot = 1;
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
        
    else
        RIFlagVolume = lensRIVolume < -1;
        
        lensRIVolume(lensRIVolume == metaData.interconeValue) = 1.47;
    end
    
    lensRIVolume(RIFlagVolume) = -lensRIVolume(RIFlagVolume);

    % make meshes for each layer
    meshAngles = (0:2:360)/180*pi;

    %%% Note, border cones won't be defined using this process
        % Add an error if they are present (rotated volume won't work anyway)
    
    % For cornea - flat
    % single Z at end
    tempZ = (metaData.coneLengthToUse + metaData.outerCorneaLengthToUse + metaData.epicorneaLengthToUse)* ...
        metaData.voxSize/metaData.voxSize3D;

    [xGrid, yGrid, zGrid] = meshgrid(1:volumeSize(1), 1:volumeSize(2), tempZ);
    xGrid = xGrid(:); yGrid = yGrid(:); zGrid = zGrid(:);
    
    corneaVertices = [xGrid, yGrid, zGrid];
    corneaVertices(:,3) = corneaVertices(:,3) + metaData.tipOffset*metaData.voxSize/metaData.voxSize3D;
    
    % Get 2D triangulation
    tempTriangulation = delaunayTriangulation(corneaVertices(:,1:2));

    corneaSurface.faces = tempTriangulation.ConnectivityList;
    corneaSurface.vertices = corneaVertices;

    if flipSurfaces
        corneaSurface.vertices(:,3) = -(corneaVertices(:,3)-volumeSize(3)/2)+volumeSize(3)/2;
    end
    corneaBorderVolume = polygon2voxel(corneaSurface, volumeSize, 'none');

    
    
    % for top of outer cornea - can be curved
    if metaData.displayProfiles.EpicorneaCone
        warning('needs testing')
        
        %%% needs testing before copy down - confirm against code for removing base of cone...
        
        % has a curve
        tempProfile = metaData.epicorneaProfileToUse*metaData.voxSize/metaData.voxSize3D;
        tempZ = metaData.corneaXRef*metaData.voxSize/metaData.voxSize3D;
        
        tempZ(isnan(tempProfile)) = [];
        tempProfile(isnan(tempProfile)) = [];

        % Adjust to actual length
        tempZ(1) = (metaData.coneLengthToUse + metaData.outerCorneaLengthToUse)*metaData.voxSize/metaData.voxSize3D;
        
        epicorneaVertices = zeros(length(tempProfile)*length(meshAngles),3);

        profileMax = max(tempProfile);
        
        for i = 1:length(tempProfile)
            for j = 1:length(meshAngles)
                epicorneaVertices((i-1)*length(meshAngles)+j,:) = [tempProfile(i)*sin(meshAngles(j)), ...
                    tempProfile(i)*cos(meshAngles(j)), tempZ(i)];
            end
        end
        
        % Get cap to pad back
        [xGrid, yGrid, zGrid] = meshgrid(-floor(profileMax):ceil(profileMax), -floor(profileMax):ceil(profileMax), tempZ(end));

        xGrid = xGrid(:); yGrid = yGrid(:); zGrid = zGrid(:);
        rGrid = sqrt(xGrid.^2 + yGrid.^2);

        endInds = find(rGrid > tempProfile(end));

        epicorneaVertices = [[xGrid(endInds), yGrid(endInds), zGrid(endInds)]', epicorneaVertices']';

        epicorneaVertices(:,1) = epicorneaVertices(:,1) + volumeSize(1)/2;
        epicorneaVertices(:,2) = epicorneaVertices(:,2) + volumeSize(2)/2;
        epicorneaVertices(:,3) = epicorneaVertices(:,3) + metaData.tipOffset*metaData.voxSize/metaData.voxSize3D;

        % Make another layer at base to force 3D
        [xGrid, yGrid, zGrid] = meshgrid(1:volumeSize(1), 1:volumeSize(2), tempZ(1)-1);
        xGrid = xGrid(:); yGrid = yGrid(:); zGrid = zGrid(:);
        
        epicorneaVertices = [epicorneaVertices [xGrid, yGrid, zGrid]' ]';
        forcedInds = find(epicorneaVertices(:,3) == tempZ(1)-1);
        
        tempTriangulation = delaunayTriangulation(epicorneaVertices);
        [tempFaces, tempVertices] = freeBoundary(tempTriangulation); 
        
        epicorneaSurface.faces = tempFaces;
        epicorneaSurface.vertices = tempVertices;
        
    else
       % just flat, can copy cornea with adjusted Z 
       epicorneaVertices = corneaVertices; 
       epicorneaVertices(:,3) = epicorneaVertices(:,3) - metaData.epicorneaLengthToUse*metaData.voxSize/metaData.voxSize3D;
       
       epicorneaSurface = corneaSurface;
       epicorneaSurface.vertices = epicorneaVertices;
    end
    
    if flipSurfaces
        epicorneaSurface.vertices(:,3) = -(epicorneaVertices(:,3)-volumeSize(3)/2)+volumeSize(3)/2;
    end
    epicorneaBorderVolume = polygon2voxel(epicorneaSurface, volumeSize, 'none');

    
    
    % for base outer cornea/intercone around cone -  flat
    % Make flat grid as for cornea
    tempZ = (metaData.coneLengthToUse)* metaData.voxSize/metaData.voxSize3D;

    [xGrid, yGrid, zGrid] = meshgrid((1:volumeSize(1))-volumeSize(1)/2, (1:volumeSize(2))-volumeSize(2)/2, tempZ);
    xGrid = xGrid(:); yGrid = yGrid(:); zGrid = zGrid(:);
    
    % Shrink vertices that are inside radius at base of cone
    rGrid = sqrt(xGrid.^2 + yGrid.^2);
    
    if ~isnan(metaData.coneProfileToUse(end))
        coneBaseRadius = metaData.coneProfileToUse(end)*metaData.voxSize/metaData.voxSize3D;
        verticesToRemove = find(rGrid < coneBaseRadius);
    else
        error('Need to get last good value in cone profile')
    end
    
    cInCVertices = [xGrid, yGrid, zGrid];
    
    % Step in slightly so these don't get incorperated into faces along base of cone
    cInCVertices(verticesToRemove, 1:2) = cInCVertices(verticesToRemove, 1:2)*0.9;
    
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

    if flipSurfaces
        cInCSurface.vertices(:,3) = -(cInCVertices(:,3)-volumeSize(3)/2)+volumeSize(3)/2;
    end
    cInCBorderVolume = polygon2voxel(cInCSurface, volumeSize, 'none');



    %%% Would be good to make a funciton for doing the rotation in these...
    % for cone - is curved note voxSize is actually 2D pixel size
    tempProfile = metaData.coneProfileToUse*metaData.voxSize/metaData.voxSize3D;
    tempZ = metaData.coneXRef*metaData.voxSize/metaData.voxSize3D;

    tempZ(isnan(tempProfile)) = [];
    tempProfile(isnan(tempProfile)) = [];
    
    % Adjust to actual length
    tempZ(end) = metaData.coneLengthToUse*metaData.voxSize/metaData.voxSize3D;
    
    profileMax = max(tempProfile);

    coneVertices = zeros(length(tempProfile)*length(meshAngles),3);

    for i = 1:length(tempProfile)
        for j = 1:length(meshAngles)
            if ~metaData.createCylinder
                coneVertices((i-1)*length(meshAngles)+j,:) = [tempProfile(i)*sin(meshAngles(j)), ...
                    tempProfile(i)*cos(meshAngles(j)), tempZ(i)];
            else
                coneVertices((i-1)*length(meshAngles)+j,:) = [profileMax*sin(meshAngles(j)), ...
                    profileMax*cos(meshAngles(j)), tempZ(i)];
            end
        end
    end
    
    % Get caps to pad front and back, otherwise there are degenerate faces
    [xGrid, yGrid, zGrid] = meshgrid(-floor(profileMax):ceil(profileMax), -floor(profileMax):ceil(profileMax), 1);

    xGrid = xGrid(:); yGrid = yGrid(:); zGrid = zGrid(:);
    rGrid = sqrt(xGrid.^2 + yGrid.^2);

    firstInds = find(rGrid < tempProfile(1));
    endInds = find(rGrid < tempProfile(end));
    
    if metaData.displayProfiles.CinC
        error('Need to make profile at base of cone')
        
        coneVertices = [ [xGrid(firstInds), yGrid(firstInds), tempZ(1)*zGrid(firstInds)]', coneVertices'];
    else
        coneVertices = [ [xGrid(firstInds), yGrid(firstInds), tempZ(1)*zGrid(firstInds)]', ...
            coneVertices', [xGrid(endInds), yGrid(endInds), (tempZ(end))*zGrid(endInds)]']'; % (tempZ(end)+1) to remove cap
    end
    
    % Store end cap inds
    forcedCapInds = find(coneVertices(:,3) == tempZ(end)+1);
    
    % Step in slightly so these don't get incorperated into faces along base of cone
    if ~isempty(forcedCapInds)
        coneVertices(forcedCapInds, 1:2) = coneVertices(forcedCapInds, 1:2)*0.9;
    end
    
    coneVertices(:,1) = coneVertices(:,1) + volumeSize(1)/2;
    coneVertices(:,2) = coneVertices(:,2) + volumeSize(2)/2;
    coneVertices(:,3) = coneVertices(:,3) + metaData.tipOffset*metaData.voxSize/metaData.voxSize3D;
    
    tempTriangulation = delaunayTriangulation(coneVertices);
    [tempFaces, tempVertices] = freeBoundary(tempTriangulation); % this seems to flip vertices in z
    
    if ~isempty(forcedCapInds)
        % remove forced inds at base
        % not all of the original forced inds are used in the final trinagulation
        usedForcedInds = find(tempVertices(:,3) == coneVertices(forcedCapInds(1),3));

        facesToRemove = zeros(length(tempFaces),1, 'logical');
        for i = 1:length(usedForcedInds)
            % Flag for each column
            facesToRemove(tempFaces(:,1) == usedForcedInds(i)) = 1;
            facesToRemove(tempFaces(:,2) == usedForcedInds(i)) = 1;
            facesToRemove(tempFaces(:,3) == usedForcedInds(i)) = 1;
        end

        tempFaces(facesToRemove,:) = [];
    %     tempVertices(usedForcedInds,:) = [];

        coneVertices(forcedCapInds,:) = [];
    end
    
    grinSurface.faces = tempFaces;
    grinSurface.vertices = tempVertices;

    % Get cone border volume
    if flipSurfaces
        grinSurface.vertices(:,3) = -(tempVertices(:,3)-volumeSize(3)/2)+volumeSize(3)/2;
    end
    coneBorderVolume = polygon2voxel(grinSurface, volumeSize, 'none');
    
    

    % for inner intercone - is curved
    if any(metaData.exposedInterconeProfileToUseLeft - metaData.exposedInterconeProfileToUseRight)
        error('Right and left intercone profiles arent equal - rotated volume also wont work')
    end
    
    profileToUse = find(~isnan(metaData.exposedInterconeProfileToUseLeft));
    
    if ~metaData.createCylinder
        tempProfile = (metaData.exposedInterconeProfileToUseLeft(profileToUse)+metaData.coneProfileToUse(profileToUse))*metaData.voxSize/metaData.voxSize3D;
    
        % Shrink 1st in profile to cone
        % Tried adding a step before hand, but always caused problems somewhere
        tempProfile(1) = metaData.coneProfileToUse(profileToUse(1))*metaData.voxSize/metaData.voxSize3D;
    else
        tempProfile = (metaData.exposedInterconeProfileToUseLeft(profileToUse)+max(metaData.coneProfileToUse))*metaData.voxSize/metaData.voxSize3D;
    
        tempProfile(1) = max(metaData.coneProfileToUse)*metaData.voxSize/metaData.voxSize3D;
    end
    
    tempZ = metaData.coneXRef(profileToUse)*metaData.voxSize/metaData.voxSize3D;
    
    interconeVertices = zeros(length(tempProfile)*length(meshAngles),3);

    profileMax = max(tempProfile);  
    
    for i = 1:length(tempProfile)
        for j = 1:length(meshAngles)
            interconeVertices((i-1)*length(meshAngles)+j,:) = [tempProfile(i)*sin(meshAngles(j)), ...
                tempProfile(i)*cos(meshAngles(j)), tempZ(i)];
        end
    end
    
    % Get fake cap to pad top
    [xGrid, yGrid, zGrid] = meshgrid(-floor(profileMax):ceil(profileMax), -floor(profileMax):ceil(profileMax), (tempZ(1)-1));

    xGrid = xGrid(:); yGrid = yGrid(:); zGrid = zGrid(:);
    rGrid = sqrt(xGrid.^2 + yGrid.^2);

    startInds = find(rGrid < tempProfile(1));

    interconeVertices = [ [xGrid(startInds), yGrid(startInds), zGrid(startInds)]', interconeVertices']';

    forcedCapInds = find(interconeVertices(:,3) == tempZ(1)-1);
    % Need a lower multiple for alpha shape
    interconeVertices(forcedCapInds, 1:2) = interconeVertices(forcedCapInds, 1:2)*0.75;
    
    interconeVertices(:,1) = interconeVertices(:,1) + volumeSize(1)/2;
    interconeVertices(:,2) = interconeVertices(:,2) + volumeSize(2)/2;
    interconeVertices(:,3) = interconeVertices(:,3) + metaData.tipOffset*metaData.voxSize/metaData.voxSize3D;
    
    % Now get a layer to force base
    [xGrid, yGrid, zGrid] = meshgrid(1:volumeSize(1), 1:volumeSize(2), tempZ(end)+metaData.tipOffset*metaData.voxSize/metaData.voxSize3D);
    xGrid = xGrid(:); yGrid = yGrid(:); zGrid = zGrid(:);
    
    rGrid = sqrt((xGrid-volumeSize(1)/2).^2 + (yGrid-volumeSize(2)/2).^2);
    endInds = find(rGrid > tempProfile(end));
    
    interconeVertices = [interconeVertices', [xGrid(endInds), yGrid(endInds), zGrid(endInds)]']';
    
    % Make another layer at base to force 3D
    [xGrid, yGrid, zGrid] = meshgrid(1:volumeSize(1), 1:volumeSize(2), tempZ(end)+1+metaData.tipOffset*metaData.voxSize/metaData.voxSize3D);
    xGrid = xGrid(:); yGrid = yGrid(:); zGrid = zGrid(:);

    interconeVertices = [interconeVertices', [xGrid, yGrid, zGrid]']';
    forcedGridInds = find(interconeVertices(:,3) == tempZ(end)+1+metaData.tipOffset*metaData.voxSize/metaData.voxSize3D);

%     tempTriangulation = delaunayTriangulation(interconeVertices);
%     [tempFaces, tempVertices] = freeBoundary(tempTriangulation); 
    
    % as exposed cone is concave has to be done with alpha - but leads to some flattening around border
    tempShape = alphaShape(interconeVertices, 10);
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

    if flipSurfaces
        interconeSurface.vertices(:,3) = -(tempVertices(:,3)-volumeSize(3)/2)+volumeSize(3)/2;
    end
    interconeBorderVolume = polygon2voxel(interconeSurface, volumeSize, 'none');

    
    
    % plots 
    figure;
%     plot3(coneVertices(:,1), coneVertices(:,2), coneVertices(:,3), '.')
    hold on; axis equal
    trisurf(grinSurface.faces, grinSurface.vertices(:,1), grinSurface.vertices(:,2), grinSurface.vertices(:,3), ...
       'FaceColor','g','FaceAlpha',0.8);
   
%     plot3(corneaVertices(:,1), corneaVertices(:,2), corneaVertices(:,3), '.')
    trisurf(corneaSurface.faces, corneaSurface.vertices(:,1), corneaSurface.vertices(:,2), corneaSurface.vertices(:,3), ...
       'FaceColor','m','FaceAlpha',0.8);
   
%    plot3(epicorneaVertices(:,1), epicorneaVertices(:,2), epicorneaVertices(:,3), '.')
   trisurf(epicorneaSurface.faces, epicorneaSurface.vertices(:,1), epicorneaSurface.vertices(:,2), epicorneaSurface.vertices(:,3), ...
       'FaceColor','c','FaceAlpha',0.8);
   
%    plot3(cInCVertices(:,1), cInCVertices(:,2), cInCVertices(:,3), '.')
   trisurf(cInCSurface.faces, cInCSurface.vertices(:,1), cInCSurface.vertices(:,2), cInCSurface.vertices(:,3), ...
       'FaceColor','b','FaceAlpha',0.8);
   
%    plot3(interconeVertices(:,1), interconeVertices(:,2), interconeVertices(:,3), '.')
   trisurf(interconeSurface.faces, interconeSurface.vertices(:,1), interconeSurface.vertices(:,2), interconeSurface.vertices(:,3), ...
       'FaceColor','b','FaceAlpha',0.8);
   
   figure;
   subplot(1,2,1)
   tempVolume = coneBorderVolume + corneaBorderVolume + epicorneaBorderVolume + cInCBorderVolume + interconeBorderVolume;
   imshow(permute(tempVolume(round(volumeSize(1)/2),:,:), [2 3 1])');
   
   subplot(1,2,2)
   imshow(permute((lensRIVolume(round(volumeSize(1)/2),:,:)-1.45)/(1.54-1.45), [2 3 1])');

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
    % Connect GRIN vertexes to an exterior RI
    vertexIndexes = sub2ind(volumeSize, round(grinSurface.vertices(:,1)), round(grinSurface.vertices(:,2)), round(grinSurface.vertices(:,3)));

    
    tempRIVolume = lensRIVolume;

    tempRIVolume(RIFlagVolume) = 0;

    tempRIVolume = imdilate(tempRIVolume, strel('Sphere',1));

    vertexExteriorRI = tempRIVolume(vertexIndexes);

    % If any are missing loop until they are caught
    while any(vertexExteriorRI == 0)
        tempRIVolume = imdilate(tempRIVolume, strel('Sphere',1));

        vertexExteriorRI(vertexExteriorRI == 0) = tempRIVolume(vertexIndexes(vertexExteriorRI == 0));
    end
    
    % Adjust vertices
    grinSurface.vertices = grinSurface.vertices*voxelSize;

    corneaSurface.vertices = corneaSurface.vertices*voxelSize;
     
    epicorneaSurface.vertices = epicorneaSurface.vertices*voxelSize;
      
    cInCSurface.vertices = cInCSurface.vertices*voxelSize;
       
    interconeSurface.vertices = interconeSurface.vertices*voxelSize;
    
    % Get normals for all
    grinNormals = meshFaceNormals(grinSurface.vertices, grinSurface.faces);
    
    corneaNormals = meshFaceNormals(corneaSurface.vertices, corneaSurface.faces);
    
    epicorneaNormals = meshFaceNormals(epicorneaSurface.vertices, epicorneaSurface.faces);
    
    cInCNormals = meshFaceNormals(cInCSurface.vertices, cInCSurface.faces);
    
    interconeNormals = meshFaceNormals(interconeSurface.vertices, interconeSurface.faces);
    
    % It seems that some vertices can end up outside of the border volume..
    % Could try add border subscripts as well?
    coneBorderVolume = imdilate(coneBorderVolume, strel('sphere',1));
    
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

zSteps = (1:volumeSize(3))*voxelSize;

if useRealData
    %%% Replace with meshgrid later on
    xStartPoints = 1:xSpacing:volumeSize(1); 
    
    if ~isempty(testOrigin)
        xStartPoints = xStartPoints(testOrigin);
    end
else
    if ~plotReferenceFiber | createLunebergLens
        xStartPoints = 1:xSpacing:volumeSize(1); 

    else
        xStartPoints = volumeSize(1)/2+radius/voxelSize*0.5;
    end

end

rayOrigins = [xStartPoints(:), ones(numel(xStartPoints),1)*volumeSize(2)/2,  ...
    ones(numel(xStartPoints),1)]*voxelSize; 

nOrigins = size(rayOrigins,1);

startRayT = [sin(incidenceAngle/180*pi), 0, cos(incidenceAngle/180*pi)];

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

rayPathArray = zeros(3, volumeSize(3), nOrigins)*NaN;

rayPathLengthArray = zeros(volumeSize(3), nOrigins)*NaN;

rayMap = zeros(volumeSize);

if testPlot
    figure; hold on; axis equal
    view(0, 0)
    if useRealData
        plot3(grinSurface.vertices(plotInds,1), grinSurface.vertices(plotInds,2),...
            grinSurface.vertices(plotInds,3), '.')
    elseif useTestData
        plot3(surfaceX(plotInds)*voxelSize, surfaceY(plotInds)*voxelSize,...
            surfaceZ(plotInds)*voxelSize, '.')
    end
    
%     plot3(0, 0, zInterfaceRI,  'cx') 

    
%     pH = patch(grinSurface);
%     
%     pH.FaceColor = 'g';
%     pH.EdgeColor = 'none';
end

% Turn of warning on singular matrix, otherwise backslash will slow down a lot
warning('off', 'MATLAB:nearlySingularMatrix')

minDeltaS = ones(nOrigins, 1);

firstIntersect = zeros(nOrigins, 3);
finalIntersect = zeros(nOrigins, 3);

finalPathLength = zeros(nOrigins, 1);
finalRayT = zeros(nOrigins, 3);
finalRayTRefract = zeros(nOrigins, 3);

if createGradedFiber
   % Period is for parabolic profile, gaussian period is different and independent of entry.
   %%% This is approx as specific value varies depending on entry RI
   halfPeriodInFiber = 2*pi*n0/beta*0.5;   
   
   numPeriods = floor(fiberLength/halfPeriodInFiber);
   
   periodSpotsArray = zeros(3, numPeriods, nOrigins);
   
end

for iOrigin = 1:nOrigins
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
    
    propogateFinalRay = 1;
    
    numberOfExits = 0;
    
    currentPeriod = 1; % only used for graded fiber

    while go
        
        % If outside, extend ray to grin area
        if ~inGraded
            x0 = rayX;
            t0 = rayT;
        
            % entering, need to find intersect to graded volume
            if useRealData
                
                % X and T combined on input, T needs to point back
                lineDef = [rayX rayT];
                
                % Check intersect against all surfaces
                [intersectPointsGrin, intersectDistanceGrin, intersectFacesGrin] = intersectLineMesh3d(lineDef, grinSurface.vertices, grinSurface.faces);
                
                [intersectPointsCornea, intersectDistanceCornea, intersectFacesCornea] = intersectLineMesh3d(lineDef, corneaSurface.vertices, corneaSurface.faces);
                
                [intersectPointsEpicornea, intersectDistanceEpicornea, intersectFacesEpicornea] = intersectLineMesh3d(lineDef, epicorneaSurface.vertices, epicorneaSurface.faces);
                
                [intersectPointsCinC, intersectDistanceCinC, intersectFacesCinC] = intersectLineMesh3d(lineDef, cInCSurface.vertices, cInCSurface.faces);
                
                [intersectPointsIntercone, intersectDistanceIntercone, intersectFacesIntercone] = intersectLineMesh3d(lineDef, interconeSurface.vertices, interconeSurface.faces);
                
                % put miniminum intersection distances into array - could be better way to do this
                if ~isempty(intersectDistanceGrin)
                    % Remove points with dist less than zero - should modify function to do this...
                    intersectPointsGrin(intersectDistanceGrin <= epsilon, :) = [];
                    intersectFacesGrin(intersectDistanceGrin <= epsilon) = [];
                    intersectDistanceGrin(intersectDistanceGrin <= epsilon) = [];
                end
                if ~isempty(intersectDistanceGrin); minDistArray = min(intersectDistanceGrin); else; minDistArray = Inf; end
                
                if ~isempty(intersectDistanceCornea);  
                    intersectPointsCornea(intersectDistanceCornea <= epsilon, :) = [];
                    intersectFacesCornea(intersectDistanceCornea <= epsilon) = [];
                    intersectDistanceCornea(intersectDistanceCornea <= epsilon) = [];
                end
                if ~isempty(intersectDistanceCornea); minDistArray = [minDistArray min(intersectDistanceCornea)]; else; minDistArray = [minDistArray Inf]; end
                
                if ~isempty(intersectDistanceEpicornea);  
                    intersectPointsEpicornea(intersectDistanceEpicornea <= epsilon, :) = [];
                    intersectFacesEpicornea(intersectDistanceEpicornea <= epsilon) = [];
                    intersectDistanceEpicornea(intersectDistanceEpicornea <= epsilon) = [];
                end
                if ~isempty(intersectDistanceEpicornea); minDistArray = [minDistArray min(intersectDistanceEpicornea)]; else; minDistArray = [minDistArray Inf]; end
                
                if ~isempty(intersectDistanceCinC);  
                    intersectPointsCinC(intersectDistanceCinC <= epsilon, :) = [];
                    intersectFacesCinC(intersectDistanceCinC <= epsilon) = [];
                    intersectDistanceCinC(intersectDistanceCinC <= epsilon) = [];
                end
                if ~isempty(intersectDistanceCinC);   minDistArray = [minDistArray min(intersectDistanceCinC)]; else; minDistArray = [minDistArray Inf]; end
                
                if ~isempty(intersectDistanceIntercone);  
                    intersectPointsIntercone(intersectDistanceIntercone <= epsilon, :) = [];
                    intersectFacesIntercone(intersectDistanceIntercone <= epsilon) = [];
                    intersectDistanceIntercone(intersectDistanceIntercone <= epsilon) = [];
                end
                if ~isempty(intersectDistanceIntercone); minDistArray = [minDistArray min(intersectDistanceIntercone)]; else; minDistArray = [minDistArray Inf]; end
                
                if any(~isinf(minDistArray))
                    [~, minIntersect] = min(minDistArray);

                    switch minIntersect

                        case 1 % Cone is first

                        %%% Check this, could it be simpler    
                            
                        % Sort out intersections
%                         if length(intersectDistance) > 1
%                             backInds = find(intersectDistance < 0);
%                             frontInds = find(intersectDistance > 0);
% 
%                             [~, closestBackInd] = min(abs(intersectDistance(backInds)));
%                             [~, closestFrontInd] = min(intersectDistance(frontInds));
% 
%                             if isempty(backInds) | isempty(frontInds)
%                                 %%% Have a treatment for ray exit
%                                 error('Check treatment - all inds either in front or behind')
%                             end
% 
%                             % Expect back ind to be closer
%                             if abs(intersectDistance(backInds(closestBackInd))) < intersectDistance(frontInds(closestFrontInd))
%                                 nearestInd = backInds(closestBackInd);
% 
%                             else
%                                error('Check treatment - front ind closer') 
%                             end
% 
%                         elseif length(intersectDistance) == 1
%                                 nearestInd = 1;
% 
%                         elseif isempty(intersectDistance)
%                                 error('No intersect')
%                         end
% 
%                         rayX = intersectPoints(nearestInd, :);
% 
%                         faceIndices = intersectFaces(nearestInd);

                            inds2Use = find(intersectDistanceGrin == min(intersectDistanceGrin));

                            rayX = mean(intersectPointsGrin(inds2Use,:),1);

                            faceIndices = intersectFacesGrin(inds2Use);

                            intersectResult = 1;

                        case 2 % cornea is first
                            inds2Use = find(intersectDistanceCornea == min(intersectDistanceCornea));

                            rayX = mean(intersectPointsCornea(inds2Use,:),1);

                            faceIndices = intersectFacesCornea(inds2Use);

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
                    end
                    
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
                    else
                       go = 0; 
                    end
                   
                else
                    go = 0;
                end
            end
            
            if go
                firstIntersect(iOrigin,:) = rayX;

                voxelX = round(rayX/voxelSize); 
                
                % extend ray path up to intersection
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
    
                % Entering into graded region
                if useTestData | (useRealData & minIntersect == 1)
                   
                    if numberOfExits > 0 & blockMultipleExits
                       plot3(rayX(1), rayX(2), rayX(3), 'ro', 'markersize', 8)
                        
                       error('ray is reentering') 
                    end
                    
                    % Get interior RI at entry point
                    [~, rOut] = numerical_dT_dt(rayX, volCoords(RIFlagInds,:), RIFlagInds, lensRIVolume, []);

                    if interfaceRefraction
                        % Get surface normal and exterior RI at entry point
                        if useRealData
                            rIn = vertexExteriorRI(grinSurface.faces(faceIndices,:));
                            rIn = mean(rIn(:));    
                            
                            surfaceNormal = mean(grinNormals(faceIndices,:),1);

                        elseif useTestData
                            if createLunebergLens
                                surfaceNormal = rayX - volumeSize/2*voxelSize;

                                surfaceNormal = surfaceNormal/norm(surfaceNormal);

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

                        % Test from surface tracer
        %                 sqrt((grinNormals(faceIndex(1),1)-rayT(1))^2+(grinNormals(faceIndex(1),2)-rayT(2))^2+...
        %                         (grinNormals(faceIndex(1),3)-rayT(3))^2) < sqrt(2)

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

                            plot3(rayX(1), rayX(2), rayX(3), 'go')
                        else
                            % TIR
                            rayT = rayT-2*dot(surfaceNormal,rayT)*surfaceNormal;

                            inGraded = 0; % ??? 

                            plot3(rayX(1), rayX(2), rayX(3), 'g*')
                        end

                    else
                        inGraded = 1; 

                        plot3(rayX(1), rayX(2), rayX(3), 'gd')
                    end

                    % now multiplied by initial RI
                    rayT = rayT/norm(rayT)*rOut;
                    
                elseif useRealData & minIntersect > 1
                    
                    if interfaceRefraction
                        % do refraction between layers
                        
                        %%% currently sets RI assuming ray comes from base
                            % Could test direction to normal and invert if negative?
                        switch minIntersect
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

                                rOut = metaData.interconeValue;
                                
                                surfaceNormal = mean(cInCNormals(faceIndices,:),1);
                                
                                if limitToConeBase
                                
                                    propogateFinalRay = 0;
                                    
                                    go = 0;
                                end

                            case 5 % intercone 
                                rIn = metaData.interconeValue;

                                rOut = metaData.innerValue;
                                
                                surfaceNormal = mean(interconeNormals(faceIndices,:),1);

                        end
                        
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

                            plot3(rayX(1), rayX(2), rayX(3), 'ko')
                        else
                            % TIR
                            rayT = rayT-2*dot(surfaceNormal,rayT)*surfaceNormal;

                            plot3(rayX(1), rayX(2), rayX(3), 'k*')
                        end
                        
                    else
                        plot3(rayX(1), rayX(2), rayX(3), 'kd')
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

                %%% If all voxels equal in patch can flag to skip calc and maintain rayT

                % only nearest 128 used, but seems to need some extra
                [~, nearInds] = sort(sqrt((testCoords(:,1)-rayX(1)).^2+(testCoords(:,2)-rayX(2)).^2+...
                    (testCoords(:,3)-rayX(3)).^2));

                if length(nearInds) > 200
                    nearInds = nearInds(1:200);
                end

                testInds = testInds(nearInds);

                testCoords = testCoords(nearInds,:);

                %%% If doing direction updating, calculate RI for current direction

                [rayX rayT deltaS*10^6 length(testInds)]
            end
            
            x0 = rayX;
            t0 = rayT;
            
            % Calc is fixed to isotropic RI
            [x, t, deltaS] = ray_interpolation(interpType, 'iso', x0', t0', deltaS, testCoords, ...
                testInds, lensRIVolume, tolerance);
            
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

                    plot3(periodSpotsArray(1,currentPeriod,iOrigin), periodSpotsArray(2,currentPeriod,iOrigin), periodSpotsArray(3,currentPeriod,iOrigin), 'xr');

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

                %%% If all voxels equal in patch can flag to skip calc and maintain rayT

                % Sort forwards for closest 128
                [~, nearInds] = sort(sqrt((testCoords(:,1)-x0(1)).^2+(testCoords(:,2)-x0(2)).^2+...
                    (testCoords(:,3)-x0(3)).^2));

                if length(nearInds) > 200
                    nearInds = nearInds(1:200);
                end

                testInds_forwards = testInds(nearInds);
                testCoords_forwards = testCoords(nearInds,:);

                lambda0 = lambdaFn(rayX, x0);

                if iterativeFinal
                    % Get backwards sorted as well
                    [~, nearInds] = sort(sqrt((testCoords(:,1)-rayX(1)).^2+(testCoords(:,2)-rayX(2)).^2+...
                        (testCoords(:,3)-rayX(3)).^2));

                    if length(nearInds) > 200
                        nearInds = nearInds(1:200);
                    end

                    testInds_backwards = testInds(nearInds);
                    testCoords_backwards = testCoords(nearInds,:);

                    % Test intersection, leaving could just be a rounding point           
                    deltaS0_final = deltaS;

                    its = 1;

                    while deltaS0_final*norm(t0) > epsilon && abs(lambda0-1) > 10^-6 
                        %%% using 10^-3 as test for zero in tests above and below

                        SF = 1; %A from paper - saftey factor

                        % included in loop on paper seems to never get hit
                        if lambda0 < 10^-6
                           rayT = t0;
                           warning('Not sure on this');
                           break
                        end

                        %this loop steps back from overshoot
                        while ~intersectionFn(rayX, x0)
                            if (SF >= 0.1)
                               SF = SF - 0.1; 
                            else
                               SF = 0.1*SF;
                            end

                            deltaS_final = SF*lambda0*deltaS0_final;

                            %%% step ray backward with RK5 with fixed step
                            [x, t] = ray_interpolation('5RKN', 'iso', x0', t0', deltaS_final, testCoords_backwards, ...
                                testInds_backwards, lensRIVolume, tolerance);

                            rayX = x'; rayT = t';
                        end

                        deltaS_final = (1-SF)*lambda0*deltaS0_final;

                        % This loop steps forward to find new contact position just before overshoot.
                        while intersectionFn(rayX, x0)
                            x0 = rayX; %think this is the case, we are stepping forward
                            t0 = rayT;

                            [x, t] = ray_interpolation('5RKN', 'iso', x0', t0', deltaS_final, testCoords_forwards, ...
                                testInds_forwards, lensRIVolume, tolerance);

                            rayX = x'; rayT = t';

                            pathLength = pathLength + sqrt((x0(1) - rayX(1))^2 + (x0(2) - rayX(2))^2 + (x0(3) - rayX(3))^2);
                        end

                        %update parameters for next loop
                        lambda0 = lambdaFn(rayX, x0);

                        deltaS0_final = deltaS_final;
                    end

                else
                    [x, t] = ray_interpolation(interpType, 'iso', x0', t0', deltaS*lambda0, testCoords_forwards, ...
                        testInds_forwards, lensRIVolume, tolerance);

                   rayX = x';
                   rayT = t';

                   pathLength = pathLength + sqrt((x0(1) - rayX(1))^2 + (x0(2) - rayX(2))^2 + (x0(3) - rayX(3))^2);
                end


                % Ray might not exit if reflected, but these will get over written when it does
                finalIntersect(iOrigin,:) = rayX;

                finalPathLength(iOrigin,:) = pathLength;

                finalRayT(iOrigin,:) = rayT;

                % Get surface normal and RI at exit point
                [~, rIn] = numerical_dT_dt(rayX, volCoords(RIFlagInds,:), RIFlagInds, lensRIVolume, []); 
                
                % In test case normalization seems to be very similar to dividing by RI at exit
    %           rayT = rayT/norm(rayT);
                rayT = rayT/rIn;
                
                % Do refraction at border
                if interfaceRefraction

                    %%% Check this, could it be simpler    
                    
                    if useRealData
                        lineDef = [rayX rayT];

                        [intersectPoints, intersectDistance, intersectFaces] = intersectLineMesh3d(lineDef, grinSurface.vertices, grinSurface.faces);

                        if length(intersectDistance) > 1
                            % Sort out intersections
                            backInds = find(intersectDistance < 0);
                            frontInds = find(intersectDistance > 0);

                            if isempty(backInds) | isempty(frontInds)
                                % Inds just on one side
                                [~, sortInds] = sort(abs(intersectDistance));

                                nearestInd = sortInds(1);
                                furtherInd = sortInds(2);

                            elseif ~isempty(backInds) & ~isempty(frontInds)
                                % Inds front and back
                                [~, closestBackInd] = min(abs(intersectDistance(backInds)));
                                [~, closestFrontInd] = min(intersectDistance(frontInds));

                                % swap to nearest and furthest
                                if abs(intersectDistance(backInds(closestBackInd))) < intersectDistance(frontInds(closestFrontInd))
                                    nearestInd = backInds(closestBackInd);
                                    furtherInd = frontInds(closestFrontInd);
                                else
                                    furtherInd = backInds(closestBackInd);
                                    nearestInd = frontInds(closestFrontInd);
                                end
                            end

                            if abs(intersectDistance(nearestInd))*10 > abs(intersectDistance(furtherInd))
                                error('Check treatment - both inds are very close')
                            end
                        elseif length(intersectDistance) == 1
                            nearestInd = 1;

                        elseif isempty(intersectDistance)
                                error('No intersect')
                        end

                        faceIndex = intersectFaces(nearestInd);

                        surfaceNormal = grinNormals(faceIndex,:);

                        rOut = mean(vertexExteriorRI(grinSurface.faces(faceIndex,:)));

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
                       % Assuming all refracted, none reflected
                       rayT = nRatio*rayT + (nRatio*cosI-cosT)*surfaceNormal;

                       plot3(rayX(1), rayX(2), rayX(3), 'mo') 

                       inGraded = 0;

                       numberOfExits = numberOfExits + 1;
                       
                       rayT = rayT/norm(rayT);

                       finalRayTRefract(iOrigin,:) = rayT;
                    else
                        % TIR
                        rayT = rayT-2*dot(surfaceNormal,rayT)*surfaceNormal;

                        rayT = rayT/norm(rayT)*rIn;
                        
                        plot3(rayX(1), rayX(2), rayX(3), 'm*') 

                        % Does not exit be default
                        inGraded = 1;

                    end
                else 
                   plot3(rayX(1), rayX(2), rayX(3), 'mo')

                   inGraded = 0;

                   numberOfExits = numberOfExits + 1;
                end
            end
            
            if rayT(3) < 0
               % kill ray if it reverses
               %%% Should flag for final plotting
               go = 0; 
               
               propogateFinalRay = 0;
               
               plot3(rayX(1), rayX(2), rayX(3), 'rd', 'markersize', 8)
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

                if testPlot
                    plot3(rayX(1), rayX(2), rayX(3), '.', 'color', rayCols(iOrigin,:));
                end
            end
        end
        
        voxelX = round(rayX/voxelSize);
        
%         rayMap(voxelX(1), voxelX(2), voxelX(3)) = 1;
        
        % Check ray has exitied volume
        if any(voxelX > volumeSize) | any(voxelX < 1)
           go = 0; 
        end
    end
    
    % extend ray path to end of z steps
    if propogateFinalRay
        if currentStep < length(zSteps)
            rayT = rayT/rayT(3);
            for i = currentStep+1:length(zSteps)
                rayPathArray(:, i, iOrigin) = rayX + rayT.*(zSteps(i) - rayX(3));
            end

            if testPlot
                plot3(rayPathArray(1, currentStep+1:length(zSteps), iOrigin), rayPathArray(2, currentStep+1:length(zSteps), iOrigin), ...
                    rayPathArray(3, currentStep+1:length(zSteps), iOrigin), '.', 'color', rayCols(iOrigin,:));
            end
        end
    end
    
    pause(0.01)
end

warning('on', 'MATLAB:nearlySingularMatrix')

tolerance
minDeltaS(minDeltaS < 1)*10^6
mean(minDeltaS(minDeltaS < 1))*10^6

%% Plot results
% close all

if useRealData

    figure; %subplot(1,2,1); 
    hold on; axis equal;
    view(0, 0)
    for iOrigin = 1:nOrigins
        rayPath = permute(rayPathArray(:, :, iOrigin), [2 1]);

        plot3(rayPath(:,1), rayPath(:,2), rayPath(:,3), 'color', rayCols(iOrigin,:));
        
%         plot3(firstIntersect(iOrigin,1), firstIntersect(iOrigin,2), firstIntersect(iOrigin,3), 'o', 'color', rayCols(iOrigin,:))
%         plot3(finalIntersect(iOrigin,1), finalIntersect(iOrigin,2), finalIntersect(iOrigin,3), 'o', 'color', rayCols(iOrigin,:))
        
        extendedRay = rayPath(end-1,:) + extendRayLength*finalRayTRefract(iOrigin,:);
        
        line([rayPath(end-1,1) extendedRay(1)], [rayPath(end-1,2) extendedRay(2)], [rayPath(end-1,3) extendedRay(3)], 'color', rayCols(iOrigin,:) )
    end

    xlim([0 volumeSize(1)*voxelSize])
    
    zlim([0 2])
    
    title(sprintf('%i deg',incidenceAngle))
    
    % Plot border
    inds = find(permute(coneBorderVolume(:, round(rayOrigins(1,2)/voxelSize), :),[1 3 2]));
    
    [tempX, tempZ] = ind2sub(volumeSize([1, 3]), inds);
    
    plot3(tempX*voxelSize, ones(length(tempZ),1)*rayOrigins(1,2), tempZ*voxelSize, 'k.');
    
    % Plot ray map through central slice
%     subplot(1,2,2); hold on; axis equal;
%     imshow(flipud( permute( rayMap(:, round(rayOrigins(1,2)/voxelSize), :), [3 1 2])))
    
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
                %%% Can gernalize to X component and off-axis rays from Eq 2 in Babayigit ... Turduev 2019
                
                error('Analytic solution only set up for Y on midline and parallel rays')
            end
                       
            % Check it's within fiber radius
            if abs(rayOrigins(iOrigin,1)-lensCentre(1)) < radius

                %Find first pixel in lens for ray
                startZ = find(RIFlagVolume(round(rayOrigins(iOrigin,1)/voxelSize),...
                    round(rayOrigins(iOrigin,2)/voxelSize), :));

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

                % Put actual start and end points into recorded path
                rayPath(startZPerOrigin(iOrigin),:) = firstIntersect(iOrigin,:);
                rayPath(topZ,:) = finalIntersect(iOrigin,:);
                
                plot((rayPath(:,1) - tempRayPath(:,1))*10^6, zSteps, 'color', rayCols(iOrigin,:))
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
            if finalIntersect(iOrigin,3) ~= 0
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
                [~, rStart] = numerical_dT_dt(firstCoord, volCoords(RIFlagInds,:), RIFlagInds, lensRIVolume, []);  

                tempRayPath(1, bottomZ:topZ) = tempRayPath(1,bottomZ)*cos(tempRayPath(3,bottomZ:topZ)*beta/rStart);
                
                % Analytic expresion for ray after exit
                % Get final Coord and RI
                finalX = tempRayPath(1,bottomZ)*cos(fiberLength*beta/rStart);
                finalCoord = [finalX + fiberCentre(1), fiberCentre(2), fiberCentre(3) + fiberLength/2];
                
                [~, rEnd] = numerical_dT_dt(finalCoord, volCoords(RIFlagInds,:), RIFlagInds, lensRIVolume, []);  

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

                rayPath(bottomZ,:) = firstIntersect(iOrigin,:);
                rayPath(topZ,:) = finalIntersect(iOrigin,:);
                
                tempRayPath = permute(trueRayPathArray(:, :, iOrigin), [2 1]);

                plot((tempRayPath(:,1) - rayPath(:,1))*10^6, zSteps, 'color', rayCols(iOrigin,:))
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
    
   if ~( any(voxelX > size(volumeFull)) | any(voxelX < 1) )
        % Check if on border voxel and fine test required
        if volumeBorder(voxelX(1), voxelX(2), voxelX(3))
            % inpolyhedron is really slow!
%             intersect = inpolyhedron(surface, x);

            % Faster but also causes problems as intersectLineMesh3d tends to shift intersect backwards...
                % Note sure what the above comment meant?
                
            lineDef = [x (x-x0)];  

            %%% Check this, could it be simpler    
            
            [~, intersectDistance] = intersectLineMesh3d(lineDef, surface.vertices, surface.faces);  
            
            backInds = find(intersectDistance < 0);
            frontInds = find(intersectDistance > 0);

            if length(intersectDistance) > 1
                if isempty(backInds) | isempty(frontInds)
                    % all inds on one side, so not in volume
                    intersect = 0;
                elseif ~isempty(backInds) & ~isempty(frontInds)
                   % must be in volume as intersects on both sides
                   intersect = 1;
                else
                   error('Check treatment') 
                end

            elseif length(intersectDistance) == 1
                % get distance from x0
                lineDef_x0 = [x0 (x-x0)];  

                [~, intersectDistance_x0] = intersectLineMesh3d(lineDef_x0, surface.vertices, surface.faces);  

                if intersectDistance_x0 > 0 & intersectDistance < 0
                    % has stepped out
                    intersect = 0;
                elseif intersectDistance_x0 < 0 & intersectDistance > 0
                    % has stepped in ? 
                    intersect = 0;
                else
                   error('Check treatment') 
                end

            elseif isempty(intersectDistance)
               % Not intersect can't be in volume 
               intersect = 0;
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

    lineDef = [x0 (x1-x0)];  

    [intersectPoints, intersectDistance] = intersectLineMesh3d(lineDef, surface.vertices, surface.faces);
    
    % Sort out intersections
    if length(intersectDistance) > 1
        backInds = find(intersectDistance < 0);
        frontInds = find(intersectDistance > 0);

        [~, closestBackInd] = min(abs(intersectDistance(backInds)));
        [~, closestFrontInd] = min(intersectDistance(frontInds));

        if isempty(backInds) | isempty(frontInds)
            % This should be called while x0 is inside mesh and x1 is outside
            error('Check treatment - all inds either in front or behind')
            % It can be that deltaS is too large
        end

        nearestInd = frontInds(closestFrontInd);

    elseif length(intersectDistance) == 1
            nearestInd = 1;

    elseif isempty(intersectDistance)
            error('No intersect')
    end
    
    % Same as distance value in this context
    lambda = norm(intersectPoints(nearestInd,:)-x0)/norm(x1-x0);
    
    if lambda > 1
        lambda;
        error('Large lambda step, intersect was triggered too early')
    elseif lambda < 0
       error('Check this') 
    end
end
