% Get radius profiles for cones

%%% For new data technique...

clear
clc
close all

voxSize = 2.18;

numCones = 6; % full segmented, extra 4 not included...

% Load in angle data
% dataFolder = '/Users/gavintaylor/Documents/Shared VM Folder/Amira/CT images for full cone profile/Labels';

dataFolder = '/Users/gavintaylor/Documents/Shared VM Folder/Amira/CT images for full cone profile/Labels_w_CinC';

% for original
% coneExposedValue = 5; %%% exposed w/o CinC is 6
% internalInterconeValue = 7;
% interconeExposedValue = 8;
% coneInConeValue = 9; 
% corneaExteriorValue = 10;

% for revised
coneExposedValue = 6; 
coneExposedWOCinCValue = 7;
coneExposedNoCleanValue = 15;
interconeOfConeExposedValue = 14;

internalInterconeValue = 8;
interconeExposedValue = 9;

cInCtoConeValue = 10; 
cInCtoInterconeValue = 11; 

epicorneaInnerValue = 12;
epicorneaOuterValue = 13;

coneData = loadtiffstack(dataFolder, 0);

%%% Note, should add rescaling to all profiles on plots here

volumeSize = size(coneData);

% Seperate cones and CinCs
tempVolume = logical(coneData == coneExposedValue);
coneExposedParam = regionprops3(tempVolume, 'VoxelList', 'VoxelIdxList', 'Volume');

tempVolume = logical(coneData == interconeOfConeExposedValue);
interconeOfConeExposedParam = regionprops3(tempVolume, 'VoxelList', 'VoxelIdxList', 'Volume');

tempVolume = logical(coneData == cInCtoConeValue);
cInCtoConeParam = regionprops3(tempVolume, 'VoxelList', 'VoxelIdxList', 'Volume');

clear tempVolume 
% Note voxel lists from regionprops have flipped x y coords

% Get subs for others
tempInds = find(coneData == interconeExposedValue);
[tempX, tempY, tempZ] = ind2sub(volumeSize, tempInds);
interconeExposedSubs = [tempX, tempY, tempZ];

% Take intersect between exposed cone and exposed intercone
% Order of intersections puts internal intercone between non-cleaned exposed cone and exposed intercone 
    % correct for that by growing non-cleaned exposed cone label into internal intercone
tempVolume = zeros(volumeSize,'logical');  
tempVolume(coneData == coneExposedNoCleanValue) = 1;
tempVolume = imdilate(tempVolume, strel('cube',3));
coneData(coneData == internalInterconeValue & tempVolume) = coneExposedNoCleanValue;

tempVolume(:) = 0;
tempVolume(coneData == coneExposedValue | coneData == coneExposedWOCinCValue | coneData == coneExposedNoCleanValue) = 1;
tempVolume = imdilate(tempVolume, strel('cube',3));

coneRimInds = find((interconeExposedValue == coneData) & tempVolume);
[~, coneRimInds] = intersect(tempInds, coneRimInds);

clear tempVolume

tempInds = find(coneData == internalInterconeValue);
[tempX, tempY, tempZ] = ind2sub(volumeSize, tempInds);
interconeInternalSubs = [tempX, tempY, tempZ];

tempInds = find(coneData == cInCtoInterconeValue);
[tempX, tempY, tempZ] = ind2sub(volumeSize, tempInds); 
cInCtoInterconeSubs = [tempX, tempY, tempZ];

tempInds = find(coneData == epicorneaInnerValue);
[tempX, tempY, tempZ] = ind2sub(volumeSize, tempInds); 
epicorneaInnerSubs = [tempX, tempY, tempZ];

tempInds = find(coneData == epicorneaOuterValue);
[tempX, tempY, tempZ] = ind2sub(volumeSize, tempInds); 
epicorneaOuterSubs = [tempX, tempY, tempZ];

% Tidy up found regions
goodCones = find(coneExposedParam.Volume > 1000);
if length(goodCones) == numCones
    coneExposedParam = coneExposedParam(goodCones,:);
else
   error('Wrong number of cones found') 
end

goodConeInterior = find(interconeOfConeExposedParam.Volume > 1000);
if length(goodConeInterior) == numCones
    interconeOfConeExposedParam = interconeOfConeExposedParam(goodConeInterior,:);
else
   error('Wrong number of intercones found') 
end

goodCinC = find(cInCtoConeParam.Volume > 1000);
if length(goodCinC) == numCones
    cInCtoConeParam = cInCtoConeParam(goodCinC,:);
else
   error('Wrong number of CinC found') 
end


%% check plots 

coneCols = lines(numCones);

fCone = figure; 

% Plot to check match
for i = 1:numCones
    subplot(1,2,1); hold on; axis equal
    tempVoxels = coneExposedParam(i,:).VoxelList{1};
    plot3(tempVoxels(:,2), tempVoxels(:,1), tempVoxels(:,3), '.', 'color', coneCols(i,:))

    tempVoxels = interconeOfConeExposedParam(i,:).VoxelList{1};
    plot3(tempVoxels(:,2), tempVoxels(:,1), tempVoxels(:,3), '.', 'color', coneCols(i,:))

    subplot(1,2,2); hold on; axis equal
    tempVoxels = cInCtoConeParam(i,:).VoxelList{1};
    plot3(tempVoxels(:,2), tempVoxels(:,1), tempVoxels(:,3), '.', 'color', coneCols(i,:))
    
end

display('Check colours correctly match for cones and CinC and intercones')

display('Check vectors for cones and CinC both point forward')

coneCenter = zeros(numCones,3);
coneAxes = zeros(numCones,3);

% Turn off for backslash in ellipsoid fit
warning('off', 'MATLAB:nearlySingularMatrix')

for i = 1:numCones
    % Get cone axis
    subplot(1,2,1); hold on; axis equal

    tempVoxels = [(coneExposedParam(i,:).VoxelList{1})' (interconeOfConeExposedParam(i,:).VoxelList{1})']';

    tempTipDirection = mean(coneExposedParam(i,:).VoxelList{1}) - mean(interconeOfConeExposedParam(i,:).VoxelList{1});
    tempTipDirection = tempTipDirection/norm(tempTipDirection);
    
    coneCenter(i,:) = mean(tempVoxels(:,[2 1 3]));
    
%   PCA axis one for cone because it is long    
%     tempAxes = pca(tempVoxels(:,[2 1 3]) - coneCenter(i,:));
%     line([0 100]*tempAxes(1,1) + coneCenter(i,1), [0 100]*tempAxes(2,1) + coneCenter(i,2),...
%         [0 100]*tempAxes(3,1) + coneCenter(i,3), 'color', coneCols(i,:))
%     coneAxes(i,:) = tempAxes(:,1);

    % Fit ellipsoid for non-closed/tilted cone
    [ center, radii, evecs, v, chi2 ] = ellipsoid_fit( tempVoxels );
    
    % hyperboloid has two negative radii
    if sum(radii < 0) == 2
        hyperboloid = 1;
    elseif sum(radii < 0) == 0   
        hyperboloid = 0;
    else
       error('Check what shape this is') 
    end
    
    [~, useR] = max(radii);
   
    %%% draw fit to test
%     mind = min( tempVoxels );
%     maxd = max( tempVoxels );
%     nsteps = 50;
%     step = ( maxd - mind ) / nsteps;
%     [ x, y, z ] = meshgrid( linspace( mind(1) - step(1), maxd(1) + step(1), nsteps ), linspace( mind(2) - step(2), maxd(2) + step(2), nsteps ), linspace( mind(3) - step(3), maxd(3) + step(3), nsteps ) );
% 
%     Ellipsoid = v(1) *x.*x +   v(2) * y.*y + v(3) * z.*z + ...
%           2*v(4) *x.*y + 2*v(5)*x.*z + 2*v(6) * y.*z + ...
%           2*v(7) *x    + 2*v(8)*y    + 2*v(9) * z;
%     p = patch( isosurface( y, x, z, Ellipsoid, -v(10) ) );
     
    coneAxes(i,:) = evecs([2 1 3],useR);

    % flip angle if needed - should point towards tip
    if sqrt((coneAxes(i,1)-tempTipDirection(1))^2 + (coneAxes(i,2)-tempTipDirection(2))^2 + ...
            (coneAxes(i,3)-tempTipDirection(3))^2) > sqrt(2)
       coneAxes(i,:) = -coneAxes(i,:); 
    end
    
    line([0 300]*coneAxes(i,1) + coneCenter(i,1), [0 300]*coneAxes(i,2) + coneCenter(i,2),...
        [0 300]*coneAxes(i,3) + coneCenter(i,3), 'color', coneCols(i,:))
   
 end

 warning('on', 'MATLAB:nearlySingularMatrix')

%% get tip and rotate for cones

grid3D.nx = volumeSize(1);
grid3D.ny = volumeSize(2);
grid3D.nz = volumeSize(3);
grid3D.minBound = [1 1 1]';
grid3D.maxBound = [volumeSize(1), volumeSize(2), volumeSize(3)]';

depthTests = (-10:1000);
depth0Ind = find(depthTests == 0);

% Storing profile data
% Assume external epicornea and CinC to intercone is flat
coneAverage = zeros(numCones,length(depthTests))*NaN;
coneStandard = zeros(numCones,length(depthTests))*NaN;

cInCAverage = zeros(numCones,length(depthTests))*NaN;
cInCStandard = zeros(numCones,length(depthTests))*NaN;

exposedInterconeAverage = zeros(numCones,length(depthTests))*NaN;
exposedInterconeStandard = zeros(numCones,length(depthTests))*NaN;

interconeRatio = zeros(numCones,length(depthTests))*NaN;

% Note this will probably have step out at bottom - will need to remove
internalInterconeAverage = zeros(numCones,length(depthTests))*NaN;
internalInterconeStandard = zeros(numCones,length(depthTests))*NaN;

epicorneaInnerAverage = zeros(numCones,length(depthTests))*NaN;
epicorneaInnerStandard = zeros(numCones,length(depthTests))*NaN;

embededConesLinked = cell(numCones,2);

% From tip to CinC, epiCornea, cornea.
lengthsToIntersects = zeros(numCones,3);
% From tip to top intercone, intercone, CinC, inner epicornea, outer epicornea
lengthsToPlanes = zeros(numCones,4);

radiusMult = 1.5;
minAngleRange = pi*3/2;
planeSepLimit = 0;

fExamples = figure;
fSummary = figure; hold on; 
fIms = figure;

for i = 1:numCones
    % Get intersect for cone
    % Create vectors
    [xCords, yCords, zCords] = amanatideswooalgorithm_efficient(coneCenter(i,:), coneAxes(i,:), grid3D, ...
        0, [], [], 0); 
    frontInds = sub2ind(volumeSize, xCords, yCords, zCords);
    
    [xCords, yCords, zCords] = amanatideswooalgorithm_efficient(coneCenter(i,:), -coneAxes(i,:), grid3D, ...
        0, [], [], 0); 
    backInds = sub2ind(volumeSize, xCords, yCords, zCords);
    
    % check intersects 
    frontIntersects = find(coneData(frontInds) == coneExposedValue);
    cInCIntersects = find(coneData(backInds) == cInCtoConeValue);
    epicorneaIntersects = find(coneData(backInds) == epicorneaInnerValue);
    corneaIntersects = find(coneData(backInds) == epicorneaOuterValue);
    
    if ~isempty(frontIntersects) 
        tipInd = frontInds(max(frontIntersects));
    else
        error('No exposed cone intersect')
    end
    
    if ~isempty(cInCIntersects)
        coneInConeInd = backInds(min(cInCIntersects));
    else
        error('No cone in cone intersect')
    end
    
    if ~isempty(epicorneaIntersects)
        epicorneaInd = backInds(min(epicorneaIntersects));
    else
        error('No epicornea intersect')
    end
    
    if ~isempty(corneaIntersects)
        corneaInd = backInds(min(corneaIntersects));
    else
        error('No cornea intersect')
    end
    
    % get original cone
    tempConeInds = [(coneExposedParam(i,:).VoxelIdxList{1})' (interconeOfConeExposedParam(i,:).VoxelIdxList{1})']';
    
    tempConeVoxels = [(coneExposedParam(i,:).VoxelList{1})' (interconeOfConeExposedParam(i,:).VoxelList{1})']';
    tempConeVoxels = tempConeVoxels(:, [2 1 3]);
    
    % get tip
    tipIndex = find(tempConeInds == tipInd);
    originalConeTip = tempConeVoxels(tipIndex,:);
    
    % get lengths to each intersection
    [tempX, tempY, tempZ] = ind2sub(volumeSize, coneInConeInd);
    lengthsToIntersects(i,1) = sqrt((originalConeTip(1) - tempX)^2 + (originalConeTip(2) - tempY)^2 + ...
        (originalConeTip(3) - tempZ)^2);
    
    [tempX, tempY, tempZ] = ind2sub(volumeSize, epicorneaInd);
    lengthsToIntersects(i,2) = sqrt((originalConeTip(1) - tempX)^2 + (originalConeTip(2) - tempY)^2 + ...
        (originalConeTip(3) - tempZ)^2);
    
    [tempX, tempY, tempZ] = ind2sub(volumeSize, corneaInd);
    lengthsToIntersects(i,3) = sqrt((originalConeTip(1) - tempX)^2 + (originalConeTip(2) - tempY)^2 + ...
        (originalConeTip(3) - tempZ)^2);
    
    % Get rotation vector
    rotVec = matrix2rotatevectors([0,0,-1], coneAxes(i,:));

    % Rotate cone around tip
    %%% Note, this doesn't guarantee that rotational axes (ie. yaw) is kept constant between cones
    rotatedConeExposed = (tempConeVoxels - originalConeTip)*rotVec;
    
    % Round height and take radius
    coneExposedRadius = sqrt(rotatedConeExposed(:,1).^2 + rotatedConeExposed(:,2).^2);
    coneExposedAngles = atan2(rotatedConeExposed(:,1), rotatedConeExposed(:,2))+pi;
        
    radiusLimit = radiusMult*max(coneExposedRadius);

    % Do for other surfaces.
    rotatedInterconeExposed = (interconeExposedSubs - originalConeTip)*rotVec;
    interconeExposedRadius = sqrt(rotatedInterconeExposed(:,1).^2 + rotatedInterconeExposed(:,2).^2);
    interconeExposedAngles = atan2(rotatedInterconeExposed(:,1), rotatedInterconeExposed(:,2))+pi;
    
    % Merge main internal intercone with borders of labelled cones.
    tempList = 1:numCones; tempList(i) = [];
    rotatedInterconeInternal = interconeInternalSubs;

    for j = tempList
        temp = (interconeOfConeExposedParam(j,:).VoxelList{1});
        rotatedInterconeInternal = [rotatedInterconeInternal' temp(:,[2 1 3])']';
    end
    
    rotatedInterconeInternal = (rotatedInterconeInternal - originalConeTip)*rotVec;
    interconeInternalRadius = sqrt(rotatedInterconeInternal(:,1).^2 + rotatedInterconeInternal(:,2).^2);
    interconeInternalAngles = atan2(rotatedInterconeInternal(:,1), rotatedInterconeInternal(:,2))+pi;

    % CinC is already matched pair
    rotatedCinCToCone = cInCtoConeParam(i,:).VoxelList{1}; 
    rotatedCinCToCone = (rotatedCinCToCone(:, [2 1 3]) - originalConeTip)*rotVec;
    cInCtoConeRadius = sqrt(rotatedCinCToCone(:,1).^2 + rotatedCinCToCone(:,2).^2);
    cInCtoConeAngles = atan2(rotatedCinCToCone(:,1), rotatedCinCToCone(:,2))+pi;

    rotatedCInCtoIntercone = (cInCtoInterconeSubs - originalConeTip)*rotVec;
    cInCtoInterconeRadius = sqrt(rotatedCInCtoIntercone(:,1).^2 + rotatedCInCtoIntercone(:,2).^2);
    cInCtoInterconeAngles = atan2(rotatedCInCtoIntercone(:,1), rotatedCInCtoIntercone(:,2))+pi;
    
    rotatedEpicorneaInner = (epicorneaInnerSubs - originalConeTip)*rotVec;
    epicorneaInnerRadius = sqrt(rotatedEpicorneaInner(:,1).^2 + rotatedEpicorneaInner(:,2).^2);
    epicorneaInnerAngles = atan2(rotatedEpicorneaInner(:,1), rotatedEpicorneaInner(:,2))+pi;
    
    rotatedEpicorneaOuter = (epicorneaOuterSubs - originalConeTip)*rotVec;
    epicorneaOuterRadius = sqrt(rotatedEpicorneaOuter(:,1).^2 + rotatedEpicorneaOuter(:,2).^2);
    epicorneaOuterAngles = atan2(rotatedEpicorneaOuter(:,1), rotatedEpicorneaOuter(:,2))+pi;

    % Remove radius and angles from each that are beyond limit

    tempFlag = zeros(length(interconeExposedRadius),1);
    tempFlag(coneRimInds) = 1;
    tempFlag(interconeExposedRadius > radiusLimit) = [];
    clippedConeRimInds = find(tempFlag);

    interconeExposedAngles(interconeExposedRadius > radiusLimit) = [];
    rotatedInterconeExposed(interconeExposedRadius > radiusLimit, :) = [];
    interconeExposedRadius(interconeExposedRadius > radiusLimit) = [];

    interconeInternalAngles(interconeInternalRadius > radiusLimit) = [];
    rotatedInterconeInternal(interconeInternalRadius > radiusLimit, :) = [];
    interconeInternalRadius(interconeInternalRadius > radiusLimit) = [];
    
    cInCtoInterconeAngles(cInCtoInterconeRadius > radiusLimit) = [];
    rotatedCInCtoIntercone(cInCtoInterconeRadius > radiusLimit, :) = [];
    cInCtoInterconeRadius(cInCtoInterconeRadius > radiusLimit) = [];
    
    epicorneaInnerAngles(epicorneaInnerRadius > radiusLimit) = [];
    rotatedEpicorneaInner(epicorneaInnerRadius > radiusLimit, :) = [];
    epicorneaInnerRadius(epicorneaInnerRadius > radiusLimit) = [];

    epicorneaOuterAngles(epicorneaOuterRadius > radiusLimit) = [];
    rotatedEpicorneaOuter(epicorneaOuterRadius > radiusLimit, :) = [];
    epicorneaOuterRadius(epicorneaOuterRadius > radiusLimit) = [];
    
     % Test plot
    figure(fExamples) 
    subplot(2,numCones,i); hold on; axis equal; view(0,0)
    
    cols = lines(8);

    plot3(rotatedConeExposed(:,1), rotatedConeExposed(:,2), rotatedConeExposed(:,3), '.', 'color', cols(1,:))

    plot3(rotatedCinCToCone(:,1), rotatedCinCToCone(:,2), rotatedCinCToCone(:,3), '.', 'color', cols(2,:))
    plot3(rotatedCInCtoIntercone(:,1), rotatedCInCtoIntercone(:,2), rotatedCInCtoIntercone(:,3), '.', 'color', cols(3,:));

    plot3(rotatedInterconeInternal(:,1), rotatedInterconeInternal(:,2), rotatedInterconeInternal(:,3), '.', 'color', cols(4,:));
    
    plot3(rotatedInterconeExposed(:,1), rotatedInterconeExposed(:,2), rotatedInterconeExposed(:,3), '.', 'color', cols(5,:));
%     plot3(rotatedInterconeExposed(clippedConeRimInds,1), rotatedInterconeExposed(clippedConeRimInds,2), rotatedInterconeExposed(clippedConeRimInds,3), 'rx')

    plot3(rotatedEpicorneaInner(:,1), rotatedEpicorneaInner(:,2), rotatedEpicorneaInner(:,3), '.', 'color', cols(6,:));
    plot3(rotatedEpicorneaOuter(:,1), rotatedEpicorneaOuter(:,2), rotatedEpicorneaOuter(:,3), '.', 'color', cols(7,:));

    line([0 0],[0 0],[0 300], 'color', coneCols(i,:))

    % Fit planes to each surface (was previously doing rotations)
    fExposedInterconeTop = fit( rotatedInterconeExposed(clippedConeRimInds,1:2), rotatedInterconeExposed(clippedConeRimInds,3), 'poly11' );
    tempSubs = [rotatedCInCtoIntercone' rotatedCinCToCone']';
    fCInCBoth = fit( tempSubs(:,1:2), tempSubs(:,3), 'poly11' );
    fEpicorneaInner = fit( rotatedEpicorneaInner(:,1:2), rotatedEpicorneaInner(:,3), 'poly11' );    
    fEpicorneaOuter = fit( rotatedEpicorneaOuter(:,1:2), rotatedEpicorneaOuter(:,3), 'poly11' );

%     plot(fExposedInterconeTop, rotatedInterconeExposed(clippedConeRimInds,1:2), rotatedInterconeExposed(clippedConeRimInds,3));
%     plot(fCInCBoth, tempSubs(:,1:2), tempSubs(:,3));
%     plot(fEpicorneaInner, rotatedEpicorneaInner(:,1:2), rotatedEpicorneaInner(:,3));
%     plot(fEpicorneaOuter, rotatedEpicorneaOuter(:,1:2), rotatedEpicorneaOuter(:,3));

    % Get normal vectors - can take ange to vertical axis as well
        % They should all be roughly co-planer
    normExposedInterconeTop = cross([1, 0, fExposedInterconeTop.p10], [0, 1, fExposedInterconeTop.p01]);
    normCInCBoth = cross([1, 0, fCInCBoth.p10], [0, 1, fCInCBoth.p01]);
    normEpicorneaInner = cross([1, 0, fEpicorneaInner.p10], [0, 1, fEpicorneaInner.p01]);
    normEpicorneaOuter = cross([1, 0, fEpicorneaOuter.p10], [0, 1, fEpicorneaOuter.p01]);

    % Flip direction so they all point down
    if normExposedInterconeTop(3) > 0; normExposedInterconeTop = -normExposedInterconeTop; end
    if normCInCBoth(3) > 0; normCInCBoth = -normCInCBoth; end
    if normEpicorneaInner(3) > 0; normEpicorneaInner = -normEpicorneaInner; end
    if normEpicorneaOuter(3) > 0; normEpicorneaOuter = -normEpicorneaOuter; end

    % Distance from cone tip are just z values at centre
    lengthsToPlanes(i,1) = fExposedInterconeTop.p00;
    lengthsToPlanes(i,2) = fCInCBoth.p00;
    lengthsToPlanes(i,3) = fEpicorneaInner.p00;
    lengthsToPlanes(i,4) = fEpicorneaOuter.p00;

    % rotate each CinC and epicornea inner and outer to parallel around there own centres
    tempRotVec = matrix2rotatevectors([0,0,-1], normCInCBoth);
    rotatedCInCtoIntercone = (rotatedCInCtoIntercone - [0 0 fCInCBoth.p00])*tempRotVec + [0 0 fCInCBoth.p00];
    rotatedCinCToCone = (rotatedCinCToCone - [0 0 fCInCBoth.p00])*tempRotVec + [0 0 fCInCBoth.p00];

    tempRotVec = matrix2rotatevectors([0,0,-1], normEpicorneaInner);
    rotatedEpicorneaInner = (rotatedEpicorneaInner - [0 0 fEpicorneaInner.p00])*tempRotVec + [0 0 fEpicorneaInner.p00];

    tempRotVec = matrix2rotatevectors([0,0,-1], normEpicorneaOuter);
    rotatedEpicorneaOuter = (rotatedEpicorneaOuter - [0 0 fEpicorneaOuter.p00])*tempRotVec + [0 0 fEpicorneaOuter.p00];

    subplot(2,numCones,numCones+i);hold on; axis equal; view(0,0)
    
    plot3(rotatedConeExposed(:,1), rotatedConeExposed(:,2), rotatedConeExposed(:,3), '.', 'color', cols(1,:))

    plot3(rotatedEpicorneaInner(:,1), rotatedEpicorneaInner(:,2), rotatedEpicorneaInner(:,3), '.', 'color', cols(6,:));
    plot3(rotatedEpicorneaOuter(:,1), rotatedEpicorneaOuter(:,2), rotatedEpicorneaOuter(:,3), '.', 'color', cols(7,:));
    plot3(rotatedCinCToCone(:,1), rotatedCinCToCone(:,2), rotatedCinCToCone(:,3), '.', 'color', cols(2,:))
    plot3(rotatedCInCtoIntercone(:,1), rotatedCInCtoIntercone(:,2), rotatedCInCtoIntercone(:,3), '.', 'color', cols(3,:));

    line([0 0],[0 0],[0 300], 'color', coneCols(i,:))

    % rotating intercone is more complicated as I want to rotate around base and then apply a shear in z 
    % average top and bottom vectors
    tempRotVec = matrix2rotatevectors([0,0,-1], (normExposedInterconeTop+normCInCBoth)/2);
    
    % rotate top centre around base centre
    newTopCentre = ([0 0 fExposedInterconeTop.p00] - [0 0 fCInCBoth.p00])*tempRotVec;
    % Then get shear components
    shearX = -newTopCentre(1)/newTopCentre(3);
    shearY = -newTopCentre(2)/newTopCentre(3);

    newTopCentre(1) = newTopCentre(1)+shearX*newTopCentre(3);
    newTopCentre(2) = newTopCentre(2)+shearY*newTopCentre(3);
    newTopCentre = newTopCentre + [0 0 fCInCBoth.p00];

    % rotate main intercone
    rotatedInterconeInternal = (rotatedInterconeInternal - [0 0 fCInCBoth.p00])*tempRotVec;
    % apply shear in z
    rotatedInterconeInternal(:,1) = rotatedInterconeInternal(:,1) + shearX*rotatedInterconeInternal(:,3);
    rotatedInterconeInternal(:,2) = rotatedInterconeInternal(:,2) + shearY*rotatedInterconeInternal(:,3);
    % Shift back
    rotatedInterconeInternal = rotatedInterconeInternal + [0 0 fCInCBoth.p00];

    plot3(rotatedInterconeInternal(:,1), rotatedInterconeInternal(:,2), rotatedInterconeInternal(:,3), '.', 'color', cols(4,:));

    % Try applying to exposed intercone
    rotatedInterconeExposed = (rotatedInterconeExposed - [0 0 fCInCBoth.p00])*tempRotVec;
    % apply shear in z
    rotatedInterconeExposed(:,1) = rotatedInterconeExposed(:,1) + shearX*rotatedInterconeExposed(:,3);
    rotatedInterconeExposed(:,2) = rotatedInterconeExposed(:,2) + shearY*rotatedInterconeExposed(:,3);
    % Shift back
    rotatedInterconeExposed = rotatedInterconeExposed + [0 0 fCInCBoth.p00];
    plot3(rotatedInterconeExposed(:,1), rotatedInterconeExposed(:,2), rotatedInterconeExposed(:,3), '.', 'color', cols(5,:));

    % Clip based on heights to planes...
    cInCtoConeRadius(rotatedCinCToCone(:,3) > fCInCBoth.p00) = [];
    cInCtoConeAngles(rotatedCinCToCone(:,3) > fCInCBoth.p00 ) = [];
    rotatedCinCToCone(rotatedCinCToCone(:,3) > fCInCBoth.p00 , :) = [];

    interconeInternalAngles(rotatedInterconeInternal(:,3) > fCInCBoth.p00 - planeSepLimit) = [];
    interconeInternalRadius(rotatedInterconeInternal(:,3) > fCInCBoth.p00 - planeSepLimit) = [];
    rotatedInterconeInternal(rotatedInterconeInternal(:,3) > fCInCBoth.p00 - planeSepLimit, :) = [];

    epicorneaInnerRadius(rotatedEpicorneaInner(:,3) < fEpicorneaInner.p00) = [];
    epicorneaInnerAngles(rotatedEpicorneaInner(:,3) < fEpicorneaInner.p00) = [];
    rotatedEpicorneaInner(rotatedEpicorneaInner(:,3) < fEpicorneaInner.p00, :) = [];

    % Do z rounding 
    rotatedConeExposed(:,3) = round(rotatedConeExposed(:,3));
    rotatedCinCToCone(:,3) = round(rotatedCinCToCone(:,3));
    rotatedInterconeInternal(:,3) = round(rotatedInterconeInternal(:,3));
    rotatedInterconeExposed(:,3) = round(rotatedInterconeExposed(:,3));
    rotatedEpicorneaInner(:,3) = round(rotatedEpicorneaInner(:,3));
    % These two are not used
%     rotatedCInCtoIntercone(:,3) = round(rotatedCInCtoIntercone(:,3));
%     rotatedEpicorneaOuter(:,3) = round(rotatedEpicorneaOuter(:,3));

    % Split internal intercone points by cone
    % First split to indvidual cones - w/ unkown number it's easiest to do from image rather than clustering...
    tempRange = ceil(max(rotatedInterconeExposed(clippedConeRimInds,1:2)) - min(rotatedInterconeExposed(clippedConeRimInds,1:2))+1);
    tempImage = zeros(tempRange, 'logical');
    tempSubs = round(rotatedInterconeExposed(clippedConeRimInds,1:2) - min(rotatedInterconeExposed(clippedConeRimInds,1:2))+1);
    tempInd = sub2ind(tempRange, tempSubs(:,1), tempSubs(:,2));
    tempImage(tempInd) = 1;
    % Need to check, could merge cones...
    tempImage= imdilate(tempImage, strel('disk',1));
    
    figure(fIms); subplot(2,numCones,i); 
    imshow(tempImage); hold on

    % Get regions
    coneRimsParam = regionprops(tempImage, 'Centroid', 'PixelIdxList', 'PixelList');

    % Clean of small regions, adjust centroid and remove closest to zero
    regionLengths = zeros(length(coneRimsParam),1);
    centroidDist = zeros(length(coneRimsParam),1);
    for j = 1:length(coneRimsParam)
        regionLengths(j) = length(coneRimsParam(j).PixelIdxList);

        coneRimsParam(j).Centroid = coneRimsParam(j).Centroid + min(rotatedInterconeExposed(clippedConeRimInds,1:2)) -1;
        centroidDist(j) = sqrt( sum(coneRimsParam(j).Centroid.^2));
    end
    centroidDist(regionLengths < 20) = [];
    coneRimsParam(regionLengths < 20) = [];
    
    [~, minI] = min(centroidDist);
    coneRimsParam(minI) = [];

    numExtraCones = length(coneRimsParam);
    extraConeCols = winter(numExtraCones);

    set(gca, 'clipping', 'off')

    centroids = zeros(numExtraCones,2);
    for j = 1:numExtraCones
        % Fit circle to border
        fitParam = CircleFitByPratt(coneRimsParam(j).PixelList);

        viscircles(fitParam(1:2),fitParam(3), 'color', extraConeCols(j,:));
        plot(coneRimsParam(j).PixelList(:,1), coneRimsParam(j).PixelList(:,2),'.', 'color', extraConeCols(j,:))
        plot(fitParam(1), fitParam(2), 'x', 'color', extraConeCols(j,:))

        centroids(j,:) = fliplr(fitParam(1:2) + min(rotatedInterconeExposed(clippedConeRimInds,1:2)) -1);
    end

    % Get closest point at each depth level on each cone
    internalInterconeID = zeros(length(depthTests),numExtraCones)*NaN;
    sideProfileSubs = zeros(length(depthTests),numExtraCones,2)*NaN;

    for j = depthTests
        tempInds = find(rotatedInterconeInternal(:,3) == j);

        if ~isempty(tempInds)
            for k = 1:numExtraCones
                % Get distances to line
                dists = abs( (centroids(k,1)-0)*(0-rotatedInterconeInternal(tempInds,2)) - ...
                    (0-rotatedInterconeInternal(tempInds,1))*(centroids(k,2)-0))/sqrt(centroids(k,1)^2+centroids(k,2)^2);
                % Also get angle between pts and centroid
                dirs = acos((centroids(k,1)*rotatedInterconeInternal(tempInds,1) + centroids(k,2)*rotatedInterconeInternal(tempInds,2)) ./ ...
                    (sqrt(centroids(k,1)^2+centroids(k,2)^2).*sqrt(rotatedInterconeInternal(tempInds,1).^2+rotatedInterconeInternal(tempInds,2).^2)) );

                % Get rid of points w/ high angles
                dists(dirs > pi/2) = 100;

                [minD, minI] = min(dists); 

                if minD < 1
                    internalInterconeID(j-min(depthTests)+1,k) = minI;

                    sideProfileSubs(j-min(depthTests)+1,k,:) = rotatedInterconeInternal(tempInds(minI),1:2);

                end
            end
        end
    end
    
    figure(fIms);
    subplot(2,numCones,i+numCones); hold on

    %%% If radius is large can have points on both sides of cone that oscialate
        % Mostly because of height discritization?

    % check for continuity in side profiles
    for j = 1:numExtraCones
        tempSubs = [permute(sideProfileSubs(:,j,:), [1 3 2]), depthTests']; 

        % get distance between each and point above that isn't nan
        points = find(~isnan(tempSubs(:,1)));
        dists = zeros(length(depthTests),1);

        for k = 2:length(points)

            temp = sqrt((tempSubs(points(k),1) - tempSubs(points(k-1),1)).^2 + ...
                (tempSubs(points(k),2) - tempSubs(points(k-1),2)).^2 + (tempSubs(points(k),3) - tempSubs(points(k-1),3)).^2);

            dists(points(k)) = min(temp);
        end

        steps = find(dists > 2.5);
%         plot(dists, '-', 'color', extraConeCols(j,:))
%         plot(steps, dists(steps),'o', 'color', extraConeCols(j,:))

        if ~isempty(steps)
            pointsTOCut =[];
            for k = 1:length(steps)
                % if below center cut all in front
                if steps(k) < mean(points)
                    pointsTOCut = [pointsTOCut 1:steps(k)-1];
                %if above centre cut all behind    
                elseif steps(k) > mean(points)
                    pointsTOCut = [pointsTOCut steps(k):length(dists)];
                end
%                 plot(pointsTOCut, dists(pointsTOCut), 'x', 'color', extraConeCols(j,:))
            end
        end

        plot3(tempSubs(:,1), tempSubs(:,2), tempSubs(:,3), '.', 'color', extraConeCols(j,:));

        plot3(tempSubs(pointsTOCut,1), tempSubs(pointsTOCut,2), tempSubs(pointsTOCut,3), 'rx');

        internalInterconeID(pointsTOCut, j) = NaN;
        sideProfileSubs(pointsTOCut, j, 1) = NaN;
    end

    % plot cleaned sides 
    figure((fExamples))
    subplot(2,numCones,numCones+i);
    for j = 1:numExtraCones
        plot3(sideProfileSubs(:,j,1), sideProfileSubs(:,j,2), depthTests, ...
            'x', 'color', extraConeCols(j,:));
    end


    % Get profile of all features along and past cone
    for j = depthTests
        % for cone
%         tempInds = find(rotatedConeExposed(:,3) == j);
%         
%         if ~isempty(tempInds)
%             tempRadius = coneExposedRadius(tempInds);
%             
%             tempAngles = coneExposedAngles(tempInds);
%             
%             % Check angle range - aim for 3/4 of circle
%             if max(tempAngles) - min(tempAngles) > pi*3/2
% 
%                 coneAverage(i, j-min(depthTests)+1) = mean(tempRadius);
% 
%                 coneStandard(i, j-min(depthTests)+1) = std(tempRadius);
% 
%                 interconeRatio(i, j-min(depthTests)+1) = sum(coneData(tempConeInds(tempInds)) == interconeOfConeExposedValue)/length(tempInds);
%             end
%         end
        
        [coneAverage(i, j-min(depthTests)+1), coneStandard(i, j-min(depthTests)+1), tempInds] = getProfileValue(rotatedConeExposed(:,3), j, ...
                coneExposedRadius, coneExposedAngles, minAngleRange);

        if ~isempty(tempInds)
            interconeRatio(i, j-min(depthTests)+1) = sum(coneData(tempConeInds(tempInds)) == interconeOfConeExposedValue)/length(tempInds);
        end

        [cInCAverage(i, j-min(depthTests)+1), cInCStandard(i, j-min(depthTests)+1)] = getProfileValue(rotatedCinCToCone(:,3), j, ...
                cInCtoConeRadius, cInCtoConeAngles, minAngleRange);

        %%% This will be different as will be gap in values
        [exposedInterconeAverage(i, j-min(depthTests)+1), exposedInterconeStandard(i, j-min(depthTests)+1)] = getProfileValue(rotatedInterconeExposed(:,3), j, ...
                interconeExposedRadius, interconeExposedAngles, minAngleRange);

        % for internal intercone just average closest point on side profiles
        tempInds = find(rotatedInterconeInternal(:,3) == j);
        profileInds = internalInterconeID(j-min(depthTests)+1,:); profileInds(isnan(profileInds)) = [];

        if ~isempty(profileInds)
            internalInterconeAverage(i, j-min(depthTests)+1) = mean( interconeInternalRadius( tempInds( profileInds))); 
            internalInterconeStandard(i, j-min(depthTests)+1) = std( interconeInternalRadius( tempInds( profileInds)));
        end

        [epicorneaInnerAverage(i, j-min(depthTests)+1), epicorneaInnerStandard(i, j-min(depthTests)+1)] = getProfileValue(rotatedEpicorneaInner(:,3), j, ...
                epicorneaInnerRadius, epicorneaInnerAngles, minAngleRange);
    end

    % Adjust offset to narrow point
    [narrowVal] = min(coneAverage(i,1:depth0Ind));
    
    % Add just for duplicates
    narrowInd = find(coneAverage(i,1:depth0Ind) == narrowVal);
    narrowInd = narrowInd(end);
    
    % Shift all back to make start point 
    coneAverage(i, depth0Ind:end) = coneAverage(i, narrowInd:end-(depth0Ind-narrowInd));
    coneStandard(i, depth0Ind:end) = coneStandard(i, narrowInd:end-(depth0Ind-narrowInd));
    cInCAverage(i, depth0Ind:end) = cInCAverage(i, narrowInd:end-(depth0Ind-narrowInd));
    cInCStandard(i, depth0Ind:end) = cInCStandard(i, narrowInd:end-(depth0Ind-narrowInd));
    exposedInterconeAverage(i, depth0Ind:end) = exposedInterconeAverage(i, narrowInd:end-(depth0Ind-narrowInd));
    exposedInterconeStandard(i, depth0Ind:end) = exposedInterconeStandard(i, narrowInd:end-(depth0Ind-narrowInd));
    interconeRatio(i, depth0Ind:end) = interconeRatio(i, narrowInd:end-(depth0Ind-narrowInd));
    internalInterconeAverage(i, depth0Ind:end) = internalInterconeAverage(i, narrowInd:end-(depth0Ind-narrowInd));
    internalInterconeStandard(i, depth0Ind:end) = internalInterconeStandard(i, narrowInd:end-(depth0Ind-narrowInd));
    epicorneaInnerAverage(i, depth0Ind:end) = epicorneaInnerAverage(i, narrowInd:end-(depth0Ind-narrowInd));
    epicorneaInnerStandard(i, depth0Ind:end) = epicorneaInnerStandard(i, narrowInd:end-(depth0Ind-narrowInd));

    % Remove part above narrow point
    coneAverage(i, 1:depth0Ind-1) = NaN;
    coneStandard(i, 1:depth0Ind-1) = NaN;
    cInCAverage(i, 1:depth0Ind-1) = NaN;
    cInCStandard(i, 1:depth0Ind-1) = NaN;
    exposedInterconeAverage(i, 1:depth0Ind-1) = NaN;
    exposedInterconeStandard(i, 1:depth0Ind-1) = NaN;
    interconeRatio(i, 1:depth0Ind-1) = NaN;
    internalInterconeAverage(i, 1:depth0Ind-1) = NaN;
    internalInterconeStandard(i, 1:depth0Ind-1) = NaN;
    epicorneaInnerAverage(i, 1:depth0Ind-1) = NaN;
    epicorneaInnerStandard(i, 1:depth0Ind-1) = NaN;

    % plot indvidual
    figure(fSummary)
    subplot(2,3,1); hold on;
%     plot(depthTests*voxSize, coneAverage(i,:)*voxSize)
    errorbar(depthTests*voxSize, coneAverage(i,:)*voxSize, coneStandard(i,:)*voxSize)

    subplot(2,3,2); hold on;
%     plot(depthTests*voxSize, cInCAverage(i,:)*voxSize)
    errorbar(depthTests*voxSize, cInCAverage(i,:)*voxSize, cInCStandard(i,:)*voxSize)

    subplot(2,3,3); hold on
    plot(depthTests*voxSize, interconeRatio(i,:)*100)

    subplot(2,3,4); hold on;
%     plot(depthTests*voxSize, exposedInterconeAverage(i,:)*voxSize)
    errorbar(depthTests*voxSize, exposedInterconeAverage(i,:)*voxSize, exposedInterconeStandard(i,:)*voxSize)

    subplot(2,3,5); hold on;
%     plot(depthTests*voxSize, internalInterconeAverage(i,:)*voxSize)
    errorbar(depthTests*voxSize, internalInterconeAverage(i,:)*voxSize, internalInterconeStandard(i,:)*voxSize)
   
    subplot(2,3,6); hold on
%     plot(depthTests*voxSize, epicorneaInnerAverage(i,:)*voxSize)
    errorbar(depthTests*voxSize, epicorneaInnerAverage(i,:)*voxSize, epicorneaInnerStandard(i,:)*voxSize)

end

for i = 1:6
    subplot(2,3,i)
    line([1 1]*mean(lengthsToPlanes(:,3))*voxSize, [0 200]); 
    line([1 1]*mean(lengthsToPlanes(:,4))*voxSize, [0 200]); 

    if i ~= 6
        ylim([0 200]); xlim([0 600])
        line([1 1]*mean(lengthsToPlanes(:,1))*voxSize, [0 200]); 
        line([1 1]*mean(lengthsToPlanes(:,2))*voxSize, [0 200]); 
    end
end
coneAverageMean = nanmean(coneAverage);
coneAverageStd = nanstd(coneAverage);

aaExportCone = [depthTests', coneAverage', coneAverageMean', coneAverageStd'];

interconeRatioMean = nanmean(interconeRatio);
interconeRatioStd = nanstd(interconeRatio);

aaExportNumber = [depthTests', interconeRatio', interconeRatioMean', interconeRatioStd'];

%%% Looks like average is making a bit of a mess of this because one cone is a bit wider, but general trend is to inflect in at end. 
    %%% Add radius based normalization as well as height based.

figure(fSummary)

subplot(1,3,1);
plot(depthTests*voxSize, coneAverageMean*voxSize, 'r-', 'linewidth', 2)

plot(depthTests*voxSize, (coneAverageMean+coneAverageStd)*voxSize, 'r:', 'linewidth', 2)
plot(depthTests*voxSize, (coneAverageMean-coneAverageStd)*voxSize, 'r:', 'linewidth', 2)

title('Cone radius profile')
xlabel('Distance from tip')
ylabel('Radius')

subplot(1,3,3);
plot(depthTests, interconeRatioMean, 'r-', 'linewidth', 2)

plot(depthTests, (interconeRatioMean+interconeRatioStd), 'r:', 'linewidth', 2)
plot(depthTests, (interconeRatioMean-interconeRatioStd), 'r:', 'linewidth', 2)
ylim([0 1])

title('Intercone embedded')
xlabel('Distance from tip')
ylabel('% Embeded')

%% Now do for cornea

depthTests = -10:200;
depth0Ind = find(depthTests == 0);

corneaAverage = zeros(numCones,length(depthTests))*NaN;
corneaStandard = zeros(numCones,length(depthTests))*NaN;

for i = 1:numCones
    % Get intersect for cornea
    % Create vectors
    [xCords, yCords, zCords] = amanatideswooalgorithm_efficient(corneaCenter(i,:), corneaAxes(i,:), grid3D, ...
        0, [], [], 0); 
    outInds = sub2ind(volumeSize, xCords, yCords, zCords);
    
    [xCords, yCords, zCords] = amanatideswooalgorithm_efficient(corneaCenter(i,:), -corneaAxes(i,:), grid3D, ...
        0, [], [], 0); 
    reverseInds = sub2ind(volumeSize, xCords, yCords, zCords);
    
    % check intersects 
    outIntersects = find(coneData(outInds) == 4);
    
    reverseIntersects = find(coneData(reverseInds) == 4);
    
    if ~isempty(outIntersects) 
        finalInd = outInds(max(outIntersects));
    elseif ~isempty(reverseIntersects)
        finalInd = reverseInds(max(reverseIntersects));
    end

    tempVoxels = corneaParam(i,:).VoxelList{1};
    tempVoxels = tempVoxels(:, [2 1 3]);
    
    tipIndex = find(corneaParam(i,:).VoxelIdxList{1} == finalInd);
    
    % Get rotation vector
    rotVec = matrix2rotatevectors([0,0,1], corneaAxes(i,:));
    
    % Rotate around center - note, this doesn't guarantee that rotational
    % axis is kept constant between corneas
    rotatedVoxels = (tempVoxels-corneaCenter(i,:))*rotVec;
    
    % Offset Tip to zero
    tipOffset = rotatedVoxels(tipIndex,:);
    rotatedVoxels = rotatedVoxels - tipOffset;
    
%     plot3(rotatedVoxels(:,1), rotatedVoxels(:,2), rotatedVoxels(:,3), 'x', 'color', cols(i,:))
    rotatedVoxels(:,3) = round(rotatedVoxels(:,3));

    corneaNum = zeros(length(depthTests),1);
    
    % Get profile along cornea
    for j = depthTests
        tempInds = find(rotatedVoxels(:,3) == j);
        
        corneaNum(j-min(depthTests)+1) = length(tempInds);
        
        if ~isempty(tempInds)
            tempDists = sqrt(rotatedVoxels(tempInds,1).^2 + rotatedVoxels(tempInds,2).^2);
            
            corneaAverage(i, j-min(depthTests)+1) = mean(tempDists);
            
            corneaStandard(i, j-min(depthTests)+1) = std(tempDists);
        end
    end
    
    figure(fTemp); 
    subplot(1,2,2); hold on
    plot(corneaNum, 'color', cols(i,:));
    
    [~, maxI] = max(corneaNum);
    plot(maxI, corneaNum(maxI), 'x', 'color', cols(i,:));
    
    % Adjust offset to narrow point
    [narrowVal] = min(corneaAverage(i,1:depth0Ind));
    
    % Add just for duplicates
    narrowInd = find(corneaAverage(i,1:depth0Ind) == narrowVal);
    narrowInd = narrowInd(end);
    
    % Shift back to make start point 
    corneaAverage(i, depth0Ind:end) = corneaAverage(i, narrowInd:end-(depth0Ind-narrowInd));
    corneaStandard(i, depth0Ind:end) = corneaStandard(i, narrowInd:end-(depth0Ind-narrowInd));

    % Remove part above narrow point
    corneaAverage(i, 1:depth0Ind-1) = NaN;
    corneaStandard(i, 1:depth0Ind-1) = NaN;

    
    figure(fSummary)
    
    subplot(1,3,2); hold on;
    plot(depthTests, corneaAverage(i,:))
    
%     subplot(1,2,2); hold on
%     plot(depthTests, corneaStandard(i,:))

    figure(fCone)
    
    subplot(1,2,2)
    
    plot3(tempVoxels(tipIndex,1), tempVoxels(tipIndex,2), tempVoxels(tipIndex,3), 'x', 'color', cols(i,:))
end

clear tempVolume

averageCornea = nanmean(corneaAverage);
averageCorneaStd = nanstd(corneaAverage);

aaExportCornea = [depthTests', corneaAverage', averageCornea', averageCorneaStd'];

figure(fSummary)

plot(depthTests, averageCornea, 'r-', 'linewidth', 2)

plot(depthTests, averageCornea+averageCorneaStd, 'r:', 'linewidth', 2)
plot(depthTests, averageCornea-averageCorneaStd, 'r:', 'linewidth', 2)

title('Cornea radius profile')
xlabel('Distance from tip')
ylabel('Radius')