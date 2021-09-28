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
% coneExposedValue = 5; %%% exposed w/o intercone is 6
% internalInterconeValue = 7;
% interconeExposedValue = 8;
% coneInConeValue = 9; 
% corneaExteriorValue = 10;

% for revised
coneExposedValue = 6; %%% exposed w/o intercone is 7
interconeOfConeExposedValue = 14;

% internalInterconeValue = 8;
interconeExposedValue = 9;

cInCtoConeValue = 10; 
cInCtoInterconeValue = 11; 

epicorneaInnerValue = 12;
epicorneaOuterValue = 13;

coneData = loadtiffstack(dataFolder, 0);

%%% Note, should add rescaling to all profiles on plots here

volumeSize = size(coneData);

% Seperate cones
%%% Removed closes as they caused problems with intercone overlap
tempVolume = logical(coneData == coneExposedValue);
% tempVolume = imclose(tempVolume, strel('sphere',1));
coneExposedParam = regionprops3(tempVolume, 'VoxelList', 'VoxelIdxList', 'Volume');

tempVolume = logical(coneData == interconeOfConeExposedValue);
interconeOfConeExposedParam = regionprops3(tempVolume, 'VoxelList', 'VoxelIdxList', 'Volume');

clear tempVolume 
% Note voxel lists from regionprops have flipped x y coords

% Get subs for others
tempInds = find(coneData == interconeExposedValue);
[tempX, tempY, tempZ] = ind2sub(volumeSize, tempInds);
interconeExposedSubs = [tempX, tempY, tempZ];

tempInds = find(coneData == cInCtoConeValue);
[tempX, tempY, tempZ] = ind2sub(volumeSize, tempInds); 
cInCtoConeSubs = [tempX, tempY, tempZ];

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
   error('Wrong number of corneas found') 
end

%% check plots 

cols = lines(numCones);

fCone = figure; 

% Plot to check match
subplot(1,2,1); hold on; axis equal
for i = 1:numCones
    tempVoxels = coneExposedParam(i,:).VoxelList{1};
    plot3(tempVoxels(:,2), tempVoxels(:,1), tempVoxels(:,3), '.', 'color', cols(i,:))

    tempVoxels = interconeOfConeExposedParam(i,:).VoxelList{1};
    plot3(tempVoxels(:,2), tempVoxels(:,1), tempVoxels(:,3), 'x', 'color', cols(i,:))
end

display('Check colours correctly match for cones and CinC')

display('Check vectors for cones and CinC both point forward')

coneCenter = zeros(numCones,3);
coneAxes = zeros(numCones,3);

% angleDiff = zeros(numCones,1);

% Turn off for backslash in ellipsoid fit
warning('off', 'MATLAB:nearlySingularMatrix')

% Plot to check axes
% subplot(1,2,2);  hold on; axis equal
for i = 1:numCones

    tempVoxels = [(coneExposedParam(i,:).VoxelList{1})' (interconeOfConeExposedParam(i,:).VoxelList{1})']';

    tempTipDirection = mean(coneExposedParam(i,:).VoxelList{1}) - mean(interconeOfConeExposedParam(i,:).VoxelList{1});
    tempTipDirection = tempTipDirection/norm(tempTipDirection);
    
    coneCenter(i,:) = mean(tempVoxels(:,[2 1 3]));
    
% PCA axis one for cone because it is long    
%     tempAxes = pca(tempVoxels(:,[2 1 3]) - coneCenter(i,:));
%     line([0 100]*tempAxes(1,1) + coneCenter(i,1), [0 100]*tempAxes(2,1) + coneCenter(i,2),...
%         [0 100]*tempAxes(3,1) + coneCenter(i,3), 'color', cols(i,:))
%     coneExposedAxes(i,:) = tempAxes(:,1);

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
    
    [~, maxR] = max(radii);
   
    %%% draw fit to test
    mind = min( tempVoxels );
    maxd = max( tempVoxels );
    nsteps = 50;
    step = ( maxd - mind ) / nsteps;
    [ x, y, z ] = meshgrid( linspace( mind(1) - step(1), maxd(1) + step(1), nsteps ), linspace( mind(2) - step(2), maxd(2) + step(2), nsteps ), linspace( mind(3) - step(3), maxd(3) + step(3), nsteps ) );

    Ellipsoid = v(1) *x.*x +   v(2) * y.*y + v(3) * z.*z + ...
          2*v(4) *x.*y + 2*v(5)*x.*z + 2*v(6) * y.*z + ...
          2*v(7) *x    + 2*v(8)*y    + 2*v(9) * z;
    p = patch( isosurface( y, x, z, Ellipsoid, -v(10) ) );
     
%   coneCenter(i,:) = center;
    coneAxes(i,:) = evecs([2 1 3],maxR);

    % flip angle if needed - should point towards tip
    if sqrt((coneAxes(i,1)-tempTipDirection(1))^2 + (coneAxes(i,2)-tempTipDirection(2))^2 + ...
            (coneAxes(i,3)-tempTipDirection(3))^2) > sqrt(2)
       coneAxes(i,:) = -coneAxes(i,:); 
    end
    
    line([0 300]*coneAxes(i,1) + coneCenter(i,1), [0 300]*coneAxes(i,2) + coneCenter(i,2),...
        [0 300]*coneAxes(i,3) + coneCenter(i,3), 'color', cols(i,:))
    
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

coneAverage = zeros(numCones,length(depthTests))*NaN;
coneStandard = zeros(numCones,length(depthTests))*NaN;
interconeRatio = zeros(numCones,length(depthTests))*NaN;

embededConesLinked = cell(numCones,2);

coneToCICLength = zeros(numCones,1);

fSummary = figure; hold on; 

fTest = figure; hold on

for i = 1:numCones
    % Get intersect for cone
    % Create vectors
    [xCords, yCords, zCords] = amanatideswooalgorithm_efficient(coneCenter(i,:), coneAxes(i,:), grid3D, ...
        0, [], [], 0); 
    outInds = sub2ind(volumeSize, xCords, yCords, zCords);
    
    [xCords, yCords, zCords] = amanatideswooalgorithm_efficient(coneCenter(i,:), -coneAxes(i,:), grid3D, ...
        0, [], [], 0); 
    reverseInds = sub2ind(volumeSize, xCords, yCords, zCords);
    
    % check intersects 
    outIntersects = find(coneData(outInds) == coneExposedValue);
    
    reverseIntersects = find(coneData(reverseInds) == cInCtoConeValue);
    
    if ~isempty(outIntersects) 
        tipInd = outInds(max(outIntersects));
    else
        error('No exposed cone intersect')
    end
    
    if ~isempty(reverseIntersects)
        coneInConeInd = reverseInds(min(reverseIntersects));
    else
        error('No cone in cone intersect')
    end
    
    tempConeInds = [(coneExposedParam(i,:).VoxelIdxList{1})' (interconeOfConeExposedParam(i,:).VoxelIdxList{1})']';
    
    tempConeVoxels = [(coneExposedParam(i,:).VoxelList{1})' (interconeOfConeExposedParam(i,:).VoxelList{1})']';
    tempConeVoxels = tempConeVoxels(:, [2 1 3]);
    
    tipIndex = find(tempConeInds == tipInd);
    [CinCX, CinCY, CinCZ] = ind2sub(volumeSize, coneInConeInd);
    
    originalConeTip = tempConeVoxels(tipIndex,:);
    
    coneToCICLength(i) = sqrt((originalConeTip(1) - CinCX)^2 + (originalConeTip(2) - CinCY)^2 + ...
        (originalConeTip(3) - CinCZ)^2);
    
    % Get rotation vector
    rotVec = matrix2rotatevectors([0,0,-1], coneAxes(i,:));
    
    % Rotate cone. 
    %%% Note, this doesn't guarantee that rotational axes (ie. yaw) is kept constant between cones
    
    % Rotate around center 
%     rotatedVoxels = (tempConeVoxels-coneCenter(i,:))*rotVec;

    % Rotate around tip
    rotatedConeExposed = (tempConeVoxels - originalConeTip)*rotVec;
    
    % Removing unnecersary
%     rotatedInternalIntercone = (internalInterconeSubs - originalConeTip)*rotVec;
%     rotatedInterconeExposed = (interconeExposedSubs - originalConeTip)*rotVec;
    
    % Offset Tip to zero
    %%% No longer neccersary
%     tipOffset = rotatedConeExposed(tipIndex,:)
%     rotatedConeExposed = rotatedConeExposed - tipOffset;

    rotatedConeExposed(:,3) = round(rotatedConeExposed(:,3));
    
%     rotatedInternalIntercone(:,3) = round(rotatedInternalIntercone(:,3));
%     rotatedInterconeExposed(:,3) = round(rotatedInterconeExposed(:,3));
    
    coneExposedRadius = sqrt(rotatedConeExposed(:,1).^2 + rotatedConeExposed(:,2).^2);
    
%     internalInterconeRadius = sqrt(rotatedInternalIntercone(:,1).^2 + rotatedInternalIntercone(:,2).^2);
%     interconeExposedRadius = sqrt(rotatedInterconeExposed(:,1).^2 + rotatedInterconeExposed(:,2).^2);
    
    % Get profile along cone
    coneExposedAngles = atan2(rotatedConeExposed(:,1), rotatedConeExposed(:,2))+pi;
%     embededConeAngles = atan2(rotatedEmbededCone(:,1), rotatedEmbededCone(:,2))+pi;
    
    for j = depthTests
        tempInds = find(rotatedConeExposed(:,3) == j);
        
        if ~isempty(tempInds)
            tempRadius = coneExposedRadius(tempInds);
            
            tempAngles = coneExposedAngles(tempInds);
        else
            tempRadius = [];
            tempAngles = [];
        end
        
        if ~isempty(tempRadius)
            
            % Check angle range aim for 3/4 of circle
            if max(tempAngles) - min(tempAngles) > pi*3/2

                coneAverage(i, j-min(depthTests)+1) = mean(tempRadius);

                coneStandard(i, j-min(depthTests)+1) = std(tempRadius);

                interconeRatio(i, j-min(depthTests)+1) = sum(coneData(tempConeInds(tempInds)) == interconeOfConeExposedValue)/length(tempInds);
            end
        end
    end

    % Adjust offset to narrow point
    [narrowVal] = min(coneAverage(i,1:depth0Ind));
    
    % Add just for duplicates
    narrowInd = find(coneAverage(i,1:depth0Ind) == narrowVal);
    narrowInd = narrowInd(end);
    
    % Shift back to make start point 
    coneAverage(i, depth0Ind:end) = coneAverage(i, narrowInd:end-(depth0Ind-narrowInd));
    coneStandard(i, depth0Ind:end) = coneStandard(i, narrowInd:end-(depth0Ind-narrowInd));
    interconeRatio(i, depth0Ind:end) = interconeRatio(i, narrowInd:end-(depth0Ind-narrowInd));
    
    % Remove part above narrow point
    coneAverage(i, 1:depth0Ind-1) = NaN;
    coneStandard(i, 1:depth0Ind-1) = NaN;
    interconeRatio(i, 1:depth0Ind-1) = NaN;
    
    figure(fSummary)
    
    subplot(1,3,1); hold on;
    plot(depthTests*voxSize, coneAverage(i,:)*voxSize)
    
    subplot(1,3,3); hold on
    plot(depthTests, interconeRatio(i,:))
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