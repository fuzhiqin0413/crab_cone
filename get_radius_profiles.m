% Get radius profiles for cones

%%% For new data technique...

clear
clc
close all

voxSize = 2.18;

numCones = 6; % full segmented, extra 4 not included...

% Load in angle data
dataFolder = '/Users/gavintaylor/Documents/Shared VM Folder/Amira/CT images for full cone profile/Labels';

coneExposedValue = 5; %%% exposed w/o intercone is 6
coneEmbededValue = 7;
interconeExposedValue = 8;
coneInConeValue = 9; 
corneaExteriorValue = 10;

coneData = loadtiffstack(dataFolder, 0);

%%% Note, should add rescaling to all profiles on plots here

%%

volumeSize = size(coneData);

% Seperate cones
%%% Removed closes as they caused problems with intercone overlap
tempVolume = logical(coneData == coneExposedValue);
% tempVolume = imclose(tempVolume, strel('sphere',1));
coneExposedParam = regionprops3(tempVolume, 'VoxelList', 'VoxelIdxList', 'Volume');

tempVolume = logical(coneData == coneEmbededValue);
% tempVolume = imclose(tempVolume, strel('sphere',1));
coneEmbeddedParam = regionprops3(tempVolume, 'VoxelList', 'VoxelIdxList', 'Volume');

tempVolume = logical(coneData == interconeExposedValue);
% tempVolume = imclose(tempVolume, strel('sphere',1));
interconeExposedParam = regionprops3(tempVolume, 'VoxelList', 'VoxelIdxList', 'Volume');

tempVolume = logical(coneData == coneInConeValue);
% tempVolume = imclose(tempVolume, strel('sphere',1));
coneInConeParam = regionprops3(tempVolume, 'VoxelList', 'VoxelIdxList', 'Volume');

tempVolume = logical(coneData == corneaExteriorValue);
% tempVolume = imclose(tempVolume, strel('sphere',1));
corneaExteriorParam = regionprops3(tempVolume, 'VoxelList', 'VoxelIdxList', 'Volume');

% Note voxel lists here have flipped x y coords

clear tempVolume

% Tidy up found regions
goodCones = find(coneExposedParam.Volume > 1000);
if length(goodCones) == numCones
    coneExposedParam = coneExposedParam(goodCones,:);
else
   error('Wrong number of cones found') 
end

goodCinC = find(coneInConeParam.Volume > 1000);
if length(goodCinC) == numCones
    coneInConeParam = coneInConeParam(goodCinC,:);
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

    tempVoxels = coneInConeParam(i,:).VoxelList{1};
    plot3(tempVoxels(:,2), tempVoxels(:,1), tempVoxels(:,3), 'x', 'color', cols(i,:))
end

display('Check colours correctly match for cones and CinC')

display('Check vectors for cones and CinC both point forward')

coneInConeCenter = zeros(numCones,3);
coneInConeAxes = zeros(numCones,3);

coneExposedCenter = zeros(numCones,3);
coneExposedAxes = zeros(numCones,3);

angleDiff = zeros(numCones,1);

% Turn off for backslash in ellipsoid fit
warning('off', 'MATLAB:nearlySingularMatrix')

% Plot to check axes
% subplot(1,2,2);  hold on; axis equal
for i = 1:numCones

    tempVoxels = coneExposedParam(i,:).VoxelList{1};

    coneExposedCenter(i,:) = mean(tempVoxels(:,[2 1 3]));
    
% PCA xis one for cone because it is long    
%     tempAxes = pca(tempVoxels(:,[2 1 3]));
%     line([0 100]*tempAxes(1,1) + coneExposedCenter(i,1), [0 100]*tempAxes(2,1) + coneExposedCenter(i,2),...
%         [0 100]*tempAxes(3,1) + coneExposedCenter(i,3), 'color', cols(i,:))
%     coneExposedAxes(i,:) = tempAxes(:,1);

    % Fit ellipsoid for non-closed cone
    
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
     
%   coneExposedCenter(i,:) = center;
    coneExposedAxes(i,:) = evecs([2 1 3],maxR);

    % PCA axis three for CinC because it is flat
    tempVoxels = coneInConeParam(i,:).VoxelList{1};
    coneInConeCenter(i,:) = mean(tempVoxels(:,[2 1 3]));
    
    tempAxes = pca(tempVoxels(:,[2 1 3]));
    coneInConeAxes(i,:) = -tempAxes(:,3); % flip up as always points down
    
    % flip angle if needed - CinC always points up but cone varies
    if sqrt((coneExposedAxes(i,1)-coneInConeAxes(i,1))^2 + (coneExposedAxes(i,2)-coneInConeAxes(i,2))^2 + ...
            (coneExposedAxes(i,3)-coneInConeAxes(i,3))^2) > sqrt(2)
       coneExposedAxes(i,:) = -coneExposedAxes(i,:); 
    end
    
    line([0 100]*coneExposedAxes(i,1) + coneExposedCenter(i,1), [0 100]*coneExposedAxes(i,2) + coneExposedCenter(i,2),...
        [0 100]*coneExposedAxes(i,3) + coneExposedCenter(i,3), 'color', cols(i,:))
    
    line([0 100]*coneInConeAxes(i,1) + coneInConeCenter(i,1), [0 100]*coneInConeAxes(i,2) + coneInConeCenter(i,2),...
        [0 100]*coneInConeAxes(i,3) + coneInConeCenter(i,3), 'linestyle', ':', 'color', cols(i,:))
    
    angleDiff(i) = acos(dot(coneExposedAxes(i,:), coneInConeAxes(i,:))/(norm(coneExposedAxes(i,:)) * norm(coneInConeAxes(i,:))) )/pi*180;
end

warning('on', 'MATLAB:nearlySingularMatrix')

angleDiff

%% get tip and rotate for cones

grid3D.nx = volumeSize(1);
grid3D.ny = volumeSize(2);
grid3D.nz = volumeSize(3);
grid3D.minBound = [1 1 1]';
grid3D.maxBound = [volumeSize(1), volumeSize(2), volumeSize(3)]';

depthTests = -10:1000;
depth0Ind = find(depthTests == 0);

coneAverage = zeros(numCones,length(depthTests))*NaN;
coneStandard = zeros(numCones,length(depthTests))*NaN;
interconeRatio = zeros(numCones,length(depthTests))*NaN;

coneToCICLength = zeros(numCones,1);

fSummary = figure; hold on; 

for i = 1:numCones
    % Get intersect for cone
    % Create vectors
    [xCords, yCords, zCords] = amanatideswooalgorithm_efficient(coneExposedCenter(i,:), coneExposedAxes(i,:), grid3D, ...
        0, [], [], 0); 
    outInds = sub2ind(volumeSize, xCords, yCords, zCords);
    
    [xCords, yCords, zCords] = amanatideswooalgorithm_efficient(coneExposedCenter(i,:), -coneExposedAxes(i,:), grid3D, ...
        0, [], [], 0); 
    reverseInds = sub2ind(volumeSize, xCords, yCords, zCords);
    
    % check intersects 
    outIntersects = find(coneData(outInds) == coneExposedValue);
    
    reverseIntersects = find(coneData(reverseInds) == coneInConeValue);
    
    if ~isempty(outIntersects) 
        internalInd = outInds(max(outIntersects));
    else
        error('No exposed cone intersect')
    end
    
    if ~isempty(reverseIntersects)
        coneInConeInd = reverseInds(min(reverseIntersects));
    else
        error('No cone in cone intersect')
    end
    
    tempInds = coneExposedParam(i,:).VoxelIdxList{1};
    tempVoxels = coneExposedParam(i,:).VoxelList{1};
    tempVoxels = tempVoxels(:, [2 1 3]);
    
    tipIndex = find(tempInds == internalInd);
    
    %%% WORK FROM HERE %%%
    
    [corneaX, corneaY, corneaZ] = ind2sub(volumeSize, corneaInd);
    
    interConeInds = find(coneData(tempInds) == 5);
    
    interConePoints = zeros(length(tempInds),1);
    interConePoints(interConeInds) = 1;
    
    coneToCICLength(i) = sqrt((tempVoxels(tipIndex,1) - corneaX)^2 + (tempVoxels(tipIndex,2) - corneaY)^2 + ...
        (tempVoxels(tipIndex,3) - corneaZ)^2);
    
    % Get rotation vector
    rotVec = matrix2rotatevectors([0,0,1], coneAxes(i,:));
    
    warning('Shift to tip for rotation center')
    
    % Rotate around center - note, this doesn't guarantee that rotational
    % axis is kept constant between cones
    rotatedVoxels = (tempVoxels-coneCenter(i,:))*rotVec;
    
    % Offset Tip to zero
    tipOffset = rotatedVoxels(tipIndex,:);
    rotatedVoxels = rotatedVoxels - tipOffset;
    
    rotatedVoxels(:,3) = round(rotatedVoxels(:,3));

%     figure; hold on; axis equal
%     plot3(rotatedVoxels(:,1), rotatedVoxels(:,2), rotatedVoxels(:,3), '.', 'color', cols(i,:))
%     plot3(rotatedVoxels(interConeInds,1), rotatedVoxels(interConeInds,2), rotatedVoxels(interConeInds,3), 'x', 'color', 'm')

    coneFullNumber = zeros(length(depthTests),1);

    % Get profile along cone
    for j = depthTests
        tempInds = find(rotatedVoxels(:,3) == j);
        
        if ~isempty(tempInds)
            tempDists = sqrt(rotatedVoxels(tempInds,1).^2 + rotatedVoxels(tempInds,2).^2);
            
            coneAverage(i, j-min(depthTests)+1) = mean(tempDists);
            
            coneStandard(i, j-min(depthTests)+1) = std(tempDists);
            
            interconeRatio(i, j-min(depthTests)+1) = sum(interConePoints(tempInds))/length(tempInds)*100;
            
            coneFullNumber(j-min(depthTests)+1) = sum(interConePoints(tempInds));
        end
    end
    
    figure(fTemp);
    subplot(1,2,1); hold on
    plot(coneFullNumber, 'color', cols(i,:));
    
    [~, maxI] = max(coneFullNumber);
    plot(maxI, coneFullNumber(maxI), 'x', 'color', cols(i,:));
    
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
    
    figure(fCone)
    
    subplot(1,2,2)
    
    plot3(corneaX, corneaY, corneaZ, 'o', 'color', cols(i,:))
    
    
    
    figure(fSummary)
    
    subplot(1,3,1); hold on;
    plot(depthTests, coneAverage(i,:))
    
    subplot(1,3,3); hold on
    plot(depthTests, interconeRatio(i,:))
end

averageCone = nanmean(coneAverage);
averageConeStd = nanstd(coneAverage);

aaExportCone = [depthTests', coneAverage', averageCone', averageConeStd'];

averageinterconeRatio = nanmean(interconeRatio);
averageinterconeRatioStd = nanstd(interconeRatio);

aaExportNumber = [depthTests', interconeRatio', averageinterconeRatio', averageinterconeRatioStd'];

figure(fSummary)

subplot(1,3,1);
plot(depthTests, averageCone, 'r-', 'linewidth', 2)

plot(depthTests, averageCone+averageConeStd, 'r:', 'linewidth', 2)
plot(depthTests, averageCone-averageConeStd, 'r:', 'linewidth', 2)

title('Cone radius profile')
xlabel('Distance from tip')
ylabel('Radius')

subplot(1,3,3);
plot(depthTests, averageinterconeRatio, 'r-', 'linewidth', 2)

plot(depthTests, averageinterconeRatio+averageinterconeRatioStd, 'r:', 'linewidth', 2)
plot(depthTests, averageinterconeRatio-averageinterconeRatioStd, 'r:', 'linewidth', 2)
ylim([0 100])

title('Intercone embedded')
xlabel('Distance from tip')
ylabel('% Embedded')

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