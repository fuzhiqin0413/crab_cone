% Get radius profiles for cones

clear
clc
close all

voxSize = 2;

numCones = 7;

% Load in angle data
dataFolder = '/Users/gavintaylor/Documents/Shared VM Folder/ForProfiles_1_allCone_Borders';

coneData = loadtiffstack(dataFolder, 0);


%%

volumeSize = size(coneData);

% Seperate cones
% Removed closes as they caused problems with intercone overlap
tempVolume = logical(coneData == 5 | coneData == 6);
% tempVolume = imclose(tempVolume, strel('sphere',1));
coneParam = regionprops3(tempVolume, 'VoxelList', 'VoxelIdxList', 'Volume');

tempVolume = logical(coneData == 4);
% tempVolume = imclose(tempVolume, strel('sphere',1));
corneaParam = regionprops3(tempVolume, 'VoxelList', 'VoxelIdxList', 'Volume');

ringVolume = logical(coneData == 7);

% Note voxel lists here have flipped x y coords

clear tempVolume

% Tidy up found regions
goodCones = find(coneParam.Volume > 1000);
if length(goodCones) == 7
    coneParam = coneParam(goodCones,:);
else
   error('Wrong number of cones found') 
end

goodCorneas = find(corneaParam.Volume > 1000);
if length(goodCones) == 7
    corneaParam = corneaParam(goodCorneas,:);
else
   error('Wrong number of corneas found') 
end

%% check plots 

cols = winter(7);

figure; 

% Plot to check match
subplot(1,2,1); hold on; axis equal
for i = 1:7
    
    tempVoxels = coneParam(i,:).VoxelList{1};
    plot3(tempVoxels(:,2), tempVoxels(:,1), tempVoxels(:,3), '.') % , 'color', cols(i,:)

%     tempVoxels = corneaParam(i,:).VoxelList{1};
%     plot3(tempVoxels(:,2), tempVoxels(:,1), tempVoxels(:,3), 'x', 'color', cols(i,:))
end


display('Check colours correctly match for cones and corneas')

coneAxes = zeros(7,3);
corneaAxes = zeros(7,3);

coneCenter = zeros(7,3);
corneaCenter = zeros(7,3);

angleDiff = zeros(7,1);

% Plot to check axes
subplot(1,2,2);  hold on; axis equal
for i = 1:7
    
    % Axis one for cone because it is long
    tempVoxels = coneParam(i,:).VoxelList{1};
    coneCenter(i,:) = mean(tempVoxels(:,[2 1 3]));
    
    tempAxes = pca(tempVoxels(:,[2 1 3]));
    line([0 100]*tempAxes(1,1) + coneCenter(i,1), [0 100]*tempAxes(2,1) + coneCenter(i,2),...
        [0 100]*tempAxes(3,1) + coneCenter(i,3), 'color', cols(i,:))
    coneAxes(i,:) = tempAxes(:,1);
    
    %Axis three for cornea because it is flat
    tempVoxels = corneaParam(i,:).VoxelList{1};
    corneaCenter(i,:) = mean(tempVoxels(:,[2 1 3]));
    
    tempAxes = pca(tempVoxels(:,[2 1 3]));
    line([0 100]*tempAxes(1,3) + corneaCenter(i,1), [0 100]*tempAxes(2,3) + corneaCenter(i,2),...
        [0 100]*tempAxes(3,3) + corneaCenter(i,3), 'color', cols(i,:), 'linestyle',':')
    corneaAxes(i,:) = tempAxes(:,3);
    
    angleDiff(i) = acos(dot(coneAxes(i,:), corneaAxes(i,:))/(norm(corneaAxes(i,:)) * norm(coneAxes(i,:))) )/pi*180;
end

angleDiff


%% get tip and rotate for cones
grid3D.nx = volumeSize(1);
grid3D.ny = volumeSize(2);
grid3D.nz = volumeSize(3);
grid3D.minBound = [1 1 1]';
grid3D.maxBound = [volumeSize(1), volumeSize(2), volumeSize(3)]';

depthTests = -10:1000;
depth0Ind = find(depthTests == 0);

coneAverage = zeros(7,length(depthTests))*NaN;
coneStandard = zeros(7,length(depthTests))*NaN;
coneNumber = zeros(7,length(depthTests))*NaN;

coneLength = zeros(7,1);

figure; hold on; 
for i = 1:7
    % Get intersect for cone
    % Create vectors
    [xCords, yCords, zCords] = amanatideswooalgorithm_efficient(coneCenter(i,:), coneAxes(i,:), grid3D, ...
        0, [], [], 0); 
    outInds = sub2ind(volumeSize, xCords, yCords, zCords);
    
    [xCords, yCords, zCords] = amanatideswooalgorithm_efficient(coneCenter(i,:), -coneAxes(i,:), grid3D, ...
        0, [], [], 0); 
    reverseInds = sub2ind(volumeSize, xCords, yCords, zCords);
    
    % check intersects 
    outIntersects = find(coneData(outInds) == 5 | coneData(outInds) == 6);
    
    reverseIntersects = find(coneData(reverseInds) == 5 | coneData(reverseInds) == 6);
    
    if ~isempty(outIntersects) 
        finalInd = outInds(max(outIntersects));
        
        corneaIntersect = find(coneData(reverseInds) == 4);
        corneaInd = reverseInds(max(corneaIntersect));
        
    elseif ~isempty(reverseIntersects)
        finalInd = reverseInds(max(reverseIntersects));
        
        corneaIntersect = find(coneData(outInds) == 4);
        corneaInd = outInds(max(corneaIntersect));
    end
    
    tempInds = coneParam(i,:).VoxelIdxList{1};
    tempVoxels = coneParam(i,:).VoxelList{1};
    tempVoxels = tempVoxels(:, [2 1 3]);
    
    tipIndex = find(tempInds == finalInd);
    [corneaX, corneaY, corneaZ] = ind2sub(volumeSize, corneaInd);
    
    interConeInds = find(coneData(tempInds) == 5);
    
    interConePoints = zeros(length(tempInds),1);
    interConePoints(interConeInds) = 1;
    
    coneLength(i) = sqrt((tempVoxels(tipIndex,1) - corneaX)^2 + (tempVoxels(tipIndex,2) - corneaY)^2 + ...
        (tempVoxels(tipIndex,3) - corneaZ)^2);
    
    % Get rotation vector
    rotVec = matrix2rotatevectors([0,0,1], coneAxes(i,:));
    
    % Rotate around center - note, this doesn't guarantee that rotational
    % axis is kept constant between cones
    rotatedVoxels = (tempVoxels-coneCenter(i,:))*rotVec;
    
    % Offset Tip to zero
    rotatedVoxels = rotatedVoxels - rotatedVoxels(tipIndex,:);
    
    rotatedVoxels(:,3) = round(rotatedVoxels(:,3));

%     figure; hold on; axis equal
%     plot3(rotatedVoxels(:,1), rotatedVoxels(:,2), rotatedVoxels(:,3), '.', 'color', cols(i,:))
%     plot3(rotatedVoxels(interConeInds,1), rotatedVoxels(interConeInds,2), rotatedVoxels(interConeInds,3), 'x', 'color', 'm')
    
    % Get profile along cone
    for j = depthTests
        tempInds = find(rotatedVoxels(:,3) == j);
        
        if ~isempty(tempInds)
            tempDists = sqrt(rotatedVoxels(tempInds,1).^2 + rotatedVoxels(tempInds,2).^2);
            
            coneAverage(i, j-min(depthTests)+1) = mean(tempDists);
            
            coneStandard(i, j-min(depthTests)+1) = std(tempDists);
            
            coneNumber(i, j-min(depthTests)+1) = sum(interConePoints(tempInds))/length(tempInds)*100;
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
    coneNumber(i, depth0Ind:end) = coneNumber(i, narrowInd:end-(depth0Ind-narrowInd));
    
    % Remove part above narrow point
    coneAverage(i, 1:depth0Ind-1) = NaN;
    coneStandard(i, 1:depth0Ind-1) = NaN;
    coneNumber(i, 1:depth0Ind-1) = NaN;
    
    subplot(1,3,1); hold on;
    plot(depthTests, coneAverage(i,:))
    
    subplot(1,3,3); hold on
    plot(depthTests, coneNumber(i,:))
end

averageCone = nanmean(coneAverage);
averageConeStd = nanstd(coneAverage);

aaExportCone = [depthTests', coneAverage', averageCone', averageConeStd'];

averageConeNumber = nanmean(coneNumber);
averageConeNumberStd = nanstd(coneNumber);

aaExportNumber = [depthTests', coneNumber', averageConeNumber', averageConeNumberStd'];

subplot(1,3,1);
plot(depthTests, averageCone, 'r-', 'linewidth', 2)

plot(depthTests, averageCone+averageConeStd, 'r:', 'linewidth', 2)
plot(depthTests, averageCone-averageConeStd, 'r:', 'linewidth', 2)

title('Cone radius profile')
xlabel('Distance from tip')
ylabel('Radius')

subplot(1,3,3);
plot(depthTests, averageConeNumber, 'r-', 'linewidth', 2)

plot(depthTests, averageConeNumber+averageConeNumberStd, 'r:', 'linewidth', 2)
plot(depthTests, averageConeNumber-averageConeNumberStd, 'r:', 'linewidth', 2)
ylim([0 100])

title('Intercone embedded')
xlabel('Distance from tip')
ylabel('% Embedded')

coneLength

%% Now do for cornea

depthTests = -10:200;
depth0Ind = find(depthTests == 0);

corneaAverage = zeros(7,length(depthTests))*NaN;
corneaStandard = zeros(7,length(depthTests))*NaN;

ringMean = zeros(7,1);
ringStandard = zeros(7,1);

for i = 1:7
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
    rotatedVoxels = rotatedVoxels - rotatedVoxels(tipIndex,:);
    
%     plot3(rotatedVoxels(:,1), rotatedVoxels(:,2), rotatedVoxels(:,3), 'x', 'color', cols(i,:))
    rotatedVoxels(:,3) = round(rotatedVoxels(:,3));

    % Get profile along cornea
    for j = depthTests
        tempInds = find(rotatedVoxels(:,3) == j);
        
        if ~isempty(tempInds)
            tempDists = sqrt(rotatedVoxels(tempInds,1).^2 + rotatedVoxels(tempInds,2).^2);
            
            corneaAverage(i, j-min(depthTests)+1) = mean(tempDists);
            
            corneaStandard(i, j-min(depthTests)+1) = std(tempDists);
        end
    end
    
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
    
    subplot(1,3,2); hold on;
    plot(depthTests, corneaAverage(i,:))
    
%     subplot(1,2,2); hold on
%     plot(depthTests, corneaStandard(i,:))

    % Get ring points associated to this cornea
    tempVolume = zeros(size(coneData), 'logical');
    tempVolume(corneaParam(i,:).VoxelIdxList{1}) = 1;
    tempVolume = imdilate(tempVolume, strel('sphere',1));
    
    ringInds = find(ringVolume & tempVolume);
    ringVoxels = zeros(length(ringInds),3);
    [ringVoxels(:,1), ringVoxels(:,2), ringVoxels(:,3)] = ind2sub(volumeSize, ringInds);
    
    % Transform
    rotatedRing = (ringVoxels-corneaCenter(i,:))*rotVec;
    
    ringRadius = sqrt(rotatedRing(:,1).^2 + rotatedRing(:,2).^2);
    
    ringMean(i) = mean(ringRadius);
    ringStandard(i) = std(ringRadius);
    
%     figure; hold on
%     plot3(rotatedVoxels(:,1), rotatedVoxels(:,2), rotatedVoxels(:,3), 'r.')
%     plot3(rotatedRing(:,1), rotatedRing(:,2), rotatedRing(:,3), 'bx')
    
end

[ringMean, ringStandard]

clear tempVolume

averageCornea = nanmean(corneaAverage);
averageCorneaStd = nanstd(corneaAverage);

aaExportCornea = [depthTests', corneaAverage', averageCornea', averageCorneaStd'];

plot(depthTests, averageCornea, 'r-', 'linewidth', 2)

plot(depthTests, averageCornea+averageCorneaStd, 'r:', 'linewidth', 2)
plot(depthTests, averageCornea-averageCorneaStd, 'r:', 'linewidth', 2)

title('Cornea radius profile')
xlabel('Distance from tip')
ylabel('Radius')