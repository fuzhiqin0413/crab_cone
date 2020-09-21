clear
clc
close all

% Load in angle data
dataFolder = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/Data/liquid angle csv/800nm_thick';

convertAngleConvention = 0;

% Convert to mm
voxelSize = 800/10^6;

scaleSteps = 800/voxelSize/10^6;

% Load phi angles
thetaVolume = loadcsvstack(sprintf('%s/theta', dataFolder), 0);

% Load theta angles
phiVolume = loadcsvstack(sprintf('%s/phi', dataFolder), 0);

if convertAngleConvention
    
    temp = thetaVolume;

    thetaVolume = 90 - phiVolume;

    phiVolume = temp;

    clear temp;
end

volumeSize = size(phiVolume);

% Convert to radians
thetaVolume = thetaVolume/180*pi;

phiVolume = phiVolume/180*pi;

% Get voxel coordinates
goodVolume = ~isnan(phiVolume) & ~isnan(thetaVolume);

goodInds = find(goodVolume);

[tempX, tempY, tempZ] = meshgrid(1:volumeSize(1), 1:volumeSize(2), 1:volumeSize(3));

volInds = sub2ind(volumeSize, tempX(:), tempY(:), tempZ(:));

volCoords = ([tempX(:), tempY(:), tempZ(:)])*voxelSize;

volInds = volInds(goodInds);

volCoords = volCoords(goodInds,:);

%% Plot fibers
% wrap angles to positive - does not make a big difference...
thetaVolume(thetaVolume < 0) = thetaVolume(thetaVolume < 0) + pi;
% 
phiVolume(phiVolume > pi) = phiVolume(phiVolume > pi) - pi;

% Get intersections for plotting
[tempX, tempY, tempZ] = meshgrid(1:1:volumeSize(1), 1:1:volumeSize(2), 1:20:volumeSize(3));

% Get octant
tempAngles = atan2(tempY(:)-volumeSize(2)/2, tempX(:)-volumeSize(1)/2);

tempInds = find(tempAngles > 0 &tempAngles < pi/8);

[~, tempInds] = intersect(volCoords, ...
    ([tempX(tempInds), tempY(tempInds), tempZ(tempInds)])*voxelSize, 'rows');

plotInds = volInds(tempInds);

plotCoords = volCoords(tempInds,:);

inds2Use = 1:length(plotInds);

%tempRadius = sqrt((volCoords(:,1) - volumeSize(1)/2).^2 + (volCoords(:,2) - volumeSize(2)/2).^2);

%figure; axis equal
% plot3(volCoords(:,1)/voxelSize, volCoords(:,2)/voxelSize, volCoords(:,3)/voxelSize, '.')

inds2Use = find(plotCoords(:,3) == min(plotCoords(:,3)))';

[length(plotInds) length(inds2Use)]

figure; 

% subplot(1,3,1);
% 
% imshow((phiVolume(:,:,1))/pi/2); 
% 
% subplot(1,3,2);
% 
% imshow((thetaVolume(:,:,1) + pi)/pi/2)
% 
% subplot(1,3,3); 
hold on; axis equal

set(gca, 'Clipping', 'off'); axis off

for i = inds2Use
    
    coord = plotCoords(i,:);
    
    phi = phiVolume(plotInds(i));
    
    theta = thetaVolume(plotInds(i)); 

    % calculate fiber vector
    vecX = cos(phi)*cos(theta);
    
    vecY = sin(phi)*cos(theta);

    vecZ = sin(theta);
    
    line( [0 voxelSize*vecX]+coord(1), [0 voxelSize*vecY]+coord(2),...
        [0 voxelSize*vecZ]+coord(3) )
end