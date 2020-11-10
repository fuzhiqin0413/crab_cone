% Set up to do both tracing on data and have testing mode on either
% lundaberg lens or graded fibre, as in nishdate 2011 papers
    
    % Need to determine correct exterior point treatment for boundary condition
        % Could be none, less or equal
 
    % Also need to add interface refraction - not required for test cases but
    % for normal tracing. Should also test this for my patch method.
        % This also requires stepping towards the border
    
    % Should test on smoothed RI map calculated from 200 nm maps
        % Have exterior RI loading and processing function
        % May also want to smooth

    % Start with fixed on axis RI but later update to recalc with angle dependence
    
    % Note grin is basically only in trasition area, so can skip calc if
    % all values in region are equal
        % Should test this speed up

clear; clc; close all

%% Set parameters
voxelSize = 0.065; % mm

interpType = '5';

useRealData = 0;
%%% Add parameters here

useTestData = 1;

radius = 1; %mm - used for both lens and fiber

    createLunebergLens = 1;
        
    createCreateGradedFiber = 0;
        fiberLength = 10;

%% Load or create test data     
if useRealData
    
elseif useTestData
    if createLunebergLens
        volumeSize = 1.5*radius/voxelSize*[1 1];
        
        lensRIVolume = createluneberglens(radius/voxelSize, volumeSize);
    
    elseif createCreateGradedFiber
        volumeSize = [radius radius fiberLength]*1.5/voxelSize;
        
        lensRIVolume = createGradedFiber(radius/voxelSize, fiberLength/voxelSize, volumeSize);
    end
    
    goodVolume = ~isnan(lensRIVolume);
    
    lensRIVolume(isnan(lensRIVolume)) = 1;
end

%% Do general set up
goodInds = find(goodVolume);

% Get surface voxels
tempVolume = imdilate(goodVolume, strel('Sphere',1)) - goodVolume;

surfaceInds = find(tempVolume);

[surfaceX, surfaceY, surfaceZ] = ind2sub(volumeSize, surfaceInds);

% Get voxel coordinates
[tempX, tempY, tempZ] = meshgrid(1:volumeSize(1), 1:volumeSize(2), 1:volumeSize(3));

volInds = sub2ind(volumeSize, tempX(:), tempY(:), tempZ(:));

volCoords = ([tempX(:), tempY(:), tempZ(:)])*voxelSize;

zSteps = (1:volumeSize(3))*voxelSize;

%%% Replace with meshgrid later on
xStartPoints = 1:10*scaleSteps:volumeSize(1);

rayOrigins = [xStartPoints(:), ones(numel(xSteps),1)*round(volumeSize(2)/2), ...
    ones(numel(xSteps),1)]*voxelSize; 

%% Run ray tracer


%% Plot results

if useTestData
   %%% Calcualte results. 
end