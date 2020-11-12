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

useRealData = 0;
%%% Add parameters here

useTestData = 1;

radius = 1; %mm - used for both lens and fiber

    createLunebergLens = 1;
        
    createCreateGradedFiber = 0;
        fiberLength = 10;

% 0 will search all voxels        
limitToLens = 1;
    %Only used if above is zero
    exteriorRI = 1;

% Options, 1, 4, 5, 45 - could add 3        
interpType = '1';        

initialDeltaS = 10^-4;
        
testPlot = 1;
%% Load or create test data     
if useRealData
    
    %%% If not doing direction updating, calcualte RI for all
    
elseif useTestData
    if createLunebergLens
        volumeSize = ceil(1.5*radius*2/voxelSize*[1 1 1]);
        
        lensRIVolume = createluneberglens(radius/voxelSize, volumeSize);
    
    elseif createCreateGradedFiber
        volumeSize = ceil([radius*2 radius*2 fiberLength]*1.5/voxelSize);
        
        lensRIVolume = createGradedFiber(radius/voxelSize, fiberLength/voxelSize, volumeSize);
    end
    
    goodVolume = ~isnan(lensRIVolume);
    
    lensRIVolume(isnan(lensRIVolume)) = exteriorRI;
end

%% Do general set up
% Search is limited to good in points
if limitToLens
    goodInds = find(goodVolume);
else
    goodInds = find(~isnan(lensRIVolume));
    
    if any(isnan(lensRIVolume(:)))
        error('NaNs in lensRIVolume?')
    end
end

% Get surface voxels
tempVolume = imdilate(goodVolume, strel('Sphere',1)) - goodVolume;

surfaceInds = find(tempVolume);

[surfaceX, surfaceY, surfaceZ] = ind2sub(volumeSize, surfaceInds);

plotInds = 1:length(surfaceX);

% Get voxel coordinates
[tempX, tempY, tempZ] = meshgrid(1:volumeSize(1), 1:volumeSize(2), 1:volumeSize(3));

volInds = sub2ind(volumeSize, tempX(:), tempY(:), tempZ(:));

volCoords = ([tempX(:), tempY(:), tempZ(:)])*voxelSize;

zSteps = (1:volumeSize(3))*voxelSize;

%%% Replace with meshgrid later on
xStartPoints = 1:5:volumeSize(1);

rayOrigins = [xStartPoints(:), ones(numel(xStartPoints),1)*round(volumeSize(2)/2), ...
    ones(numel(xStartPoints),1)]*voxelSize; 

nOrigins = size(rayOrigins,2);

%% Run ray tracer

rayPathArray = zeros(3, volumeSize(3), nOrigins)*NaN;

if testPlot
    figure; hold on; axis equal;
    
    plot3(surfaceX(plotInds)*voxelSize, surfaceY(plotInds)*voxelSize,...
        surfaceZ(plotInds)*voxelSize, '.')
end

for iOrigin = 1:nOrigins
    %%% Could set up to do for series of ray angles
    rayT = [0, 0, 1]; %3x1 : ray will move from bottom to top

    % In phyisical coordiantes
    %%%% Was in voxel coordinates
    rayX = rayOrigins(iOrigin,:);
    
    %%% Note rayX and x (voxelX) swapped!!! 

    inGraded = 0;
    
    %%%% Was ray x in phyiscial coordinates
    voxelX = round(rayOrigins(iOrigin,:)/voxelSize);
    
    lastNearbyX = rayX;
    
    go = 1;

    % Get initial values for ray-tracer to use (it takes closest 128)
    % For unit matrix, radius 4 will take ~268 (*2?) points
    [tempX, tempY, tempZ] = meshgrid((voxelX(1)-4:voxelX(1)+4)', (voxelX(2)-4:voxelX(2)+4)',...
        (voxelX(3)-4:voxelX(3)+4)');

    testCoords= [tempX(:), tempY(:), tempZ(:)];
    
    tempInd = find(testCoords(:,1) < 1 | testCoords(:,2) < 1 | testCoords(:,3) < 1 | ...
            testCoords(:,1) > volumeSize(1) | testCoords(:,2) > volumeSize(2) | testCoords(:,3) > volumeSize(3));

    testCoords(tempInd, :) = [];

    testInds = sub2ind(volumeSize, testCoords(:,1), testCoords(:,2), testCoords(:,3));

    [testInds, tempInds] = intersect(testInds, goodInds);
    
    testCoords = testCoords(tempInds,:);
    
    %%% If doing direction updating, calculate RI for initial direction
    
    deltaS = initialDeltaS;

    currentStep = 1;

    rayPathArray(:, currentStep, iOrigin) = rayX;
    
    while go
        % test if nearest voxel is 
        if goodVolume(voxelX(1), voxelX(2), voxelX(3)) == 0
            %indicate where border touched
            if inGraded & testPlot
                plot3(rayX(1), rayX(2), rayX(3), 'mx')
            end
            
            inGraded = 0;
        else
            if ~inGraded & testPlot
                plot3(rayX(1), rayX(2), rayX(3), 'cx')
            end
            
            inGraded = 1;
        end
    
        %check if x has moved to include new voxels
        if (norm(rayX-lastNearbyX) > voxelSize) & inGraded
            %if so, reselect nearby point

            lastNearbyX = rayX;

            [tempX, tempY, tempZ] = meshgrid((voxelX(1)-4:voxelX(1)+4)', (voxelX(2)-4:voxelX(2)+4)', ...
                (voxelX(3)-4:voxelX(3)+4)');

            testCoords= [tempX(:), tempY(:), tempZ(:)];

            tempInd = find(testCoords(:,1) < 1 | testCoords(:,2) < 1 | testCoords(:,3) < 1 | ...
                    testCoords(:,1) > volumeSize(1) | testCoords(:,2) > volumeSize(2) | testCoords(:,3) > volumeSize(3));
                
            testCoords(tempInd, :) = [];
            
            testInds = sub2ind(volumeSize, testCoords(:,1), testCoords(:,2), testCoords(:,3));
            
            testInds = intersect(testInds, goodInds);
            
            testCoords = testCoords(tempInds,:);
            
            %%% If doing direction updating, calculate RI for current direction
            
            [rayX rayT deltaS*10^3]
        end
        
        %propogate ray
        if ~inGraded
            rayX = rayX + rayT*deltaS; 
        else

            x0 = rayX;
            t0 = rayT;

            % Calc is fixed to isotropic RI
            [x, t, deltaS] = ray_interpolation(interpType, 'iso', x0', t0', deltaS, testCoords*voxelSize, ...
                testInds, lensRIVolume);

            rayX = x';
            
            rayT = t';
        end 
        
        voxelX = round(rayX/voxelSize);
        
        % At each voxel step, record path for plotting later 
        if rayX(3) >= zSteps(currentStep+1)
            if currentStep + 1 < length(zSteps)
                currentStep = currentStep + 1;
                
                rayPathArray(:, currentStep, iOrigin) = rayX;
                
            else
                % Break ray if it has reach end of z steps
                go = 0;
            end
            
            if testPlot
                plot3(rayX(1), rayX(2), rayX(3), '.b');
                
                pause(0.01)
            end
        end
        
        % Check ray has exitied volume
        if any(voxelX > volumeSize) | any(voxelX < 1)
           go = 0; 
        end
    end
end
%% Plot results

% Plot ray path in object
figure; hold on; axis equal;
for iOrigin = 1:nOrigins
    rayPath = permute(rayPathArray(:, :, iOrigin), [2 1]);

    plot(rayPath(:,1), rayPath(:,2));
end

%%% Plot spot diagram
figure; hold on

%%% Plot ray error - if in fiber
figure; hold on

if useTestData
   %%% Calcualte results. 
end