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

    %%% To do after debug:
        % Revise above
        % Do step in to get entry and exit

clear; 
clc; close all

%% Set parameters
voxelSize = 0.065; % mm

useRealData = 0;
%%% Add parameters here

useTestData = 1;

radius = 1; %mm - used for both lens and fiber

    createLunebergLens = 0;
        
    createCreateGradedFiber = 1;
        fiberLength = 10;
        n0 = sqrt(2.5);
        alpha = 1/sqrt(2.5);
        
% 0 will search all voxels        
limitToLens = 1;
    %Only used if above is zero
    exteriorRI = 1;

% Options, 1, 4(S B F), 5(B F), 45       
interpType = '45';        

% Manual testing, a bit crap...
iRayTest = 3;

testResults = [];

if iRayTest == 1
    testTypes = cell(3,1);
    testResults = cell(3,1);
end

initialDeltaS = 10^-3;
        
testPlot = 1;
%% Load or create test data     
if useRealData
    
    %%% If not doing direction updating, calcualte RI for all
    
elseif useTestData
    if ~xor(createLunebergLens, createCreateGradedFiber)
        error('Either or neither lens/fiber flags set. Check.')
    end
    
    if createLunebergLens
        volumeSize = ceil(radius*2*1.5/voxelSize*[1 1 1]);
        
        lensRIVolume = createluneberglens(radius/voxelSize, volumeSize);
    
    elseif createCreateGradedFiber
        volumeSize = ceil([radius*2 radius*2 fiberLength]*1.5/voxelSize);
        
        lensRIVolume = createGradedFiber(radius/voxelSize, fiberLength/voxelSize, volumeSize, ...
            n0, alpha, voxelSize);
    end
    
    goodVolume = ~isnan(lensRIVolume);
    
    lensRIVolume(isnan(lensRIVolume)) = exteriorRI;

    % Get surface voxels
    tempVolume = imdilate(~goodVolume, strel('Sphere',1)) & goodVolume;

    surfaceInds = find(tempVolume);

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
    
    centreLine = find(goodVolume(:,round(volumeSize(2)/2),bottomZ));
    tempRadius = sqrt((centreLine - volumeSize(1)/2).^2 + ...
        (round(volumeSize(2)/2) - volumeSize(2)/2)^2 );
    
    plot(centreLine, sqrt(2.5-(tempRadius*voxelSize/radius).^2), 'rx');
    
    plot(centreLine, sqrt(n0^2*(1-alpha^2*(tempRadius*voxelSize).^2)), 'g');
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

plotInds = 1:length(surfaceX);

% Get voxel coordinates
[tempX, tempY, tempZ] = meshgrid(1:volumeSize(1), 1:volumeSize(2), 1:volumeSize(3));

volInds = sub2ind(volumeSize, tempX(:), tempY(:), tempZ(:));

volCoords = ([tempX(:), tempY(:), tempZ(:)])*voxelSize;

zSteps = (1:volumeSize(3))*voxelSize;

%%% Replace with meshgrid later on
xStartPoints = volumeSize(1)*0.25; %1:5:volumeSize(1);

rayOrigins = [xStartPoints(:), ones(numel(xStartPoints),1)*volumeSize(2)/2,  ...
    ones(numel(xStartPoints),1)]*voxelSize; 

nOrigins = size(rayOrigins,1);

startRayT = [0, 0, 1];

rayCols = lines(nOrigins);

%% Run ray tracer

rayPathArray = zeros(3, volumeSize(3), nOrigins)*NaN;

enterancePoints = zeros(3, nOrigins)*NaN;

exitPoints = zeros(3, nOrigins)*NaN;

if testPlot
    figure; hold on; axis equal;
    view(0, 0)
    plot3(surfaceX(plotInds)*voxelSize, surfaceY(plotInds)*voxelSize,...
        surfaceZ(plotInds)*voxelSize, '.')
end

for iOrigin = 1% 1:nOrigins
    %%% Could set up to do for series of ray angles
    rayT = startRayT; %3x1 : ray will move from bottom to top

    % In phyisical coordiantes
    %%%% Was in voxel coordinates
    rayX = rayOrigins(iOrigin,:);
    
    inGraded = 0;
    
    %%%% Was ray x in phyiscial coordinates
    voxelX = round(rayOrigins(iOrigin,:)/voxelSize);
    
    lastNearbyX = [rayX(1), rayX(2), -1];
    
    go = 1;

    %%% just do in main loop, remove from here.
    
    % Get initial values for ray-tracer to use (it takes closest 128)
    % For unit matrix, radius 4 will take ~26 8 (*2?) points
%     [tempX, tempY, tempZ] = meshgrid((voxelX(1)-4:voxelX(1)+4)', (voxelX(2)-4:voxelX(2)+4)',...
%         (voxelX(3)-4:voxelX(3)+4)');
% 
%     testCoords= [tempX(:), tempY(:), tempZ(:)];
%     
%     tempInd = find(testCoords(:,1) < 1 | testCoords(:,2) < 1 | testCoords(:,3) < 1 | ...
%             testCoords(:,1) > volumeSize(1) | testCoords(:,2) > volumeSize(2) | testCoords(:,3) > volumeSize(3));
% 
%     testCoords(tempInd, :) = [];
% 
%     testInds = sub2ind(volumeSize, testCoords(:,1), testCoords(:,2), testCoords(:,3));
% 
%     [testInds, tempInds] = intersect(testInds, goodInds);
%     
%     testCoords = testCoords(tempInds,:);
    
    %%% If doing direction updating, calculate RI for initial direction
    
    deltaS = initialDeltaS;

    currentStep = 1;

    rayPathArray(:, currentStep, iOrigin) = rayX;
    
    while go
        % test if nearest voxel is graded or not
        if goodVolume(voxelX(1), voxelX(2), voxelX(3)) == 0
            %indicate if border transition
            if inGraded & testPlot
                plot3(rayX(1), rayX(2), rayX(3), 'mx')
                
                %%% Set up better intersect
                exitPoints(:, iOrigin) = rayX;
            end
            
            inGraded = 0;
        else
            if ~inGraded & testPlot
                plot3(rayX(1), rayX(2), rayX(3), 'cx')
                
                %Set up better intersect
                enterancePoints(:, iOrigin) = rayX;
            end
            
            inGraded = 1;
        end
    
        %%% Need to add step out at interafece and refraction in and out
        
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
            
            [testInds, tempInds] = intersect(testInds, goodInds);
            
            testCoords = testCoords(tempInds,:);
            
%             figure; plot3(testCoords(:,1), testCoords(:,2), testCoords(:,3), '.')
%             hold on; plot3(voxelX(1), voxelX(2), voxelX(3), 'rx')
            
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

%             [x, t, deltaS] = ray_interpolation(interpType, 'iso', x0', t0', deltaS, volCoords(goodInds,:), ...
%                 goodInds, lensRIVolume);
            
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
                plot3(rayX(1), rayX(2), rayX(3), '.', 'color', rayCols(iOrigin,:));
                
            end
        end
        
        % Check ray has exitied volume
        if any(voxelX > volumeSize) | any(voxelX < 1)
           go = 0; 
        end
    end

    pause(0.01)
end

if ~isempty(testResults)
    testResults{iRayTest} = rayPathArray;
    testTypes{iRayTest} = interpType
end
%% Plot results
% close all

if ~useTestData

else
    
    if createLunebergLens
        % Plot ray path in object
        figure; subplot(2,2,1); hold on; axis equal;
        view(0, 0)
        for iOrigin = 1:nOrigins
            rayPath = permute(rayPathArray(:, :, iOrigin), [2 1]);

            plot3(rayPath(:,1), rayPath(:,2), rayPath(:,3), 'color', rayCols(iOrigin,:));
        end
        %%% Consider replacing lines with circle drawing...
        
        line([1 volumeSize(1)]*voxelSize, [1 1]*volumeSize(2)/2*voxelSize, ...
            [1 1]*zSteps(bottomZ), 'color', 'b')
        
        line([1 volumeSize(1)]*voxelSize, [1 1]*volumeSize(2)/2*voxelSize, ...
            [1 1]*zSteps(topZ), 'color', 'b')
        
        title('Ray path')
        
        % Plot analytic ray path - for on axis rays along central y axis
        subplot(2,2,2); hold on; axis equal;
        view(0, 0)
        trueRayPathArray = zeros(3, volumeSize(3), nOrigins)*NaN;
        
        lensCentre = volumeSize/2*voxelSize;
        
        for iOrigin = 1:nOrigins
            tempRayPath = zeros(3, volumeSize(3))*NaN;
            
            if rayOrigins(iOrigin,2) ~= lensCentre(2) | ... 
                    any(startRayT - [0, 0, 1])
                %%% Can gernalize to X component and off-axis rays from
                %%% Eq 2 in Babayigit ... Turduev 2019
                
                error('Analytic solution only set up for Y on midline and parallel rays')
            end
            
            % Check it's within fiber radius
            if abs(rayOrigins(iOrigin,1)-lensCentre(1)) < radius

                %%% Later, use exact points intersection point for start and end
                
                %Find first pixel in lens for ray
                startZ = find(goodVolume(round(rayOrigins(iOrigin,1)/voxelSize),...
                    round(rayOrigins(iOrigin,2)/voxelSize), :));

                startZ = startZ(1);
                
                % Just a line up to enterance
                tempRayPath(1,1:startZ-1) = rayOrigins(iOrigin,1)-lensCentre(1); % offset to fiber centre
                tempRayPath(2,:) = rayOrigins(iOrigin,2);
                tempRayPath(3,:) = zSteps - lensCentre(3);            

                %%% Doesn't account for any y offset in ray as fiber does
                tempRayPath(1, startZ:topZ) = tempRayPath(1,1)*(tempRayPath(3,startZ)*tempRayPath(3,startZ:topZ) + ...
                    radius*sqrt(radius^2 + tempRayPath(3,startZ)^2 - tempRayPath(3,startZ:topZ).^2))/...
                    (tempRayPath(3,startZ)^2 + radius^2);

                %%% Doesn't plot ray after exit
                tempRayPath(:, topZ+1:end) = NaN;
                
                % Add fiber centre back for X
                tempRayPath(1,:) = tempRayPath(1,:) + lensCentre(1);
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
        % plot lines at top and bottom
        line([1 volumeSize(1)]*voxelSize, [1 1]*volumeSize(2)/2*voxelSize, ...
            [1 1]*zSteps(bottomZ), 'color', 'b')
        
        line([1 volumeSize(1)]*voxelSize, [1 1]*volumeSize(2)/2*voxelSize, ...
            [1 1]*zSteps(topZ), 'color', 'b')
        
        title('Plot true ray path')

        % Plot X error along Z
        subplot(2,2,3); hold on;
        for iOrigin = 1:nOrigins
            rayPath = permute(rayPathArray(:, :, iOrigin), [2 1]);
            
            tempRayPath = permute(trueRayPathArray(:, :, iOrigin), [2 1]);
            
            plot((rayPath(:,1) - tempRayPath(:,1))*10^3, zSteps, 'color', rayCols(iOrigin,:))
        end
        line([-20 20], [1 1]*zSteps(bottomZ), 'color', 'b')
        
        line([-20 20], [1 1]*zSteps(topZ), 'color', 'b')
        
        
        xlabel('Error in micron')
        title('Ray error')
        
        % plot spot diagram
        subplot(2,2,4); hold on; axis equal;
        for iOrigin = 1:nOrigins
%             rayPath = permute(rayPathArray(:, :, iOrigin), [2 1]);
% 
%             plot3(rayPath(:,1), rayPath(:,2), rayPath(:,3), 'color', rayCols(iOrigin,:));
        end
        title('Spot diagram')
    else
        % Plot ray path in object
        figure; subplot(1,3,1); hold on; axis equal;
        view(0, 0)
        
        if ~isempty(testResults)
           for iResult = 1:length(testResults)
               
                rayPathArray = testResults{iResult};
        
                for iOrigin = 1:nOrigins
                    rayPath = permute(rayPathArray(:, :, iOrigin), [2 1]);
                    if ~isempty(rayPath)

                        plot3(rayPath(:,1), rayPath(:,2), rayPath(:,3));
                    end
                end

           end
           
           legend(testTypes)
        else
            for iOrigin = 1:nOrigins
                rayPath = permute(rayPathArray(:, :, iOrigin), [2 1]);

                plot3(rayPath(:,1), rayPath(:,2), rayPath(:,3), 'color', rayCols(iOrigin,:));
            end
        end
        title('Plot ray path')
        
        % plot lines at top and bottom
        line([1 volumeSize(1)]*voxelSize, [1 1]*volumeSize(2)/2*voxelSize, ...
            [1 1]*zSteps(bottomZ), 'color', 'b')
        
        line([1 volumeSize(1)]*voxelSize, [1 1]*volumeSize(2)/2*voxelSize, ...
            [1 1]*zSteps(topZ), 'color', 'b')
        
        % Plot analytic ray path - for on axis rays along central y axis
        subplot(1,3,2); hold on; axis equal;
        view(0, 0)
        
        trueRayPathArray = zeros(3, volumeSize(3), nOrigins)*NaN;
        
        fiberCentre = volumeSize(1:2)/2*voxelSize;
        
        %%% Test compared to paper       
        for iOrigin = 1:nOrigins
            tempRayPath = zeros(3, volumeSize(3))*NaN;
            
            if rayOrigins(iOrigin,2) ~= fiberCentre(2) | ... 
                    any(startRayT - [0, 0, 1])
                %%% Can gernalize to X component and off-axis rays from
                %%% Section 5.3 in Merchland 1978
                
                error('Analytic solution only set up for Y on midline and parallel rays')
            end
            
            % Check it's within fiber radius
            if abs(rayOrigins(iOrigin,1)-fiberCentre(1)) < radius
            
                % Just a line up to enterance
                tempRayPath(1,1:bottomZ-1) = rayOrigins(iOrigin,1)-fiberCentre(1); % offset to fiber centre
                tempRayPath(2,:) = rayOrigins(iOrigin,2)-fiberCentre(2);
                tempRayPath(3,:) = zSteps;

                % Get initial RI
                initalRI = lensRIVolume(round(rayOrigins(iOrigin,1)/voxelSize), ...
                    round(rayOrigins(iOrigin,2)/voxelSize), bottomZ);

                % Ray path in fiber - Nishidate 2011, eqn 31
                %%% This seems incorrect
%                 tempRayPath(1, bottomZ:topZ) = tempRayPath(1,bottomZ-1)*cos((zSteps(bottomZ:topZ) - ...
%                     zSteps(bottomZ))/initalRI);

                % Gives correct result From Merchland 5.41 and other refs,
                tempRayPath(1, bottomZ:topZ) = tempRayPath(1,bottomZ-1)*cos((zSteps(bottomZ:topZ) - ...
                    zSteps(bottomZ))*alpha*n0);
                
                % Need to get analytic expression after exit
%  
%                  % Add fiber centre back for X
                 tempRayPath(1,:) = tempRayPath(1,:) + fiberCentre(1);
                 tempRayPath(2,:) = tempRayPath(2,:) + fiberCentre(2);
                 
                 plotVec = 0;
            else
                % outside fiber, just propogate along z
                
                tempRayPath(1,:) = rayOrigins(iOrigin,1); 
                tempRayPath(2,:) = rayOrigins(iOrigin,2);
                tempRayPath(3,:) = zSteps;
                
                plotVec = 0;
            end
            
            if plotVec
                % plot vectors for comparison
                line([0 5]*rVec(1) + tempRayPath(1, topZ), [0 5]*rVec(2) + tempRayPath(2, topZ), ...
                    [0 5]*rVec(3) + tempRayPath(3, topZ), 'color', 'm')
                
                line([0 4]*vec(1) + tempRayPath(1, topZ), [0 4]*vec(2) + tempRayPath(2, topZ), ...
                    [0 4]*vec(3) + tempRayPath(3, topZ), 'color', 'g')
                
                % Point on one period
                %plot3(tempRayPath(1,1), tempRayPath(2,1), tempRayPath(3,bottomZ) + rayPeriod, 'rx');
                
                % Point at second center after some dispersion
                plot3(fiberCentre(1), tempRayPath(2,1), tempRayPath(3,bottomZ) + rayPeriod*3/4, 'rx');
            end
            
            plot3(tempRayPath(1,:), tempRayPath(2,:), tempRayPath(3,:), 'color', rayCols(iOrigin,:)); 
            
            trueRayPathArray(:,:,iOrigin) = tempRayPath;
        end
        % plot lines at top and bottom
        line([1 volumeSize(1)]*voxelSize, [1 1]*volumeSize(2)/2*voxelSize, ...
            [1 1]*zSteps(bottomZ), 'color', 'b')
        
        line([1 volumeSize(1)]*voxelSize, [1 1]*volumeSize(2)/2*voxelSize, ...
            [1 1]*zSteps(topZ), 'color', 'b')
        
        title('Plot true ray path - check all')
        
        % Plot X error along Z
        subplot(1,3,3); hold on; %axis equal;
        if ~isempty(testResults)
           for iResult = 1:length(testResults)
               
               rayPathArray = testResults{iResult};
               
               for iOrigin = 1:nOrigins
                    rayPath = permute(rayPathArray(:, :, iOrigin), [2 1]);

                    tempRayPath = permute(trueRayPathArray(:, :, iOrigin), [2 1]);
                    if ~isempty(rayPath)

                        plot(rayPath(:,1) - tempRayPath(:,1), zSteps)
                    end
               end
           end
       else
            for iOrigin = 1:nOrigins
                rayPath = permute(rayPathArray(:, :, iOrigin), [2 1]);

                tempRayPath = permute(trueRayPathArray(:, :, iOrigin), [2 1]);

                plot(rayPath(:,1) - tempRayPath(:,1), zSteps, 'color', rayCols(iOrigin,:))
           end
       end
       
        line([-2 2]*voxelSize, [1 1]*zSteps(bottomZ), 'color', 'b')
        
        line([-2 2]*voxelSize, [1 1]*zSteps(topZ), 'color', 'b')
        
        xlim([-voxelSize voxelSize]*2)
        title('Ray error')
    end
end