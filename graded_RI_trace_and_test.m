% Set up to do both tracing on data and have testing mode on either
% lundaberg lens or graded fibre, as in nishdate 2011 papers
    
    % Also need to add interface refraction - not required for test cases but
    % for normal tracing. Should also do accuracy test for my patch method.
    
    % Should test on smoothed RI map calculated from 200 nm maps
        % Have exterior RI loading and processing function
        % May also want to smooth

    % Start with fixed on axis RI but later update to recalc with angle dependence
    
    % Note grin is basically only in trasition area, so can skip calc if
    % all values in region are equal
        % Should test this speed up
        
    % Create 3d view of rays and spot diagram (for fibre do on phase and half phase)
        % Add fiber test case of inclined rays
                 
clear; 
clc; close all

%% Set parameters
voxelSize = 0.065; % mm

useRealData = 0;
%%% Add parameters here

useTestData = 1;

exteriorRI = 1;

radius = 1; %mm - used for both lens and fiber

    createLunebergLens = 1;
        
    createGradedFiber = 0;
        fiberLength = 10;
        n0 = sqrt(2.5);
        alpha = 1/(2.5); sqrt(2.5);
        
% Options, 1, 4S, 4RKN, 5RKN, 45RKN 
% Difference 4 to 5 is quite small
interpType = '45RKN';
    %This doesn't have big effect on error < 10^-9, ok even at 10^-6 but then collapses
        % May effect peripheral rays more
    tolerance = 10^-9; % Nishidate uses 10^-12, seems a bit too tight
    
    % Seems good to start a bit lower than the general deltaS
        % in fixed step, 10^-3 vs. 10^-4 gives rough order of magnitude error
    initialDeltaS = 10^-3;
    
    % This keeps spot size surpsingly small, even on loose tolerance
        % (I guess it ends up stepping further at loose tolerance...)
    iterativeFinal = 1;

% Includes extra points around GRIN region in gradient calc for entry and exit points
    % on entry effects inital notch and error along whole path, doesn't effect end much with or without iterative
    % This is probably wrong as increases border RI
    %%% should remove!!!
growGoodInds = 0;    
    growOnFinal = 0;

% For fiber - 0.25 is border, 0.375 is 25% and 0.5 is center.    
startPos = 0.375;    
    
% Manual testing, a bit crap...
iRayTest = 3;

testResults = [];

if iRayTest == 1
    testTypes = cell(3,1);
    testResults = cell(3,1);
end
        
testPlot = 1;
%% Load or create test data     
if useRealData
    
    %%% If not doing direction updating, calcualte RI for all
    
    if growGoodInds
        growGoodInds = 0;
        warning('Turned of growing good inds - not handled for real data')
    end
    
elseif useTestData
    if ~xor(createLunebergLens, createGradedFiber)
        error('Both or neither lens/fiber flags set. Check.')
    end
    
    if createLunebergLens
        volumeSize = ceil(radius*2*1.5/voxelSize*[1 1 1]);
        
        lensRIVolume = createluneberglens(radius/voxelSize, volumeSize);
    
        if growGoodInds & exteriorRI ~= 1
           error('Growing good inds will create error when RI not equal to one') 
        end
        
    elseif createGradedFiber
        volumeSize = ceil([radius*2 radius*2 fiberLength]*1.5/voxelSize);
        
        lensRIVolume = createGradedFiber(radius/voxelSize, fiberLength/voxelSize, volumeSize, ...
            n0, alpha, voxelSize);
        
        if growGoodInds
            growGoodInds = 0;
            warning('Turned of growing good inds - not handled for fiber')
            %%% Probably also need to grow RI of GRIN region
        end
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
    
    % Nishidate eqn
    if createGradedFiber
        plot(centreLine, sqrt(2.5-(tempRadius*voxelSize/radius).^2), 'rx');

        plot(centreLine, sqrt(n0^2*(1-alpha^2*(tempRadius*voxelSize).^2)), 'g');
    end
end

%% Do general set up

% Search is limited to good in points
tempVolume = imdilate(goodVolume, strel('sphere',1));

extraInds = find(tempVolume);

goodInds = find(goodVolume);


plotInds = 1:length(surfaceX);

% Get voxel coordinates
[tempX, tempY, tempZ] = meshgrid(1:volumeSize(1), 1:volumeSize(2), 1:volumeSize(3));

volInds = sub2ind(volumeSize, tempX(:), tempY(:), tempZ(:));

volCoords = ([tempX(:), tempY(:), tempZ(:)])*voxelSize;

zSteps = (1:volumeSize(3))*voxelSize;

%%% Replace with meshgrid later on
xStartPoints = 1:5:volumeSize(1); volumeSize(1)*startPos; %

% xStartPoints = xStartPoints(3:5);

rayOrigins = [xStartPoints(:), ones(numel(xStartPoints),1)*volumeSize(2)/2,  ...
    ones(numel(xStartPoints),1)]*voxelSize; 

nOrigins = size(rayOrigins,1);

startRayT = [0, 0, 1];

rayCols = lines(nOrigins);

%% Run ray tracer

% Get border test value
if useRealData

elseif useTestData
    if createLunebergLens
        % Get sphere intersection

       intersectionFn = @(x1)(sqrt(sum((x1 - volumeSize/2*voxelSize).^2)) > ...
           radius);

       lambdaFn = @(x1, x0)((sqrt(sum((x1 - volumeSize/2*voxelSize).^2))-radius)/...
           (sqrt(sum((x1 - volumeSize/2*voxelSize).^2))-sqrt(sum((x0 - volumeSize/2*voxelSize).^2))));

    elseif createGradedFiber
       referenceDistance = zSteps(topZ+1);

       intersectionFn = @(x1)(x1(3) > referenceDistance);

       lambdaFn = @(x1, x0)((x1(3)-referenceDistance)/(x1(3)-x0(3)));
    end
end

rayPathArray = zeros(3, volumeSize(3), nOrigins)*NaN;

rayPathLengthArray = zeros(volumeSize(3), nOrigins)*NaN;

if testPlot
    figure; hold on; axis equal;
    view(0, 0)
    plot3(surfaceX(plotInds)*voxelSize, surfaceY(plotInds)*voxelSize,...
        surfaceZ(plotInds)*voxelSize, '.')
end

% Turn of warning on singular matrix, otherwise backslash will slow down a lot
warning('off', 'MATLAB:nearlySingularMatrix')

minDeltaS = ones(nOrigins, 1);

firstIntersect = zeros(nOrigins, 3);
finalIntersect = zeros(nOrigins, 3);

finalPathLength = zeros(nOrigins, 1);

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
    
    intersectResult = 1;
    
    pathLength = 0; % in grin area
    
    exitedFlag = 0;
    extraIndsUsed = 0;
    
    useExtraInds = 0;
    
    while go

        % Test if entering a GRIN area
        if ~intersectResult & ~inGraded & ~exitedFlag

            % entering, need to find intersect to graded volume
            if useRealData
                
            elseif useTestData
                if createLunebergLens                 
                    % Get sphere intersection
                    origin = rayX - volumeSize/2*voxelSize;
                    centre = 0;
                    originalDistance = centre - origin;
                    
                    t_ca = dot(originalDistance,rayT);
                    d = sqrt(dot(originalDistance,originalDistance)-t_ca^2);
                    t_hc = sqrt(radius^2-d^2);
                    deltaTemp = t_ca - t_hc;

                    % Steps ray back to intersection
                    rayX = rayX + rayT*deltaTemp; %3x1 coordinates
                    
                elseif createGradedFiber
                    % Get z distance past bottom of gradient section
                    deltaTemp = (zSteps(bottomZ) - rayX(3))/rayT(3);
                    
                    % Step back to algin exactly
                    rayX = rayX + rayT*deltaTemp;
                end
            end
            
            firstIntersect(iOrigin,:) = rayX;
            
            voxelX = round(rayX/voxelSize); 
            
            % Get surface normal at entry point
            if useRealData
                
            elseif useTestData
                if createLunebergLens
                    % Not really important, shouldn't have refraction in air
                    surfaceNormal = rayX - volumeSize/2*voxelSize;

                    surfaceNormal = surfaceNormal/norm(surfaceNormal);
 
                    rIn = 1;
                elseif createGradedFiber
                    % Just points backwards
                    surfaceNormal = [0 0 -1];
                    
                    [~, rIn] = numerical_dT_dt(rayX, volCoords(goodInds,:), goodInds, lensRIVolume);   
                end
            end
            
            rOut = exteriorRI;
            
            nRatio = rIn/rOut;
            cosI = -dot(surfaceNormal, rayT);
            sinT2 = nRatio^2*(1-cosI^2);
            cosT = sqrt(1-sinT2);
            
            % Assuming all refracted, non reflected
            rayT = nRatio*rayT + (nRatio*cosI-cosT)*surfaceNormal;
            
            inGraded = 1;
        end
        
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
            [testInds, tempInds] = intersect(testInds, goodInds);
            
            testCoords = testCoords(tempInds,:)*voxelSize;
            
            % only nearest 128 used, but seems to need some extra
            [~, nearInds] = sort(sqrt((testCoords(:,1)-rayX(1)).^2+(testCoords(:,2)-rayX(2)).^2+...
                (testCoords(:,3)-rayX(3)).^2));
            
            if length(nearInds) > 200
                nearInds = nearInds(1:200);
            end
            
            testInds = testInds(nearInds);
            
            testCoords = testCoords(nearInds,:);
            
%             figure; plot3(testCoords(:,1), testCoords(:,2), testCoords(:,3), '.')
%             hold on; plot3(voxelX(1), voxelX(2), voxelX(3), 'rx')
            
            %%% If doing direction updating, calculate RI for current direction
            
            [rayX rayT deltaS*10^6 length(testInds)]
        end
        
        %propogate ray
        if ~inGraded
            rayX = rayX + rayT*voxelSize/3;
            
            %%% Could also draw line until intersect reached, as in original surface ray-tracer
        else

            x0 = rayX;
            t0 = rayT;

            % Calc is fixed to isotropic RI
            if ~useExtraInds
                [x, t, deltaS] = ray_interpolation(interpType, 'iso', x0', t0', deltaS, testCoords, ...
                    testInds, lensRIVolume, tolerance);
            else
                [x, t, deltaS] = ray_interpolation(interpType, 'iso', x0', t0', deltaS, volCoords(extraInds,:), ...
                    extraInds, lensRIVolume, tolerance);
                
                extraIndsUsed = 1;
                
                useExtraInds = 0;
            end
            
            rayX = x';
            
            rayT = t';
            
            pathLength = pathLength + sqrt((x0(1) - rayX(1))^2 + (x0(2) - rayX(2))^2 + (x0(3) - rayX(3))^2);
            
            if deltaS < minDeltaS(iOrigin)
                minDeltaS(iOrigin) = deltaS;
            end
            
        end 
        
        intersectResult = intersectionFn(rayX);
        
        % Test if leaving a GRIN area
        if intersectResult & inGraded & ~exitedFlag

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
            if growOnFinal & growGoodInds
                [testInds, tempInds] = intersect(testInds, extraInds);
            else
                [testInds, tempInds] = intersect(testInds, goodInds);
            end
            
            testCoords = testCoords(tempInds,:)*voxelSize;
            
            % only nearest 128 used, but seems to need some extra
            [~, nearInds] = sort(sqrt((testCoords(:,1)-x0(1)).^2+(testCoords(:,2)-x0(2)).^2+...
                (testCoords(:,3)-x0(3)).^2));
            
            if length(nearInds) > 200
                nearInds = nearInds(1:200);
            end
            
            testInds_forwards = testInds(nearInds);
            testCoords_forwards = testCoords(nearInds,:);
            
            % Get backwards sorted as well
            [~, nearInds] = sort(sqrt((testCoords(:,1)-rayX(1)).^2+(testCoords(:,2)-rayX(2)).^2+...
                (testCoords(:,3)-rayX(3)).^2));
            
            if length(nearInds) > 200
                nearInds = nearInds(1:200);
            end
            
            testInds_backwards = testInds(nearInds);
            testCoords_backwards = testCoords(nearInds,:);
            
            if any(goodVolume(testInds) == 0)
                extraIndsUsed = 1;
            end
            
            if iterativeFinal
                % Test intersection, leaving could just be a rounding point           
                deltaS0_final = deltaS;
                lambda0 = lambdaFn(rayX, x0);

                %error term, sqrt of machine precision, from Nishidate paper sqrt(1.49e-8)
                    % Take geometric mean of single and double precision, double takes ages to evaluate
                epsilon = sqrt(geomean([1.19e-7 2.22e-16])); 

                its = 1;

                if createLunebergLens 
                    difference_init = (radius - sqrt(sum((x0 - volumeSize/2*voxelSize).^2)))*10^6;
                elseif createGradedFiber
                    difference_init = (referenceDistance - x0(3))*10^6;
                end

                while deltaS0_final*norm(t0) > epsilon && abs(lambda0-1) > 10^-5 
                    %%% using 10^-3 as test for zero in tests above and below

                    SF = 1; %A from paper - saftey factor

                    % included in loop on paper but never hit
                    if lambda0 < 10^-5
                       rayT = t0;
                       warning('Not sure on this');
                       break
                    end

                    %this loop steps back from overshoot
                    while intersectionFn(rayX)
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
                    while ~intersectionFn(rayX)
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

                    its = its + 1;
                end

                if createLunebergLens 
                    difference_end = (radius - sqrt(sum((x0 - volumeSize/2*voxelSize).^2)))*10^6;
                elseif createGradedFiber
                    difference_end = (referenceDistance - x0(3))*10^6;
                end

                if useTestData
                    [its-1 difference_init difference_end]
                end

                % Get surface normal at exit point
                if useRealData

                elseif useTestData
                    if createLunebergLens
                        % Not really important, no refraction
                        surfaceNormal = -(rayX - volumeSize/2*voxelSize);

                        surfaceNormal = surfaceNormal/norm(surfaceNormal);

                    elseif createGradedFiber

                        % Just points backwards
                        surfaceNormal = [0 0 -1];
                    end
                end

            else
               lambda0 = lambdaFn(rayX, x0);
               
                [x, t] = ray_interpolation(interpType, 'iso', x0', t0', deltaS*lambda0, testCoords_forwards, ...
                    testInds_forwards, lensRIVolume, tolerance);
            
               rayX = x';

               rayT = t';
               
               pathLength = pathLength + sqrt((x0(1) - rayX(1))^2 + (x0(2) - rayX(2))^2 + (x0(3) - rayX(3))^2);
            end
            
            finalIntersect(iOrigin,:) = rayX;

            finalPathLength(iOrigin,:) = pathLength;
            
            exitedFlag = 1;
            
            % Do refraction at border
            %%% What case is this for
            
%             nExterior = lensRIVolume(voxelX(1), voxelX(2), voxelX(3));
%             % Final RI (x0 will be effectively equal x, but use first...)
%             [~, nInterior] = numerical_dT_dt(x0, testCoords*voxelSize, testInds, lensRIVolume);
% 
%             % Calculate refraction (normalize first)
%                 % Note surface normal should point into GRIN region
%             rayT = rayT/norm(rayT);
% 
%             nRatio = (nInterior/nExterior);
%             cosI = -dot(surfaceNormal, rayT);
%             sinT2 = nRatio^2*(1-cosI^2);
%             cosT = sqrt(1-sinT2);
% 
%             if sinT2 < 1
%                 % normal refraction
%                 rayT = nRatio*rayT + (nRatio*cosI-cosT)*surfaceNormal;
%             else
%                % total internal reflection
% 
%                % haven't really delt decide how to deal with this yet
%                asin(nExterior/nInterior)/pi*180 %critical angle
%                error('TIR')
% 
%                % reflection
%                rayT = rayT-2*dot(surfaceNormal,rayT)*surfaceNormal;
%             end
%             rayT - t0/norm(t0)

            % Normalize ray trajectory 
            rayT = rayT/norm(rayT);

            inGraded = 0;
            
        end
        
        voxelX = round(rayX/voxelSize);
        
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
                
                if extraIndsUsed
                    plot3(rayX(1), rayX(2), rayX(3), 'rx');
                end
                
            end
            
            extraIndsUsed = 0;
        end
        
        % Check ray has exitied volume
        if any(voxelX > volumeSize) | any(voxelX < 1)
           go = 0; 
        end
    end
    
    pause(0.01)
end

warning('on', 'MATLAB:nearlySingularMatrix')

if ~isempty(testResults)
    testResults{iRayTest} = rayPathArray;
    testTypes{iRayTest} = interpType
end

tolerance
minDeltaS(minDeltaS < 1)*10^6
mean(minDeltaS(minDeltaS < 1))*10^6

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
                %%% Can gernalize to X component and off-axis rays from
                %%% Eq 2 in Babayigit ... Turduev 2019
                
                error('Analytic solution only set up for Y on midline and parallel rays')
            end
            
            % Check it's within fiber radius
            if abs(rayOrigins(iOrigin,1)-lensCentre(1)) < radius

                %Find first pixel in lens for ray
                startZ = find(goodVolume(round(rayOrigins(iOrigin,1)/voxelSize),...
                    round(rayOrigins(iOrigin,2)/voxelSize), :));

                startZ = startZ(1);
                startZPerOrigin(iOrigin) = startZ;
                
                % Just a line up to enterance
                tempRayPath(1,1:startZ-1) = rayOrigins(iOrigin,1)-lensCentre(1); % offset to sphere centre
                tempRayPath(2,:) = rayOrigins(iOrigin,2) - lensCentre(2);
%                 tempRayPath(3,:) = zSteps - lensCentre(3);
                % Use actual z values
                tempRayPath(3,:) = rayPathArray(3, :, iOrigin) - lensCentre(3);
                
                % get non stepped z values for start and end
                tempRayPath(3,startZ) = firstIntersect(iOrigin,3)- lensCentre(3);
                tempRayPath(3,topZ) = finalIntersect(iOrigin,3)- lensCentre(3);
                
                %%% Would be good to have solution as function of optical path length
                
                % From eqn 3 in Turduev cloaking paper (eqn 2 for non parallel rays)
                tempRayPath(1, startZ:topZ) = tempRayPath(1,1)*(tempRayPath(3,startZ)*tempRayPath(3,startZ:topZ) + ...
                    radius*sqrt(radius^2 + tempRayPath(3,startZ)^2 - tempRayPath(3,startZ:topZ).^2))/...
                    (tempRayPath(3,startZ)^2 + radius^2);

                %%% Doesn't plot ray after exit
                tempRayPath(:, topZ+1:end) = NaN;
                
                % Add fiber centre back for X
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
                
                plot((rayPath(topZ,1) - tempRayPath(topZ,1))*10^6, zSteps(topZ), 'o', 'color', rayCols(iOrigin,:))
            end
        end
%         plot(circleX, circleY, 'b')
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
                % Use actual z value after step.
                tempRayPath(3,:) = rayPathArray(3, :, iOrigin) - zSteps(bottomZ);

                % Get initial RI
                initalRI = lensRIVolume(round(rayOrigins(iOrigin,1)/voxelSize), ...
                    round(rayOrigins(iOrigin,2)/voxelSize), bottomZ);

                %%% Can format in terms of path length?
                
                % Ray path in fiber - Nishidate 2011, eqn 31
                %%% This seems incorrect
%                 tempRayPath(1, bottomZ:topZ) = tempRayPath(1,bottomZ-1)*cos((tempRayPath(3,bottomZ:topZ) - ...
%                     zSteps(bottomZ))/1.5);
                
%                 tempRayPath(1, bottomZ:topZ) = tempRayPath(1,bottomZ-1)*cos((tempRayPath(3,bottomZ:topZ) - ...
%                     zSteps(bottomZ))/initalRI);

                % Gives correct result From Merchland 5.41 and other refs,
                tempRayPath(1, bottomZ:topZ) = tempRayPath(1,bottomZ-1)*cos(tempRayPath(3,bottomZ:topZ)*alpha*n0);
                
                % Need to get analytic expression after exit
%  
%                  % Add fiber centre back for X
                 tempRayPath(1,:) = tempRayPath(1,:) + fiberCentre(1);
                 tempRayPath(2,:) = tempRayPath(2,:) + fiberCentre(2);
                 tempRayPath(3,:) = tempRayPath(3,:) + zSteps(bottomZ);
                 
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
                
                worstError = max(abs(rayPath(:,1) - tempRayPath(:,1))) 
                
                endError = abs(rayPath(topZ,1) - tempRayPath(topZ,1))
           end
       end
       
        line([-2 2]*voxelSize, [1 1]*zSteps(bottomZ), 'color', 'b')
        
        line([-2 2]*voxelSize, [1 1]*zSteps(topZ), 'color', 'b')
        
        xlim([-1 1]*10^-5)
        title('Ray error')
    end
end