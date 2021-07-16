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

useRealData = 1;
if useRealData
    dataFile = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/Data/Matlab RI Volumes/Test1_SD0_vox2.mat';
    voxelSize = 2*10^-3; % mm

    createGradedFiber = 0;
    createLunebergLens = 0;
    plotReferenceFiber = 0;
    useTestData = 0;
    
    exteriorRI = 1.33;
else
    useTestData = 1;
    
    radius = 1; %mm - used for both lens and fiber
    voxelSize = 0.065; % mm

    exteriorRI = 1;

    createLunebergLens = 1;
        
    createGradedFiber = 0;
        fiberLength = 10;
        n0 = sqrt(2.5);
        %%% 1/sqrt(2.5) matches refractive index profile described, 1/(2.5) gives curve similar to result shown
        alpha = 1/(2.5); 
        plotReferenceFiber = 0; % will plot data from Nishidate and only trace 0.5
end

% Options, 1, 4S, 4RKN, 5RKN, 45RKN 
% Difference 4 to 5 is quite small
interpType = '4S';
    %This doesn't have big effect on error < 10^-9, ok even at 10^-6 but then collapses
        % May effect peripheral rays more
    tolerance = 10^-9; % Nishidate uses 10^-12, seems a bit too tight
    
    % Seems good to start a bit lower than the general deltaS
        % in fixed step, 10^-3 vs. 10^-4 gives rough order of magnitude error
    initialDeltaS = 10^-3;
    
    % This keeps spot size surpsingly small, even on loose tolerance
        % (I guess it ends up stepping further at loose tolerance...)
    iterativeFinal = 0;  
    
testPlot = 1;
%% Load or create test data     
if useRealData
    temp = load(dataFile);
    lensRIVolume = temp.volume;
    clear temp
    
    warning('Flip volume to correct axis in original data (then image will show upside down)')
    
    % Flip to change z direction
    lensRIVolume = permute(flipud(permute(lensRIVolume,[3 1 2])), [2 3 1]);
    
    volumeSize = size(lensRIVolume);
    
    % Using negative coding for GRIN
    goodVolume = lensRIVolume < 0;
    
    lensRIVolume(goodVolume) = -lensRIVolume(goodVolume);
    
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
    
    goodVolume = ~isnan(lensRIVolume);

    lensRIVolume(isnan(lensRIVolume)) = exteriorRI;
end

%% Do general set up

% Get surface voxels
borderVolume = imdilate(~goodVolume, strel('Sphere',1)) & goodVolume;

surfaceInds = find(borderVolume);

[surfaceX, surfaceY, surfaceZ] = ind2sub(volumeSize, surfaceInds);

bottomZ = min(surfaceZ);

topZ = max(surfaceZ);

% Get border volume to test against
borderVolume = imdilate(borderVolume, strel('Sphere',1));

% Compute trigulation
%%% Move to other file later - replace surface points with exact values...
grinSurface = isosurface(goodVolume, 0.5);
grinSurface.vertices = grinSurface.vertices*voxelSize;

if useTestData
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

% Search is limited to good in points
goodInds = find(goodVolume);

if useRealData
    plotInds = 1:50:length(surfaceX);
    
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
    xStartPoints = 1:20:volumeSize(1); 
else
    if ~plotReferenceFiber
        xStartPoints = 1:5:volumeSize(1); 

    else
        xStartPoints = volumeSize(1)/2+radius/voxelSize*0.5;
    end

end

rayOrigins = [xStartPoints(:), ones(numel(xStartPoints),1)*volumeSize(2)/2,  ...
    ones(numel(xStartPoints),1)]*voxelSize; 

nOrigins = size(rayOrigins,1);

startRayT = [0, 0, 1];

rayCols = lines(nOrigins);

%% Run ray tracer

volumeCenter = volumeSize/2*voxelSize;

% Get border test values

%%% Note intersection functions are 1 when outside of GRIN, 0 inside
if useRealData
    % Calculating inpolyhedron is very slow, so test on volume first
    intersectionFn = @(x1)(~surfaceIntersectFunction(borderVolume, x1,voxelSize, grinSurface));
       
    lambdaFn = @(x1, x0)(0.5);
    
    warning('Improve this and add lambda')
    
elseif useTestData
    if createLunebergLens
        % Get sphere intersection
       intersectionFn = @(x1)(sqrt(sum((x1 - volumeCenter).^2)) > ...
           radius);

       % Calculate as radius of ratios.
       lambdaFn = @(x1, x0)((sqrt(sum((x1 - volumeCenter).^2))-radius)/...
           (sqrt(sum((x1 - volumeCenter).^2))-sqrt(sum((x0 - volumeCenter).^2))));

    elseif createGradedFiber
       % Just checks X position, not Y
       intersectionFn = @(x1)(~(x1(3) < volumeCenter(3) + fiberLength/2 & ...
            x1(3) > volumeCenter(3) - fiberLength/2  & ...
            x1(1) > volumeCenter(1) - radius & x1(1) < volumeCenter(1) + radius));

       lambdaFn = @(x1, x0)((x1(3)- (volumeCenter(3) + fiberLength/2))/(x1(3)-x0(3)));
    end
end

rayPathArray = zeros(3, volumeSize(3), nOrigins)*NaN;

rayPathLengthArray = zeros(volumeSize(3), nOrigins)*NaN;

if testPlot
    figure; hold on; axis equal;
    view(0, 0)
    plot3(surfaceX(plotInds)*voxelSize, surfaceY(plotInds)*voxelSize,...
        surfaceZ(plotInds)*voxelSize, '.')
    
    pH = patch(grinSurface);
    
    pH.FaceColor = 'g';
    pH.EdgeColor = 'none';
end

% Turn of warning on singular matrix, otherwise backslash will slow down a lot
warning('off', 'MATLAB:nearlySingularMatrix')

minDeltaS = ones(nOrigins, 1);

firstIntersect = zeros(nOrigins, 3);
finalIntersect = zeros(nOrigins, 3);

finalPathLength = zeros(nOrigins, 1);
finalRayT = zeros(nOrigins, 3);

if createGradedFiber
   % Change to half period as focus points are actually at 1/4 and 3/4 points
   halfPeriodInFiber = 2*pi/(n0*alpha)*0.5;   
   
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
    
    intersectResult = 1;
    
    pathLength = 0; % in grin area
    
    exitedFlag = 0;
    
    currentPeriod = 1; % only used for graded fiber
    
    warning('Need to calculate refraction between layers')
    
    while go

        % Test if entering a GRIN area
        if ~intersectResult & ~inGraded & ~exitedFlag

            % entering, need to find intersect to graded volume
            if useRealData
                warning('Need to find correct enterance intersect')
                
                % X and T combined on input, T needs to point back
                lineDef = [rayX -rayT];
                
                [intersect, position] = intersectLineMesh3d(lineDef, grinSurface.vertices, grinSurface.faces);
                
                if length(positions) > 1
                    error('Multiple intersects - not treated')
                end
                
                error('Check from here')
                
                rayX = intersect;
                
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
                    deltaTemp = ((volumeCenter(3) - fiberLength/2) - rayX(3))/rayT(3);
                    
                    rayX = rayX + rayT*deltaTemp;
                end
            end
            
            firstIntersect(iOrigin,:) = rayX;
            
            voxelX = round(rayX/voxelSize); 
            
            % Get surface normal and RI at entry point
            [~, rOut] = numerical_dT_dt(rayX, volCoords(goodInds,:), goodInds, lensRIVolume, []);
            
            if useRealData
                warning('Need to find correct enterance normal and rIn')
                
                surfaceNormal = [0 0 -1];
                
                rIn = exteriorRI;
                
            elseif useTestData
                if createLunebergLens
                    surfaceNormal = rayX - volumeSize/2*voxelSize;

                    surfaceNormal = surfaceNormal/norm(surfaceNormal);
 
                    % Should be very close to 1, but allowing this removes notch in error profile
                        % However, it does increase error spot slightly
                    %rOut = 1;
                    
                elseif createGradedFiber
                    % Just points backwards
                    surfaceNormal = [0 0 -1];
                    
                end
                
                rIn = exteriorRI;
            end
            
            % Calculate entry refraction
            nRatio = rIn/rOut;
            cosI = -dot(surfaceNormal, rayT);
            sinT2 = nRatio^2*(1-cosI^2);
            cosT = sqrt(1-sinT2);
            
            % Assuming all refracted, non reflected
            rayT = nRatio*rayT + (nRatio*cosI-cosT)*surfaceNormal;
            rayT = rayT/norm(rayT);
            
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
        
        %propogate ray
        if ~inGraded
            rayX = rayX + rayT*voxelSize/3;
            
            %%% Could also draw line until intersect reached, as in original surface ray-tracer
        else

            x0 = rayX;
            t0 = rayT;

            % Calc is fixed to isotropic RI
            [x, t, deltaS] = ray_interpolation(interpType, 'iso', x0', t0', deltaS, testCoords, ...
                testInds, lensRIVolume, tolerance, []);
            
            rayX = x';
            
            rayT = t';
            
            pathLength = pathLength + sqrt((x0(1) - rayX(1))^2 + (x0(2) - rayX(2))^2 + (x0(3) - rayX(3))^2);
            
            if deltaS < minDeltaS(iOrigin)
                minDeltaS(iOrigin) = deltaS;
            end
            
            if createGradedFiber
                % test if focal point is passed and capture if so
                if rayX(3)-(volumeSize(3)/2*voxelSize - fiberLength/2) > halfPeriodInFiber*(currentPeriod-0.5)

                    % Step back
                    deltaTemp = (halfPeriodInFiber*(currentPeriod-0.5) - (rayX(3)-(volumeSize(3)/2*voxelSize - fiberLength/2)))/rayT(3);

                    periodSpotsArray(:,currentPeriod,iOrigin) = rayX + deltaTemp*rayT; 

                    plot3(periodSpotsArray(1,currentPeriod,iOrigin), periodSpotsArray(2,currentPeriod,iOrigin), periodSpotsArray(3,currentPeriod,iOrigin), 'xr');

                    currentPeriod = currentPeriod + 1;
                end
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
            [testInds, tempInds] = intersect(testInds, goodInds);
            
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
                lambda0 = lambdaFn(rayX, x0);

                %error term, sqrt of machine precision, from Nishidate paper sqrt(1.49e-8)
                    % Take geometric mean of single and double precision, double takes ages to evaluate
                epsilon = sqrt(geomean([1.19e-7 2.22e-16])); 

                its = 1;

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
                            testInds_backwards, lensRIVolume, tolerance, []);

                        rayX = x'; rayT = t';
                    end

                    deltaS_final = (1-SF)*lambda0*deltaS0_final;

                    % This loop steps forward to find new contact position just before overshoot.
                    while ~intersectionFn(rayX)
                        x0 = rayX; %think this is the case, we are stepping forward
                        t0 = rayT;

                        [x, t] = ray_interpolation('5RKN', 'iso', x0', t0', deltaS_final, testCoords_forwards, ...
                            testInds_forwards, lensRIVolume, tolerance, []);

                        rayX = x'; rayT = t';
                        
                        pathLength = pathLength + sqrt((x0(1) - rayX(1))^2 + (x0(2) - rayX(2))^2 + (x0(3) - rayX(3))^2);
                    end

                    %update parameters for next loop
                    lambda0 = lambdaFn(rayX, x0);

                    deltaS0_final = deltaS_final;
                end

            else
               lambda0 = lambdaFn(rayX, x0);
               
                [x, t] = ray_interpolation(interpType, 'iso', x0', t0', deltaS*lambda0, testCoords_forwards, ...
                    testInds_forwards, lensRIVolume, tolerance, []);
            
               rayX = x';
               rayT = t';
               
               pathLength = pathLength + sqrt((x0(1) - rayX(1))^2 + (x0(2) - rayX(2))^2 + (x0(3) - rayX(3))^2);
            end
            
            % Normalize no before refraction
            rayT = rayT/norm(rayT);
            
            finalIntersect(iOrigin,:) = rayX;

            finalPathLength(iOrigin,:) = pathLength;
            
            finalRayT(iOrigin,:) = rayT;
            
            % Do refraction at border

            % Get surface normal and RI at exit point
            [~, rIn] = numerical_dT_dt(rayX, volCoords(goodInds,:), goodInds, lensRIVolume, []); 
                                 
            if useRealData
                warning('Need to calculate normal correctly and rOut')
                
                surfaceNormal = [0 0 -1]; 
                
                rOut = exteriorRI;
                            
            elseif useTestData
                if createLunebergLens
                    surfaceNormal = -(rayX - volumeSize/2*voxelSize);

                    surfaceNormal = surfaceNormal/norm(surfaceNormal);
                    
                    % Usually not exactly 1, this allows refraction but actually decreases error on exit trajectory
%                     rIn = 1;
  
                elseif createGradedFiber
                    % Just points backwards
                    surfaceNormal = [0 0 -1];   
                end
                
                rOut = exteriorRI;
            end
            
            % Calculate exit refraction
            nRatio = rIn/rOut;
            cosI = -dot(surfaceNormal, rayT);
            sinT2 = nRatio^2*(1-cosI^2);
            cosT = sqrt(1-sinT2);
            
            % Assuming all refracted, non reflected
            rayT = nRatio*rayT + (nRatio*cosI-cosT)*surfaceNormal;
            rayT = rayT/norm(rayT);

%             finalRayT(iOrigin,:)
%             rayT
            
            exitedFlag = 1;
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
            end
        end
        
        % Check ray has exitied volume
        if any(voxelX > volumeSize) | any(voxelX < 1)
           go = 0; 
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
                startZ = find(goodVolume(round(rayOrigins(iOrigin,1)/voxelSize),...
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

%                 finalRayT(iOrigin,:)/finalRayT(iOrigin,3)
                
                % Propogate ray
                tempRayPath(1, topZ+1:end) = (tempRayPath(3, topZ+1:end)-finalR)*finalT(1) + finalX;
                
%                                 tempRayPath(1, topZ+1:end) = NaN;
                
                
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
                    any(startRayT - [0, 0, 1])
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
                    %%% Seems incorrect, or not general
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
                
                % Gives correct result From Merchland 5.41 and other refs, 
                % I think this should be sqrt(alpha) from book, but only works if alpha = 1/n0^2
                beta = alpha*n0;
                tempRayPath(1, bottomZ:topZ) = tempRayPath(1,bottomZ)*cos(tempRayPath(3,bottomZ:topZ)*beta);
                
                % Analytic expresion for ray after exit
                % Get final Coord and RI
                finalX = tempRayPath(1,bottomZ)*cos(fiberLength*beta);
                finalCoord = [finalX + fiberCentre(1), fiberCentre(2), fiberCentre(3) + fiberLength/2];
                
                [~, rIn] = numerical_dT_dt(finalCoord, volCoords(goodInds,:), goodInds, lensRIVolume, []);  

                % Get final vector and normalize
                finalP = -beta*tempRayPath(1,bottomZ)*sin(fiberLength*beta);
                finalT = [finalP, 0, 1];
                finalT = finalT/norm(finalT);
                
                surfaceNormal = [0 0 -1];
                rOut = exteriorRI;
                
                % Refract
                nRatio = rIn/rOut;
                cosI = -dot(surfaceNormal, finalT);
                sinT2 = nRatio^2*(1-cosI^2);
                cosT = sqrt(1-sinT2);

                % Assuming all refracted, non reflected
                refractT = nRatio*finalT + (nRatio*cosI-cosT)*surfaceNormal;
                
                % No need to normalize as rescaled
                refractT = refractT/refractT(3);
                
                % Propogate ray
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

                plot((rayPath(:,1) - tempRayPath(:,1))*10^6, zSteps, 'color', rayCols(iOrigin,:))
            end
        end
       
       if plotReferenceFiber
          plot(nishidatError(:,2), nishidatError(:,1)+volumeSize(3)/2*voxelSize-fiberLength/2, 'r') 
          
          % Recalcualte error for traced vs calculated
          % First interp to same zsteps
          interpNishidatPath = interp1(nishidatPath(:,2), nishidatPath(:,1), tempRayPath(:,3), 'linear', 0);
          
          interpNishidatPath(interpNishidatPath == 0) = NaN;
          
          interpNishidatPath = interpNishidatPath - volumeSize(1)/2*voxelSize;
          
          plot((interpNishidatPath - nishidateSolution')*100, zSteps, 'r--')

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

function  intersect = surfaceIntersectFunction(volume, x, scale, surface) 
    
    if volume(round(x(1)/scale),round(x(2)/scale),round(x(3)/scale)) 
        intersect = inpolyhedron(surface, x);
    else
        intersect = 0;
    end
end