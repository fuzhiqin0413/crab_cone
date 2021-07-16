clear
clc
close all

% Load in angle data
dataFolder = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/Data/liquid angle csv/200nm_thick';

convertAngleConvention = 0;

testingLunda = 0;

% Convert to mm
voxelSize = 800/10^6;

scaleSteps = 800/voxelSize/10^6;

% Note that epsilon (permitivitty) is n (RI) squared
nPerp = 1.53; % Short axis

nPara = 1.47; %(1.53-0.0024); % Long axis

nDif = nPerp - nPara;

%%% note angles are converted from Oliver's convention to Hwang & Rey 2005
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

% Calculate vector components
fiberXVolume = cos(phiVolume).*cos(thetaVolume);

fiberYVolume = sin(phiVolume).*cos(thetaVolume);

fiberZVolume = sin(thetaVolume);

%% Create lundberg lens for testingLunda
if ~testingLunda
    % Get good volume
    goodVolume = ~isnan(phiVolume) & ~isnan(thetaVolume);
else
    test_volume = ones(volumeSize);

    radius = volumeSize(3)*voxelSize/2;

    for i = 1:volumeSize

        for j = 1:volumeSize

            for k = 1:volumeSize

                x = (i-volumeSize(1)/2)*voxelSize;

                y = (j-volumeSize(2)/2)*voxelSize;

                z = (k-volumeSize(3)/2)*voxelSize;

                r = sqrt(x.^2+y.^2+z.^2);

                if r < radius
                    test_volume(i,j,k) = sqrt(2-(r/radius)^2);
                end
            end
        end
    end

    goodVolume = ~isnan(test_volume);
    
    test_volume(isnan(test_volume)) = 1;
end

%% Continue set up
goodInds = find(goodVolume);

% Get surface voxels
tempVolume = imdilate(goodVolume, strel('Sphere',1)) - goodVolume;

surfaceInds = find(tempVolume);

[surfaceX, surfaceY, surfaceZ] = ind2sub(volumeSize, surfaceInds);

% Get voxel coordinates
[tempX, tempY, tempZ] = meshgrid(1:volumeSize(1), 1:volumeSize(2), 1:volumeSize(3));

volInds = sub2ind(volumeSize, tempX(:), tempY(:), tempZ(:));

volCoords = ([tempX(:), tempY(:), tempZ(:)])*voxelSize;

% volInds = volInds(goodInds);
% 
% volCoords = volCoords(goodInds,:);

zSteps = (1:volumeSize(3))*voxelSize;

% set up initial rays.
xSteps = 1:10*scaleSteps:volumeSize(1);

%[xSteps, ySteps] = meshgrid(xSteps, xSteps);

%rayOrigins = [xSteps(:), ySteps(:), ones(numel(xSteps),1)]*voxelSize; 

rayOrigins = [xSteps(:),  ones(numel(xSteps),1)*round(volumeSize(2)/2), ones(numel(xSteps),1)]*voxelSize; 

% Just get origins on base of cone
rayOrigins = intersect(rayOrigins, volCoords(goodInds,:), 'rows');

%%
% test plot
figure; hold on, axis equal

plot3(volCoords(goodInds,1), volCoords(goodInds,2), volCoords(goodInds,3), '.'); 

plot3(rayOrigins(:,1), rayOrigins(:,2), rayOrigins(:,3), 'rx')

topInds = find(volCoords(goodInds,3) == min(volCoords(:,3)));

plot3(volCoords(goodInds(topInds),1), volCoords(goodInds(topInds),2), volCoords(goodInds(topInds),3), 'g.'); 

bottomInds = find(volCoords(goodInds,3) == max(volCoords(:,3)));

plot3(volCoords(goodInds(bottomInds),1), volCoords(goodInds(bottomInds),2), volCoords(goodInds(bottomInds),3), 'g.'); 

% output plot
figure; subplot(1,3,1); hold on, axis equal

plot3(rayOrigins(:,1), rayOrigins(:,2), rayOrigins(:,3), 'rx')

plot3(volCoords(goodInds(topInds),1), volCoords(goodInds(topInds),2), volCoords(goodInds(topInds),3), 'g.'); 

plot3(volCoords(goodInds(bottomInds),1), volCoords(goodInds(bottomInds),2), volCoords(goodInds(bottomInds),3), 'g.'); 

interpType = '5';

pathDifference = zeros(size(rayOrigins,1), 3);

% Exterior will be nPara (lowest)
RIVolume = ones(volumeSize)*1*nPara;

% Set for all
warning('update angle to take half period not quarter')
    %%% And below
angle2Fiber = acos( (fiberXVolume(goodInds)*0 + fiberYVolume(goodInds)*0 + fiberZVolume(goodInds)*1) ./ ...
        (sqrt(fiberXVolume(goodInds).^2 + fiberYVolume(goodInds).^2 + fiberZVolume(goodInds).^2) * norm([0 0 1])) );

RIVolume(goodInds) = nPara + nDif*sin(angle2Fiber);

% figure;
% imshow((RIVolume(:,:,round(volumeSize(3)*0.2))-1.3)/0.3)
% 
% figure;
% imshow((flipud(permute(RIVolume(round(volumeSize(1)/2),:,:), [3 2 1]))-1.3)/0.3)
% line([volumeSize(1)/2 volumeSize(1)/2], [1 volumeSize(3)], 'color', 'g')


dsP = 50;

figure; hold on; axis equal
toPlotInds = randperm(length(surfaceX));
toPlotInds = toPlotInds(1:dsP:length(surfaceX));

plot3(surfaceX(toPlotInds)*voxelSize, surfaceY(toPlotInds)*voxelSize, surfaceZ(toPlotInds)*voxelSize, 'g.');
line([volumeSize(1)/2 volumeSize(1)/2]*voxelSize, [volumeSize(1)/2 volumeSize(1)/2]*voxelSize, ...
    [-50 volumeSize(3)+50]*voxelSize, 'color', 'g')

axis off; set(gcf, 'Color','k', 'InvertHardcopy', 'off');

for i = 1:size(rayOrigins,1); %[5 6 12 13] % for 19
    rayT = [0, 0, 1]; %3x1 : ray will move from bottom to top
    
    rayX = round(rayOrigins(i,:)/voxelSize);
    
    inGraded = 0;
    
    x = rayOrigins(i,:);
    
    go = 1;

    lastNearbyX = x;

    % Calculate RI volume for all to start
    % Based on angle between ray and fiber vector
        %%% update angle to take half period not quarter
%     angle2Fiber = acos( (fiberXVolume(goodInds)*rayT(1) + fiberYVolume(goodInds)*rayT(2) + fiberZVolume(goodInds)*rayT(3)) ./ ...
%         (sqrt(fiberXVolume(goodInds).^2 + fiberYVolume(goodInds).^2 + fiberZVolume(goodInds).^2) * norm(rayT)) );
% 
%     RIVolume(goodInds) = nPara + nDif*sin(angle2Fiber);
    
    % Get values to test
    % For unit matrix, radius 4 will take ~268 (*2?) points
    [tempX, tempY, tempZ] = meshgrid((rayX(1)-4:rayX(1)+4)', (rayX(2)-4:rayX(2)+4)', (rayX(3)-4:rayX(3)+4)');

    testCoords= [tempX(:), tempY(:), tempZ(:)];
    
    tempInd = find(testCoords(:,1) < 1 | testCoords(:,2) < 1 | testCoords(:,3) < 1 | ...
            testCoords(:,1) > volumeSize(1) | testCoords(:,2) > volumeSize(2) | testCoords(:,3) > volumeSize(3));

    testCoords(tempInd, :) = [];

    testInds = sub2ind(volumeSize, testCoords(:,1), testCoords(:,2), testCoords(:,3));
    
    deltaS = 10^-4;
    
    rayPath = zeros(volumeSize(3), 3)*NaN;
    
    currentStep = 1;
    
    exitStep = 0;  
    
    rayPath(currentStep,:) = x;
    
    its = 1;
    
    while go
        % test if nearest voxel is 
        if goodVolume(rayX(1), rayX(2), rayX(3)) == 0
            %indicate where border touched
            if inGraded
                %plot3(x(1), x(2), x(3), 'mx')
                
                exitStep = currentStep - 1;
            end
            
            inGraded = 0;
        else
            if ~inGraded
                %plot3(x(1), x(2), x(3), 'cx')
            end
            
            inGraded = 1;
        end
        
        %check if x has moved to include new voxels
        if (norm(x-lastNearbyX) > voxelSize) & inGraded
            %if so, reselect nearby point

            lastNearbyX = x;

            [tempX, tempY, tempZ] = meshgrid((rayX(1)-4:rayX(1)+4)', (rayX(2)-4:rayX(2)+4)', (rayX(3)-4:rayX(3)+4)');

            testCoords= [tempX(:), tempY(:), tempZ(:)];

            tempInd = find(testCoords(:,1) < 1 | testCoords(:,2) < 1 | testCoords(:,3) < 1 | ...
                    testCoords(:,1) > volumeSize(1) | testCoords(:,2) > volumeSize(2) | testCoords(:,3) > volumeSize(3));
                
            testCoords(tempInd, :) = [];
            
            testInds = sub2ind(volumeSize, testCoords(:,1), testCoords(:,2), testCoords(:,3));
            
            % update RI volume
            %%% Need to update so calc is just on good inds.
                %%% update angle to take half period not quarter
%             angle2Fiber = acos( (fiberXVolume(volInds(nearbyInds))*rayT(1) + fiberYVolume(volInds(nearbyInds))*rayT(2) + fiberZVolume(volInds(nearbyInds))*rayT(3) ) ./ ...
%                 (sqrt(fiberXVolume(volInds(nearbyInds)).^2 + fiberYVolume(volInds(nearbyInds)).^2 + fiberZVolume(volInds(nearbyInds)).^2) * norm(rayT)) );
% 
%             RIVolume(volInds(nearbyInds)) = nPara + nDif*sin(angle2Fiber);
            

            [x rayT deltaS*10^3]
        end
        
        if ~inGraded | length(testInds) < 128
            x = x + rayT*deltaS; 
        else

            x0 = x;
            t0 = rayT;

            if ~testingLunda
                [x, t, deltaS] = ray_interpolation(interpType, 'iso', x0', t0', deltaS, testCoords*voxelSize, ...
                    testInds, RIVolume);
            else
                [x, t, deltaS] = ray_interpolation(interpType, 'iso', x0', t0', deltaS, testCoords*voxelSize, ...
                    testInds, test_volume);
            end

            x = x';
            
            rayT = t';
            
        end 
        
        rayX = round(x/voxelSize);
        
        its = its + 1;
        
        % record path for plotting later
        if x(3) >= zSteps(currentStep+1)
            if currentStep + 1 < length(zSteps)
                currentStep = currentStep + 1;
            end
            
            rayPath(currentStep,:) = x;
            
            if inGraded
%                 plot3(x(1), x(2), x(3), '.b');
% 
%                 pause(0.01)
            end
        end
        
        if any(rayX > volumeSize) | any(rayX < 1)
           go = 0; 
        end
    end
    
    rayPath(isnan(rayPath(:,1)),:) = [];
    
    if exitStep == 0
        exitStep = currentStep;
    end
    
    if exitStep > 1
        plot3(rayPath(1:exitStep,1), rayPath(1:exitStep,2), rayPath(1:exitStep,3), 'r', 'linewidth', 1.5);

        plot3(rayPath(exitStep+1:end,1), rayPath(exitStep+1:end,2), rayPath(exitStep+1:end,3), 'w');
        
        line([1 1]*rayPath(1,1), [1 1]*rayPath(1,2), ...
            [-30 0]*voxelSize, 'color', 'w');  
        
        line([0 rayT(1)*30*voxelSize]+rayPath(end,1), [0 rayT(2)*30*voxelSize]+rayPath(end,2), ...
            [0 rayT(3)*30*voxelSize]+rayPath(end,3), 'color', 'w')    
    
    end

    %line(rayPath([1 end],1), rayPath([1 end],2), rayPath([end],3)*[1 1], 'color', 'r')
    
    pathDifference(i,1) = sqrt((rayPath(1,1)-rayPath(end,1))^2 + (rayPath(1,2)-rayPath(end,2))^2);
    
    pathDifference(i,2) = rayPath(end,1) - rayPath(1,1);
    
    pathDifference(i,3) = rayPath(end,2) - rayPath(1,2);
    
    pause(0.01)
end

set(gca, 'Clipping', 'off'); axis off


subplot(1,3,2);

plot((rayOrigins(:,1)-mean(volCoords(topInds,1)))*1000, pathDifference(:,1)*1000);

xlabel('Radial distance from slab centre (um)')

ylabel('Radial change in ray position  (um)')

line(min(volCoords(topInds,1) - mean(volCoords(topInds,1)))*1000*[1 1], [0 0.5], 'color', 'g')

line(max(volCoords(topInds,1) - mean(volCoords(topInds,1)))*1000*[1 1], [0 0.55], 'color', 'g')

% get psuedointerference pattern
subplot(1,3,3);

xRef = (1+5:0.01:volumeSize(1)-5)*voxelSize;

finalPoints = interp1(rayOrigins(:,1), (pathDifference(:,2)+rayOrigins(:,1)), xRef, 'linear', 100);

xRef = (1+5:1:volumeSize(1)-5)*voxelSize;

[n,x] = hist(finalPoints, xRef);

n(1) = 100; n(end) = 100;

plot((xRef-mean(volCoords(topInds,1)))*1000, n/100);

line(min(volCoords(topInds,1) - mean(volCoords(topInds,1)))*1000*[1 1], [0 1], 'color', 'g')

line(max(volCoords(topInds,1) - mean(volCoords(topInds,1)))*1000*[1 1], [0 1], 'color', 'g')

ylabel('Intensity')

xlabel('Radial distance from slab centre (um)')


figure; 

subplot(2,2,1); plot(pathDifference(:,1));

subplot(2,2,2); plot(rayOrigins(:,1), pathDifference(:,1));

subplot(2,2,3); plot(rayOrigins(:,1), pathDifference(:,2));

subplot(2,2,4); plot(rayOrigins(:,1), pathDifference(:,3));


min(pathDifference(:,2))