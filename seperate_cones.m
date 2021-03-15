clear
clc
close all

voxSize = 6.17;

% Load in ct data
dataFolder = '/Users/gavintaylor/Documents/Shared VM Folder/CT labels for Gavin as tif/LP2_matlab/Labels';
dataVolume = loadtiffstack(dataFolder, 1);

volumeSize = size(dataVolume);

% Split Control point location (@3) and Absolute coordinate location (@6)
% from ascii path set. Actual Path is not very important as we just want
% the points we clicked

points = load('/Users/gavintaylor/Documents/Shared VM Folder/CT labels for Gavin as tif/LP2_matlab/ControlPoints.pathset');
pointsCoords = load('/Users/gavintaylor/Documents/Shared VM Folder/CT labels for Gavin as tif/LP2_matlab/Coordinates.pathset');

%%% Flip coords 1 and 2 to match loaded volume

pointsCoords = pointsCoords(:, [2 1 3]);

if size(points, 1) ~= size(pointsCoords,1)
    error('Length of path files dont match')
end

tipInds = find(points);
numTips = length(tipInds);

tipCoords = pointsCoords(tipInds,:);

figure; hold on; axis equal
plot3(pointsCoords(:,1), pointsCoords(:,2), pointsCoords(:,3), '.')
plot3(tipCoords(:,1), tipCoords(:,2), tipCoords(:,3), 'x')

% Get cones exterior voxels
coneExterior = logical(~dataVolume);

% Exterior is outer volume, and grown into other volume.
STREL_6_CONNECTED = strel('sphere', 1); 
temp = STREL_6_CONNECTED.Neighborhood; 
temp(:,:,2) = 1; temp([1 3],2,:) = 1; temp(2,[1 3],:) = 1;
STREL_18_CONNECTED = strel('arbitrary', temp); 

coneExterior = imdilate(coneExterior, STREL_18_CONNECTED);

coneExterior = coneExterior & dataVolume;

% Get surface coords
surfaceIndexList = find(coneExterior);

nIndex = length(surfaceIndexList); coneSurfaceCoords = zeros(nIndex, 3);

[coneSurfaceCoords(:,1), coneSurfaceCoords(:,2), coneSurfaceCoords(:,3)] = ...
    ind2sub(volumeSize, surfaceIndexList);

plot3(coneSurfaceCoords(:,1), coneSurfaceCoords(:,2), coneSurfaceCoords(:,3), '.')

%% Flood fill growth approach
% Two volumes, one with labels, one with availible
coneExteriorLabels = coneExterior*0;
coneExteriorAvailible = logical(coneExterior);

figure; hold on; axis equal
plot3(tipCoords(:,1), tipCoords(:,2), tipCoords(:,3), 'gx')

% Place label on nearest voxel of volume for each cone - some coords aren't
% exactly on voxels, so need to search
for iTip = 1:numTips
    tempDists = sqrt((coneSurfaceCoords(:,1) - tipCoords(iTip, 1)).^2 + (coneSurfaceCoords(:,2) - tipCoords(iTip, 2)).^2 + ...
        (coneSurfaceCoords(:,3) - tipCoords(iTip, 3)).^2);
    
    [~, closestInd] = min(tempDists);
    
    tipCoords(iTip, :) = coneSurfaceCoords(closestInd, :);
    
    coneExteriorLabels(coneSurfaceCoords(closestInd,1), coneSurfaceCoords(closestInd,2), ...
        coneSurfaceCoords(closestInd,3)) = iTip;
    
    coneExteriorAvailible(coneSurfaceCoords(closestInd,1), coneSurfaceCoords(closestInd,2), ...
        coneSurfaceCoords(closestInd,3)) = 0;
    
    plot3(coneSurfaceCoords(closestInd,1), coneSurfaceCoords(closestInd,2), ...
        coneSurfaceCoords(closestInd,3), 'ro')
end

% Grow in steps - set until good connection - 50 to far, 10 to short, 30 is a bit too much
% May need to clean volume more
extendLength = 15;
for iStep = 1:extendLength
    % And to exterior limits growth to surface voxels
    tic
    enlargedLabels = imdilate(coneExteriorLabels, STREL_18_CONNECTED);
    
    % Fast if done with logical and indexes
    enlargedSurface =  logical(enlargedLabels) & coneExteriorAvailible;
    
    newSurfInds = find(enlargedSurface);
    
    coneExteriorLabels(newSurfInds) = enlargedLabels(newSurfInds);
    
    coneExteriorAvailible(newSurfInds) = 0;
    toc
end

labelIndexList = find(coneExteriorLabels);

nIndex = length(labelIndexList); labelSurfaceCoords = zeros(nIndex, 3);

[labelSurfaceCoords(:,1), labelSurfaceCoords(:,2), labelSurfaceCoords(:,3)] = ...
    ind2sub(volumeSize, labelIndexList);

% plot3(labelSurfaceCoords(:,1), labelSurfaceCoords(:,2), labelSurfaceCoords(:,3), '.')

coneLabels = coneExteriorLabels(labelIndexList);

for iTip = 1:numTips
    inds = find(coneLabels == iTip);
    
    plot3(labelSurfaceCoords(inds,1), labelSurfaceCoords(inds,2), labelSurfaceCoords(inds,3), '.')
end

%% manipulate cones
coneCoordsCells = cell(numTips, 1);

fC = figure; hold on; axis equal

%%% may need to do something extra to tidy up border shapes
    %%% Clean volume, fill gaps cut through on border.
    %%% Make cut around edge that constrains shape.
    
%%% Messing around in this section to get reliable centre line
    % PCA isn't great, neither is pixel centre to tip.
    % Centre line from Skeleton is better
    % Fitting ellipse seems to be quite good, unless there is a lot of
    % material hanging of side, in which case a ball is made...
        % Should be easy to flag
        
    %%% Could also try point could registration done on iffy ones against avg. shape
        % Allow scale to accomodate size variance, but then just adjust orientation
    
grid3D.nx = volumeSize(1); grid3D.ny = volumeSize(2); grid3D.nz = volumeSize(3);
grid3D.minBound = [1 1 1]';
grid3D.maxBound = volumeSize';
                
% Turn off for backslash in ellipsoid fit
warning('off', 'MATLAB:nearlySingularMatrix')

xStep = 0;

yStep = 0;

offSet = 30;

radiiHist = zeros(numTips,1);

goodConeShapes = cell(numTips,1);

coneTipsCoords = zeros(numTips,3);

coneTipsAngles = zeros(numTips,2);

for iTip = 1:numTips;

    inds = find(coneLabels == iTip);
    
    coneCoords = labelSurfaceCoords(inds,:);
    
    %Fitting a an ellipsoid
    
    [ center, radii, evecs, v, chi2 ] = ellipsoid_fit( coneCoords );
    
    % hyperboloid has two negative radii
    if sum(radii < 0) == 2
        hyperboloid = 1;
    elseif sum(radii < 0) == 0   
        hyperboloid = 0;
    else
       error('Check what shape this is') 
    end
    
    [~, maxR] = max(radii);
    
    mainVec = evecs(:,maxR);
    
    otherRadii = [1 2 3];
    
    otherRadii(maxR) = [];
    
    radiiHist(iTip) = radii(maxR)/mean(radii(otherRadii));
    
    
    % Main vec can point either direction, 
    % check angle between vector from center vs line from centre to tip
    tipVec = tipCoords(iTip,:) - center';
    tipVec = tipVec/norm(tipVec);
    
    if sqrt((mainVec(1)-tipVec(1))^2 + (mainVec(2)-tipVec(2))^2 + (mainVec(3)-tipVec(3))^2) > ...
            sqrt(2)
       mainVec = -mainVec; 
    end
    
    % Take intersection to surface in volume
        % Potential issues if there are multiple intersections/fold backs..
    [intersectX, intersectY, intersectZ] = amanatideswooalgorithm_efficient(center, mainVec, grid3D, 0,...
                    [], [], 0);
    
    indexList = sub2ind(volumeSize, intersectX, intersectY, intersectZ);
    
    intersectInds = find(coneExteriorLabels(indexList) == iTip);
    
    if ~isempty(intersectInds)
        % Take correct intersect and flip vector if needed
        if hyperboloid
            mainVec = -mainVec;

            % hyperboloid centre is outside, so first intersect will be tip
            intersectInds = intersectInds(1);
            
            col = 'b';
            
            % These basically always seem to be well oriented
        else
            % ellipsoid centre is inside, so last intersect will be tip
            intersectInds = intersectInds(end);
            
            % Groups seem to be above 2.5 is good, above 2 is ok
            
            if abs(radiiHist(iTip)) >2.5
                col = 'b';
            elseif radiiHist(iTip) > 2
                col = 'g';
            else
                col = 'm';
            end
            
        end

        intersectCoord = [intersectX(intersectInds), intersectY(intersectInds), intersectZ(intersectInds)];

        % rotate cone
        coneCoords = coneCoords - intersectCoord;

        % Limit to those quite close
    %     coneCoords(sqrt(coneCoords(:,1).^2 + coneCoords(:,2).^2 + coneCoords(:,3).^2) > 15, :) = [];

        elevationAngle = acos(mainVec(3)/norm(mainVec));

        azimuthAngle = atan2(mainVec(2), mainVec(1));

        coneCoords = coneCoords*vrrotvec2mat([0 0 1 azimuthAngle])*vrrotvec2mat([0 1 0 elevationAngle]);

        %figure(fC)
        subplot(1,2,1); hold on; axis equal
        plot3(coneCoords(:,1) + xStep*offSet, coneCoords(:,2), coneCoords(:,3) + yStep*offSet, '.', 'color', col)
        
        subplot(1,2,2); hold on
        if col == 'b'
            plot3(coneCoords(:,1)*voxSize, coneCoords(:,2)*voxSize, coneCoords(:,3)*voxSize, '.')
            
            goodConeShapes{iTip} = coneCoords;
            
            coneTipsCoords(iTip,:) = intersectCoord;

            coneTipsAngles(iTip,:) = [azimuthAngle, elevationAngle];

        end
    else
        %%% Happens on very short cone
        coneCoords = coneCoords - tipCoords(iTip,:);
        
        subplot(1,2,1); hold on; axis equal
        plot3(coneCoords(:,1) + xStep*offSet, coneCoords(:,2), coneCoords(:,3) + yStep*offSet, '.', 'color', 'r')
        
        % To plot
%         figure; axis equal; hold on
% 
%         %Test in original coords
%         plot3(coneCoords(:,1), coneCoords(:,2), coneCoords(:,3), '.')
% 
%         plot3(tipCoords(iTip,1), tipCoords(iTip,2), tipCoords(iTip,3), 'go')
% 
%         plot3(intersectCoord(1), intersectCoord(2), intersectCoord(3), 'rx')
% 
%         %Test ellipsoid center and vector
%         plot3(center(1), center(2), center(3), 'rx')
% 
%         line([0 mainVec(1)*100] + center(1), [0 mainVec(2)*100] + center(2),...
%             [0 mainVec(3)*100] + center(3))
% 
%           % Plut full ellipsoid
%         mind = min( coneCoords );
%         maxd = max( coneCoords );
%         nsteps = 50;
%         step = ( maxd - mind ) / nsteps;
%         [ x, y, z ] = meshgrid( linspace( mind(1) - step(1), maxd(1) + step(1), nsteps ), linspace( mind(2) - step(2), maxd(2) + step(2), nsteps ), linspace( mind(3) - step(3), maxd(3) + step(3), nsteps ) );
% 
%         Ellipsoid = v(1) *x.*x +   v(2) * y.*y + v(3) * z.*z + ...
%                   2*v(4) *x.*y + 2*v(5)*x.*z + 2*v(6) * y.*z + ...
%                   2*v(7) *x    + 2*v(8)*y    + 2*v(9) * z;
%         p = patch( isosurface( x, y, z, Ellipsoid, -v(10) ) ); 
    end

    xStep = xStep + 1;
        
    if xStep > 30
        xStep = 0;

        yStep = yStep + 1;
    end

end

warning('on', 'MATLAB:nearlySingularMatrix')

figure;
hist(radiiHist, 400)

%% Do shape analysis - 3D

distStep = 1;
distStepNumber = 10;

angleStepNumber = 10;
angleStep = 360/angleStepNumber;

interpolatedConeRadius = zeros(numTips, distStepNumber*angleStepNumber)*NaN;

% distVals = (1:distStepNumber)*distStep;
% Make finer at start
distVals = [0.25 0.5 1 2 3 4 6 8 10 12];

[distInterp, angleInterp] = meshgrid(distVals, (1:angleStepNumber)*angleStep);

distInterp = -distInterp(:);
angleInterp = angleInterp(:)/180*pi - pi;

figure; axis equal

toRemove = [];

for iTip = 1:numTips
    
    tempCone = goodConeShapes{iTip};
    if ~isempty(tempCone)
        % convert x/y to radial cords
        
%         distValue = sqrt(tempCone(:,1).^2 + tempCone(:,2).^2 + tempCone(:,3).^2);
        
        radialValue = sqrt(tempCone(:,1).^2 + tempCone(:,2).^2);
        
        angleValue = atan2(tempCone(:,2), tempCone(:,1));
        
        % Wont cross over +/- pi wrap around, could adjust for this.
        
%         % make interpolant - works quite well, but kind of misses top   
        interpRadius = scatteredInterpolant(tempCone(:,3), angleValue, radialValue, 'linear', 'nearest');
        
        tempRadius = interpRadius(distInterp, angleInterp);

%       plot3(tempCone(:,1)*voxSize, tempCone(:,2)*voxSize, tempCone(:,3)*voxSize,'b.')
        hold on
        
        plot3(tempRadius.*cos(angleInterp)*voxSize, tempRadius.*sin(angleInterp)*voxSize, distInterp*voxSize, 'ro') 

        interpolatedConeRadius(iTip,:) = tempRadius;
    else
       toRemove = [toRemove iTip]; 
    end
end
 
interpolatedConeRadius(toRemove, :) = [];
%% Plot modes - 3D
[radiusModes,~,~,~,explained] = pca(interpolatedConeRadius);

meanRadius = mean(interpolatedConeRadius);
stdRadius = std(interpolatedConeRadius);

figure;
subplot(4,2,1); axis equal; 
plot3(meanRadius'.*cos(angleInterp)*voxSize, meanRadius'.*sin(angleInterp)*voxSize, distInterp*voxSize, '.g')
hold on;

plusRadius = meanRadius + stdRadius;
minusRadius = meanRadius - stdRadius;

plot3(minusRadius'.*cos(angleInterp)*voxSize, minusRadius'.*sin(angleInterp)*voxSize, distInterp*voxSize, '.b')
plot3(plusRadius'.*cos(angleInterp)*voxSize, plusRadius'.*sin(angleInterp)*voxSize, distInterp*voxSize, '.r')
    
title('Mean and SD')
legend('Mean shape', 'Negative variation', 'Positive variation')

subplot(4,2,2)
plot(cumsum(explained), 'x')
xlabel('Mode')
ylabel('Total variance explained (%)')
ylim([0 100])

% 1st mode
for iMode = 1:6
    subplot(4,2,2 + iMode); axis equal; 
    plusRadius = meanRadius + radiusModes(:, iMode)'*10;
    minusRadius = meanRadius - radiusModes(:, iMode)'*10;

    plot3(meanRadius'.*cos(angleInterp)*voxSize, meanRadius'.*sin(angleInterp)*voxSize, distInterp*voxSize, '.g')
    hold on
    plot3(minusRadius'.*cos(angleInterp)*voxSize, minusRadius'.*sin(angleInterp)*voxSize, distInterp*voxSize, '.b')
    plot3(plusRadius'.*cos(angleInterp)*voxSize, plusRadius'.*sin(angleInterp)*voxSize, distInterp*voxSize, '.r')

    title(sprintf('Mode %i', iMode))
end

%% Cone shape analysis - 2D

% Rough approach. Round X and take all voxels zero, then rotate

angleStep_2D = 30;
% Check value gives clean steps
angleInterp_2D = 0:angleStep_2D:(180-angleStep_2D)

for iTip = 1:numTips
    
    %%% Try doing on surface image directly
    if all(coneTipsCoords(iTip,:) > 0)
        % set new bounds
        minCoords = coneTipsCoords(iTip,:)-(extendLength+10);
        minCoords(minCoords < 1) = 1;

        maxCoords = coneTipsCoords(iTip,:)+(extendLength+10);
        maxCoords(maxCoords > volumeSize) = volumeSize(maxCoords > volumeSize);

        % Take subvolume
        coneSubvolume = coneExteriorLabels(minCoords(1):maxCoords(1), minCoords(2):maxCoords(2), ...
            minCoords(3):maxCoords(3));

        % Remove other stuff
        coneSubvolume(coneSubvolume ~= iTip) = 0;

        coneSubvolume = logical(coneSubvolume);
        
        % Generally centered at tip 
        if any(size(coneSubvolume) ~= (extendLength+10)*2+1)
            % Catch if not
            error('Adjust centre')
        end
        temp = zeros(4,4); temp(4,4) = 1;
        
        M2 = temp; M2(1:3,1:3) = vrrotvec2mat([0 0 1 coneTipsAngles(iTip,1)]);

        M3 = temp; M3(1:3,1:3) = vrrotvec2mat([0 1 0 coneTipsAngles(iTip,2)]);
        
        coneSubvolume = logical( affine_transform_full(single(coneSubvolume), M2*M3, 5));

        % Note small holes appear after affine transform, close fixes these. 
        % Probably bug with nearest neighbour interp I added to affine transform c file
            %%% Note close may not fix surface well
            % could dilate twice, transform, close, then erode twice
%         coneSubvolume = imclose(coneSubvolume, STREL_6_CONNECTED);
        
        baseCone = goodConeShapes{iTip};
        
        figure; 
        
        %plot3(baseCone(:,1), baseCone(:,2), baseCone(:,3), 'rx')
        
        for jAngle = 1:length(angleInterp_2D)
           % rotate volume
           M = temp; M(1:3,1:3) = vrrotvec2mat([0 0 1 angleInterp_2D(jAngle)/180*pi]);
           
           rotSubvolume = logical( affine_transform_full(single(coneSubvolume), M, 5));
           
           % Maybe best to do close now, otherwise it's done twice 
           rotSubvolume = imclose(rotSubvolume, STREL_6_CONNECTED);
           
           % get coords
           rotInds = find(rotSubvolume);
           
           rotCoords = zeros(length(rotInds), 3);
           
           [rotCoords(:,1), rotCoords(:,2), rotCoords(:,3)] = ...
                ind2sub(size(coneSubvolume), rotInds);
            
           xAxisInds = find(rotCoords(:,2) == 26);
           
           %plot3(rotCoords(:,1)-26, rotCoords(:,2)-26, rotCoords(:,3)-26, '.') 
           [~, sInds] = sort(rotCoords(xAxisInds,1));
           
           subplot(2,3,jAngle); 
           imshow(permute(rotSubvolume(:,26,:), [1 3 2])); hold on
           plot(26, 26, 'rx')
           plot(rotCoords(xAxisInds(sInds),3), rotCoords(xAxisInds(sInds),1), '-')
        end
        %%% Take points on zero axis, then rotate and do again --- how does
        %%% it compare?
    end
    
    
    %%% Rotated coordiantes
%     baseCone = goodConeShapes{iTip};
%     
%     figure; 
%     
%     for jAngle = 1:length(angleInterp_2D)
%         rotCone = baseCone*vrrotvec2mat([0 0 1 angleInterp_2D(jAngle)/180*pi]);
%         
%         % Round Y and get zero values (on X axis)     
%         xAxisInds = find(round(rotCone(:, 2)) == 0);
%         % Seems this leaves a lot of oversampling - try 1D interp
%         
%         subplot(1,3,1); hold on
%         plot3(baseCone(xAxisInds,1), baseCone(xAxisInds,2), baseCone(xAxisInds,3), 'o');
%         
%         subplot(1,3,2); hold on
%         [~, sInds] = sort(rotCone(xAxisInds,1));
%         plot(rotCone(xAxisInds(sInds),1), rotCone(xAxisInds(sInds),3))
%         
%         % Test interpolate for smoother sampling
%         % 1d can't have duplicate data points
%             %%% Suffers from disc artifacts...
% %         uniqueX = unique(rotCone(xAxisInds,1));
% %         Zvals = zeros(length(uniqueX),1);
% %         
% %         for kX = 1:length(uniqueX)
% %             vals2Use = rotCone(xAxisInds,1) == uniqueX(kX);
% %             Zvals(kX) = mean( rotCone(xAxisInds(vals2Use),3));
% %         end
% %         
% %         interpZ = interp1(uniqueX, Zvals, -5:5, 'linear');
% 
%         % Try scattered on 2D (w/ non rounded y vals - 1d Interp essentially uses all Y values at 0)
%             %%% This is not better, probably worse...
%         tempZInterp = scatteredInterpolant(rotCone(xAxisInds,1), rotCone(xAxisInds,2), rotCone(xAxisInds,3), 'linear', 'nearest');
% 
%         tempZValues = tempZInterp(-5:5, zeros(1,11));
%         
%         subplot(1,3,3); hold on
%         plot(-5:5, tempZValues)
% 
%     end
%     
%     subplot(1,3,1); hold on
%     plot3(baseCone(:,1), baseCone(:,2), baseCone(:,3), '.');
        
end