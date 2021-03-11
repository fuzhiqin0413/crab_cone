clear
clc
close all

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
for iStep = 1:15
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

% figure; hold on; axis equal

%%% may need to do something extra to tidy up border shapes
    %%% Clean volume, fill gaps cut through on border.
    %%% Make cut around edge that constrains shape.
    
%%% Messing around in this section to get reliable centre line
    % PCA isn't great, neither is pixel centre to tip.
    % Centre line from Skeleton is better
    % Fitting ellipse seems to be quite good, unless there is a lot of
    % material hanging of side, in which case a ball is made...
        % Should be easy to flag
    
grid3D.nx = volumeSize(1); grid3D.ny = volumeSize(2); grid3D.nz = volumeSize(3);
grid3D.minBound = [1 1 1]';
grid3D.maxBound = volumeSize';
                
% Turn off for backslash in ellipsoid fit
warning('off', 'MATLAB:nearlySingularMatrix')

for iTip = 400:20:600; %1:numTips;

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
    
    % Take correct intersect and flip vector if needed
    if hyperboloid
        mainVec = -mainVec;
        
        % hyperboloid centre is outside, so first intersect will be tip
        intersectInds = intersectInds(1);
    else
        % ellipsoid centre is inside, so last intersect will be tip
        intersectInds = intersectInds(end);
    end
    
    intersectCoord = [intersectX(intersectInds), intersectY(intersectInds), intersectZ(intersectInds)];

    % To plot
    figure; axis equal; hold on

    % Test in original coords
%     plot3(coneCoords(:,1), coneCoords(:,2), coneCoords(:,3), '.')
% 
%     plot3(tipCoords(iTip,1), tipCoords(iTip,2), tipCoords(iTip,3), 'go')
%     
%     plot3(intersectCoord(1), intersectCoord(2), intersectCoord(3), 'rx')
    
    % Test ellipsoid center and vector
%     plot3(center(1), center(2), center(3), 'rx')
%     
%     line([0 mainVec(1)*100] + center(1), [0 mainVec(2)*100] + center(2),...
%         [0 mainVec(3)*100] + center(3))
    
      % Plut full ellipsoid
%     mind = min( coneCoords );
%     maxd = max( coneCoords );
%     nsteps = 50;
%     step = ( maxd - mind ) / nsteps;
%     [ x, y, z ] = meshgrid( linspace( mind(1) - step(1), maxd(1) + step(1), nsteps ), linspace( mind(2) - step(2), maxd(2) + step(2), nsteps ), linspace( mind(3) - step(3), maxd(3) + step(3), nsteps ) );
% 
%     Ellipsoid = v(1) *x.*x +   v(2) * y.*y + v(3) * z.*z + ...
%               2*v(4) *x.*y + 2*v(5)*x.*z + 2*v(6) * y.*z + ...
%               2*v(7) *x    + 2*v(8)*y    + 2*v(9) * z;
%     p = patch( isosurface( x, y, z, Ellipsoid, -v(10) ) );

    % rotate cone
    coneCoords = coneCoords - intersectCoord;
    
    % Limit to those quite close
%     coneCoords(sqrt(coneCoords(:,1).^2 + coneCoords(:,2).^2 + coneCoords(:,3).^2) > 15, :) = [];

    elevationAngle = acos(mainVec(3)/norm(mainVec));
    
    azimuthAngle = atan2(mainVec(2), mainVec(1));
    
    coneCoords = coneCoords*vrrotvec2mat([0 0 1 azimuthAngle])*vrrotvec2mat([0 1 0 elevationAngle]);
    
    plot3(coneCoords(:,1), coneCoords(:,2), coneCoords(:,3), '.')

    % Looking at taking centre line on subvolume
    %%% Branches could be a problem as border of neighboring cone often included. 
    %%% Also tends to wander at tip rather than staying straight
    
%     coneBounds = [min(coneCoords(:,1))-1, min(coneCoords(:,2))-1, min(coneCoords(:,3))-1; ...
%                   max(coneCoords(:,1))+1, max(coneCoords(:,2))+1, max(coneCoords(:,3))+1];
%     
%     coneVolume = dataVolume(coneBounds(1,1):coneBounds(2,1), coneBounds(1,2):coneBounds(2,2), ...
%         coneBounds(1,3):coneBounds(2,3));          
%     
%     tempIndex = find(coneVolume);
% 
%     nIndex = length(tempIndex); tempCoords = zeros(nIndex, 3);
% 
%     [tempCoords(:,1), tempCoords(:,2), tempCoords(:,3)] = ...
%         ind2sub(size(coneVolume), tempIndex);
%     
% %     figure; axis equal
% %     plot3(tempCoords(:,1), tempCoords(:,2), tempCoords(:,3), '.')
% 
%     S = skeleton(coneVolume);
%     
%     figure,
%   FV = isosurface(coneVolume,0.5);
%   patch(FV,'facecolor',[1 0 0],'facealpha',0.3,'edgecolor','none');
%   view(3)
%   camlight
% % Display the skeleton
%   hold on;
%   for i=1:length(S)
%     L=S{i};
%     plot3(L(:,2),L(:,1),L(:,3),'-','Color',rand(1,3));
%   end

    
  % Playing with PCA and volume center
%     coneCoords = labelSurfaceCoords(inds,:) - tipCoords(iTip,:);
%     
%     % Limit to those quite close
%     %coneCoords(sqrt(coneCoords(:,1).^2 + coneCoords(:,2).^2 + coneCoords(:,3).^2) > 15, :) = [];
%     
%     coneCenter = mean(coneCoords);
% 
%     %%% Tried taking centreLine with PCA but didn't work well for open shape
%     
%     elevationAngle = acos(coneCenter(3)/norm(coneCenter));
%     
%     azimuthAngle = atan2(coneCenter(2), coneCenter(1));
%     
%     coneCoords = coneCoords*vrrotvec2mat([0 0 1 azimuthAngle])*vrrotvec2mat([0 1 0 elevationAngle]);
%     
%     coneCenter = coneCenter*vrrotvec2mat([0 0 1 azimuthAngle])*vrrotvec2mat([0 1 0 elevationAngle]);
%     
%     %plot3(coneCoords(:,1), coneCoords(:,2), coneCoords(:,3), '.')
% 
% %     line([0 coneCenter(1)], [0 coneCenter(2)],... 
% %         [0 coneCenter(3)], 'color', 'r')
end

warning('on', 'MATLAB:nearlySingularMatrix')

%% Distance based approach didn't work well as caught borders of adjacent cones
if 0
    %% Take average distance to nearest tip neighbor

    tipDists = zeros(numTips,1);

    for iTip = 1:numTips
        tempDists = sqrt((tipCoords(:,1) - tipCoords(iTip, 1)).^2 + (tipCoords(:,2) - tipCoords(iTip, 2)).^2 + ...
            (tipCoords(:,3) - tipCoords(iTip, 3)).^2);

        % Remove self
        tempDists(tempDists == 0) = [];

        tipDists(iTip) = min(tempDists);
    end

    %% Get points out to neighbor
    %%% Should test for connectivity later

    figure; hold on; axis equal
    plot3(tipCoords(:,1), tipCoords(:,2), tipCoords(:,3), 'x')

    for iTip = 1:numTips
        tempDists = sqrt((coneSurfaceCoords(:,1) - tipCoords(iTip, 1)).^2 + (coneSurfaceCoords(:,2) - tipCoords(iTip, 2)).^2 + ...
            (coneSurfaceCoords(:,3) - tipCoords(iTip, 3)).^2);

        closeInds = find(tempDists < tipDists(iTip)/2);

        plot3(tipCoords(iTip,1), tipCoords(iTip,2), tipCoords(iTip,3), 'or')
        plot3(coneSurfaceCoords(closeInds,1), coneSurfaceCoords(closeInds,2), coneSurfaceCoords(closeInds,3), '.')
    end
end