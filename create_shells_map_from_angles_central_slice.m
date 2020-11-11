%%% Full map is too large to work on 200, 
%%% so split to seperate file to save slices for FDTD

%%% Now adding processing based on single file at a time.

clear; clc; close all

%% Set up parameters
dataFolder = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/Data/liquid angle csv/200nm_full';

convertAngleConvention = 0;

% Check conversion is set correctly
if isempty(strfind(dataFolder, 'thick')) & ~convertAngleConvention
    warning('Check if convertAngleConvention should be set')  

elseif ~isempty(strfind(dataFolder, 'thick')) & convertAngleConvention
    warning('Check if convertAngleConvention should be unset')
    
end

voxelSize = 200/10^6; %mm

nPerp = 1.53; % Short axis

nPara = 1.51;

nDif = nPerp - nPara;

%% Set up file loading
% Pulled ou of loadcsvstack.m as needs to be done singularly.

currentDirectory = pwd;          

thetaDirectory = sprintf('%s/theta', dataFolder);

phiDirectory = sprintf('%s/phi', dataFolder);

% Find which files to load - for theta
cd(thetaDirectory); directoryContents = dir; thetaIncludedFiles = [];

% natural langueg sort then remove from cells
thetaDirectoryNames = natsortfiles({directoryContents.name});

for i = 1:length(directoryContents)
    if strfind(thetaDirectoryNames{i}, 'csv')   

        thetaIncludedFiles = [ thetaIncludedFiles i ];
    end
end

% Find which files to load - for theta
cd(phiDirectory); directoryContents = dir; phiIncludedFiles = [];

% natural langueg sort then remove from cells
phiDirectoryNames = natsortfiles({directoryContents.name});

for i = 1:length(directoryContents)
    if strfind(phiDirectoryNames{i}, 'csv')   

        phiIncludedFiles = [ phiIncludedFiles i ];
    end
end

if length(thetaIncludedFiles) ~= length(phiIncludedFiles)
    error('Different number of files for theta and phi')
end

nImages = length(thetaIncludedFiles);

%% Run everything on single sliced

% Load test image
[tempImage] = readmatrix(thetaDirectoryNames{thetaIncludedFiles(1)}, 'TreatAsMissing', '-');
            
% Make square
size2Use = max(size(tempImage));
            
volumeSize = [size2Use, size2Use, nImages];

%%% Need to alocate arrays for saving data here.
phiCentralSlice = zeros(size2Use, nImages);

thetaCentralSlice = zeros(size2Use, nImages);

angle2FiberSlice = zeros(size2Use, nImages);

directRISlice = zeros(size2Use, nImages);

ringSlice = zeros(size2Use, nImages);

averagedRISlice = zeros(size2Use, nImages);

interpolatedRISlice = zeros(size2Use, nImages);

for iImage = 1:nImages
    % Load theta
    thetaSlice = zeros(size2Use, size2Use)*NaN;
    
    cd(thetaDirectory);
    
    tempImage = readmatrix(thetaDirectoryNames{thetaIncludedFiles(iImage)}, 'TreatAsMissing', '-');
    
    thetaSlice(size2Use-size(tempImage,1)+1:end,size2Use-size(tempImage,2)+1:end) = tempImage;
    
    % Load phi
    phiSlice = zeros(size2Use, size2Use)*NaN;
    
    cd(phiDirectory);
    
    tempImage = readmatrix(phiDirectoryNames{phiIncludedFiles(iImage)}, 'TreatAsMissing', '-');
    
    phiSlice(size2Use-size(tempImage,1)+1:end,size2Use-size(tempImage,2)+1:end) = tempImage;
    
    if convertAngleConvention 
        warning('May need to rethink later angle corrections if corrected convenction used')

        temp = phiSlice;

        thetaSlice = phiSlice;

        phiSlice = temp;

        clear temp;
    end
    
    % some Nan on theta border where phi has value, replace with 0
    thetaSlice(isnan(thetaSlice) & phiSlice == 0) = 0;

    % convert to radians
    thetaSlice = thetaSlice/180*pi;

    phiSlice = phiSlice/180*pi;

    if iImage == 1
        figure; 
        subplot(1,2,1);
        imshow(phiSlice/2/pi); title('phi');

        subplot(1,2,2);
        imshow((thetaSlice+pi/2)/pi); title('theta');
    end
    
    % Store central line values
    thetaCentralSlice(:,iImage) = thetaSlice(volumeSize(1)/2,:);

    phiCentralSlice(:,iImage) = phiSlice(volumeSize(1)/2,:);
    
    % Convert map directly to RI (Given rays along optical axis)
    outSideValue = NaN;

    goodSlice = ~isnan(phiSlice); 

    goodInds = find(goodSlice);

    % Calculate vector components
    fiberXSlice = cos(phiSlice).*cos(thetaSlice);

    fiberYSlice = sin(phiSlice).*cos(thetaSlice);

    fiberZSlice = sin(thetaSlice);

    % Flip z if pointing down
    downZInds = find(fiberZSlice < 0);

    fiberXSlice(downZInds) = -fiberXSlice(downZInds); 
    fiberYSlice(downZInds) = -fiberYSlice(downZInds);
    fiberZSlice(downZInds) = -fiberZSlice(downZInds); 

    directRI = ones(size2Use, size2Use)*outSideValue;

    % For on axis light.
    angle2Fiber = acos( (fiberXSlice(goodInds)*0 + fiberYSlice(goodInds)*0 + fiberZSlice(goodInds)*1) ./ ...
            (sqrt(fiberXSlice(goodInds).^2 + fiberYSlice(goodInds).^2 + fiberZSlice(goodInds).^2) * norm([0 0 1])) );

    % Adjusted range compared to original
    directRI(goodInds) = nPara + sin(angle2Fiber*2-pi/2)*nDif/2 + nDif/2;
    
    % Store central line values
    tempSlice = zeros(size2Use, size2Use);
    
    tempSlice(goodInds) = angle2Fiber;
    
    angle2FiberSlice(:,iImage) = tempSlice(volumeSize(1)/2,:);
    
    directRISlice(:,iImage) = directRI(volumeSize(1)/2,:);
    
    % Split into rings
    phiMask = phiSlice;
    
    phiMask = isnan(phiMask);

    ringLayer = zeros(volumeSize(1:2));

    minVal = 0;

    ringNum = 1;

    % Dilate in while some phi not equal to zero
    while any(phiMask(:) == 0)
        % Intersection of mask with phi values 
        tempSlice = imdilate(phiMask, strel('Disk', 1)) & ~isnan(phiSlice);

        % Avoid accidentally adding pieces from next ring
        if minVal > 0
          tempSlice = tempSlice & phiSlice >= minVal;
        end

        tempInds = find(tempSlice);  

        if ~isempty(tempInds)
          ringVals = unique(phiSlice(tempInds));

          if min(ringVals) > 0 & minVal == 0
            minVal = min(ringVals);
          end

          phiMask(tempInds) = 1;

          phiSlice(tempInds) = NaN;

          ringLayer(tempInds) = ringNum;

          [tempX, tempY] = ind2sub(volumeSize(1:2), tempInds);

        else
         % Click to next ring 
         minVal = 0;

         ringNum = ringNum + 1;
        end
    end

    % Put central row into slice collector
    ringSlice(:,iImage) = ringLayer(volumeSize(1)/2,:);

    % Average over RI rings on slice
    tempRILine = directRISlice(:, iImage);
    
    tempRILine(tempRILine == -1) = NaN;

    % Set up interpolate
    numRings = length( unique( ringSlice(ringSlice(:,iImage) > 0,iImage)));

    % Store radial post and intesnity for interpolation
    interpPoints = zeros(numRings*2-1,2);

    for jRing = 1:numRings
        % 1d, so inds are also y pos
        yInds = find(ringSlice(:, iImage) == jRing);

        % Check if ring doesn't cross middle
        if ~any(yInds == volumeSize(2)/2)
            % put points on both sides
            indsLeft = yInds(yInds < volumeSize(2)/2);

            interpPoints(jRing, :) = [mean(indsLeft) nanmean(tempRILine(indsLeft))];

            indsRight = yInds(yInds > volumeSize(2)/2);

            interpPoints(numRings*2 - 1 - jRing + 1, :) = [mean(indsRight) nanmean(tempRILine(indsRight))];
        else
            % no break, put in center
            interpPoints(numRings, :) = [mean(yInds) nanmean(tempRILine(yInds))];
        end

        tempRILine(yInds) = nanmean(tempRILine(yInds));
    end

    averagedRISlice(:, iImage) = tempRILine;

    % make interpolated slice
    yInds = find(ringSlice(:, iImage) > 0);

    interpolatedRISlice(yInds, iImage) = interp1(interpPoints(:,1) , interpPoints(:,2), yInds, 'linear', interpPoints(1,2));
end  
%% Plot and save    

% Plotting - phi and theta
figure; 
subplot(1,2,1);
imshow(flipud(phiCentralSlice')/2/pi); title('phi')

subplot(1,2,2);
imshow((flipud(thetaCentralSlice')+pi/2)/pi); title ('theta')

% Potting - fiber angles and direct RI conversion
% Test angle range
goodInds = find(directRISlice > 0);  

[tempVal, tempInd] = unique(directRISlice(goodInds));

figure; 
plot(angle2FiberSlice(goodInds(tempInd))/pi*180, tempVal, 'o'); hold on
plot(0:5:90, nPara+sin(((0:5:90)/180*pi)*2-pi/2)*nDif/2 + nDif/2, 'x')
title('Compare RI vs. angle range')

figure;
imshow(flipud(directRISlice'-1.50)/0.03)
line([volumeSize(1)/2 volumeSize(1)/2], [1 volumeSize(3)], 'color', 'g')
title('RI by voxel')

% Plot - coloured rings
cols = lines(100);

maxRings = max(ringSlice(:));

figure; subplot(1,2,1);
imshow(flipud(ringSlice')/maxRings)

subplot(1,2,2); hold on
imshow(flipud(phiCentralSlice')/2/pi)

for iRing = 1:maxRings
    inds = find(ringSlice == iRing);

    [tempX, tempY] = ind2sub(volumeSize(1:2), inds);

    plot(tempX, -tempY+volumeSize(3), '.', 'color', cols(iRing, :))
end

% Plot - Averaged and interpolated volumes
figure;
subplot(1,2,1); hold on
imshow(flipud(averagedRISlice'-1.51)/0.02)
line([volumeSize(2)/2 volumeSize(2)/2], [1 volumeSize(3)], 'color', 'g')

subplot(1,2,2); hold on
imshow(flipud(interpolatedRISlice'-1.51)/0.02)
line([volumeSize(2)/2 volumeSize(2)/2], [1 volumeSize(3)], 'color', 'g')

% Saving
currentDirectory = pwd; 
cd('/Users/gavintaylor/Desktop')

warning('check names are correct')

directRISlice(isnan(directRISlice)) = -1;

writematrix(directRISlice,'200_cone_full_section_lamellae.csv') 

averagedRISlice(isnan(averagedRISlice)) = -1;

writematrix(averagedRISlice,'200_cone_full_section_averaged.csv') 

interpolatedRISlice(isnan(interpolatedRISlice)) = -1;

writematrix(interpolatedRISlice,'200_cone_full_section_interpolated.csv') 
cd(pwd)