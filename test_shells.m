clear
clc
close all

% Load in angle data
dataFolder = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/Data/liquid angle csv/400nm_thick';

convertAngleConvention = 0;

% Convert to mm
voxelSize = 800/10^6;

scaleSteps = 800/voxelSize/10^6;

% Load theta angles
thetaVolume = loadcsvstack(sprintf('%s/theta', dataFolder), 0);

% Load phi angles
phiVolume = loadcsvstack(sprintf('%s/phi', dataFolder), 0);

if convertAngleConvention
    
    temp = thetaVolume;

    thetaVolume = phiVolume;

    phiVolume = temp;

    clear temp;
end

volumeSize = size(phiVolume);

phiVals = unique(phiVolume(~isnan(phiVolume)));

figure; 
subplot(1,2,1);
imshow(phiVolume(:,:,10)/360)

subplot(1,2,2);
imshow((thetaVolume(:,:,10)+90)/180)

%% Break up into shells by layer
close all

cols = lines(100);

ringVolume = zeros(volumeSize);

for iSlice = 1:volumeSize(3)
   phiSlice = phiVolume(:,:,iSlice);
   
   phiMask = phiSlice;
   phiMask = isnan(phiMask);
   
   ringLayer = zeros(volumeSize(1:2));
   
   minVal = 0;
   
   ringNum = 1;
   
%    figure; 
%    imshow(phiVolume(:,:,iSlice)/360); hold on
   
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

          %plot(tempX, tempY, '.', 'color', cols(ringNum, :))
      else
         % Click to next ring 
         minVal = 0;
         
         ringNum = ringNum + 1;
      end
   end
   
   ringVolume(:,:,iSlice) = ringLayer;
end

maxRings = max(ringVolume(:));

figure; subplot(1,3,1);
imshow(ringVolume(:,:,1)/maxRings)
subplot(1,3,2);
imshow(ringVolume(:,:,end)/maxRings)
subplot(1,3,3);
imshow( permute( ringVolume(volumeSize(1)/2,:,:), [2 3 1])/maxRings)

% Look at distribution of angles in ring
if 0
    for iSlice = 1; %:volumeSize(3)

        figure; hold on

        tempRingSlice = ringVolume(:,:,iSlice);

        tempThetaSlice = thetaVolume(:,:,iSlice);

        tempPhiSlice = phiVolume(:,:,iSlice);

        numRings = max(tempRingSlice(:));

        % get points in each ring
        for jRing = 1:numRings
           ringInds = find(tempRingSlice == jRing); 

           % get orientation within each ring
           for kPhi = 1:length(phiVals)

               phiInds = find(tempPhiSlice(ringInds) == phiVals(kPhi));

               [n, x] = hist(tempThetaSlice(ringInds(phiInds)), -90:10:90);

               if kPhi < 10
                   subplot(2,5,kPhi); hold on
                   plot(x, n, 'color', cols(jRing, :))
               end

           end
        end
    end
end
%% For plotting inclination of all pixels
if 0
    for iSlice = 1 %:volumeSize(3)
        figure; 
        subplot(1,2,1);
        imshow(phiVolume(:,:,iSlice)/360); hold on

        subplot(1,2,2);
        imshow((thetaVolume(:,:,iSlice)+90)/180); hold on

        tempRingSlice = ringVolume(:,:,iSlice);

        tempThetaSlice = thetaVolume(:,:,iSlice);

        tempPhiSlice = phiVolume(:,:,iSlice);

        numRings = max(tempRingSlice(:));

        % get points in each ring
        for jRing = 1 %:numRings
           ringInds = find(tempRingSlice == jRing); 

           [ringX, ringY] = ind2sub(volumeSize(1:2), ringInds);

           % Get zero orientation within each ring
           zeroInds = find(tempPhiSlice(ringInds) == phiVals(1));

           % Get orientation of each zero ind
           for kInd = 1:50:length(zeroInds)
               % firstly get nearby within thickness
                    % Less than phi vals should always get entire band...?
               nearbyInds = find(sqrt((ringX - ringX(zeroInds(kInd))).^2 + ...
                   (ringY - ringY(zeroInds(kInd))).^2) < length(phiVals));

               for mNearby = 1:length(nearbyInds)
                    tempPhiPoint = tempPhiSlice(ringInds(nearbyInds(mNearby)))/180*pi;

                    tempThetaPoint = tempThetaSlice(ringInds(nearbyInds(mNearby)))/180*pi;

                    % Border with zero phi has nan for theta.
                    if isnan(tempThetaPoint)
                       tempThetaPoint = 0; 
                    end

                    vecX = cos(tempPhiPoint)*cos(tempThetaPoint);

                    vecY = sin(tempPhiPoint)*cos(tempThetaPoint);

                    vecZ = sin(tempThetaPoint);

                    % if pointing down in Z, flip
                    if vecZ < 0
                        vecX = -vecX; vecY = -vecY; vecZ = -vecZ; 
                    end

                    % Check angle if vector is pointing towards center 
                    % (based on shorter distance at top of vector than base)

                    % Now use for flat as well
                    if 1; %vecZ ~= 0 
                        basedDist = sqrt((ringX(nearbyInds(mNearby)) - volumeSize(1)/2)^2 + ...
                            (ringY(nearbyInds(mNearby)) - volumeSize(2)/2)^2);

                        topDist = sqrt((ringX(nearbyInds(mNearby)) + vecX - volumeSize(1)/2)^2 + ...
                            (ringY(nearbyInds(mNearby)) + vecY - volumeSize(2)/2)^2);

                        if topDist > basedDist
                           % Roate by 180 degrees 

                           vecX = -vecX;
                           vecY = -vecY;

                           tempCols = 'r';
                        else
                            tempCols = 'b';
                        end
                    else
                       tempCols = cols(jRing, :);
                    end

                    subplot(1,2,1)
                    line( [0 4*vecX]+ringX(nearbyInds(mNearby)), [0 4*vecY]+ringY(nearbyInds(mNearby)),...
                        [0 4*vecZ], 'color', tempCols)

                    subplot(1,2,2)
                    line( [0 4*vecX]+ringX(nearbyInds(mNearby)), [0 4*vecY]+ringY(nearbyInds(mNearby)),...
                        [0 4*vecZ], 'color', tempCols) 

               end
           end
        end
    end
end

%% Get max inclination per band and plot
if 1
    volumeVectorX = zeros(volumeSize)*NaN;
    volumeVectorY = zeros(volumeSize)*NaN;   
    volumeVectorZ = zeros(volumeSize)*NaN;
    
    for iSlice = 1:volumeSize(3)

%         figure; 
%         imshow((thetaVolume(:,:,iSlice)+90)/180); hold on

        tempRingSlice = ringVolume(:,:,iSlice);

        tempThetaSlice = thetaVolume(:,:,iSlice);

        tempPhiSlice = phiVolume(:,:,iSlice);

        numRings = max(tempRingSlice(:));

        % get points in each ring
        for jRing = 1:numRings
           ringInds = find(tempRingSlice == jRing); 

           [ringX, ringY] = ind2sub(volumeSize(1:2), ringInds);

           % Get zero inclination within each ring (should be two bands)
           %zeroInds = find(tempThetaSlice(ringInds) == 0 | isnan(tempThetaSlice(ringInds)));
           zeroInds = find((tempPhiSlice(ringInds) == 0 | tempPhiSlice(ringInds) == 180 | ...
               tempPhiSlice(ringInds) == 360) & ...
               (ringX == volumeSize(1)/2 |ringY == volumeSize(2)/2));
           %%% Note cutting along middle of volume
           
           %plot(ringX(zeroInds), ringY(zeroInds), '.')
           
           % Get orientation of each zero ind
           for kInd = 1:length(zeroInds)
               % get nearby within half ring thickness, and that are closer
               % to center
               
               distToCenter = sqrt((ringX - volumeSize(1)/2).^2 + ...
                   (ringY - volumeSize(2)/2).^2);
               
               distToPoint = sqrt((ringX - ringX(zeroInds(kInd))).^2 + ...
                   (ringY - ringY(zeroInds(kInd))).^2);
               
               nearbyInds = find(distToPoint < length(phiVals)/2 & ...
                   distToCenter < distToCenter(zeroInds(kInd)));

               if ~isempty(nearbyInds)
                   % plot(ringX(nearbyInds), ringY(nearbyInds), '.', 'color', cols(jRing, :))

                   % Take point with largest theta to use
                   [~, pointToUse] = max( abs( tempThetaSlice(ringInds(nearbyInds)) ));

                   tempPhiPoint = tempPhiSlice(ringInds(nearbyInds(pointToUse)))/180*pi;

                   tempThetaPoint = tempThetaSlice(ringInds(nearbyInds(pointToUse)))/180*pi;

                   % Border with zero phi has nan for theta.
                    if isnan(tempThetaPoint)
                       tempThetaPoint = 0; 
                    end

                    % Will not work if completely vertical..
                    if tempThetaPoint == pi/2
                        tempThetaPoint = 85/180*pi;
                    elseif tempThetaPoint == -pi/2
                        tempThetaPoint = -85/180*pi;
                    end

                    vecX = cos(tempPhiPoint)*cos(tempThetaPoint);

                    vecY = sin(tempPhiPoint)*cos(tempThetaPoint);

                    vecZ = sin(tempThetaPoint);

                    % if pointing down in Z, flip
                    if vecZ < 0
                        vecX = -vecX; vecY = -vecY; vecZ = -vecZ; 
                    end

                    % Check angle if vector is pointing towards center 
                    % (based on shorter distance at top of vector than base)

                    % Now use for flat as well
                    if 1; %vecZ ~= 0
                        basedDist = sqrt((ringX(nearbyInds(pointToUse)) - volumeSize(1)/2)^2 + ...
                            (ringY(nearbyInds(pointToUse)) - volumeSize(2)/2)^2);

                        topDist = sqrt((ringX(nearbyInds(pointToUse)) + vecX - volumeSize(1)/2)^2 + ...
                            (ringY(nearbyInds(pointToUse)) + vecY - volumeSize(2)/2)^2);

                        if topDist > basedDist
                           % Rotate by 180 degrees 

                           vecX = -vecX;
                           vecY = -vecY;
                        end
                    end

    %                line( [0 4*vecX]+ringX(zeroInds(kInd)), [0 4*vecY]+ringY(zeroInds(kInd)),...
    %                         [0 4*vecZ], 'color', cols(jRing, :))

                   volumeVectorX(ringX(zeroInds(kInd)), ringY(zeroInds(kInd)), iSlice) = vecX;
                   volumeVectorY(ringX(zeroInds(kInd)), ringY(zeroInds(kInd)), iSlice) = vecY;
                   volumeVectorZ(ringX(zeroInds(kInd)), ringY(zeroInds(kInd)), iSlice) = vecZ;
               end
           end
        end
    end
end

%% Plot vectors along 2D slice

% Center of X
ringsSlice = permute(ringVolume(volumeSize(1)/2,:,:),[2 3 1]);
sliceVectorY = permute(volumeVectorY(volumeSize(1)/2,:,:),[2 3 1]);
sliceVectorZ = permute(volumeVectorZ(volumeSize(1)/2,:,:),[2 3 1]);

figure; hold on

indsToPlot = find(~isnan(sliceVectorY));

[plotX, plotY] = ind2sub(volumeSize(2:3), indsToPlot);

for iVec = 1:length(indsToPlot)
    ang = atan2(sliceVectorZ(indsToPlot(iVec)), sliceVectorY(indsToPlot(iVec)))/pi*180;
    
    %if ang > 20 & ang < 160
        line( [0 4*sliceVectorY(indsToPlot(iVec))]+plotX(iVec), [0 4*sliceVectorZ(indsToPlot(iVec))]+plotY(iVec),...
            'color', cols(ringsSlice(indsToPlot(iVec)), :))
    %end
end

% Center of Y
ringsSlice = permute(ringVolume(:,volumeSize(2)/2,:),[1 3 2]);
sliceVectorX = permute(volumeVectorX(:,volumeSize(2)/2,:),[1 3 2]);
sliceVectorY = permute(volumeVectorY(:,volumeSize(2)/2,:),[1 3 2]);
sliceVectorZ = permute(volumeVectorZ(:,volumeSize(2)/2,:),[1 3 2]);

figure; hold on

indsToPlot = find(~isnan(sliceVectorX));

[plotX, plotY] = ind2sub(volumeSize([1, 2]), indsToPlot);

for iVec = 1:length(indsToPlot)
    ang = atan2(sliceVectorZ(indsToPlot(iVec)), sliceVectorX(indsToPlot(iVec)))/pi*180;
    
    %if ang > 20 & ang < 160
        line( [0 4*sliceVectorX(indsToPlot(iVec))]+plotX(iVec), [0 4*sliceVectorZ(indsToPlot(iVec))]+plotY(iVec), ...
            'color', cols(ringsSlice(indsToPlot(iVec)), :))
    %end
end

%% Get some shells - old
if 0

    %for phiI = 1:length(phiVals)

        %shellInds = find(phiVolume == phiVals(phiI));

        shellInds = find(thetaVolume == 0);

        [shellX, shellY, shellZ] = ind2sub(volumeSize, shellInds);

        % figure; hold on;
        % plot3(shellX, shellY, shellZ, '.')

        %% Try to break up individual shells
        shellVolume = zeros(volumeSize, 'logical');

        shellVolume(shellInds) = 1;

        connectedComponents = bwconncomp(shellVolume,26);

        shellProps = regionprops(connectedComponents, 'PixelList');

         figure; hold on

        maxValues = zeros(length(shellProps), 3);

        for i = 1:length(shellProps)
             plot3(shellProps(i).PixelList(:,1), shellProps(i).PixelList(:,2), ...
                shellProps(i).PixelList(:,3), '.');

            maxValues(i,:) = [max(shellProps(i).PixelList(:,1) - volumeSize(1)/2), ...
                max(shellProps(i).PixelList(:,3)), length(shellProps(i).PixelList(:,1))];
        end

        layers2Use = find(maxValues(:,2) > 5 & maxValues(:,3) > 100);

        layerCols = lines(max(layers2Use));

        %% Plot as 2D vectors
        figure; 

        for i = layers2Use'
            % plot for phi
            subplot(1,2,1); hold on; axis equal
            tempInds = find(shellProps(i).PixelList(:,2) == volumeSize(2)/2);

            testX = shellProps(i).PixelList(tempInds,1);
            testZ = shellProps(i).PixelList(tempInds,3);

            % Get theta values
            for j = 1:length(testX)
                tempTheta = thetaVolume(testX(j), volumeSize(2)/2, testZ(j))/180*pi;

                if testX(j) < volumeSize(2)/2
                    tempTheta = -tempTheta;
                end

                line(testX(j) + [0 2*cos(tempTheta)], testZ(j) + [0 2*sin(tempTheta)], 'color', layerCols(i,:));
            end

            % plot for theta
            subplot(1,2,2); hold on; axis equal
            tempInds = find(shellProps(i).PixelList(:,3) == 1);

            testX = shellProps(i).PixelList(tempInds,1);
            testY = shellProps(i).PixelList(tempInds,2);

            % Get theta values
            for j = 1:length(testX)
                tempPhi = phiVolume(testX(j), testY(j), 1)/180*pi;

    %             if testX(j) < volumeSize(2)/2
    %                 tempPhi = -tempPhi;
    %                 
    %             end

                line(testX(j) + [0 2*cos(tempPhi)], testY(j) + [0 2*sin(tempPhi)], 'color', layerCols(i,:));
            end
        end
    %end

    %% Now fit curve to indvidual shells
    % Just take 2D to start with. 8th order poly seemed to be required to
    % capture shape of large shell.

    %%% Should fit to exterior of all included poly to get exit threshold

    figure; hold on; axis equal

    for i = layers2Use'
        tempInds = find(shellProps(i).PixelList(:,2) == volumeSize(2)/2);

        testX = shellProps(i).PixelList(tempInds,1);
        testZ = shellProps(i).PixelList(tempInds,3);

        plot(testX, testZ, '.')    

        fittedPoly = fit(testX,testZ,'poly8');

        plot(fittedPoly, testX, testZ);
    end    
    ylim([0 300])    
end
