clear
clc
close all

% Load in angle data
dataFolder = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/Data/liquid angle csv/400nm';

convertAngleConvention = 1;

shrinkForRISlice = 00;

% Convert to mm
voxelSize = 200/10^6;

nPerp = 1.53; % Short axis

nPara = 1.51;

nDif = nPerp - nPara;

% Load theta angles
thetaVolume = loadcsvstack(sprintf('%s/theta', dataFolder), 0);

% Load phi angles
phiVolume = loadcsvstack(sprintf('%s/phi', dataFolder), 0);

%%% May need to rethink later angle corrections if corrected convenction used

if convertAngleConvention
    
    temp = thetaVolume;

    thetaVolume = phiVolume;

    phiVolume = temp;

    clear temp;
end

if shrinkForRISlice
    phiVolume = single(phiVolume);

    thetaVolume = single(thetaVolume);
end

volumeSize = size(phiVolume);

phiVals = unique(phiVolume(~isnan(phiVolume)));

thetaVals = unique(thetaVolume(~isnan(thetaVolume)));

% some Nan on theta border where phi has value, replace with 0
thetaVolume(isnan(thetaVolume) & phiVolume == 0) = 0;

%Convert to radians
thetaVolume = thetaVolume/180*pi;

phiVolume = phiVolume/180*pi;


%Get good inds
if ~shrinkForRISlice
    goodVolume = ~isnan(phiVolume); 

    goodInds = find(goodVolume);
end

figure; 
subplot(1,2,1);
imshow(phiVolume(:,:,10)/2/pi)

subplot(1,2,2);
imshow((thetaVolume(:,:,10)+pi/2)/pi)

%% Convert map directly to RI (Given rays along optical axis)

%%% Compare to vector corrected image
%%% Phi will be different, but on-axis RI should just depend on theta

outSideValue = NaN;

if shrinkForRISlice
    phiVolumeTemp = double( permute( phiVolume(volumeSize(1)/2,:,:),[2 3 1]));

    thetaVolumeTemp = double( permute( thetaVolume(volumeSize(1)/2,:,:),[2 3 1]));

    goodSlice = ~isnan(phiVolumeTemp); 

    goodInds = find(goodSlice);
else
    phiVolumeTemp = phiVolume;
    
    thetaVolumeTemp = thetaVolume;
end

% Calculate vector components
fiberXVolume = cos(phiVolumeTemp).*cos(thetaVolumeTemp);

fiberYVolume = sin(phiVolumeTemp).*cos(thetaVolumeTemp);

fiberZVolume = sin(thetaVolumeTemp);

% Flip z if pointing down
downZInds = find(fiberZVolume < 0);

fiberXVolume(downZInds) = -fiberXVolume(downZInds); 
fiberYVolume(downZInds) = -fiberYVolume(downZInds);
fiberZVolume(downZInds) = -fiberZVolume(downZInds); 

if shrinkForRISlice
    directRIVolume = ones(size(phiVolumeTemp))*outSideValue;
else
    directRIVolume = ones(volumeSize)*outSideValue;
end

% For on axis light.
angle2Fiber = acos( (fiberXVolume(goodInds)*0 + fiberYVolume(goodInds)*0 + fiberZVolume(goodInds)*1) ./ ...
        (sqrt(fiberXVolume(goodInds).^2 + fiberYVolume(goodInds).^2 + fiberZVolume(goodInds).^2) * norm([0 0 1])) );
    
% Adjusted range compared to original
directRIVolume(goodInds) = nPara + sin(angle2Fiber*2-pi/2)*nDif/2 + nDif/2;

% Test angle range
[tempVal, tempInd] = unique(directRIVolume(goodInds));

figure; 
plot(angle2Fiber(tempInd)/pi*180, tempVal); hold on
plot(0:5:90, nPara+sin(((0:5:90)/180*pi)*2-pi/2)*nDif/2 + nDif/2, 'x')

figure;

if ~shrinkForRISlice
    subplot(1,2,1)
    imshow((directRIVolume(:,:,round(volumeSize(3)*0.2))-1.5)/0.1)

    subplot(1,2,2)
    imshow((flipud(permute(directRIVolume(round(volumeSize(1)/2),:,:), [3 2 1]))-1.5)/0.1)
    line([volumeSize(1)/2 volumeSize(1)/2], [1 volumeSize(3)], 'color', 'g')
else
    imshow(flipud(directRIVolume'-1.50)/0.03)
    line([volumeSize(1)/2 volumeSize(1)/2], [1 volumeSize(3)], 'color', 'g')

    directRIVolume(isnan(directRIVolume)) = -1;
    writematrix(directRIVolume,'200_cone_long_section.csv') 
end
%% Split into rings
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
   
   ringVolume(:,:,iSlice) = ringLayer;
end

maxRings = max(ringVolume(:));

figure; subplot(1,3,1);
imshow(ringVolume(:,:,1)/maxRings)

subplot(1,3,2); hold on
imshow(phiVolume(:,:,1)/360)

for iRing = 1:maxRings
    inds = find(ringVolume(:,:,1) == iRing);
    
    [tempX, tempY] = ind2sub(volumeSize(1:2), inds);
    
    plot(tempX, tempY, '.', 'color', cols(iRing, :))
end
    
subplot(1,3,3);
imshow( permute( ringVolume(volumeSize(1)/2,:,:), [2 3 1])/maxRings)

%% For each band, identify vectors and RI in shells
%%% Had to add ring sampled/averaged aniso RI
%%% Also have to add interpolate at finer steps (say 2/4 values between points)
%%% Convert ring sampled back to full (also consider exterior if not)

% If set 90 theta inclined to 85
inclineVertical = 1;

calcFullVec = 1;

% Where to place RI (shell is always set to edge)
riPlacment = 'centre'; % alternative is 'edge' which is at phi 0/180

% Sparse allocation is being stored in full matrix
    % Probably more efficient to determine points and then use lists...
%%% Shell has out of plane discontinuities, need to check.
warning('point to check')
shellVolumeX = zeros(volumeSize)*NaN;
shellVolumeY = zeros(volumeSize)*NaN; 
shellVolumeZ = zeros(volumeSize)*NaN;

%For isotropic
ringSampledRIVolume = zeros(volumeSize)*NaN;
    
if calcFullVec
    vectorRIVolume = zeros(volumeSize)*NaN;

    fullVolumeX = zeros(volumeSize)*NaN;
    fullVolumeY = zeros(volumeSize)*NaN; 
    fullVolumeZ = zeros(volumeSize)*NaN;
end

for iSlice = 1 %1:volumeSize(3)

    figure; 
    subplot(1,2,1)
    imshow((thetaVolume(:,:,iSlice) + pi/2)/pi); hold on
    subplot(1,2,2)
    imshow(phiVolume(:,:,iSlice)/pi/2); hold on

    tempRingSlice = ringVolume(:,:,iSlice);

    tempThetaSlice = thetaVolume(:,:,iSlice);

    tempPhiSlice = phiVolume(:,:,iSlice);

    numRings = max(tempRingSlice(:));

    % get points in each ring
    for jRing = 1:numRings
       ringInds = find(tempRingSlice == jRing); 

       [ringX, ringY] = ind2sub(volumeSize(1:2), ringInds);

       % Get zero inclination within each ring (should be two bands)
       zeroInds = find((tempPhiSlice(ringInds) == 0 | tempPhiSlice(ringInds) == pi | ...
           tempPhiSlice(ringInds) == 2*pi));

       % Get orientation of each zero ind
       for kInd = 1:length(zeroInds)
           % Get vectors for zero point
           tempPhiPoint = tempPhiSlice(ringInds(zeroInds(kInd)));

           tempThetaPoint = tempThetaSlice(ringInds(zeroInds(kInd)));
           
           [vecXZero, vecYZero, vecZZero] = calculateshellvectorsfromangles(tempPhiPoint, tempThetaPoint, ...
               inclineVertical, volumeSize, ringX(zeroInds(kInd)), ringY(zeroInds(kInd)));
           
           % Get angle to zero fibre    
           angle2Fiber = acos( (vecXZero*0 + vecYZero*0 + vecZZero*1) ./ ...
                (sqrt(vecXZero.^2 + vecYZero.^2 + vecZZero.^2) * norm([0 0 1])) );
    
           % Adjusted range compared to original
           zeroRI = nPara + sin(angle2Fiber*2-pi/2)*nDif/2 + nDif/2;
           
           % get nearby within half ring thickness, and that are closer to center
           distToCenter = sqrt((ringX - volumeSize(1)/2).^2 + ...
               (ringY - volumeSize(2)/2).^2);

           distToPoint = sqrt((ringX - ringX(zeroInds(kInd))).^2 + ...
               (ringY - ringY(zeroInds(kInd))).^2);

           nearbyInds = find(distToPoint < length(phiVals)/2 & ...
               distToCenter < distToCenter(zeroInds(kInd)));

           if ~isempty(nearbyInds)
               % Take point with largest theta to use
               [maxV, pointToUse] = max( abs( tempThetaSlice(ringInds(nearbyInds)) ));

               % just did for 1
               if 1
               
    %                subplot(1,2,1)
    %                plot(ringX(zeroInds(kInd)), ringY(zeroInds(kInd)), 'x', 'color', cols(jRing, :))
    %                plot(ringX(nearbyInds(tempInds)), ringY(nearbyInds(tempInds)), '.', 'color', cols(jRing, :))
    %                
    %                subplot(1,2,2)
    %                plot(ringX(zeroInds(kInd)), ringY(zeroInds(kInd)), 'x', 'color', cols(jRing, :))
    %                plot(ringX(nearbyInds(tempInds)), ringY(nearbyInds(tempInds)), '.', 'color', cols(jRing, :))

                   %tempPhiSlice(ringInds(nearbyInds(tempInds)))

                   tempPhiPoint = tempPhiSlice(ringInds(nearbyInds(pointToUse)));

                   tempThetaPoint = tempThetaSlice(ringInds(nearbyInds(pointToUse)));

                   [vecX, vecY, vecZ] = calculateshellvectorsfromangles(tempPhiPoint, tempThetaPoint, ...
                        inclineVertical, volumeSize, ringX(nearbyInds(pointToUse)), ringY(nearbyInds(pointToUse)));

               else
                   maxInds = find(abs(tempThetaSlice(ringInds(nearbyInds))) == maxV);


                   %%% Test getting max angle between vectors (given theta is already max)
                   tempVecs = zeros(length(maxInds), 3);

                   angDifs = zeros(length(maxInds), 1);

                   for lInd = 1:length(maxInds)
                       tempPhiPoint = tempPhiSlice(ringInds(nearbyInds(maxInds(lInd))));

                       tempThetaPoint = tempThetaSlice(ringInds(nearbyInds(maxInds(lInd))));

                       [tempVecs(lInd,1), tempVecs(lInd,2), tempVecs(lInd,3)] = calculateshellvectorsfromangles(tempPhiPoint, tempThetaPoint, ...
                            inclineVertical, volumeSize, ringX(nearbyInds(maxInds(lInd))), ringY(nearbyInds(maxInds(lInd))));

                       angDifs(lInd) = acos((vecXZero*tempVecs(lInd,1)+ vecYZero*tempVecs(lInd,2)) / ...
                                (norm([vecXZero vecYZero]) * norm(tempVecs(lInd,1:2)) ));
                   end

                   [~, pointToUse] = max(angDifs);

                   %[angDifs, tempPhiSlice(ringInds(nearbyInds(maxInds)))]

                   vecX = tempVecs(pointToUse,1); vecY = tempVecs(pointToUse,2); vecZ = tempVecs(pointToUse,3);
               end
               
               % Place vector data in shells 
               shellVolumeX(ringX(zeroInds(kInd)), ringY(zeroInds(kInd)), iSlice) = vecX;
               shellVolumeY(ringX(zeroInds(kInd)), ringY(zeroInds(kInd)), iSlice) = vecY;
               shellVolumeZ(ringX(zeroInds(kInd)), ringY(zeroInds(kInd)), iSlice) = vecZ;
               
               % Place averaged GRIN data
               if strcmp(riPlacment, 'centre')
                   riX = round( mean([ringX(nearbyInds(pointToUse)) ringX(zeroInds(kInd))]));
                   riY = round( mean([ringY(nearbyInds(pointToUse)) ringY(zeroInds(kInd))]));
                   
               elseif strcmp(riPlacment, 'edge')
                   riX = ringX(zeroInds(kInd));
                   riY = ringY(zeroInds(kInd));
                   
               else
                  error('Unknown RI placment') 
               end
               
               % Get angle to other fiber
               angle2Fiber = acos( (vecX*0 + vecY*0 + vecZ*1) ./ ...
                    (sqrt(vecX.^2 + vecY.^2 + vecZ.^2) * norm([0 0 1])) );
    
               % Adjusted range compared to original
               otherRI = nPara + sin(angle2Fiber*2-pi/2)*nDif/2 + nDif/2;
               
               % Average both RIs
               ringSampledRIVolume(riX, riY, iSlice) = (zeroRI + otherRI)/2;
           end
       end
       
       if calcFullVec
           % Also calc RI diretly on each voxel for comparison
           for kInd = 1:length(ringInds)
               % Get vectors for zero point
               tempPhiPoint = tempPhiSlice(ringInds(kInd));

               tempThetaPoint = tempThetaSlice(ringInds(kInd));

               [vecX, vecY, vecZ] = calculateshellvectorsfromangles(tempPhiPoint, tempThetaPoint, ...
                   inclineVertical, volumeSize, ringX(kInd), ringY(kInd) );

               if isnan(vecX)
                    b = 1;
               end
               fullVolumeX(ringX(kInd), ringY(kInd), iSlice) = vecX;
               fullVolumeY(ringX(kInd), ringY(kInd), iSlice) = vecY;
               fullVolumeZ(ringX(kInd), ringY(kInd), iSlice) = vecZ;

               % Get angle to zero fibre    
               angle2Fiber = acos( (vecX*0 + vecY*0 + vecZ*1) ./ ...
                    (sqrt(vecX.^2 + vecY.^2 + vecZ.^2) * norm([0 0 1])) );

               % Adjusted range compared to original
               vectorRIVolume(ringX(kInd), ringY(kInd), iSlice) = nPara + sin(angle2Fiber*2-pi/2)*nDif/2 + nDif/2;

           end
       end
    end
end

% Plot to check
figure; 
subplot(2,2,1);
imshow((directRIVolume(:,:,iSlice)-1.5)/0.1)

temp = vectorRIVolume(:,:,iSlice)-directRIVolume(:,:,iSlice);
max(abs(temp(~isnan(temp))))

subplot(2,2,2);
imshow(temp)

subplot(2,2,3);
imshow((ringSampledRIVolume(:,:,iSlice)-1.5)/0.1)

subplot(2,2,4);
imshow((vectorRIVolume(:,:,iSlice)-1.5)/0.1)

% Plot hists to check
figure; hold on
temp = ringSampledRIVolume(:,:,iSlice);
inds2Use = find(~isnan(temp));
[n, x] = hist(temp(inds2Use), 1.5:0.001:1.55);
plot(x, n)

temp = directRIVolume(:,:,iSlice);
[n, x] = hist(temp(inds2Use), 1.5:0.001:1.55);
plot(x,n);

% Plot vectors on slice to check
figure; 

% Plot for full X cross.section
subplot(2,2,1); hold on
ringsSlice = permute(ringVolume(volumeSize(1)/2,:,:),[2 3 1]);
sliceVectorX = permute(fullVolumeX(volumeSize(1)/2,:,:),[2 3 1]);
sliceVectorY = permute(fullVolumeY(volumeSize(1)/2,:,:),[2 3 1]);
sliceVectorZ = permute(fullVolumeZ(volumeSize(1)/2,:,:),[2 3 1]);

indsToPlot = find(~isnan(sliceVectorZ));
[plotX, plotY] = ind2sub(volumeSize(2:3), indsToPlot);

for iVec = 1:length(indsToPlot)
    line( [0 4*sliceVectorY(indsToPlot(iVec))]+plotX(iVec), [0 4*sliceVectorZ(indsToPlot(iVec))]+plotY(iVec),...
         [0 4*sliceVectorX(indsToPlot(iVec))], 'color', cols(ringsSlice(indsToPlot(iVec)), :))
end

% Plot for rings X cross.section
subplot(2,2,2); hold on
sliceVectorX = permute(shellVolumeX(volumeSize(1)/2,:,:),[2 3 1]);
sliceVectorY = permute(shellVolumeY(volumeSize(1)/2,:,:),[2 3 1]);
sliceVectorZ = permute(shellVolumeZ(volumeSize(1)/2,:,:),[2 3 1]);

indsToPlot = find(~isnan(sliceVectorZ));
[plotX, plotY] = ind2sub(volumeSize(2:3), indsToPlot);

for iVec = 1:length(indsToPlot)
    line( [0 4*sliceVectorY(indsToPlot(iVec))]+plotX(iVec), [0 4*sliceVectorZ(indsToPlot(iVec))]+plotY(iVec),...
        [0 4*sliceVectorX(indsToPlot(iVec))], 'color', cols(ringsSlice(indsToPlot(iVec)), :))
end

% Plot for full Y Cross.section
%%% Shows discontinuities in mirror angles, probably not good?
subplot(2,2,3); hold on
ringsSlice = permute(ringVolume(:,volumeSize(2)/2,:),[1 3 2]);
sliceVectorX = permute(fullVolumeX(:,volumeSize(2)/2,:),[1 3 2]);
sliceVectorY = permute(fullVolumeY(:,volumeSize(2)/2,:),[1 3 2]);
sliceVectorZ = permute(fullVolumeZ(:,volumeSize(2)/2,:),[1 3 2]);

indsToPlot = find(~isnan(sliceVectorZ));
[plotX, plotY] = ind2sub(volumeSize([1, 2]), indsToPlot);

for iVec = 1:length(indsToPlot)
    line( [0 4*sliceVectorX(indsToPlot(iVec))]+plotX(iVec), [0 4*sliceVectorZ(indsToPlot(iVec))]+plotY(iVec), ...
        [0 4*sliceVectorY(indsToPlot(iVec))], 'color', cols(ringsSlice(indsToPlot(iVec)), :))
end

% Plot for full Y Cross.section
subplot(2,2,4); hold on
sliceVectorX = permute(shellVolumeX(:,volumeSize(2)/2,:),[1 3 2]);
sliceVectorY = permute(shellVolumeY(:,volumeSize(2)/2,:),[1 3 2]);
sliceVectorZ = permute(shellVolumeZ(:,volumeSize(2)/2,:),[1 3 2]);

indsToPlot = find(~isnan(sliceVectorZ));
[plotX, plotY] = ind2sub(volumeSize([1, 2]), indsToPlot);

for iVec = 1:length(indsToPlot)
    line( [0 4*sliceVectorX(indsToPlot(iVec))]+plotX(iVec), [0 4*sliceVectorZ(indsToPlot(iVec))]+plotY(iVec), ...
        [0 4*sliceVectorY(indsToPlot(iVec))], 'color', cols(ringsSlice(indsToPlot(iVec)), :))
end

%% Exerpiment with 2D vector data
figure; hold on; axis equal; set(gca, 'Clipping', 'off')
% Take bottom slice of each
ringsSlice = ringVolume(:,:,1);
sliceVectorX = shellVolumeX(:,:,1);
sliceVectorY = shellVolumeY(:,:,1);
sliceVectorZ = shellVolumeZ(:,:,1);

for iX = 1:volumeSize(1)
    
    for jY = 1:volumeSize(2)
        % check if nan
        if ~isnan(sliceVectorX(iX, jY))
            % Plot of on center or diagonal
            toPlot = 0;
            
            if iX == volumeSize(1)/2
                toPlot = 1;
            elseif jY == volumeSize(2)/2
                toPlot = 1;
            elseif iX == jY 
                toPlot = 1;
            elseif iX == (-1*(jY - volumeSize(2)/2) + volumeSize(2)/2)
                toPlot = 1;
            end
            
            if toPlot
                line( [0 10*sliceVectorX(iX, jY)]+iX, [0 10*sliceVectorY(iX, jY)]+jY, ...
                    [0 10*sliceVectorZ(iX, jY)], 'color', cols(ringsSlice(iX, jY), :))
                
                plot(iX, jY, 'x')
            end
            
        end
    end
end
    
