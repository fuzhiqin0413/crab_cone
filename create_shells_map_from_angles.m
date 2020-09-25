clear
clc
close all

% Load in angle data
dataFolder = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/Data/liquid angle csv/400nm';

convertAngleConvention = 1;

% Convert to mm
voxelSize = 400/10^6;

nPerp = 1.53; % Short axis

nPara = 1.51;

nDif = nPerp - nPara;

% Load theta angles
thetaVolume = loadcsvstack(sprintf('%s/theta', dataFolder), 0);

% Load phi angles
phiVolume = loadcsvstack(sprintf('%s/phi', dataFolder), 0);

%%% May need to rethink later angle corrections if corrected convenction
%%% used

if convertAngleConvention
    
    temp = thetaVolume;

    thetaVolume = phiVolume;

    phiVolume = temp;

    clear temp;
end

volumeSize = size(phiVolume);

phiVals = unique(phiVolume(~isnan(phiVolume)));

thetaVals = unique(thetaVolume(~isnan(thetaVolume)));

% some Nan on theta border where phi has value, replace with 0
thetaVolume(isnan(thetaVolume) & phiVolume == 0) = 0;

% Convert to radians
thetaVolume = thetaVolume/180*pi;

phiVolume = phiVolume/180*pi;

%Get good inds
goodVolume = ~isnan(phiVolume); 

goodInds = find(goodVolume);

figure; 
subplot(1,2,1);
imshow(phiVolume(:,:,10)/pi)

subplot(1,2,2);
imshow((thetaVolume(:,:,10)+pi/2)/pi)

%% Convert map directly to RI (Given rays along optical axis)

%%% Compare to vector corrected image
%%% Phi will be different, but on-axis RI should just depend on theta

outSideValue = NaN;

% Calculate vector components
fiberXVolume = cos(phiVolume).*cos(thetaVolume);

fiberYVolume = sin(phiVolume).*cos(thetaVolume);

fiberZVolume = sin(thetaVolume);

% Flip z if pointing down
downZInds = find(fiberZVolume < 0);

fiberXVolume(downZInds) = -fiberXVolume(downZInds); 
fiberYVolume(downZInds) = -fiberYVolume(downZInds);
fiberZVolume(downZInds) = -fiberZVolume(downZInds); 

directRIVolume = ones(volumeSize)*outSideValue;

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
subplot(1,2,1)
imshow((directRIVolume(:,:,round(volumeSize(3)*0.2))-1.5)/0.1)

subplot(1,2,2)
imshow((flipud(permute(directRIVolume(round(volumeSize(1)/2),:,:), [3 2 1]))-1.5)/0.1)
line([volumeSize(1)/2 volumeSize(1)/2], [1 volumeSize(3)], 'color', 'g')

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
%%% %%% Needs testing!!!

% If set 90 theta inclined to 85
inclineVertical = 1;

% Where to place RI (shell is always set to edge)
riPlacment = 'centre'; % alternative is 'edge' which is at phi 0/180

%%% Sparse allocation is being stored in full matrix
    %%% Probably more efficient to determine points and then use lists...
shellVolumeX = zeros(volumeSize)*NaN;
shellVolumeY = zeros(volumeSize)*NaN; 
shellVolumeZ = zeros(volumeSize)*NaN;

%For isotropic
%%% Check how this matches data for direct method
ringSampledRIVolume = zeros(volumeSize)*NaN;
    
for iSlice = 1:volumeSize(3)

    %figure; 
    %imshow((thetaVolume(:,:,iSlice) + pi/2)/pi); hold on

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
           tempPhiPoint = tempPhiSlice(zeroInds(kInd));

           tempThetaPoint = tempThetaSlice(zeroInds(kInd));
           
           [vecXZero, vecYZero, vecZZero] = calculateshellvectorsfromangles(tempPhiPoint, tempThetaPoint);
           
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
               %plot(ringX(nearbyInds), ringY(nearbyInds), '.', 'color', cols(jRing, :))

               % Take point with largest theta to use
               [~, pointToUse] = max( abs( tempThetaSlice(ringInds(nearbyInds)) ));

               %%% Need to do calculation for both zero ind and max ind
               
               tempPhiPoint = tempPhiSlice(ringInds(nearbyInds(pointToUse)));

               tempThetaPoint = tempThetaSlice(ringInds(nearbyInds(pointToUse)));

               [vecX, vecY, vecZ] = calculateshellvectorsfromangles(tempPhiPoint, tempThetaPoint);

               % Place vector data in shells 
               shellVectorX(ringX(zeroInds(kInd)), ringY(zeroInds(kInd)), iSlice) = vecX;
               shellVectorY(ringX(zeroInds(kInd)), ringY(zeroInds(kInd)), iSlice) = vecY;
               shellVectorZ(ringX(zeroInds(kInd)), ringY(zeroInds(kInd)), iSlice) = vecZ;
               
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
               
               ringSampledRIVolume(riX, riY, iSlice) = (zeroRI + otherRI)/2;
               
               %%% Check how this matches data from direct method
           end
       end
    end
end