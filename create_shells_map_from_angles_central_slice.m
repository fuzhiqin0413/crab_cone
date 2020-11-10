%%% Full map is too large to work on 200, 
%%% so split to seperate file to save slices for FDTD

clear
clc
close all

% Load in angle data
dataFolder = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/Data/liquid angle csv/200nm_thick';

convertAngleConvention = 0;

% Check conversion is set correctly
if isempty(strfind(dataFolder, 'thick')) & ~convertAngleConvention
    warning('Check if convertAngleConvention should be set')

elseif ~isempty(strfind(dataFolder, 'thick')) & convertAngleConvention
    warning('Check if convertAngleConvention should be unset')
    
end

shrinkForRISlice = 0;

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

% Shrink main volume size a bit
phiVolume = single(phiVolume);

thetaVolume = single(thetaVolume);
    
volumeSize = size(phiVolume);

phiVals = unique(phiVolume(~isnan(phiVolume)));

thetaVals = unique(thetaVolume(~isnan(thetaVolume)));

% some Nan on theta border where phi has value, replace with 0
thetaVolume(isnan(thetaVolume) & phiVolume == 0) = 0;

% convert to radians
thetaVolume = thetaVolume/180*pi;

phiVolume = phiVolume/180*pi;

figure; 
subplot(1,2,1);
imshow(phiVolume(:,:,10)/2/pi)

subplot(1,2,2);
imshow((thetaVolume(:,:,10)+pi/2)/pi)

%% Convert map directly to RI (Given rays along optical axis)

%%% Compare to vector corrected image
%%% Phi will be different, but on-axis RI should just depend on theta

outSideValue = NaN;

% Shrink down to save size
phiVolumeTemp = double( permute( phiVolume(volumeSize(1)/2,:,:),[2 3 1]));

thetaVolumeTemp = double( permute( thetaVolume(volumeSize(1)/2,:,:),[2 3 1]));

goodSlice = ~isnan(phiVolumeTemp); 

goodInds = find(goodSlice);

% Calculate vector components
fiberXVolume = cos(phiVolumeTemp).*cos(thetaVolumeTemp);

fiberYVolume = sin(phiVolumeTemp).*cos(thetaVolumeTemp);

fiberZVolume = sin(thetaVolumeTemp);

% Flip z if pointing down
downZInds = find(fiberZVolume < 0);

fiberXVolume(downZInds) = -fiberXVolume(downZInds); 
fiberYVolume(downZInds) = -fiberYVolume(downZInds);
fiberZVolume(downZInds) = -fiberZVolume(downZInds); 

directRIVolume = ones(size(phiVolumeTemp))*outSideValue;

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
imshow(flipud(directRIVolume'-1.50)/0.03)
line([volumeSize(1)/2 volumeSize(1)/2], [1 volumeSize(3)], 'color', 'g')
%% Split into rings

%%% Instead of volume, write central into slice.

cols = lines(100);

ringSlice = zeros(size(phiVolumeTemp));

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
   
   %%% Put central row into slice collector
   ringSlice(:,iSlice) = ringLayer(volumeSize(1)/2,:);
end

maxRings = max(ringSlice(:));

%% 
figure; subplot(1,2,1);
imshow(flipud(ringSlice')/maxRings)

subplot(1,2,2); hold on
imshow(flipud(phiVolumeTemp')/2/pi)

for iRing = 1:maxRings
    inds = find(ringSlice == iRing);
    
    [tempX, tempY] = ind2sub(volumeSize(1:2), inds);
    
    plot(tempX, -tempY+volumeSize(3), '.', 'color', cols(iRing, :))
end
    
%%% Average over RI rings on slice
directRIVolume(directRIVolume == -1) = NaN;

averagedRIVolume = zeros(size(phiVolumeTemp))*NaN;

interpolatedRIVolume = zeros(size(phiVolumeTemp))*NaN;

for iSlice = 1:volumeSize(3)
    tempRISlice = directRIVolume(:, iSlice);
    
    % Set up interpolate
    numRings = length( unique( ringSlice(ringSlice(:, iSlice) > 0, iSlice)));
    
    % Store radial post and intesnity for interpolation
    interpPoints = zeros(numRings*2-1,2);
    
    for jRing = 1:numRings
        % 1d, so inds are also y pos
        yInds = find(ringSlice(:, iSlice) == jRing);

        % Check if ring doesn't cross middle
        if ~any(yInds == volumeSize(2)/2)
            % put points on both sides
            indsLeft = yInds(yInds < volumeSize(2)/2);
            
            interpPoints(jRing, :) = [mean(indsLeft) nanmean(tempRISlice(indsLeft))];
            
            indsRight = yInds(yInds > volumeSize(2)/2);
            
            interpPoints(numRings*2 - 1 - jRing + 1, :) = [mean(indsRight) nanmean(tempRISlice(indsRight))];
        else
            % no break, put in center
            interpPoints(numRings, :) = [mean(yInds) nanmean(tempRISlice(yInds))];
        end
        
        tempRISlice(yInds) = nanmean(tempRISlice(yInds));
    end
    
    averagedRIVolume(:, iSlice) = tempRISlice;
    
    % make interpolated slice
    yInds = find(ringSlice(:, iSlice) > 0);
    
    interpolatedRIVolume(yInds, iSlice) = interp1(interpPoints(:,1) , interpPoints(:,2), yInds, 'linear', interpPoints(1,2));
end
    
figure;
subplot(1,2,1); hold on
imshow(flipud(averagedRIVolume'-1.51)/0.02)
line([volumeSize(2)/2 volumeSize(2)/2], [1 volumeSize(3)], 'color', 'g')

subplot(1,2,2); hold on
imshow(flipud(interpolatedRIVolume'-1.51)/0.02)
line([volumeSize(2)/2 volumeSize(2)/2], [1 volumeSize(3)], 'color', 'g')

currentDirectory = pwd; 
cd('/Users/gavintaylor/Desktop')

warning('check names are correct')

directRIVolume(isnan(directRIVolume)) = -1;

writematrix(directRIVolume,'500_cone_full_section_lamellae.csv') 

averagedRIVolume(isnan(averagedRIVolume)) = -1;

writematrix(averagedRIVolume,'500_cone_full_section_averaged.csv') 

interpolatedRIVolume(isnan(interpolatedRIVolume)) = -1;

writematrix(interpolatedRIVolume,'500_cone_full_section_interpolated.csv') 
cd(pwd)