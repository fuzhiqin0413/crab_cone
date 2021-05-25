get_radius_profiles

%% 
%  close all
% Parameters from oliver
interCone = 356.79/voxSize;
interConeSD = 24.31/voxSize;
use3Dintercone = 1;

outerCornea = 104.84/voxSize;
outerCorneaSD = 6.59/voxSize;

epiCornea = 32.68/voxSize;
epiCorneaSD = 2.21/voxSize;

sliceSize = round([3*max(coneAverage(:)) max(coneRingHeightMean)*1.5]);

slice = zeros(sliceSize(1), sliceSize(2));

nPlot = 2;
tText = 'Mean'; %'-2 SD'; 'Mean'; 
SDMult = 0;
figure

tipOffset = 10; % Number of Nans at start of array

% Set up parameters
% Labels
% outerValue = 0;
% innerValue = 1;
% coneValue = 5;
% outerCorneaValue = 4;
% epiCorneaValue = 3;
% interconeValue = 2;

% RI values
outerValue = 1.33;
innerValue = 1.34;
coneValue = NaN;
outerCorneaValue = 1.5;
epiCorneaValue = 1.53;
interconeValue = 1.47;

% Change to ring height for cone length
coneLengthToUse = mean(coneRingHeightMean) + std(coneRingHeightMean)*SDMult;

innerConeLengthToUse = mean(corneaRingHeightMean) + std(corneaRingHeightMean)*SDMult;

corneaOffsetToUse = mean(coneToCICLength) + std(coneToCICLength)*SDMult;

interConeLengthToUse = interCone + interConeSD*SDMult;

outerCorneaToUse = outerCornea + outerCorneaSD*SDMult;

epiCorneaToUse = epiCornea + epiCorneaSD*SDMult;

% Trim start of buffer
bufferLength = 10;

% direct average
% coneProfileToUse =  averageCone + averageConeStd*SDMult;
% coneProfileToUse = coneProfileToUse(bufferLength+1:end);
% 
% innerConeProfileToUse =  averageCornea + averageCorneaStd*SDMult;
% innerConeProfileToUse = innerConeProfileToUse(bufferLength+1:end);

% do restretch given ring height
restretchLength_cone = 500;
coneProfilesStretch = zeros(7, restretchLength_cone);

restretchLength_cornea = 50;
corneaProfilesStretch = zeros(7,restretchLength_cornea);

depthTests = -10:1000;

for i = 1:7
   tempConeProfile = coneAverage(i,bufferLength+1:round(coneRingHeightMean(i)));
   tempDepthProfile = depthTests(bufferLength+1:round(coneRingHeightMean(i)));
   coneProfilesStretch(i,:) = interp1(tempDepthProfile, tempConeProfile, ...
       tempDepthProfile(1):(tempDepthProfile(end)/(restretchLength_cone-1)):tempDepthProfile(end), 'linear');
   
   tempCorneaProfile = corneaAverage(i,bufferLength+1:round(corneaRingHeightMean(i)));
   tempDepthProfile = depthTests(bufferLength+1:round(corneaRingHeightMean(i)));
   corneaProfilesStretch(i,:) = interp1(tempDepthProfile, tempCorneaProfile, ...
       tempDepthProfile(1):(tempDepthProfile(end)/(restretchLength_cornea-1)):tempDepthProfile(end), 'linear');
end

% Take average and then stretch to current length
meanStretchedCone = mean(coneProfilesStretch);
meanStretchedCone = interp1((0:restretchLength_cone-1), meanStretchedCone, ...
    0:(restretchLength_cone-1)/(round(coneLengthToUse)-1):(restretchLength_cone-1), 'linear');

stdStretchedCone = std(coneProfilesStretch);
stdStretchedCone = interp1((0:restretchLength_cone-1), stdStretchedCone, ...
    0:(restretchLength_cone-1)/(round(coneLengthToUse)-1):(restretchLength_cone-1), 'linear');

meanStrechedCornea = mean(corneaProfilesStretch);
meanStrechedCornea = interp1((0:restretchLength_cornea-1), meanStrechedCornea, ...
    0:(restretchLength_cornea-1)/(round(innerConeLengthToUse)-1):(restretchLength_cornea-1), 'linear');

stdStrechedCornea = std(corneaProfilesStretch);
stdStrechedCornea = interp1((0:restretchLength_cornea-1), stdStrechedCornea, ...
    0:(restretchLength_cornea-1)/(round(innerConeLengthToUse)-1):(restretchLength_cornea-1), 'linear');

coneProfileToUse =  meanStretchedCone + stdStretchedCone*SDMult;

innerConeProfileToUse =  meanStrechedCornea + stdStrechedCornea*SDMult;

% Take from full cover of center cone
if use3Dintercone
    coveringProfileToUse = interconeRatio(4,:);
    coveringProfileToUse = coveringProfileToUse(bufferLength+1:end);

    firstProfileFull = find(coveringProfileToUse == 100);
    % Use scaling to represent re-streching
    interConeLengthToUse = round(firstProfileFull(1)/round(coneRingHeightMean(4))*round(coneLengthToUse));
end

interConeCover = round(coneLengthToUse) - round(interConeLengthToUse);

slice(:) = outerValue;

slice(:,1:tipOffset) = innerValue;

% Plot cone up to end of cone in cone length
for i = 1:round(coneLengthToUse)
    
    yPos = tipOffset + i;
    
    if isnan(coneProfileToUse(i))
       coneProfileToUse(i) = coneProfileToUse(i-1); 
    end
    
    xPosTop = round(sliceSize(1)/2 + coneProfileToUse(i));
    xPosBottom = round(sliceSize(1)/2 - coneProfileToUse(i));

    % Fill between top and bottom
    if ~isnan(coneValue)
        slice(xPosBottom:xPosTop, yPos) = coneValue;
    else
        % From oliver
        % rescale relative diameter to be 80
        tempX = ((xPosBottom:xPosTop)-sliceSize(1)/2)/coneProfileToUse(i)*80;
        
        slice(xPosBottom:xPosTop, yPos) = 1.52-0.000004914*tempX.^2;
    end
    
    
    % Fill intercone
    if i >= interConeCover
        slice(xPosTop+1:end, yPos) = interconeValue;
        
        slice(1:xPosBottom-1, yPos) = interconeValue;
    else
        slice(xPosTop+1:end, yPos) = innerValue;
        
        slice(1:xPosBottom-1, yPos) = innerValue;
    end
end

% Then plot cone in cone
for i = 1:round(innerConeLengthToUse)
    
    yPos = tipOffset + round(corneaOffsetToUse) + i;
    
    if isnan(innerConeProfileToUse(i))
       innerConeProfileToUse(i) = innerConeProfileToUse(i-1); 
    end
    
    xPosTop = round(sliceSize(1)/2 + innerConeProfileToUse(i));
    xPosBottom = round(sliceSize(1)/2 - innerConeProfileToUse(i));
    
    slice(xPosBottom:xPosTop, yPos) = outerCorneaValue;
end

[round(corneaOffsetToUse) + round(innerConeLengthToUse), round(coneLengthToUse)]

% Fill up to outer cornea if cone in cone is to short
    % Seems to come from a rounding problem
    %%% Note could also be case that cone in cone is too long but not noticed
if round(corneaOffsetToUse) + round(innerConeLengthToUse) < round(coneLengthToUse)
   for yPos = round(corneaOffsetToUse) + round(innerConeLengthToUse)+1:round(coneLengthToUse)
       slice(xPosBottom:xPosTop, yPos+ tipOffset) = outerCorneaValue;
   end
end

% Place outer cornea
bottomOuterCornea = tipOffset + round(coneLengthToUse)+1;
topOuterCornea = tipOffset + round(coneLengthToUse) + round(outerCorneaToUse);
slice(:, bottomOuterCornea:topOuterCornea) = outerCorneaValue;

% Place epi cornea
bottomEpiCornea = tipOffset + round(coneLengthToUse) + round(outerCorneaToUse) + 1;
topEpiCornea = tipOffset + round(coneLengthToUse) + round(outerCorneaToUse) + round(epiCorneaToUse);
slice(:, bottomEpiCornea:topEpiCornea) = epiCorneaValue;

subplot(1,3,nPlot)
if ~isnan(coneValue)
    imshow(slice'/(max(slice(:))+1))
else
    imshow((slice'-1.45)/(1.54-1.45))
end
    title(tText);