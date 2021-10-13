% get_radius_profiles

%% 
%  close all

voxScale = 1;
writeLarge = 1;

error('Check if voxel sizing is used correctly')

% Parameters from oliver
interCone = 356.79/voxSize;
interConeSD = 24.31/voxSize;
use3Dintercone = 1;

outerCornea = 104.84/voxSize;
outerCorneaSD = 6.59/voxSize;

epiCornea = 32.68/voxSize;
epiCorneaSD = 2.21/voxSize;

nPlot = 2;
tText = 'Mean'; %'-2 SD'; '+'; 
SDMult = 0;
figure

tipOffset = 30; % Number of Nans at start of arrays

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

% Tweaking profiles - set for 0 SD
    warning('Tweaking for smooth cone top')

    % Cut first three points from top of cone and lengthen end
    coneProfileToUse(1:3) = [];
    coneProfileToUse(end:end+3) = coneProfileToUse(end);

    % Tweak radius on first three points of intercone
    innerConeProfileToUse(1:3) = innerConeProfileToUse(1:3)./[3 1.25 1.1];

% Take from full cover of center cone
if use3Dintercone
    coveringProfileToUse = interconeRatio(4,:);
    coveringProfileToUse = coveringProfileToUse(bufferLength+1:end);

    firstProfileFull = find(coveringProfileToUse == 100);
    % Use scaling to represent re-streching
    interConeLengthToUse = round(firstProfileFull(1)/round(coneRingHeightMean(4))*round(coneLengthToUse));
end

interConeCover = round(coneLengthToUse) - round(interConeLengthToUse);


% Make volume to fit dimensions
topEpiCornea = tipOffset + round(coneLengthToUse) + round(outerCorneaToUse) + round(epiCorneaToUse);

volumeSize = round([3*max(coneProfileToUse), 3*max(coneProfileToUse), topEpiCornea + tipOffset]);

volume = zeros(volumeSize(1), volumeSize(2), volumeSize(3));


volume(:) = outerValue;

volume(:,:,1:tipOffset) = innerValue;

[xSlice, ySlice] = meshgrid(1:volumeSize(1), 1:volumeSize(2));

rSlice = sqrt((xSlice(:)-volumeSize(1)/2).^2 + (ySlice(:)-volumeSize(2)/2).^2);

% Plot cone up to end of cone in cone length
for i = 1:round(coneLengthToUse)
    
    zPos = tipOffset + i;
    
    if isnan(coneProfileToUse(i))
       coneProfileToUse(i) = coneProfileToUse(i-1); 
    end
    
    % Fill between cone
    tempSlice = zeros(volumeSize(1), volumeSize(2));
    
    goodInds = find(rSlice <= coneProfileToUse(i));

    if ~isnan(coneValue)
        
        tempSlice(goodInds) = coneValue;
    else
        % From oliver
        % rescale relative diameter to be 80
        tempR = rSlice(goodInds)/coneProfileToUse(i)*80;
        
        % Using negative coding for GRIN
        tempSlice(goodInds) = -(1.52-0.000004914*tempR.^2);
    end
 
    % Fill intercone
    outerInds = find(rSlice > coneProfileToUse(i));
    
    if i >= interConeCover
        tempSlice(outerInds) = interconeValue;
    else
        tempSlice(outerInds) = innerValue;
    end
    
    % Place in volume
    volume(:, :, zPos) = tempSlice;
end

% Then plot cone in cone
for i = 1:round(innerConeLengthToUse)
    
    zPos = tipOffset + round(corneaOffsetToUse) + i;
    
    if isnan(innerConeProfileToUse(i))
       innerConeProfileToUse(i) = innerConeProfileToUse(i-1); 
    end
    
    goodInds = find(rSlice <= innerConeProfileToUse(i));
    
    tempSlice = volume(:, :, zPos);
    
    tempSlice(goodInds) = outerCorneaValue;
    
    volume(:, :, zPos) = tempSlice;
end

[round(corneaOffsetToUse) + round(innerConeLengthToUse), round(coneLengthToUse)]

% Fill up to outer cornea if cone in cone is to short
    % Seems to come from a rounding problem
    %%% Note could also be case that cone in cone is too long but haven't seen this occur
if round(corneaOffsetToUse) + round(innerConeLengthToUse) < round(coneLengthToUse)
   for zPos = round(corneaOffsetToUse) + round(innerConeLengthToUse)+1:round(coneLengthToUse)
       
       tempSlice = volume(:, :, zPos + tipOffset);
    
       tempSlice(goodInds) = outerCorneaValue;

       volume(:, :, zPos + tipOffset) = tempSlice;
   end
end

% Place outer cornea
bottomOuterCornea = tipOffset + round(coneLengthToUse)+1;
topOuterCornea = tipOffset + round(coneLengthToUse) + round(outerCorneaToUse);
volume(:, :, bottomOuterCornea:topOuterCornea) = outerCorneaValue;

% Place epi cornea
bottomEpiCornea = tipOffset + round(coneLengthToUse) + round(outerCorneaToUse) + 1;
topEpiCornea = tipOffset + round(coneLengthToUse) + round(outerCorneaToUse) + round(epiCorneaToUse);
volume(:, :, bottomEpiCornea:topEpiCornea) = epiCorneaValue;

GRINInds = find(volume < 0);

subplot(1,3,nPlot)
if ~isnan(coneValue)
    imshow(volume'/(max(volume(:))+1))
else
    volume(GRINInds) = -volume(GRINInds);
    
    subplot(1,3,1)
    imshow((volume(:,:,tipOffset + round(coneLengthToUse/2))'-1.45)/(1.54-1.45))
    
    subplot(1,3,2)
    imshow((volume(:,:,tipOffset + round(coneLengthToUse) - round(innerConeLengthToUse/2))'-1.45)/(1.54-1.45))
    
    subplot(1,3,3)
    imshow(flipud(permute(volume(round(volumeSize(1)/2),:,:), [2 3 1])'-1.45)/(1.54-1.45))
    
    volume(GRINInds) = -volume(GRINInds);
end

    
volume = imresize(volume, voxScale, 'nearest'); 

if writeLarge 
    currentDirectory = pwd; 
    cd('/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/Data/Matlab RI Volumes')

    warning('check names are correct')

    save('Test1_SD0_vox2.mat', 'volume') 
    cd(pwd)
end