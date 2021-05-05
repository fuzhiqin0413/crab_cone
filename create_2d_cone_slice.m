get_radius_profiles

%% 
%  close all
% Parameters from oliver
interCone = 356.79;
interConeSD = 24.31;

outerCornea = 104.84;
outerCorneaSD = 6.59;

epiCornea = 32.68;
epiCorneaSD = 2.21;

sliceSize = round([3*max(coneAverage(:)) max(coneLength)*1.75]);

slice = zeros(sliceSize(1), sliceSize(2));

nPlot = 2;
tText =  'Mean'; '+2 SD';
SDMult = 0;

tipOffset = 10; % Number of Nans at start of array

% Set up parameters
coneLengthToUse = mean(coneLength) + std(coneLength)*SDMult;

innerConeLengthToUse = mean(corneaLength) + std(corneaLength)*SDMult;

interConeLengthToUse = interCone + interConeSD*SDMult;

outerCorneaToUse = outerCornea + outerCorneaSD*SDMult;

epiCorneaToUse = epiCornea + epiCorneaSD*SDMult;

ringDiameterToUse = mean(ringMean) + std(ringMean)*SDMult;

% Trim start of buffer
bufferLength = 10;

coneProfileToUse =  averageCone + averageConeStd*SDMult;
coneProfileToUse = coneProfileToUse(bufferLength+1:end);

innerConeProfileToUse =  averageCornea + averageCorneaStd*SDMult;
innerConeProfileToUse = innerConeProfileToUse(bufferLength+1:end);

% Use olivers average rather than full cover from center cone
% coveringProfileToUse = coneNumber(4,:);
% coveringProfileToUse = coveringProfileToUse(bufferLength+1:end);
% 
% firstProfileFull = find(coveringProfileToUse > 1);
% firstProfileFull = firstProfileFull(1);
% (round(coneLengthToUse) + round(innerConeLengthToUse)) - firstProfileFull 

interConeCover = round(coneLengthToUse) + round(innerConeLengthToUse) - round(interConeLengthToUse);

% Plot cone up to end of cone in cone length
for i = 1:(round(coneLengthToUse) + round(innerConeLengthToUse))
    
    yPos = tipOffset + i;
    
    if isnan(coneProfileToUse(i))
       coneProfileToUse(i) = coneProfileToUse(i-1); 
    end
    
    xPosTop = round(sliceSize(1)/2 + coneProfileToUse(i));
    xPosBottom = round(sliceSize(1)/2 - coneProfileToUse(i));

    % Fill between top and bottom
    slice(xPosBottom:xPosTop, yPos) = 1;
    
    % Fill intercone
    if i >= interConeCover
        slice(xPosTop+1:end, yPos) = 4;
        
        slice(1:xPosBottom-1, yPos) = 4;
    end
 
end

% Then plot cone in cone
for i = 1:round(innerConeLengthToUse)
    
    yPos = tipOffset + round(coneLengthToUse) + i;
    
    if isnan(innerConeProfileToUse(i))
       innerConeProfileToUse(i) = innerConeProfileToUse(i-1); 
    end
    
    xPosTop = round(sliceSize(1)/2 + innerConeProfileToUse(i));
    xPosBottom = round(sliceSize(1)/2 - innerConeProfileToUse(i));
    
    slice(xPosBottom:xPosTop, yPos) = 2;
end


% Place outer cornea
bottomOuterCornea = tipOffset + round(coneLengthToUse) + round(innerConeLengthToUse)+1;
topOuterCornea = tipOffset + round(coneLengthToUse) + round(innerConeLengthToUse) + round(outerCorneaToUse);
slice(:, bottomOuterCornea:topOuterCornea) = 2;

% Place epi cornea
bottomEpiCornea = tipOffset + round(coneLengthToUse) + round(innerConeLengthToUse) + round(outerCorneaToUse) + 1
topEpiCornea = tipOffset + round(coneLengthToUse) + round(innerConeLengthToUse) +...
    round(outerCorneaToUse) + round(epiCorneaToUse);
slice(:, bottomEpiCornea:topEpiCornea) = 3;

subplot(1,3,nPlot)
imshow(slice'/5)
title(tText);