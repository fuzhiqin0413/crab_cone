get_radius_profiles

%% 
%  close all

voxScale = 4; % Scale up factor on image for FDTD
writeLarge = 0;

%%% Switch to use profiles from CT
% Parameters from oliver
% interCone = 356.79/voxSize;
% interConeSD = 24.31/voxSize;
% use3Dintercone = 1;
% 
% outerCornea = 104.84/voxSize;
% outerCorneaSD = 6.59/voxSize;
% 
% epiCornea = 32.68/voxSize;
% epiCorneaSD = 2.21/voxSize;

% For display
nPlot = 2;
tText = 'Mean'; %'-2 SD'; 'Mean'; 
%%% 1 is closest but can currently intrude as intercone defined as radius not distance
coneStepForIntercone = 1;    
SDMult = 0;
figure

tipOffset = 30; 

useSlopedIntercone = 0;

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
coneLengthToUse = mean(lengthsToPlanes(:,2)) + std(lengthsToPlanes(:,2))*SDMult;

interconeLengthToUse = mean(lengthsToPlanes(:,1)) + std(lengthsToPlanes(:,1))*SDMult;

outerCorneaLengthToUse = mean(lengthsToPlanes(:,3)) + std(lengthsToPlanes(:,3))*SDMult - coneLengthToUse;

epicorneaLengthToUse = mean(lengthsToPlanes(:,4)) + std(lengthsToPlanes(:,4))*SDMult - coneLengthToUse - outerCorneaLengthToUse;

% Trim start of buffer
bufferLength = 10;

restretchLength_cone = ceil(max(lengthsToPlanes(:,2))/100)*100;

% Should ref length be 50??? - Still uses full length???
restretchLength_cornea = ceil((max(lengthsToPlanes(:,4))-min(lengthsToPlanes(:,3)))/10)*10;

%%% Not completely sure that just mean(coneRefDiameter) is correct to use for reference radius

[meanStretchedCone, stdStretchedCone, coneXRef] = restretchProfile(coneAverage(:,bufferLength+1:end), numCones, depthTests(bufferLength+1:end), ...
        coneRefDiameter, ones(numCones,1), ceil(lengthsToPlanes(:,2)), mean(coneRefDiameter), coneLengthToUse, restretchLength_cone, 0);

% Cone in cone now stretched to full cone length
    % Doesn't guarantee that cone tips are x-aligned but should match bases well...
[meanStretchedCinC, stdStretchedCinC] = restretchProfile(cInCAverage(:,bufferLength+1:end), numCones, depthTests(bufferLength+1:end), ...
        coneRefDiameter, ones(numCones,1), ceil(lengthsToPlanes(:,2)), mean(coneRefDiameter), coneLengthToUse, restretchLength_cone, 0);

[meanStretchedExposedIntercone, stdStretchedExposedIntercone] = restretchProfile(exposedInterconeAverage(:,bufferLength+1:end), numCones, depthTests(bufferLength+1:end), ...
        coneRefDiameter, ones(numCones,1), ceil(lengthsToPlanes(:,2)), mean(coneRefDiameter), coneLengthToUse, restretchLength_cone, 0);

[meanStretchedInternalIntercone, stdStretchedInternalIntercone] = restretchProfile(internalInterconePaths(:,bufferLength+1:end,coneStepForIntercone), numCones, depthTests(bufferLength+1:end), ...
        coneRefDiameter, ones(numCones,1), ceil(lengthsToPlanes(:,2)), mean(coneRefDiameter), coneLengthToUse, restretchLength_cone, 1);

% epicornea cone stretch to its own length
[meanStretchedEpicorneaInner, stdStretchedEpicorneaInner, corneaXRef] = restretchProfile(epicorneaInnerAverage(:,bufferLength+1:end), numCones, depthTests(bufferLength+1:end), ...
        coneRefDiameter, floor(lengthsToPlanes(:,3)), ceil(lengthsToPlanes(:,4)), mean(coneRefDiameter), epicorneaLengthToUse, restretchLength_cornea, 0);

figure;
subplot(1,2,1); hold on
errorbar(coneXRef, meanStretchedCone, stdStretchedCone);
errorbar(coneXRef, meanStretchedCinC, stdStretchedCinC);
errorbar(coneXRef, meanStretchedExposedIntercone, stdStretchedExposedIntercone);
errorbar(coneXRef, meanStretchedInternalIntercone, stdStretchedInternalIntercone);

subplot(1,2,2); hold on
errorbar(corneaXRef, meanStretchedEpicorneaInner, stdStretchedEpicorneaInner);


coneProfileToUse =  meanStretchedCone + stdStretchedCone*SDMult;

CinCProfileToUse =  meanStretchedCinC + stdStretchedCinC*SDMult;

exposedInterconeProfileToUse = meanStretchedExposedIntercone + stdStretchedExposedIntercone*SDMult;

internalInterconeProfileToUse = meanStretchedInternalIntercone + stdStretchedInternalIntercone*SDMult;

epicorneaProfileToUse = meanStretchedEpicorneaInner + stdStretchedEpicorneaInner*SDMult;

%%% WTF is going on here???
if 0
    % Tweaking profiles - set for 0 SD
        display('Tweaking')
    
        % Cut first three points from top of cone and lengthen end
        coneProfileToUse(1:3) = [];
        coneProfileToUse(end:end+3) = coneProfileToUse(end);
    
        % adjust radius on first three points of intercone
        innerConeProfileToUse(1:3) = innerConeProfileToUse(1:3)./[3 1.25 1.1];
end

% Make slice to fit dimensions
topEpiCornea = tipOffset + ceil(coneLengthToUse) + ceil(outerCorneaLengthToUse) + ceil(epicorneaLengthToUse);

sliceSize = round([3*max(coneProfileToUse), topEpiCornea + tipOffset]);

slice = zeros(sliceSize(1), sliceSize(2));

%%% Continue from hear then convert to 3D

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
topOuterCornea = tipOffset + round(coneLengthToUse) + round(outerCorneaLengthToUse);
slice(:, bottomOuterCornea:topOuterCornea) = outerCorneaValue;

% Place epi cornea
bottomEpiCornea = tipOffset + round(coneLengthToUse) + round(outerCorneaLengthToUse) + 1;
topEpiCornea = tipOffset + round(coneLengthToUse) + round(outerCorneaLengthToUse) + round(epicorneaLengthToUse);
slice(:, bottomEpiCornea:topEpiCornea) = epiCorneaValue;

subplot(1,3,nPlot)
if ~isnan(coneValue)
    imshow(slice'/(max(slice(:))+1))
else
    imshow((slice'-1.45)/(1.54-1.45))
end
    title(tText);
    
largerImg = imresize(slice, voxScale,'nearest'); 

if writeLarge 
    currentDirectory = pwd; 
    cd('/Users/gavintaylor/Desktop')

    warning('check names are correct')

    writematrix(largerImg,'500_average_cone_0_sd.csv') 
    cd(pwd)

    figure;
    imshow((largerImg'-1.45)/(1.54-1.45))
end