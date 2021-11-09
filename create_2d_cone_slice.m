get_radius_profiles

%% 
%  close all
voxSize = 2.18;

scaleUpVoxels = 0;
    newVoxSize = 0.5;
    voxScale = voxSize/newVoxSize; % Scale up factor on image for FDTD

writeImage = 0;
    fileName = '500_average_cone_0_sd.csv';

displayRI = 1;
    
% For display
nPlot = 4;
tText = 'Smooth - 2.18 um'; %'-2 SD'; 'Mean';

SDMult = 0;
mainFig = figure;

%Can be from 1 - 4, can also be equal
interconeOnLeft = NaN;
interconeOnRight = NaN;

tipOffset = 30; 
bufferLength = 10;
 
reduceTipsCone = 1;
reduceTipsCinC = 1;
reduceTipsCornea = 1;

smoothProfiles = 1;
    smoothLength = 5;
    smoothLengthIntercone = 20;

tipGradientCorrection = 0;    
    correctionType = 'both'; %'distance', 'height', 'both'
    correctionScale = 1;
    
if displayRI
    outerValue = 1.33;
    innerValue = 1.34;
    coneValue = 'radial'; %'radial', 'linear' 'both'
    outerCorneaValue = 1.5;
    epicorneaValue = 1.53;
    interconeValue = 1.47;

else
% Labels
    outerValue = 0;
    innerValue = 1;
    coneValue = 5;
    outerCorneaValue = 4;
    epicorneaValue = 3;
    interconeValue = 2;
end

if reduceTipsCone
    bufferOffset = 9;
else
    bufferOffset = 10;
end

if bufferOffset > bufferLength
    error('buffer offset to large')
end

coneLengthToUse = (bufferLength - bufferOffset) + mean(lengthsToPlanes(:,2)) + std(lengthsToPlanes(:,2))*SDMult;

% interconeLengthToUse = (bufferLength - bufferOffset) + mean(lengthsToPlanes(:,1)) + std(lengthsToPlanes(:,1))*SDMult;

outerCorneaLengthToUse = (bufferLength - bufferOffset) + mean(lengthsToPlanes(:,3)) + std(lengthsToPlanes(:,3))*SDMult - coneLengthToUse;

epicorneaLengthToUse = (bufferLength - bufferOffset) + mean(lengthsToPlanes(:,4)) + std(lengthsToPlanes(:,4))*SDMult - coneLengthToUse - outerCorneaLengthToUse;

restretchLength_cone = ceil(max(lengthsToPlanes(:,2))/100)*100;

restretchLength_cornea = ceil((max(lengthsToPlanes(:,4))-min(lengthsToPlanes(:,3)))/10)*10;

[meanStretchedCone, stdStretchedCone, coneXRef] = restretchProfile(coneAverage(:,bufferOffset+1:end), numCones, depthTests(bufferLength+1:end), ...
        coneRefDiameter, ones(numCones,1), ceil(lengthsToPlanes(:,2)), mean(coneRefDiameter), coneLengthToUse, restretchLength_cone, 0, 0, 0);

% Cone in cone now stretched to full cone length
    % Doesn't guarantee that cone tips are x-aligned but should match bases well...
[meanStretchedCinC, stdStretchedCinC] = restretchProfile(cInCAverage(:,bufferLength+1:end), numCones, depthTests(bufferLength+1:end), ...
        coneRefDiameter, ones(numCones,1), ceil(lengthsToPlanes(:,2)), mean(coneRefDiameter), coneLengthToUse, restretchLength_cone, 0, 1, 1);

% normalised from top plane of exposed intercone, would be good to also normalize to base plane of exposed intercone
[meanStretchedExposedIntercone, stdStretchedExposedIntercone] = restretchProfile(exposedInterconeAverage(:,bufferLength+1:end), numCones, depthTests(bufferLength+1:end), ...
        coneRefDiameter, ones(numCones,1), ceil(lengthsToPlanes(:,2)), mean(coneRefDiameter), coneLengthToUse, restretchLength_cone, 0, 0, 1);

if isnan(interconeOnLeft); tempVal = 1; else; tempVal = interconeOnLeft; end
[meanStretchedInternalInterconeLeft, stdStretchedInternalInterconeLeft] = restretchProfile(internalInterconePaths(:,bufferLength+1:end,tempVal), numCones, depthTests(bufferLength+1:end), ...
        coneRefDiameter, ones(numCones,1), ceil(lengthsToPlanes(:,2)), mean(coneRefDiameter), coneLengthToUse, restretchLength_cone, 1, 1, 0);

if isnan(interconeOnRight); tempVal = 1; else; tempVal = interconeOnRight; end
[meanStretchedInternalInterconeRight, stdStretchedInternalInterconeRight] = restretchProfile(internalInterconePaths(:,bufferLength+1:end,tempVal), numCones, depthTests(bufferLength+1:end), ...
        coneRefDiameter, ones(numCones,1), ceil(lengthsToPlanes(:,2)), mean(coneRefDiameter), coneLengthToUse, restretchLength_cone, 1, 1, 0);
    
% epicornea cone stretch to epicornea length
[meanStretchedEpicorneaInner, stdStretchedEpicorneaInner, corneaXRef] = restretchProfile(epicorneaInnerAverage(:,bufferLength+1:end), numCones, depthTests(bufferLength+1:end), ...
        coneRefDiameter, floor(lengthsToPlanes(:,3)), ceil(lengthsToPlanes(:,4)), mean(coneRefDiameter), epicorneaLengthToUse, restretchLength_cornea, 0, -1, 1);

coneXRef = (0:(length(meanStretchedCone)+1)/length(meanStretchedCone):length(meanStretchedCone))*voxSize; 
corneaXRef = (0:(length(meanStretchedEpicorneaInner)+1)/length(meanStretchedEpicorneaInner):length(meanStretchedEpicorneaInner))*voxSize; 

figure;
subplot(1,3,1); hold on
errorbar(coneXRef, meanStretchedCone*voxSize, stdStretchedCone*voxSize);
errorbar(coneXRef, meanStretchedCinC*voxSize, stdStretchedCinC*voxSize);

if reduceTipsCone
    %%% Tried interpolation but didn't work well
    % decided to add a step with half radius at top so it's not too blunt

    % reduce cone tip
    goodConeInds = find(~isnan(meanStretchedCone));
    
    meanStretchedCone(goodConeInds(1)-1) = meanStretchedCone(goodConeInds(1))/2;
    stdStretchedCone(goodConeInds(1)-1) = stdStretchedCone(goodConeInds(1))/2;
    
    plot(coneXRef(goodConeInds(1)-1), meanStretchedCone(goodConeInds(1)-1)*voxSize, 'mx');
end

if reduceTipsCinC
    % Reduce CinC tip
    goodCinCInds = find(~isnan(meanStretchedCinC));
    
    meanStretchedCinC(goodCinCInds(1)-1) = meanStretchedCinC(goodCinCInds(1))/2;
    stdStretchedCinC(goodCinCInds(1)-1) = stdStretchedCinC(goodCinCInds(1))/2;
    
    plot(coneXRef(goodCinCInds(1)-1), meanStretchedCinC(goodCinCInds(1)-1)*voxSize, 'mx');
end

subplot(1,3,2); hold on
errorbar(coneXRef, meanStretchedExposedIntercone*voxSize, stdStretchedExposedIntercone*voxSize);
errorbar(coneXRef, meanStretchedInternalInterconeLeft*voxSize, stdStretchedInternalInterconeLeft*voxSize);
errorbar(coneXRef, meanStretchedInternalInterconeRight*voxSize, stdStretchedInternalInterconeRight*voxSize);

subplot(1,3,3); hold on
errorbar(corneaXRef, meanStretchedEpicorneaInner*voxSize, stdStretchedEpicorneaInner*voxSize);

if reduceTipsCornea
    % reduce epicornea cone tip
    goodEpicorneaInds = find(~isnan(meanStretchedEpicorneaInner));
    
    meanStretchedEpicorneaInner(goodEpicorneaInds(end)+1) = meanStretchedEpicorneaInner(goodEpicorneaInds(end))/2;
    stdStretchedEpicorneaInner(goodEpicorneaInds(end)+1) = stdStretchedEpicorneaInner(goodEpicorneaInds(end))/2;
    
    plot(corneaXRef(goodEpicorneaInds(end)+1), meanStretchedEpicorneaInner(goodEpicorneaInds(end)+1)*voxSize, 'mx');
end

% Get profiles to use
coneProfileToUse =  meanStretchedCone + stdStretchedCone*SDMult;
coneProfileToUse(coneProfileToUse < 0) = NaN;

CinCProfileToUse =  meanStretchedCinC + stdStretchedCinC*SDMult;
CinCProfileToUse(CinCProfileToUse < 0) = NaN;

exposedInterconeProfileToUse = meanStretchedExposedIntercone + stdStretchedExposedIntercone*SDMult;
exposedInterconeProfileToUse(exposedInterconeProfileToUse < 0) = 0;

internalInterconeProfileToUseLeft = meanStretchedInternalInterconeLeft + stdStretchedInternalInterconeLeft*SDMult;
internalInterconeProfileToUseLeft(internalInterconeProfileToUseLeft < 0) = 0;

internalInterconeProfileToUseRight = meanStretchedInternalInterconeRight + stdStretchedInternalInterconeRight*SDMult;
internalInterconeProfileToUseRight(internalInterconeProfileToUseRight < 0) = 0;

epicorneaProfileToUse = meanStretchedEpicorneaInner + stdStretchedEpicorneaInner*SDMult;
epicorneaProfileToUse(epicorneaProfileToUse < 0) = NaN;

% Clip top of exposed cone on both sides
exposedInterconeProfileToUseLeft = exposedInterconeProfileToUse;

if ~isnan(interconeOnLeft)
    for i = 1:length(internalInterconeProfileToUseLeft)
        if ~isnan(exposedInterconeProfileToUseLeft(i)) & ~isnan(internalInterconeProfileToUseLeft(i))

            if exposedInterconeProfileToUseLeft(i) > internalInterconeProfileToUseLeft(i)
                exposedInterconeProfileToUseLeft(i:end) = NaN;
                break
            end
        end
    end
end

exposedInterconeProfileToUseRight = exposedInterconeProfileToUse;

if ~isnan(interconeOnRight)
    for i = 1:length(internalInterconeProfileToUseRight)
        if ~isnan(exposedInterconeProfileToUseRight(i)) & ~isnan(internalInterconeProfileToUseRight(i))

            if exposedInterconeProfileToUseRight(i) > internalInterconeProfileToUseRight(i)
                exposedInterconeProfileToUseRight(i:end) = NaN;
                break
            end
        end
    end
end

if smoothProfiles
    
    temp = isnan(coneProfileToUse);
    coneProfileToUse = smooth(coneXRef, coneProfileToUse, smoothLength);
    coneProfileToUse(temp) = NaN;
    
    temp = isnan(CinCProfileToUse);
    CinCProfileToUse = smooth(coneXRef, CinCProfileToUse, smoothLength);
    CinCProfileToUse(temp) = NaN;
    
    temp = isnan(exposedInterconeProfileToUseLeft);
    exposedInterconeProfileToUseLeft = smooth(coneXRef, exposedInterconeProfileToUseLeft, smoothLength);
    exposedInterconeProfileToUseLeft(temp) = NaN;
    
    temp = isnan(exposedInterconeProfileToUseRight);
    exposedInterconeProfileToUseRight = smooth(coneXRef, exposedInterconeProfileToUseRight, smoothLength);
    exposedInterconeProfileToUseRight(temp) = NaN;
    
    temp = isnan(internalInterconeProfileToUseLeft);
    internalInterconeProfileToUseLeft = smooth(coneXRef, internalInterconeProfileToUseLeft, smoothLengthIntercone);
    internalInterconeProfileToUseLeft(temp) = NaN;
    
    temp = isnan(internalInterconeProfileToUseRight);
    internalInterconeProfileToUseRight = smooth(coneXRef, internalInterconeProfileToUseRight, smoothLengthIntercone);
    internalInterconeProfileToUseRight(temp) = NaN;
    
    temp = isnan(epicorneaProfileToUse);
    epicorneaProfileToUse = smooth(corneaXRef, epicorneaProfileToUse, smoothLength);
    epicorneaProfileToUse(temp) = NaN;
end

if scaleUpVoxels
    tempConeXRef = coneXRef(1):newVoxSize:coneXRef(end);
    tempCorenaXRef = corneaXRef(1):newVoxSize:corneaXRef(end);
    
    coneLengthToUse = floor(coneLengthToUse)*voxScale;
    outerCorneaLengthToUse = floor(outerCorneaLengthToUse)*voxScale;
    epicorneaLengthToUse = floor(epicorneaLengthToUse)*voxScale;
    
    tipOffset = round(tipOffset*voxScale)
    
    coneProfileToUse = interp1(coneXRef, coneProfileToUse, tempConeXRef,'linear')*voxScale;
    CinCProfileToUse = interp1(coneXRef, CinCProfileToUse, tempConeXRef,'linear')*voxScale;
    
    exposedInterconeProfileToUseLeft = interp1(coneXRef, exposedInterconeProfileToUseLeft, tempConeXRef,'linear')*voxScale;
    exposedInterconeProfileToUseRight = interp1(coneXRef, exposedInterconeProfileToUseRight, tempConeXRef,'linear')*voxScale;
    internalInterconeProfileToUseLeft = interp1(coneXRef, internalInterconeProfileToUseLeft, tempConeXRef,'linear')*voxScale;
    internalInterconeProfileToUseRight = interp1(coneXRef, internalInterconeProfileToUseRight, tempConeXRef,'linear')*voxScale;
    
    epicorneaProfileToUse = interp1(corneaXRef, epicorneaProfileToUse, tempCorenaXRef,'linear')*voxScale;
    
    voxSize = newVoxSize;
    
    coneXRef = tempConeXRef;
    corneaXRef = tempCorenaXRef;
end

subplot(1,3,1); hold on
plot(coneXRef, coneProfileToUse*voxSize);
plot(coneXRef, CinCProfileToUse*voxSize);

subplot(1,3,2); hold on
plot(coneXRef, exposedInterconeProfileToUseLeft*voxSize);
plot(coneXRef, exposedInterconeProfileToUseRight*voxSize);
plot(coneXRef, internalInterconeProfileToUseLeft*voxSize);
plot(coneXRef, internalInterconeProfileToUseRight*voxSize);

subplot(1,3,3); hold on
plot(corneaXRef, epicorneaProfileToUse*voxSize);

% Make slice to fit dimensions
topEpicornea = tipOffset + ceil(coneLengthToUse) + ceil(outerCorneaLengthToUse) + ceil(epicorneaLengthToUse);

sliceSize = round([3*max(coneProfileToUse), topEpicornea + tipOffset]);

slice = zeros(sliceSize(1), sliceSize(2));
coneMap = zeros(sliceSize(1), sliceSize(2));
distMap = zeros(sliceSize(1), sliceSize(2), 'logical');

%%% Continue from here then convert to 3D
    %%% Need to update profile placment (exposed, internal intercone, epicornea cone)

slice(:) = outerValue;

slice(:,1:tipOffset) = innerValue;

passedExposedConeLeft = 0;
passedExposedConeRight = 0;
passedConeTip = 0;

% Get distance map of cone shape
distMap(:,1:tipOffset) = 1;

for i = 1:sliceSize(2)
    if i > tipOffset & i <= round(coneLengthToUse) + tipOffset

        if ~isnan(coneProfileToUse(i-tipOffset))
            passedConeTip = 1;
        else
            distMap(:,i) = 1;
        end
        
        % avoids gap at end
        if passedConeTip
            if isnan(coneProfileToUse(i-tipOffset)) 
               coneProfileToUse(i-tipOffset) = coneProfileToUse(i-tipOffset-1); 
            end

            xPosTop = round(sliceSize(1)/2 + coneProfileToUse(i-tipOffset));
            xPosBottom = round(sliceSize(1)/2 - coneProfileToUse(i-tipOffset));
                        
            distMap(xPosTop+1:end, i) = 1;
            distMap(1:xPosBottom-1, i) = 1;
        end
    end
end

distMap = bwdist(distMap,'quasi-euclidean');

passedConeTip = 0;
firstGoodConeInd = Inf;
minDifDist = 0;

for i = 1:sliceSize(2)

    % Place cone profile
    if i > tipOffset & i <= round(coneLengthToUse) + tipOffset

        if ~isnan(coneProfileToUse(i-tipOffset))
            passedConeTip = 1;
            
            if i < firstGoodConeInd
                firstGoodConeInd = i;
            end
        end

        if passedConeTip
            xPosTop = round(sliceSize(1)/2 + coneProfileToUse(i-tipOffset));
            xPosBottom = round(sliceSize(1)/2 - coneProfileToUse(i-tipOffset));
        
            % Fill between top and bottom
            if isnumeric(coneValue)
                slice(xPosBottom:xPosTop, i) = coneValue;
            elseif ischar(coneValue)
                
                tempRadius = abs((xPosBottom:xPosTop)-sliceSize(1)/2);
                
                if tipGradientCorrection
                    tempDist = distMap(xPosBottom:xPosTop,i);

                    warning('Final layer seems too low relative to rest of border')
                    %%% Maybe remove the -1?
                    
                    switch correctionType
                        case 'distance'
                            tempDist = max(tempRadius) - (tempDist' - 1)*correctionScale;
                            tempRadius(tempRadius < tempDist) = tempDist(tempRadius < tempDist);
                            
                        case 'height'
                            tempDist = max(tempRadius) - (i-firstGoodConeInd)*correctionScale;
                            tempRadius(tempRadius < tempDist) = tempDist;
                            
                        case 'both'
                            tempDist = max(tempRadius) - (((tempDist' - 1) + (i-firstGoodConeInd))/2)*correctionScale;
                            tempRadius(tempRadius < tempDist) = tempDist(tempRadius < tempDist);
                    end
                    
                end
                
                % rescale relative diameter to be 80
                tempRadius = tempRadius/coneProfileToUse(i-tipOffset)*80;
                
                % rescale length relative to full cone
                tempZ = (coneLengthToUse-(i-tipOffset))/coneLengthToUse;
                
                switch coneValue
                    case 'radial'
                        % Original
                        error('switch back to original profile')
%                         tempRI = 1.52-0.000004914*tempRadius.^2;
                        
                        % For for linear contribution
                        tempRI = 1.5+(0.02)-0.00000612*tempRadius.^2;
                        
                    case 'linear'
                        tempRI = 1.5+(tempZ*0.01+0.01);
                        
                    case 'both'
                        tempRI = 1.5+(tempZ*0.01+0.01)-0.00000612*tempRadius.^2;
                        
                end
    
                slice(xPosBottom:xPosTop, i) = tempRI;
            end
    
            coneMap(xPosBottom:xPosTop, i) = (coneLengthToUse-(i-tipOffset))/coneLengthToUse;
            
            % Place inner value
            slice(xPosTop+1:end, i) = innerValue;         
            slice(1:xPosBottom-1, i) = innerValue;
        end

        % Place exposed cone and other cone profile.
        %left
        if ~isnan(exposedInterconeProfileToUseLeft(i-tipOffset))
            % Given exposed cone fill out too its level
            xPosInterconeBottom = round(exposedInterconeProfileToUseLeft(i-tipOffset));
            slice(xPosBottom-xPosInterconeBottom:xPosBottom-1, i) = interconeValue;

            % For some the profiles don't start until some distance after the planes...
            passedExposedConeLeft = 1;

        elseif passedExposedConeLeft
            if ~isnan(interconeOnLeft)
                % avoids gap at end
                if isnan(internalInterconeProfileToUseLeft(i-tipOffset))
                   internalInterconeProfileToUseLeft(i-tipOffset) = internalInterconeProfileToUseLeft(i-tipOffset-1); 
                end

                % Fill intercone out to adjacent cone level
                xPosInterconeBottom = round(internalInterconeProfileToUseLeft(i-tipOffset));

                % Add in adjacent cone
                if xPosBottom-xPosInterconeBottom-1 >= 1
                    slice(xPosBottom-xPosInterconeBottom:xPosBottom-1, i) = interconeValue;

                    if ~isnan(coneValue)      
                        slice(1:xPosBottom-xPosInterconeBottom-1, i) = coneValue;
                    else
                        % Copy in RI from current cone     
                        slice(1:xPosBottom-xPosInterconeBottom-1, i) = fliplr(tempRI(1:xPosBottom-xPosInterconeBottom-1));
                    end
                else
                    slice(1:xPosBottom-1, i) = interconeValue;
                end
            else
                slice(1:xPosBottom-1, i) = interconeValue;
            end
        end

        %right
        if ~isnan(exposedInterconeProfileToUseRight(i-tipOffset))
            % Given exposed cone fill out too its level
            xPosInterconeTop = round(exposedInterconeProfileToUseRight(i-tipOffset));
            slice(xPosTop+1:xPosTop+xPosInterconeTop, i) = interconeValue;  

            % For some the profiles don't start until some distance after the planes...
            passedExposedConeRight = 1;
        elseif passedExposedConeRight
            if ~isnan(interconeOnRight)
                % avoids gap at end
                if isnan(internalInterconeProfileToUseRight(i-tipOffset))
                   internalInterconeProfileToUseRight(i-tipOffset) = internalInterconeProfileToUseRight(i-tipOffset-1); 
                end

                % Fill intercone out to adjacent cone level
                xPosInterconeTop = round(internalInterconeProfileToUseRight(i-tipOffset));

                % Add in adjacent cone
                if xPosTop+xPosInterconeTop+1 <= sliceSize(1)
                    slice(xPosTop+1:xPosTop+xPosInterconeTop, i) = interconeValue;  

                    if ~isnan(coneValue)
                        slice(xPosTop+xPosInterconeTop+1:end, i) = coneValue;         
                    else
                        % Copy in RI from current cone
                        slice(xPosTop+xPosInterconeTop+1:end, i) = tempRI(1:sliceSize(1)-(xPosTop+xPosInterconeTop));         
                    end
                else
                    slice(xPosTop+1:end, i) = interconeValue;  
                end
            else
                slice(xPosTop+1:end, i) = interconeValue; 
            end
        end

        % Place cone in cone profile
        if ~isnan(CinCProfileToUse(i-tipOffset))
            xPosCinCTop = round(sliceSize(1)/2 + CinCProfileToUse(i-tipOffset));
            xPosCinCBottom = round(sliceSize(1)/2 - CinCProfileToUse(i-tipOffset));

            slice(xPosCinCBottom:xPosCinCTop, i) = outerCorneaValue;
            
            coneMap(xPosCinCBottom:xPosCinCTop, i) = 0;
        end
    end

    % Place outer cornea value
    if i > round(coneLengthToUse) + tipOffset & i <= round(coneLengthToUse + outerCorneaLengthToUse) + tipOffset
        slice(:, i) = outerCorneaValue;
    end    

    % Place epicornea value
    if i > round(coneLengthToUse + outerCorneaLengthToUse) + tipOffset & ...
            i <= round(coneLengthToUse + outerCorneaLengthToUse + epicorneaLengthToUse) + tipOffset
        
        slice(:, i) = epicorneaValue;

        if ~isnan(epicorneaProfileToUse(i-round(coneLengthToUse + tipOffset + outerCorneaLengthToUse)))
            xPosCorneaTop = round(sliceSize(1)/2 + epicorneaProfileToUse(i-round(coneLengthToUse + tipOffset + outerCorneaLengthToUse)));
            xPosCorneaBottom = round(sliceSize(1)/2 - epicorneaProfileToUse(i-round(coneLengthToUse + tipOffset + outerCorneaLengthToUse)));

            slice(xPosCorneaBottom:xPosCorneaTop, i) = outerCorneaValue;
        end
    end
end

figure(mainFig);

subplot(2,3,nPlot)
if ~ischar(coneValue)
    imshow(slice'/~(max(slice(:))+1))
else
    % For tips
%     imshow((slice(40:100,20:80)'-1.45)/(1.54-1.45))
    
    imshow((slice'-1.45)/(1.54-1.45))
end

title(tText);

if writeImage
    currentDirectory = pwd; 
    cd('/Users/gavintaylor/Desktop')

    writematrix(slice,fileName) 
    cd(pwd)
end
