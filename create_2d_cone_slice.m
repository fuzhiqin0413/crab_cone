get_radius_profiles

%% 
%  close all

voxScale = voxSize/0.5; % Scale up factor on image for FDTD
writeLarge = 0;

% For display
nPlot = 2;
tText = ''; %'-2 SD'; 'Mean';

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
    
% Set up parameters
% Labels
% outerValue = 0;
% innerValue = 1;
% coneValue = 5;
% outerCorneaValue = 4;
% epicorneaValue = 3;
% interconeValue = 2;

% RI values
outerValue = 1.33;
innerValue = 1.34;
coneValue = 'radial'; %'radial', 'linear' 'both'
    changeConeTip = 1;
outerCorneaValue = 1.5;
epicorneaValue = 1.53;
interconeValue = 1.47;

if reduceTipsCone
    bufferOffset = 9;
else
    bufferOffset = 10;
end

if bufferOffset > bufferLength
    error('buffer offset to large')
end

% Change to ring height for cone length
coneLengthToUse = (bufferLength - bufferOffset) + mean(lengthsToPlanes(:,2)) + std(lengthsToPlanes(:,2))*SDMult;

interconeLengthToUse = (bufferLength - bufferOffset) + mean(lengthsToPlanes(:,1)) + std(lengthsToPlanes(:,1))*SDMult;

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

figure;
subplot(1,3,1); hold on
errorbar(coneXRef, meanStretchedCone, stdStretchedCone);
errorbar(coneXRef, meanStretchedCinC, stdStretchedCinC);

if reduceTipsCone
    %%% Tried interpolation but didn't work well
    % decided to add a step with half radius at top so it's not too blunt

    % reduce cone tip
    goodConeInds = find(~isnan(meanStretchedCone));
    
    meanStretchedCone(goodConeInds(1)-1) = meanStretchedCone(goodConeInds(1))/2;
    stdStretchedCone(goodConeInds(1)-1) = stdStretchedCone(goodConeInds(1))/2;
    
    plot(coneXRef(goodConeInds(1)-1), meanStretchedCone(goodConeInds(1)-1), 'mx');
end

if reduceTipsCinC
    % Reduce CinC tip
    goodCinCInds = find(~isnan(meanStretchedCinC));
    
    meanStretchedCinC(goodCinCInds(1)-1) = meanStretchedCinC(goodCinCInds(1))/2;
    stdStretchedCinC(goodCinCInds(1)-1) = stdStretchedCinC(goodCinCInds(1))/2;
    
    plot(coneXRef(goodCinCInds(1)-1), meanStretchedCinC(goodCinCInds(1)-1), 'mx');
end

subplot(1,3,2); hold on
errorbar(coneXRef, meanStretchedExposedIntercone, stdStretchedExposedIntercone);
errorbar(coneXRef, meanStretchedInternalInterconeLeft, stdStretchedInternalInterconeLeft);
errorbar(coneXRef, meanStretchedInternalInterconeRight, stdStretchedInternalInterconeRight);

subplot(1,3,3); hold on
errorbar(corneaXRef, meanStretchedEpicorneaInner, stdStretchedEpicorneaInner);

if reduceTipsCornea
    % reduce epicornea cone tip
    goodEpicorneaInds = find(~isnan(meanStretchedEpicorneaInner));
    
    meanStretchedEpicorneaInner(goodEpicorneaInds(end)+1) = meanStretchedEpicorneaInner(goodEpicorneaInds(end))/2;
    stdStretchedEpicorneaInner(goodEpicorneaInds(end)+1) = stdStretchedEpicorneaInner(goodEpicorneaInds(end))/2;
    
    plot(corneaXRef(goodEpicorneaInds(end)+1), meanStretchedEpicorneaInner(goodEpicorneaInds(end)+1), 'mx');
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
epicorneaProfileToUse(epicorneaProfileToUse < 0) = 0;

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

subplot(1,3,1); hold on
plot(coneXRef, coneProfileToUse);
plot(coneXRef, CinCProfileToUse);

subplot(1,3,2); hold on
plot(coneXRef, exposedInterconeProfileToUseLeft);
plot(coneXRef, exposedInterconeProfileToUseRight);
plot(coneXRef, internalInterconeProfileToUseLeft);
plot(coneXRef, internalInterconeProfileToUseRight);

subplot(1,3,3); hold on
plot(corneaXRef, epicorneaProfileToUse);

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

distMap = bwdist(distMap,'chessboard');

passedConeTip = 0;

for i = 1:sliceSize(2)

    % Place cone profile
    if i > tipOffset & i <= round(coneLengthToUse) + tipOffset

        if ~isnan(coneProfileToUse(i-tipOffset))
            passedConeTip = 1;
        end

        if passedConeTip
            xPosTop = round(sliceSize(1)/2 + coneProfileToUse(i-tipOffset));
            xPosBottom = round(sliceSize(1)/2 - coneProfileToUse(i-tipOffset));
        
            % Fill between top and bottom
            if isnumeric(coneValue)
                slice(xPosBottom:xPosTop, i) = coneValue;
            elseif ischar(coneValue)
                % rescale relative diameter to be 80
%                 tempRadius = abs((xPosBottom:xPosTop)-sliceSize(1)/2)/coneProfileToUse(i-tipOffset)*80;
                
                tempRadius = distMap(xPosBottom:xPosTop,i);

                %%% Adjust this correction
                if max(tempRadius) < coneProfileToUse(i-tipOffset)/4
                    tempRadius = (1-(tempRadius/max(tempRadius))+1)/(max(tempRadius)+1);
                    
                    tempRadius = tempRadius/max(tempRadius);
                else
%                     tempRadius = tempRadius/coneProfileToUse(i-tipOffset);
%                     tempRadius = -(tempRadius - max(tempRadius));
                    
                    tempRadius = abs((xPosBottom:xPosTop)-sliceSize(1)/2)/coneProfileToUse(i-tipOffset);
                end
                
                tempRadius = tempRadius*80;
                
                % rescale length relative to full cone
                tempZ = (coneLengthToUse-(i-tipOffset))/coneLengthToUse;
                
                switch coneValue
                    case 'radial'
                        % Original
%                         tempRI = 1.52-0.000004914*tempRadius.^2;
                        
                        % For for linear contribution
                        tempRI = 1.5+(0.01)-0.00000612*tempRadius.^2;
                        
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
        if i > round(interconeLengthToUse) + tipOffset

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

if 0 %changeConeTip
    indsToPlot = find(coneMap);
    
    figure;
    plot3(distMap(indsToPlot), slice(indsToPlot), coneMap(indsToPlot), '.')
end

figure(mainFig);

subplot(1,3,nPlot)
if ~ischar(coneValue)
    imshow(slice'/~(max(slice(:))+1))
else
    imshow((slice'-1.45)/(1.54-1.45))
end

title(tText);

if writeLarge
    largerImg = imresize(slice, voxScale,'nearest'); 
    
    currentDirectory = pwd; 
    cd('/Users/gavintaylor/Desktop')

    warning('check names are correct')

    writematrix(largerImg,'500_average_cone_0_sd.csv') 
    cd(pwd)

    figure;
    imshow((largerImg'-1.45)/(1.54-1.45))
end
