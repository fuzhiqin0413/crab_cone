% Initially test on 2d using shells pulled from images in other script

clc; 
close all

% Currently only lines up well with slice through center of X
ringsSlice = permute(ringVolume(volumeSize(1)/2,:,:),[2 3 1]);
sliceVectorY = permute(volumeVectorY(volumeSize(1)/2,:,:),[2 3 1]);
sliceVectorZ = permute(volumeVectorZ(volumeSize(1)/2,:,:),[2 3 1]);

% Calculate angles and plot working volume
figure; 

subplot(4,1,2); hold on; axis equal

indsToPlot = find(~isnan(sliceVectorY));

[plotX, plotY] = ind2sub(volumeSize(2:3), indsToPlot);

anglesSlice = zeros(volumeSize(2:3))*NaN;

for iVec = 1:length(indsToPlot)
    
    % Just do for half sized volume
    if plotX(iVec) < volumeSize(2)/2
        ang = atan2(sliceVectorZ(indsToPlot(iVec)), sliceVectorY(indsToPlot(iVec)))/pi*180;

        if ang > 5 & ang < 175
            line( [0 4*sliceVectorY(indsToPlot(iVec))]+plotX(iVec), [0 4*sliceVectorZ(indsToPlot(iVec))]+plotY(iVec),...
                'color', cols(ringsSlice(indsToPlot(iVec)), :))
            
            anglesSlice(plotX(iVec), plotY(iVec)) = ang; 
        end
    end
end

% Some low angles, may be ones that are at 180
anglesSlice(anglesSlice > 90) = NaN;

%% Get angles from bodos image - very roughly from colour map

angIncBodo = [0, 50,   55,   60,  65,   70,   75,   80, 90];
reflectivityBodo = [0, 0.04, 0.07, 0.1, 0.13, 0.16, 0.22, 0.42, 1];

angIncFine = 1:90;
reflectivityFine = interp1(angIncBodo, reflectivityBodo, angIncFine, 'linear', 'extrap');

% figure; plot(anIncFine, reflectivityFine); hold on
% plot(angIncBodo, reflectivityBodo, 'rx')
%% Set up simple reflector in 2D

subplot(4,1,3); hold on; axis equal

% Set up start values in array with space
% xPos, yPos, xRay, yRay, intensity, color
xValues = min(plotX):2:volumeSize(2)/2;

startValues = zeros(round(length(xValues)*1.5),6)*NaN;

startValues(1:length(xValues),1) = xValues;

startValues(1:length(xValues),2) = 1;

startValues(1:length(xValues),3) = 0;

startValues(1:length(xValues),4) = 1;

startValues(1:length(xValues),5) = 1;

inds = sub2ind(volumeSize(2:3), xValues, ones(1,length(xValues)));

startValues(1:length(xValues),6) = ringsSlice(inds);

zSteps = 1:volumeSize(3);

intensitySlice = zeros(volumeSize(2:3));

% Now do in while loop, a bit dangerous...
while any(~isnan(startValues(:,1))) 
    
    % Get index to use
    firstInd = find(~isnan(startValues(:,1)));
    firstInd = firstInd(1);
    
    % Pull initial values out of array
    rayT = startValues(firstInd, 3:4);
    
    rayX = startValues(firstInd, 1:2);
   
    rayIntensity = startValues(firstInd, 5);
    
    colVal = startValues(firstInd,6);
    
    % Wipe array line
    startValues(firstInd,:) = NaN;
    
    rayPath = zeros(volumeSize(3), 2)*NaN;
    
    rayPath(1,:) = [rayX];
    
    currentStep = floor(rayX(2));
    
    go = 1;
    
    % Copied from other code
    while go
        
        % Check if in mirror voxel
        if ~isnan(anglesSlice(round(rayX(1)), round(rayX(2))))
           
           mirrorX = sliceVectorY(round(rayX(1)), round(rayX(2)));
              
           mirrorY = sliceVectorZ(round(rayX(1)), round(rayX(2)));
           
           % If so, check its pointing towards centre (for left side)
           if anglesSlice(round(rayX(1)), round(rayX(2))) > 0
              
               %%% May need to do differently for right side
               
               % Rotate inwards by -90 degrees to get mirror normal
               normalX = mirrorX*cos(pi/2) - mirrorY*sin(pi/2);
               
               normalY = mirrorX*sin(pi/2) - mirrorY*cos(pi/2);
           else
              error('?') 
           end
           
           ang2Normal = round(acos((normalX*rayT(1) + normalY*rayT(2)) / ...
               (norm([normalX normalY]) * norm(rayT)))/pi*180);
           
           if ang2Normal < 1
               ang2Normal = 1;
           elseif ang2Normal > 90
               % I don't really see how this happens...
               ang2Normal = 90;
           end

           if reflectivityFine(ang2Normal)*rayIntensity > 0.05
               % From Bram de Greve pdf for vectors
               cosI = -(normalX*rayT(1) + normalY*rayT(2));

%                figure; hold on; axis equal
%                line([0 mirrorX], [0 mirrorY], 'color', 'g');
%                line([0 normalX], [0 normalY], 'color', 'b'); 
%                line([0 -rayT(1)], [0 -rayT(2)], 'color', 'r'); 

               newRayT = rayT + 2*cosI*[normalX normalY];

               newRayT = newRayT/norm(newRayT);
               
               newRayX = rayX + newRayT;
               
%                line([0 rayT(1)], [0 rayT(2)], 'color', 'k'); 
               
               % Check ray does not leave area before saving
               if ~(any(newRayX > volumeSize(2:3)) | any(newRayX < 1))
               
                   % Push reflected ray into spare slot in array
                        % Effectively branched

                   firstInd = find(isnan(startValues(:,1)));
                   firstInd = firstInd(1);

                   % Stepped one point forward
                   %%% Possible some connectivity issues...
                   startValues(firstInd,1:2) = newRayX;

                   startValues(firstInd,3:4) = newRayT;

                   startValues(firstInd,5) = rayIntensity*reflectivityFine(ang2Normal);
                   
                   startValues(firstInd,6) = colVal;
               end
           end
           
           % Reduce intensity of transmited ray
           rayIntensity = rayIntensity * (1-reflectivityFine(ang2Normal));
        end
        
        rayX = rayX + rayT; 
        
        % record path for plotting later
        if rayX(2) >= zSteps(currentStep+1)
            %%% Indexing here is generally a bit weird
            
            if currentStep + 1 < length(zSteps)
                currentStep = currentStep + 1;
                
                % Also add to intensity slice
                intensitySlice(round(rayX(1)), round(rayX(2))) = ...
                    intensitySlice(round(rayX(1)), round(rayX(2))) + rayIntensity;
            end
            
            rayPath(currentStep,:) = rayX;
        end
        
        if any(rayX > volumeSize(2:3)) | any(rayX < 1) | rayIntensity < 0.05
           go = 0; 
        end
    end
    
    plot(rayPath(:,1), rayPath(:,2), 'color', cols(colVal, :))

    %text(rayPath(1,1), rayPath(1,2)-5, sprintf('%i', iVal))
end


xlim([xValues(1)-10 xValues(end)+10]);
axis off
title('Rays')

subplot(4,1,2); 
xlim([xValues(1)-10 xValues(end)+10])
axis off
title('Lamellae planes')

subplot(4,1,4);
imshow(flipud(intensitySlice'))
 xlim([xValues(1)-10 xValues(end)+10])
title('Ray power')
 
tempTheta = permute(thetaVolume(volumeSize(1)/2,:,:),[2 3 1]);
subplot(4,1,1);
imshow(flipud((tempTheta'+90)/180))
 xlim([xValues(1)-10 xValues(end)+10])
 title('Theta')