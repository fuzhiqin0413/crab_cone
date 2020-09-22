clear
clc
close all

% Load in angle data
dataFolder = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/Data/liquid angle csv/400nm';

convertAngleConvention = 1;

% Convert to mm
voxelSize = 800/10^6;

scaleSteps = 800/voxelSize/10^6;

% Load theta angles
thetaVolume = loadcsvstack(sprintf('%s/theta', dataFolder), 0);

% Load phi angles
phiVolume = loadcsvstack(sprintf('%s/phi', dataFolder), 0);

if convertAngleConvention
    
    temp = thetaVolume;

    thetaVolume = -phiVolume;

    phiVolume = temp;

    clear temp;
end

volumeSize = size(phiVolume);

phiVals = unique(phiVolume(~isnan(phiVolume)));

figure; 
subplot(1,2,1);
imshow(phiVolume(:,:,10)/360)

subplot(1,2,2);
imshow((thetaVolume(:,:,10)+90)/180)

%% Break up into shells by layer
close all

cols = lines(100);

ringVolume = zeros(volumeSize);

for iSlice = 1; %:volumeSize(3)
   phiSlice = phiVolume(:,:,iSlice);
   
   phiMask = phiSlice;
   phiMask = isnan(phiMask);
   
   ringLayer = zeros(volumeSize(1:2));
   
   minVal = 0;
   
   ringNum = 1;
   
   figure; 
   imshow(phiVolume(:,:,iSlice)/360); hold on
   
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

          plot(tempX, tempY, '.', 'color', cols(ringNum, :))
      else
         % Click to next ring 
         minVal = 0;
         
         ringNum = ringNum + 1;
      end
   end
   
   ringVolume(:,:,iSlice) = ringLayer;
end

% Look at distribution of angles in ring
for iSlice = 1; %:volumeSize(3)
    
    figure; hold on
    
    tempRingSlice = ringVolume(:,:,iSlice);
    
    tempThetaSlice = thetaVolume(:,:,iSlice);
    
    tempPhiSlice = phiVolume(:,:,iSlice);
    
    numRings = max(tempRingSlice(:));
    
    % get points in each ring
    for jRing = 1:numRings
       ringInds = find(tempRingSlice == jRing); 
       
       % get orientation within each ring
       for kPhi = 1:length(phiVals)
           
           phiInds = find(tempPhiSlice(ringInds) == phiVals(kPhi));
           
           [n, x] = hist(tempThetaSlice(ringInds(phiInds)), -90:10:90);

           subplot(2,5,kPhi); hold on
           plot(x, n, 'color', cols(jRing, :))
           
       end
    end
end

%% Get max inclinaltion and plot
for iSlice = 1; %:volumeSize(3)
    
    figure; 
    subplot(1,2,1);
    imshow(phiVolume(:,:,iSlice)/360); hold on
    
    subplot(1,2,2);
    imshow((thetaVolume(:,:,iSlice)+90)/180); hold on
    
    tempRingSlice = ringVolume(:,:,iSlice);
    
    tempThetaSlice = thetaVolume(:,:,iSlice);
    
    tempPhiSlice = phiVolume(:,:,iSlice);
    
    numRings = max(tempRingSlice(:));
    
    % get points in each ring
    for jRing = 1:numRings
       ringInds = find(tempRingSlice == jRing); 
       
       [ringX, ringY] = ind2sub(volumeSize(1:2), ringInds);
       
       % Get zero orientation within each ring
       zeroInds = find(tempPhiSlice(ringInds) == phiVals(1));
       
       % Get orientation of each zero ind
       for kInd = 1 % length(zeroInds)
           % firstly get nearby within thickness
                % Less than phi vals should always get entire band...?
           nearbyInds = find(sqrt((ringX - ringX(zeroInds(kInd))).^2 + ...
               (ringY - ringY(zeroInds(kInd))).^2) < length(phiVals));
           
           %plot(ringX(nearbyInds), ringY(nearbyInds), '.', 'color', cols(jRing, :))
           
%            unique(tempThetaSlice(ringInds(nearbyInds)))
%            
%            unique(tempPhiSlice(ringInds(nearbyInds)))
      
           
           for mNearby = 1:length(nearbyInds)
                tempPhiPoint = tempPhiSlice(ringInds(nearbyInds(mNearby)))/180*pi;
                
                tempThetaPoint = tempThetaSlice(ringInds(nearbyInds(mNearby)))/180*pi;
               
                % Border with zero phi has nan for theta.
                if isnan(tempThetaPoint)
                   tempThetaPoint = 0; 
                end
                
                vecX = cos(tempPhiPoint)*cos(tempThetaPoint);

                vecY = sin(tempPhiPoint)*cos(tempThetaPoint);

                vecZ = sin(tempThetaPoint);

                % if pointing down, flip
                if vecZ < 0
                    %vecX = -vecX; vecY = -vecY; vecZ = -vecZ; 
                end
                
                subplot(1,2,1)
                line( [0 5*vecX]+ringX(nearbyInds(mNearby)), [0 5*vecY]+ringY(nearbyInds(mNearby)),...
                    [0 5*vecZ], 'color', cols(jRing, :))

                subplot(1,2,2)
                line( [0 5*vecX]+ringX(nearbyInds(mNearby)), [0 5*vecY]+ringY(nearbyInds(mNearby)),...
                    [0 5*vecZ], 'color', cols(jRing, :) ) 
                
           end
       end
    end
end

%% Get some shells - old
if 0

    %for phiI = 1:length(phiVals)

        %shellInds = find(phiVolume == phiVals(phiI));

        shellInds = find(thetaVolume == 0);

        [shellX, shellY, shellZ] = ind2sub(volumeSize, shellInds);

        % figure; hold on;
        % plot3(shellX, shellY, shellZ, '.')

        %% Try to break up individual shells
        shellVolume = zeros(volumeSize, 'logical');

        shellVolume(shellInds) = 1;

        connectedComponents = bwconncomp(shellVolume,26);

        shellProps = regionprops(connectedComponents, 'PixelList');

         figure; hold on

        maxValues = zeros(length(shellProps), 3);

        for i = 1:length(shellProps)
             plot3(shellProps(i).PixelList(:,1), shellProps(i).PixelList(:,2), ...
                shellProps(i).PixelList(:,3), '.');

            maxValues(i,:) = [max(shellProps(i).PixelList(:,1) - volumeSize(1)/2), ...
                max(shellProps(i).PixelList(:,3)), length(shellProps(i).PixelList(:,1))];
        end

        layers2Use = find(maxValues(:,2) > 5 & maxValues(:,3) > 100);

        layerCols = lines(max(layers2Use));

        %% Plot as 2D vectors
        figure; 

        for i = layers2Use'
            % plot for phi
            subplot(1,2,1); hold on; axis equal
            tempInds = find(shellProps(i).PixelList(:,2) == volumeSize(2)/2);

            testX = shellProps(i).PixelList(tempInds,1);
            testZ = shellProps(i).PixelList(tempInds,3);

            % Get theta values
            for j = 1:length(testX)
                tempTheta = thetaVolume(testX(j), volumeSize(2)/2, testZ(j))/180*pi;

                if testX(j) < volumeSize(2)/2
                    tempTheta = -tempTheta;
                end

                line(testX(j) + [0 2*cos(tempTheta)], testZ(j) + [0 2*sin(tempTheta)], 'color', layerCols(i,:));
            end

            % plot for theta
            subplot(1,2,2); hold on; axis equal
            tempInds = find(shellProps(i).PixelList(:,3) == 1);

            testX = shellProps(i).PixelList(tempInds,1);
            testY = shellProps(i).PixelList(tempInds,2);

            % Get theta values
            for j = 1:length(testX)
                tempPhi = phiVolume(testX(j), testY(j), 1)/180*pi;

    %             if testX(j) < volumeSize(2)/2
    %                 tempPhi = -tempPhi;
    %                 
    %             end

                line(testX(j) + [0 2*cos(tempPhi)], testY(j) + [0 2*sin(tempPhi)], 'color', layerCols(i,:));
            end
        end
    %end

    %% Now fit curve to indvidual shells
    % Just take 2D to start with. 8th order poly seemed to be required to
    % capture shape of large shell.

    %%% Should fit to exterior of all included poly to get exit threshold

    figure; hold on; axis equal

    for i = layers2Use'
        tempInds = find(shellProps(i).PixelList(:,2) == volumeSize(2)/2);

        testX = shellProps(i).PixelList(tempInds,1);
        testZ = shellProps(i).PixelList(tempInds,3);

        plot(testX, testZ, '.')    

        fittedPoly = fit(testX,testZ,'poly8');

        plot(fittedPoly, testX, testZ);
    end    
    ylim([0 300])    
end
