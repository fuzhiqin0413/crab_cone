clear
clc
close all

% Load in angle data
dataFolder = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/Data/liquid angle csv/800nm_thick';

convertAngleConvention = 0;

% Convert to mm
voxelSize = 800/10^6;

scaleSteps = 800/voxelSize/10^6;

% Load theta angles
thetaVolume = loadcsvstack(sprintf('%s/theta', dataFolder), 0);

% Load phi angles
phiVolume = loadcsvstack(sprintf('%s/phi', dataFolder), 0);

if convertAngleConvention
    
    temp = thetaVolume;

    thetaVolume = 90 - phiVolume;

    phiVolume = temp;

    clear temp;
end

volumeSize = size(phiVolume);

phiVals = unique(phiVolume(~isnan(phiVolume)));

%% Get some shells
close all

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

