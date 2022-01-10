angleCols = flipud(viridis(length(incidenceAngle)));
    
outlineCol = [0.5 0.5 0.5];

spotF = figure;
rayFX = figure; 
rayFY = figure; 
tipF = figure;

acceptancePercentage = zeros(length(incidenceAngle),1);
TIRPercengtage = zeros(length(incidenceAngle),1);

focusXHeight = zeros(length(incidenceAngle),1);
focusXRadius = zeros(length(incidenceAngle),1);

focusYHeight = zeros(length(incidenceAngle),1);
focusYRadius = zeros(length(incidenceAngle),1);

colcHeight = zeros(length(incidenceAngle),1);
colcRadius = zeros(length(incidenceAngle),1);

rayReverseNum = zeros(length(incidenceAngle),1);

nPlotRows = 2; % could do dynamically, but blah...

coneTipZ = coneProfileZ(1)*voxelSize;
priorInd = find(zSteps < coneTipZ); priorInd = priorInd(end);
topConeInds = find(coneProfileZ*voxelSize > coneTipZ - exposedHeight/1000);

% Get plot points
xPlotPoints = xStartPointsOrig*0; xPlotPoints(1:plotSpacing:end) = 1;
        
xPlotPoints = [fliplr(xPlotPoints(2:end)) xPlotPoints];
 
if trace3D
    [xPlot, yPlot] = meshgrid(xPlotPoints, xPlotPoints);

    plotOrigins = xPlot(:) & yPlot(:);
end

% Ray needs to get above this height to be plotted at all
rayHeightReq = (coneProfileZ(end)+(interconeProfileZ(end)-coneProfileZ(end))/4)*voxelSize;

exposedConeInds = find(coneProfileZ*voxelSize > coneTipZ - exposedHeight/1000);

for aAngle = 1:length(incidenceAngle);

    rayPathArray = rayPathCells{aAngle}; 
    finalIntersect = finalIntersectCells{aAngle};
    finalRay = finalRayCells{aAngle};
    finalRayTRefract = finalRayTRefractCells{aAngle};
    TIRFlag = TIRFlagCells{aAngle};

    % Plot spot diagram first
    figure(spotF); 
    subplot(nPlotRows, 5, aAngle)
    hold on; axis equal
    set(gca,'TickDir','out','LineWidth', 1, 'FontSize', 20);

    % store for acceptance angle
    rayAccepted = zeros(nOrigins, 1)*NaN;
    rayForColc = zeros(nOrigins, 1)*NaN;

    baseMult = 1.5;

    for iOrigin = 1:nOrigins
        if all(~isnan(finalRay(iOrigin,:))) & all(~isnan(rayPathArray(:,:,iOrigin))) & all(~isnan(finalIntersect(iOrigin,:)))

            if finalIntersect(iOrigin,3) >= rayHeightReq
                % check if final intersect is closer than last z step    
                if abs(finalIntersect(iOrigin, 3) - coneTipZ) < abs(rayPathArray(3,priorInd,iOrigin) - coneTipZ)

                    %%% Check this final ray is used correctly

                    % if it is less than a micron away just plot directly
                    if abs(finalIntersect(iOrigin, 3) - coneTipZ) < 10^-3
                       xTip = rayPathArray(1,priorInd,iOrigin);
                       yTip = rayPathArray(2,priorInd,iOrigin);
                    else
                       % should adjust with the refracted ray? (could be intersection after...)
                       if finalIntersect(iOrigin, 3) - coneTipZ <= 0
                            % intesect is behind, can just extend refracted ray
                            zDiff = coneTipZ - finalIntersect(iOrigin, 3);
                            normRay = finalRayTRefract(iOrigin,:)/finalRayTRefract(iOrigin,3);

                            xTip = finalIntersect(iOrigin, 1) + normRay(1)*zDiff;
                            yTip = finalIntersect(iOrigin, 2) + normRay(2)*zDiff;
                       else
                            error('Need to add treatment') 
                       end
                   end
                else
                    zDiff = coneTipZ - rayPathArray(3,priorInd,iOrigin);
                    normRay = finalRay(iOrigin,:)/finalRay(iOrigin,3);

                    xTip = rayPathArray(1,priorInd,iOrigin) + normRay(1)*zDiff;
                    yTip = rayPathArray(2,priorInd,iOrigin) + normRay(2)*zDiff;
                end

                xTip = xTip - volumeSize(1)/2*voxelSize;
                yTip = yTip - volumeSize(1)/2*voxelSize;

                % Accepted if within receptor radius and angle isn't too large and left above exposed part of tip
                if sqrt(xTip.^2 + yTip.^2) <= receptorRadius/1000 & ...
                        dot(finalRay(iOrigin,:), [0 0 1])/(norm(finalRay(iOrigin,:))*norm([0 0 1])) < receptorAcceptance/180*pi & ...
                        finalIntersect(iOrigin, 3) > coneTipZ - exposedHeight/1000
                    rayAccepted(iOrigin) = 1;
                else
                    rayAccepted(iOrigin) = 0;
                end

                if plotLineCols
                    col = rayCols(iOrigin,:);
                elseif colorTIR & TIRFlag(iOrigin,1)
                    col = 'b';
                else
                   col = 'k';
                end

                if finalIntersect(iOrigin, 3) > coneTipZ - exposedHeight/1000;
                    %sqrt(xTip.^2 + yTip.^2) <= coneProfileR(1)*2.5*voxelSize
                    %if plotOrigins(iOrigin); plot(xTip*1000, yTip*1000, 'x', 'color', col); end
                    plot(xTip*1000, yTip*1000, '.', 'color', col,'markersize',6);
                    rayForColc(iOrigin) = 1;
                else
                    %if plotOrigins(iOrigin); plot(xTip*1000, yTip*1000, 'o', 'color', col); end
                    h = plot(xTip*1000, yTip*1000, 'o', 'color', col,'markersize',4); 
                    rayForColc(iOrigin) = 0;
                end

                if abs(xTip) > coneProfileR(end)*voxelSize*baseMult | ...
                       abs(yTip) > coneProfileR(end)*voxelSize*baseMult
                   % Add arrow as ray is of diagram

                   % Firstly scale to border
                   if abs(xTip) > coneProfileR(end)*voxelSize*baseMult & abs(xTip) > abs(yTip)
                       % X out of range and larger than Y
                       yTipClip = yTip*(coneProfileR(end)*voxelSize*baseMult)/abs(xTip);
                       if xTip > 0
                            xTipClip = coneProfileR(end)*voxelSize*baseMult;
                       else
                            xTipClip = -coneProfileR(end)*voxelSize*baseMult;
                       end
                   elseif abs(yTip) > coneProfileR(end)*voxelSize*baseMult & abs(xTip) < abs(yTip)
                       % Y our of range and larger than X
                       xTipClip = xTip*(coneProfileR(end)*voxelSize*baseMult)/abs(yTip);
                       if yTip > 0
                            yTipClip = coneProfileR(end)*voxelSize*baseMult;
                       else
                            yTipClip = -coneProfileR(end)*voxelSize*baseMult;
                       end
                   end

                   tipNorm = sqrt(xTipClip^2 + yTipClip^2);

                   line(xTipClip*1000+[-xTipClip/tipNorm*10 0], yTipClip*1000+[-yTipClip/tipNorm*10 0], 'linewidth', 1.5, 'color', 'k')
                end

            end
        end
    end

    if acceptanceUsingReceptor      
        acceptancePercentage(aAngle) = sum( rayAccepted(~isnan(rayAccepted)))/sum(~isnan(rayAccepted));

        viscircles([0 0],receptorRadius, 'color', 'r')
    else
       % use COLC and dont plot recetpor
       acceptancePercentage(aAngle) = sum( rayForColc(~isnan(rayForColc)))/sum(~isnan(rayForColc));
    end
       
    % Get TIR rays in COLC
    tempInds = find(rayForColc == 1);
    
    TIRPercengtage(aAngle) = sum(TIRFlag(tempInds,1)>0)/sum( rayForColc(~isnan(rayForColc)));
    
    rayForColc(isnan(rayForColc)) = 0;
    
    % Finalize spot plotting
    if ~metaData.createCylinder
        viscircles([0 0],coneProfileR(exposedConeInds(end))*1000*voxelSize, 'color', [0.5 0 0.5]);

        viscircles([0 0],coneProfileR(end)*1000*voxelSize, 'color', outlineCol);
    else
        profileMax = max(metaData.coneProfileToUse)*metaData.voxSize/metaData.voxSize3D;
        
        viscircles([0 0], profileMax*1000*voxelSize, 'color', [0.5 0 0.5]);
    end 
    
    ylim([-coneProfileR(end) coneProfileR(end)]*voxelSize*baseMult*1000)
    xlim([-coneProfileR(end) coneProfileR(end)]*voxelSize*baseMult*1000)

    set(gca, 'xtick', [-100:50:100], 'xtick', [-100:50:100])

    title(sprintf('%.1f deg',incidenceAngle(aAngle)))

    % Plot raypaths in X
    figure(rayFX); 
    subplot(nPlotRows, 5, aAngle); hold on; axis equal;

    title(sprintf('%.1f deg - X plane',incidenceAngle(aAngle)))

    % Plot borders based on original profiles
    % Cone
    if metaData.displayProfiles.CinC & ~metaData.createCylinder
        plot((volumeSize(2)/2+[-fliplr(coneProfileR) (coneProfileR)])*voxelSize,...
            [fliplr(coneProfileZ) (coneProfileZ)]*voxelSize, 'color', outlineCol, 'linewidth',2);

        plot((volumeSize(2)/2+[-coneProfileR(end) -fliplr(cInCProfileR) (cInCProfileR) coneProfileR(end)])*voxelSize,...
            [coneProfileZ(end) fliplr(cInCProfileZ) (cInCProfileZ) coneProfileZ(end)]*voxelSize, 'color', outlineCol, 'linewidth',2);
    
    elseif ~metaData.displayProfiles.CinC & ~metaData.createCylinder
        plot((volumeSize(2)/2+[-coneProfileR fliplr(coneProfileR) -coneProfileR(1)])*voxelSize,...
            [coneProfileZ fliplr(coneProfileZ) coneProfileZ(1)]*voxelSize, 'color', outlineCol, 'linewidth',2);
    
    elseif metaData.createCylinder
        plot((volumeSize(2)/2+profileMax*[-1*ones(1,length(coneProfileZ)), ones(1,length(coneProfileZ)), -1])*voxelSize,...
            [coneProfileZ fliplr(coneProfileZ) coneProfileZ(1)]*voxelSize, 'color', outlineCol, 'linewidth',2);
    end

    if ~metaData.createCylinder
        plot((volumeSize(2)/2+[-fliplr(coneProfileR(topConeInds)) (coneProfileR(topConeInds)) ])*voxelSize,...
                    [fliplr(coneProfileZ(topConeInds)) (coneProfileZ(topConeInds)) ]*voxelSize, 'color', [0.5 0 0.5], 'linewidth',3);
    else
        plot((volumeSize(2)/2+profileMax*[-1*ones(1,length(topConeInds)), ones(1,length(topConeInds))])*voxelSize,...
                    [fliplr(coneProfileZ(topConeInds)) (coneProfileZ(topConeInds)) ]*voxelSize, 'color', [0.5 0 0.5], 'linewidth',3);
    end
    
    % Others
    line([0 volumeSize(2)]*voxelSize, [1 1]*corneaZ*voxelSize, 'color', outlineCol, 'linewidth',2);

    if metaData.displayProfiles.EpicorneaCone
        plot( ([0 volumeSize(2)/2-epicorneaProfileR volumeSize(2)/2+fliplr(epicorneaProfileR) volumeSize(2)])*voxelSize,...
            [epicorneaProfileZ(1) epicorneaProfileZ fliplr(epicorneaProfileZ) epicorneaProfileZ(1)]*voxelSize, 'color', outlineCol, 'linewidth',2)
    else
        line([0 volumeSize(2)]*voxelSize, [1 1]*epicorneaZ*voxelSize, 'color', outlineCol, 'linewidth',2);
    end

    line([0 volumeSize(2)/2-coneProfileR(end)]*voxelSize, [1 1]*coneBaseZ*voxelSize, 'color', outlineCol, 'linewidth',2);
    line([volumeSize(2)/2 + coneProfileR(end) volumeSize(2)]*voxelSize, [1 1]*coneBaseZ*voxelSize, 'color', outlineCol, 'linewidth',2);

    plot( ([volumeSize(2)/2-interconeProfileR 0])*voxelSize, [interconeProfileZ interconeProfileZ(end)]*voxelSize, 'color', outlineCol, 'linewidth',2)
    plot( ([volumeSize(2)/2+interconeProfileR volumeSize(2)])*voxelSize, [interconeProfileZ interconeProfileZ(end)]*voxelSize, 'color', outlineCol, 'linewidth',2)

    set(gca, 'FontSize', 20)

    % scale bar
    if scaleBarsOnRayDiagram
        line ([0.025 0.025], [-0.1 0.1]+(interconeProfileZ(end)+coneProfileZ(end))/2*voxelSize, 'linewidth', 2, 'color', 'k')
        axis off
    end

    % Plot actual rays on top
    % limit number of lines plotted
    tempInds = find(raysOnXZPlane);
    tempRaysOnXZPlane = raysOnXZPlane*0;
    tempRaysOnXZPlane(tempInds(1:plotSpacing:end)) = 1;

    for iOrigin = 1:nOrigins
        if finalIntersect(iOrigin,3) >= rayHeightReq
            rayPath = permute(rayPathArray(:, :, iOrigin), [2 1]);

            % plot3(rayPath(:,1), rayPath(:,2), rayPath(:,3), 'color', rayCols(iOrigin,:));
            % shifted to 2 x 2d plots for assymetry

%                 if finalRay(iOrigin,3) > 0
%                     extendedRay = rayPath(end-1,:) + extendRayLength*finalRay(iOrigin,:);
%                 else
%                     extendedRay = rayPath(end-1,:);
%                 end
            %line([rayPath(end-1,1) extendedRay(1)], [rayPath(end-1,2) extendedRay(2)], [rayPath(end-1,3) extendedRay(3)], 'color', rayCols(iOrigin,:) )

            if rayForColc(iOrigin); style = '-'; else; style = ':'; end

            if plotLineCols;
                col = rayCols(iOrigin,:); 
            elseif colorTIR & TIRFlag(iOrigin,1)
                col = 'b';     
            else
                col = 'k'; 
            end


            if (~justPlotCenterRays & plotOrigins(iOrigin)) | (justPlotCenterRays & tempRaysOnXZPlane(iOrigin))
                plot(rayPath(:,1),  rayPath(:,3), 'color', col, 'linestyle', style);
%                     line([rayPath(end-1,1) extendedRay(1)], [rayPath(end-1,3) extendedRay(3)], 'color', col, 'linestyle', style)
            end
        end
    end

    if aAngle == 1 & plotRIImageOnRayDiagram
        subplot(nPlotRows, 5, nPlotRows*5); 
        imshow(fliplr(permute((lensRIVolume(round(volumeSize(1)/2),:,:)-1.45)/(1.54-1.45), [2 3 1]))');

        ylim([-(volumeSize(3)*1/(volumeSize(3)*voxelSize)-volumeSize(3)) volumeSize(3)])
    end

    % plot in Y
    figure(rayFY); 

    subplot(nPlotRows, 5, aAngle); 
    hold on; axis equal;

    if incidenceAngle(aAngle) == 0 | 1
        title(sprintf('%.1f deg - Y plane',incidenceAngle(aAngle)))
    else
        title(sprintf('%.1f deg - Y plane (Rays lean in)',incidenceAngle(aAngle))) 
    end

    % Plot borders from original profiles
    %Cone
    if metaData.displayProfiles.CinC & ~metaData.createCylinder
        plot((volumeSize(2)/2+[-fliplr(coneProfileR) (coneProfileR)])*voxelSize,...
            [fliplr(coneProfileZ) (coneProfileZ)]*voxelSize, 'color', outlineCol, 'linewidth',2);

        plot((volumeSize(2)/2+[-coneProfileR(end) -fliplr(cInCProfileR) (cInCProfileR) coneProfileR(end)])*voxelSize,...
            [coneProfileZ(end) fliplr(cInCProfileZ) (cInCProfileZ) coneProfileZ(end)]*voxelSize, 'color', outlineCol, 'linewidth',2);
    
    elseif ~metaData.displayProfiles.CinC & ~metaData.createCylinder
        plot((volumeSize(2)/2+[-coneProfileR fliplr(coneProfileR) -coneProfileR(1)])*voxelSize,...
            [coneProfileZ fliplr(coneProfileZ) coneProfileZ(1)]*voxelSize, 'color', outlineCol, 'linewidth',2);
    
    elseif metaData.createCylinder
        plot((volumeSize(2)/2+profileMax*[-1*ones(1,length(coneProfileZ)), ones(1,length(coneProfileZ)), -1])*voxelSize,...
            [coneProfileZ fliplr(coneProfileZ) coneProfileZ(1)]*voxelSize, 'color', outlineCol, 'linewidth',2);
    end

    if ~metaData.createCylinder
        plot((volumeSize(2)/2+[-fliplr(coneProfileR(topConeInds)) (coneProfileR(topConeInds)) ])*voxelSize,...
                    [fliplr(coneProfileZ(topConeInds)) (coneProfileZ(topConeInds)) ]*voxelSize, 'color', [0.5 0 0.5], 'linewidth',3);
    else
        plot((volumeSize(2)/2+profileMax*[-1*ones(1,length(topConeInds)), ones(1,length(topConeInds))])*voxelSize,...
                    [fliplr(coneProfileZ(topConeInds)) (coneProfileZ(topConeInds)) ]*voxelSize, 'color', [0.5 0 0.5], 'linewidth',3);
    end
    
    % Others
    line([0 volumeSize(1)]*voxelSize, [1 1]*corneaZ*voxelSize, 'color', outlineCol, 'linewidth',2);

    if metaData.displayProfiles.EpicorneaCone
        plot( ([0 volumeSize(2)/2-epicorneaProfileR volumeSize(2)/2+fliplr(epicorneaProfileR) volumeSize(2)])*voxelSize,...
            [epicorneaProfileZ(1) epicorneaProfileZ fliplr(epicorneaProfileZ) epicorneaProfileZ(1)]*voxelSize, 'color', outlineCol, 'linewidth',2)
    else
        line([0 volumeSize(2)]*voxelSize, [1 1]*epicorneaZ*voxelSize, 'color', outlineCol, 'linewidth',2);
    end

    line([0 volumeSize(1)/2-coneProfileR(end)]*voxelSize, [1 1]*coneBaseZ*voxelSize, 'color', outlineCol, 'linewidth',2);
    line([volumeSize(1)/2 + coneProfileR(end) volumeSize(1)]*voxelSize, [1 1]*coneBaseZ*voxelSize, 'color', outlineCol, 'linewidth',2);

    plot( ([volumeSize(1)/2-interconeProfileR 0])*voxelSize, [interconeProfileZ interconeProfileZ(end)]*voxelSize, 'color', outlineCol, 'linewidth',2)
    plot( ([volumeSize(1)/2+interconeProfileR volumeSize(1)])*voxelSize, [interconeProfileZ interconeProfileZ(end)]*voxelSize, 'color', outlineCol, 'linewidth',2)

    set(gca, 'FontSize', 20)

    % scale bar
    if scaleBarsOnRayDiagram
        line ([0.025 0.025], [-0.1 0.1]+(interconeProfileZ(end)+coneProfileZ(end))/2*voxelSize, 'linewidth', 2, 'color', 'k')
        axis off
    end

    % Plot actual rays on top
    % limit number of lines plotted
    tempInds = find(raysOnYZPlane);
    tempRaysOnYZPlane = raysOnYZPlane*0;
    tempRaysOnYZPlane(tempInds(1:plotSpacing:end)) = 1;

    for iOrigin = 1:nOrigins
        if finalIntersect(iOrigin,3) >= rayHeightReq
            rayPath = permute(rayPathArray(:, :, iOrigin), [2 1]);

            %%% Check why this is happening - it shouldn't be
%                 if finalRay(iOrigin,3) > 0
%                     extendedRay = rayPath(end-1,:) + extendRayLength*finalRay(iOrigin,:);
%                 else
%                     extendedRay = rayPath(end-1,:);
%                 end

            if rayForColc(iOrigin); style = '-'; else; style = ':'; end

            if plotLineCols;
                col = rayCols(iOrigin,:); 
            elseif colorTIR & TIRFlag(iOrigin,1)
                col = 'b';     
            else
                col = 'k'; 
            end

            if (~justPlotCenterRays & plotOrigins(iOrigin)) | (justPlotCenterRays & tempRaysOnYZPlane(iOrigin))
                plot(rayPath(:,2),  rayPath(:,3), 'color', col, 'linestyle', style);
%                     line([rayPath(end-1,2) extendedRay(2)], [rayPath(end-1,3) extendedRay(3)], 'color', col, 'linestyle', style);
            end
        end
    end

    if aAngle == 1 & plotRIImageOnRayDiagram
        subplot(nPlotRows, 5, nPlotRows*5); 
        imshow(fliplr(permute((lensRIVolume(round(volumeSize(1)/2),:,:)-1.45)/(1.54-1.45), [2 3 1]))');

        ylim([-(volumeSize(3)*1/(volumeSize(3)*voxelSize)-volumeSize(3)) volumeSize(3)])
    end

    % Get lines of least confusion
    raysToUse = find(rayForColc);

    if ~isempty(raysToUse)
        focusXRadius(aAngle) = volumeSize(2)*voxelSize;
        focusYRadius(aAngle) = volumeSize(2)*voxelSize;
        colcRadius(aAngle) = volumeSize(2)*voxelSize;

        if aAngle == 1 
            topStep = length(zSteps);
        else
            topStep = colcHeight(aAngle-1)/voxelSize;
        end
        
        for iStep = 1:topStep
            tempXRad = max(rayPathArray(1, iStep, raysToUse)) - min(rayPathArray(1, iStep, raysToUse));
            if focusXRadius(aAngle) > tempXRad
                focusXRadius(aAngle) = tempXRad;
                focusXHeight(aAngle) = zSteps(iStep);
                minXPoints = [min(rayPathArray(1, iStep, raysToUse)) max(rayPathArray(1, iStep, raysToUse))];
            end

            tempYRad = max(rayPathArray(2, iStep, raysToUse)) - min(rayPathArray(2, iStep, raysToUse));
            if focusYRadius(aAngle) > tempYRad
                focusYRadius(aAngle) = tempYRad;
                focusYHeight(aAngle) = zSteps(iStep);
                minYPoints = [min(rayPathArray(2, iStep, raysToUse)) max(rayPathArray(2, iStep, raysToUse))];
            end

            [tempColcRad, tempColcCenter] = ExactMinBoundCircle(permute(rayPathArray(1:2, iStep, raysToUse), [3 1 2]));
            if colcRadius(aAngle) > tempColcRad
                colcRadius(aAngle) = tempColcRad;
                colcHeight(aAngle) = zSteps(iStep);
                minColcCenter = tempColcCenter;
            end
        end

        figure(rayFX); 
        subplot(nPlotRows, 5, aAngle); 
%             line(minXPoints, [1 1]*focusXHeight(aAngle), 'color', angleCols(aAngle,:), 'linewidth', 4)

        line(minColcCenter(1) + colcRadius(aAngle)*[-1 1], [1 1]*colcHeight(aAngle), 'color', angleCols(aAngle,:), 'linewidth', 4, 'linestyle', '-')

        figure(rayFY); 
        subplot(nPlotRows, 5, aAngle); 
%             line(minYPoints, [1 1]*focusYHeight(aAngle), 'color', angleCols(aAngle,:), 'linewidth', 4)

        line(minColcCenter(2) + colcRadius(aAngle)*[-1 1], [1 1]*colcHeight(aAngle), 'color', angleCols(aAngle,:), 'linewidth', 4, 'linestyle', '-')

        % Plot on tip
        figure(tipF)

        if aAngle == 1

            for i = 1:3
                subplot(1,3,i); hold on; axis equal

                if ~metaData.createCylinder
                    plot((volumeSize(2)/2+[-coneProfileR fliplr(coneProfileR) -coneProfileR(1)])*voxelSize,...
                        [coneProfileZ fliplr(coneProfileZ) coneProfileZ(1)]*voxelSize, 'color', outlineCol, 'linewidth',2);

                    plot((volumeSize(2)/2+[-fliplr(coneProfileR(topConeInds)) (coneProfileR(topConeInds)) ])*voxelSize,...
                        [fliplr(coneProfileZ(topConeInds)) (coneProfileZ(topConeInds)) ]*voxelSize, 'color', [0.5 0 0.5], 'linewidth',3);
                else
                    plot((volumeSize(2)/2+profileMax*[-1*ones(1,length(coneProfileZ)), ones(1,length(coneProfileZ)), -1])*voxelSize,...
                        [coneProfileZ fliplr(coneProfileZ) coneProfileZ(1)]*voxelSize, 'color', outlineCol, 'linewidth',2);
                    
                    plot((volumeSize(2)/2+profileMax*[-1*ones(1,length(topConeInds)), ones(1,length(topConeInds))])*voxelSize,...
                        [fliplr(coneProfileZ(topConeInds)) (coneProfileZ(topConeInds)) ]*voxelSize, 'color', [0.5 0 0.5], 'linewidth',3);
                end
                
                plot( ([volumeSize(1)/2-interconeProfileR 0])*voxelSize, [interconeProfileZ interconeProfileZ(end)]*voxelSize, 'color', outlineCol, 'linewidth',2)
                plot( ([volumeSize(1)/2+interconeProfileR volumeSize(1)])*voxelSize, [interconeProfileZ interconeProfileZ(end)]*voxelSize, 'color', outlineCol, 'linewidth',2)

            end
        end

        subplot(1,3,1);
        line(mean(minColcCenter) + colcRadius(aAngle)*[-1 1], [1 1]*colcHeight(aAngle), 'color', angleCols(aAngle,:), 'linewidth', 3)

        subplot(1,3,2);
        line(minXPoints, [1 1]*focusXHeight(aAngle), 'color', angleCols(aAngle,:), 'linewidth', 3)

        subplot(1,3,3);
        line(minYPoints, [1 1]*focusYHeight(aAngle), 'color', angleCols(aAngle,:), 'linewidth', 3)
    else
        focusXRadius(aAngle) = NaN;
        focusYRadius(aAngle) = NaN;
        colcRadius(aAngle) = NaN;

        focusXHeight(aAngle) = NaN;
        focusYHeight(aAngle) = NaN;
        colcHeight(aAngle) = NaN;
    end

    rayReverseNum(aAngle) = sum(rayReverseCells{aAngle});
end

figure(rayFX);
for aAngle = 1:length(incidenceAngle)
     subplot(nPlotRows, 5, aAngle);

    xlim([0 volumeSize(2)*voxelSize])
    ylim([corneaZ*voxelSize-0.05 coneTipZ+0.15])
end

figure(rayFY);
for aAngle = 1:length(incidenceAngle)
     subplot(nPlotRows, 5, aAngle);

    xlim([0 volumeSize(2)*voxelSize])
    ylim([corneaZ*voxelSize-0.05 coneTipZ+0.15])
end

figure(tipF)
for i = 1:3
    subplot(1,3,i);
    ylim(coneTipZ+[-0.25 0.05])
    
    if ~metaData.createCylinder
        xlim(volumeSize(1)/2*voxelSize+[-0.1 0.1])
    else
        xlim(volumeSize(1)/2*voxelSize+[-0.125 0.125])
    end

    set(gca, 'YTick', coneTipZ+(-0.25:0.05:0.05), 'YTickLabel', -250:50:50, ...
        'XTick', volumeSize(1)/2*voxelSize+(-0.1:0.05:0.1), 'XTickLabel', -100:50:100, ...
        'TickDir','out', 'LineWidth', 1, 'FontSize', 20)

    colormap(angleCols);
    cH = colorbar;
    set(cH, 'TickLabels', incidenceAngle, 'Ticks', ((1:length(incidenceAngle))-1)/(length(incidenceAngle)-1),'TickDir','out')
end

subplot(1,3,1);
title('COLC by angle')

subplot(1,3,2);
title('X Focus by angle')

subplot(1,3,3);
title('Y Focus by angle')

% for some reason this didn't work on cylinder so solved directly. 
    % Should take last crossover if there is some fluctuation.
% [~, uInd] = unique(acceptancePercentage);
% uInd = sort(uInd);
% acceptanceAngle = interp1(acceptancePercentage(uInd), incidenceAngle(uInd), 0.5, 'linear');

lowInds = find(acceptancePercentage > 0.5);
if ~isempty(lowInds)
    % Take last
    lowInds = lowInds(end);
%     if  acceptancePercentage(lowInds + 1) > 0
        slope = (acceptancePercentage(lowInds + 1) - acceptancePercentage(lowInds))/...
            (incidenceAngle(lowInds + 1) - incidenceAngle(lowInds));
        acceptanceAngle = incidenceAngle(lowInds) + (0.5 - acceptancePercentage(lowInds))/slope
%     end
end

figure(spotF); set(gcf, 'position', [-1919 -149 1920 1104])
figure(rayFX); set(gcf, 'position', [-1919 -149 1920 1104])
figure(rayFY); set(gcf, 'position', [-1919 -149 1920 1104])
figure(tipF); set(gcf, 'position', [1 86 1680 869]) 

% Plot final results
sumF = figure; set(gcf, 'position', [1 86 1263 869]) 
subplot(3,3,1); hold on
plot(incidenceAngle, acceptancePercentage, 'k-', 'linewidth',2)
plot(acceptanceAngle, 0.5, 'kd','markersize',10, 'linewidth',2)
if plotColorsOnSummary
    for aAngle = 1:length(incidenceAngle)
        plot(incidenceAngle(aAngle), acceptancePercentage(aAngle), 'o', 'color', angleCols(aAngle,:))
    end
end
title('Acceptance function')
ylim([0 1]); 
if acceptanceUsingReceptor
    ylabel('% entering receptor'); 
else
    ylabel('% passing cone tip'); 
end
xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);

subplot(3,3,4); hold on
plot(incidenceAngle, TIRPercengtage, 'k-', 'linewidth',2)
if plotColorsOnSummary
    for aAngle = 1:length(incidenceAngle)
        plot(incidenceAngle(aAngle), TIRPercengtage(aAngle), 'o', 'color', angleCols(aAngle,:))
    end
end
title('TIR Inclusion')
ylim([0 1]); xlim([0 20]); ylabel('% TIR in COLC'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'YTick',[0 0.25 0.5 0.75 1],'XTick',[0 5 10 15 20]);

subplot(3,3,2); hold on
plot(incidenceAngle, (colcHeight - coneTipZ)*1000, 'k-', 'linewidth',2)
if plotColorsOnSummary
    for aAngle = 1:length(incidenceAngle)
        plot(incidenceAngle(aAngle), (colcHeight(aAngle)- coneTipZ)*1000, 'o', 'color', angleCols(aAngle,:))
    end
end
title('COLC Height from Tip')
line([0 20], [0 0], 'color', 'k', 'linewidth',2)
ylim([-250 100]); ylabel('Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', [-200 -100 0 100],'XTick',[0 5 10 15 20]);

subplot(3,3,3); hold on
plot(incidenceAngle, colcRadius*1000, 'k-', 'linewidth',2)
if plotColorsOnSummary
    for aAngle = 1:length(incidenceAngle)
        plot(incidenceAngle(aAngle), colcRadius(aAngle)*1000, 'o', 'color', angleCols(aAngle,:))
    end
end
title('COLC radius')
ylim([0 75]); ylabel('Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', [0 25 50 75],'XTick',[0 5 10 15 20]);

subplot(3,3,5); hold on
plot(incidenceAngle, (focusXHeight - coneTipZ)*1000, 'k-', 'linewidth',2)
if plotColorsOnSummary
    for aAngle = 1:length(incidenceAngle)
        plot(incidenceAngle(aAngle), (focusXHeight(aAngle) - coneTipZ)*1000, 'o', 'color', angleCols(aAngle,:))
    end
end
title('X Focus Height from Tip')
line([0 20], [0 0], 'color', 'k', 'linewidth',2)
ylim([-250 50]); ylabel('Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', -250:50:50,'XTick',[0 5 10 15 20]);

subplot(3,3,6); hold on
plot(incidenceAngle, focusXRadius*1000, 'k-', 'linewidth',2)
if plotColorsOnSummary
    for aAngle = 1:length(incidenceAngle)
        plot(incidenceAngle(aAngle), focusXRadius(aAngle)*1000, 'o', 'color', angleCols(aAngle,:))
    end
end
title('X Focus Radius')
ylim([0 75]); ylabel('Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', [0 25 50 75],'XTick',[0 5 10 15 20]);

subplot(3,3,8); hold on
plot(incidenceAngle, (focusYHeight - coneTipZ)*1000, 'k-', 'linewidth',2)
if plotColorsOnSummary
    for aAngle = 1:length(incidenceAngle)
        plot(incidenceAngle(aAngle), (focusYHeight(aAngle)- coneTipZ)*1000, 'o', 'color', angleCols(aAngle,:))
    end
end
title('Y Focus Height from Tip')
line([0 20], [0 0], 'color', 'k', 'linewidth',2)
ylim([-250 100]); ylabel('Height (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', [-200 -100 0 100],'XTick',[0 5 10 15 20]);

subplot(3,3,9); hold on
plot(incidenceAngle, focusYRadius*1000, 'k-', 'linewidth',2)
if plotColorsOnSummary
    for aAngle = 1:length(incidenceAngle)
        plot(incidenceAngle(aAngle), focusYRadius(aAngle)*1000, 'o', 'color', angleCols(aAngle,:))
    end
end
title('Y Focus Radius')
ylim([0 75]); ylabel('Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', [0 25 50 75],'XTick',[0 5 10 15 20]);

subplot(3,3,7); hold on
plot(incidenceAngle, rayReverseNum, 'k-', 'linewidth',2)
if plotColorsOnSummary
    for aAngle = 1:length(incidenceAngle)
        plot(incidenceAngle(aAngle), rayReverseNum(aAngle), 'o', 'color', angleCols(aAngle,:))
    end
end
title('Ray reversals')
ylabel('Number'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20,'XTick',[0 5 10 15 20]);