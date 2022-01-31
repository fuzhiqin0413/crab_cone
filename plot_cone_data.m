angleCols = flipud(viridis(length(incidenceAngle)));
    
outlineCol = [0.5 0.5 0.5];
baseMult = 1.5; % Size of spot diagram relative to base of cone

spotF = figure;
rayFX = figure; 
rayFY = figure; 
tipF = figure;

acceptancePercentageNight = zeros(length(incidenceAngle),1); % for night
acceptancePercentageDay = zeros(length(incidenceAngle),1);
acceptancePercentageColc = zeros(length(incidenceAngle),1);
TIRPercengtage = zeros(length(incidenceAngle),1);

focusXHeight = zeros(length(incidenceAngle),1);
focusXRadius = zeros(length(incidenceAngle),1);

focusYHeight = zeros(length(incidenceAngle),1);
focusYRadius = zeros(length(incidenceAngle),1);

colcHeight = zeros(length(incidenceAngle),1);
colcRadius = zeros(length(incidenceAngle),1);

rayReverseNum = zeros(length(incidenceAngle),1);

nPlotRows = 2; % could do dynamically, but blah...

priorInd = find(zSteps < coneTipZ); priorInd = priorInd(end);
exteriorInds = find(zSteps >= coneTipZ);
topConeInds = find(coneProfileZ*voxelSize > coneTipZ - exposedHeight/1000);

priorIndNightReceptor = find(zSteps < coneTipZ + receptorDistanceNight/1000); priorIndNightReceptor = priorIndNightReceptor(end);
priorIndDayReceptor = find(zSteps < coneTipZ + receptorDistanceDay/1000); priorIndDayReceptor = priorIndDayReceptor(end);
priorIndDayNarrowStart = find(zSteps < coneTipZ + dayNarrowStart/1000); priorIndDayNarrowStart = priorIndDayNarrowStart(end);
priorIndDayNarrowEnd = find(zSteps < coneTipZ + dayNarrowEnd/1000); priorIndDayNarrowEnd = priorIndDayNarrowEnd(end);

coneBaseZ = coneProfileZ(end)*voxelSize;
%%% Kept following ind for legacy results where firstIntersect wasn't saved
followingInd = find(zSteps < coneBaseZ); followingInd = followingInd(end);
plotFailedBase = 1;

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
    if exist('firstIntersectCells','var') == 1
        firstIntersect = firstIntersectCells{aAngle};
        noFirstIntersect = 0;
    else
        noFirstIntersect = 1;
    end
    
    % Plot spot diagram first
    figure(spotF); 
    subplot(nPlotRows, 5, aAngle)
    hold on; axis equal
    set(gca,'TickDir','out','LineWidth', 1, 'FontSize', 20);

    % store for acceptance angle
    rayAcceptedNight = zeros(nOrigins, 1)*NaN; 
    rayAcceptedDay = zeros(nOrigins, 1)*NaN; 
    rayForColc = zeros(nOrigins, 1)*NaN;

    for iOrigin = 1:nOrigins
        % Good to test, but a bit of hack to account for the fact that there is can be a missing faces along the intercone - cone base interface (seems to just occur on cylinder)
        %%% Switch to check intersect is at height of cone base using firstIntersect when saved
       
        % problem seems to be that ray sneaks in base or around side of cone, so check for change in direction at base
        if noFirstIntersect
            rayTm1 = rayPathArray(:,followingInd,iOrigin) - rayPathArray(:,followingInd-1,iOrigin);
            rayTp1 = rayPathArray(:,followingInd+1,iOrigin) - rayPathArray(:,followingInd,iOrigin);

            baseTest = any(rayTm1 - rayTp1 ~= 0);
        else
            if all(isnan(firstIntersect(iOrigin,:))) || all(~isnan(firstIntersect(iOrigin,:))) && abs(firstIntersect(iOrigin, 3) - coneBaseZ) < 1e-9
                baseTest = 1;
                
            else
                baseTest = 0;
            end
        end
        
        if metaData.displayProfiles.CinC
            % limit to near border otherwise rays going into CinC will flag either test
            %%% Can use firstIntersect when availible
            if noFirstIntersect
                if sqrt((rayPathArray(1,followingInd,iOrigin) - volumeSize(1)/2*voxelSize)^2 + ...
                        (rayPathArray(2,followingInd,iOrigin) - volumeSize(2)/2*voxelSize)^2) < coneProfileR(end)*voxelSize*0.95

                    baseTest = 1;
                end  
            else
                if sqrt((firstIntersect(1,iOrigin) - volumeSize(1)/2*voxelSize)^2 + ...
                        (firstIntersect(2,iOrigin) - volumeSize(2)/2*voxelSize)^2) < coneProfileR(end)*voxelSize*0.95

                    baseTest = 1;
                end 
            end
        end

        if all(~isnan(finalRay(iOrigin,:))) && all(all(~isnan(rayPathArray(:,:,iOrigin)))) && ...
                all(~isnan(finalIntersect(iOrigin,:))) && baseTest

            if finalIntersect(iOrigin,3) >= rayHeightReq
                % get position at tip
                % check if final intersect is closer than last z step    
                if abs(finalIntersect(iOrigin, 3) - coneTipZ) < abs(rayPathArray(3,priorInd,iOrigin) - coneTipZ)
                    % if it is less than a micron away just plot directly
                    if abs(finalIntersect(iOrigin, 3) - coneTipZ) < 10^-3
                       xTip = finalIntersect(iOrigin, 1) - volumeSize(1)/2*voxelSize;
                       yTip = finalIntersect(iOrigin, 2) - volumeSize(2)/2*voxelSize;
                    else
                       % should adjust with the refracted ray? (could be intersection after...)
                       if finalIntersect(iOrigin, 3) - coneTipZ <= 0
                            % intesect is behind, can just extend refracted ray
                            zDiff = coneTipZ - finalIntersect(iOrigin, 3);
                            normRay = finalRayTRefract(iOrigin,:)/finalRayTRefract(iOrigin,3);

                            xTip = finalIntersect(iOrigin, 1) + normRay(1)*zDiff - volumeSize(1)/2*voxelSize;
                            yTip = finalIntersect(iOrigin, 2) + normRay(2)*zDiff - volumeSize(2)/2*voxelSize;
                       else
                            error('Need to add treatment') 
                       end
                   end
                else
                    zDiff = coneTipZ - rayPathArray(3,priorInd,iOrigin);
                    normRay = finalRay(iOrigin,:)/finalRay(iOrigin,3);

                    xTip = rayPathArray(1,priorInd,iOrigin) + normRay(1)*zDiff - volumeSize(1)/2*voxelSize;
                    yTip = rayPathArray(2,priorInd,iOrigin) + normRay(2)*zDiff - volumeSize(2)/2*voxelSize;
                end
                

                % get position at night receptor
                zDiff = coneTipZ + receptorDistanceNight/1000 - rayPathArray(3,priorIndNightReceptor,iOrigin);
                normRay = finalRay(iOrigin,:)/finalRay(iOrigin,3);

                xTipNight = rayPathArray(1,priorIndNightReceptor,iOrigin) + normRay(1)*zDiff - volumeSize(1)/2*voxelSize;
                yTipNight = rayPathArray(2,priorIndNightReceptor,iOrigin) + normRay(2)*zDiff - volumeSize(2)/2*voxelSize;
                
                % get position at day receptor, and both narrows
                zDiff = coneTipZ + receptorDistanceDay/1000 - rayPathArray(3,priorIndDayReceptor,iOrigin);
                normRay = finalRay(iOrigin,:)/finalRay(iOrigin,3);

                xTipDay = rayPathArray(1,priorIndDayReceptor,iOrigin) + normRay(1)*zDiff - volumeSize(1)/2*voxelSize;
                yTipDay = rayPathArray(2,priorIndDayReceptor,iOrigin) + normRay(2)*zDiff - volumeSize(2)/2*voxelSize;
                
                zDiff = coneTipZ + dayNarrowStart/1000 - rayPathArray(3,priorIndDayNarrowStart,iOrigin);
                normRay = finalRay(iOrigin,:)/finalRay(iOrigin,3);

                xTipDayNarrowStart = rayPathArray(1,priorIndDayNarrowStart,iOrigin) + normRay(1)*zDiff - volumeSize(1)/2*voxelSize;
                yTipDayNarrowStart = rayPathArray(2,priorIndDayNarrowStart,iOrigin) + normRay(2)*zDiff - volumeSize(2)/2*voxelSize;
                
                zDiff = coneTipZ + dayNarrowEnd/1000 - rayPathArray(3,priorIndDayNarrowEnd,iOrigin);
                normRay = finalRay(iOrigin,:)/finalRay(iOrigin,3);

                xTipDayNarrowEnd = rayPathArray(1,priorIndDayNarrowEnd,iOrigin) + normRay(1)*zDiff - volumeSize(1)/2*voxelSize;
                yTipDayNarrowEnd = rayPathArray(2,priorIndDayNarrowEnd,iOrigin) + normRay(2)*zDiff - volumeSize(2)/2*voxelSize;
                
                % Night accepted if within receptor radius and left above exposed part of tip and within aperture at tip
                if sqrt(xTipNight.^2 + yTipNight.^2) <= receptorRadiusNight/1000 & ...
                        sqrt(xTip.^2 + yTip.^2) <= apertureRadiusNight/1000 & ...
                        finalIntersect(iOrigin, 3) > coneTipZ - exposedHeight/1000
                    rayAcceptedNight(iOrigin) = 1;
                else
                    rayAcceptedNight(iOrigin) = 0;
                end
                
                % Day accepted if within receptor radius and left above exposed part of tip and within aperture at both narrow  points
                if sqrt(xTipDay.^2 + yTipDay.^2) <= receptorRadiusDay/1000 & ...
                        sqrt(xTipDayNarrowStart.^2 + yTipDayNarrowStart.^2) <= apertureRadiusDay/1000 & ...
                        sqrt(xTipDayNarrowEnd.^2 + yTipDayNarrowEnd.^2) <= apertureRadiusDay/1000 & ...
                        finalIntersect(iOrigin, 3) > coneTipZ - exposedHeight/1000
                    rayAcceptedDay(iOrigin) = 1;
                else
                    rayAcceptedDay(iOrigin) = 0;
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
        elseif ~baseTest
            % Nan final intersect so it's not plotted
            finalIntersect(iOrigin,3) = NaN;
        end
        
        if ~baseTest & plotFailedBase
            temp = gcf;
            figure; hold on; axis equal
%             trisurf(grinSurface.faces, grinSurface.vertices(:,1), grinSurface.vertices(:,2), grinSurface.vertices(:,3), ...
%                        'FaceColor','g','FaceAlpha',0.1);
            % Just plot circle as surface isn't saved
            if ~metaData.createCylinder
                r = coneProfileR(end)*voxelSize;
            else
                r = profileMax*voxelSize;
            end
            
            teta=-pi:0.01:pi;
            x = r*cos(teta)+volumeSize(1)/2*voxelSize;
            y = r*sin(teta)+volumeSize(2)/2*voxelSize;
            plot3(x,y,ones(1,numel(x))*coneBaseZ)
            
            plot3(rayPathArray(1,:,iOrigin), rayPathArray(2,:,iOrigin), rayPathArray(3,:,iOrigin),'r','linewidth',1)
            plot3(rayPathArray(1,followingInd,iOrigin), rayPathArray(2,followingInd,iOrigin), rayPathArray(3,followingInd,iOrigin),'mx','linewidth',1)
            
            title(sprintf('%i %i', aAngle, iOrigin))
            figure(temp)
        end
    end

    %Get accpetance angle for all
    acceptancePercentageNight(aAngle) = sum( rayAcceptedNight(~isnan(rayAcceptedNight)))/sum(~isnan(rayAcceptedNight));
     
    acceptancePercentageDay(aAngle) = sum( rayAcceptedDay(~isnan(rayAcceptedDay)))/sum(~isnan(rayAcceptedDay));

    viscircles([0 0],receptorRadiusDay, 'color', 'r', 'linewidth', 2, 'linestyle', ':')
    
    viscircles([0 0],receptorRadiusNight, 'color', 'r', 'linewidth', 2)

    acceptancePercentageColc(aAngle) = sum( rayForColc(~isnan(rayForColc)))/sum(~isnan(rayForColc));

    
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
%                 plot3(rayPath(:,1), rayPath(:,2), rayPath(:,3), 'color', col, 'linestyle', style);
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

        %%% Had a limit that forced the colc to descend as the angle increased, was removed
        if limitCOLCToExterior
            iStep = exteriorInds;
        else
            iStep = 1:length(zSteps);
        end
        
        for iStep = iStep;
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

    if limitCOLCToExterior
       if ~metaData.createCylinder
            xlim(volumeSize(1)/2*voxelSize+[-0.075 0.075])
        else
            xlim(volumeSize(1)/2*voxelSize+[-0.125 0.125])
        end
        
        ylim(coneTipZ+[-0.025 0.05])     
       
       set(gca, 'YTick', coneTipZ+(-0.025:0.025:0.05), 'YTickLabel', -25:25:50, ...
            'XTick', volumeSize(1)/2*voxelSize+(-0.1:0.05:0.1), 'XTickLabel', -100:50:100, ...
            'TickDir','out', 'LineWidth', 1, 'FontSize', 20)
    else
        if ~metaData.createCylinder
            xlim(volumeSize(1)/2*voxelSize+[-0.1 0.1])
        else
            xlim(volumeSize(1)/2*voxelSize+[-0.125 0.125])
        end
        
        ylim(coneTipZ+[-0.25 0.05])
        set(gca, 'YTick', coneTipZ+(-0.25:0.05:0.05), 'YTickLabel', -250:50:50, ...
            'XTick', volumeSize(1)/2*voxelSize+(-0.1:0.05:0.1), 'XTickLabel', -100:50:100, ...
            'TickDir','out', 'LineWidth', 1, 'FontSize', 20)
    end

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

% Get accpetance angles
acceptancePercentageNight = acceptancePercentageNight/acceptancePercentageNight(1);
lowInds = find(acceptancePercentageNight > 0.5);
% Take last
lowInds = lowInds(end);
    slope = (acceptancePercentageNight(lowInds + 1) - acceptancePercentageNight(lowInds))/...
        (incidenceAngle(lowInds + 1) - incidenceAngle(lowInds));
    acceptanceAngleNight = incidenceAngle(lowInds) + (0.5 - acceptancePercentageNight(lowInds))/slope;

acceptancePercentageDay = acceptancePercentageDay/acceptancePercentageDay(1);
lowInds = find(acceptancePercentageDay > 0.5);
if ~isempty(lowInds)
    % Take last
    lowInds = lowInds(end);
        slope = (acceptancePercentageDay(lowInds + 1) - acceptancePercentageDay(lowInds))/...
            (incidenceAngle(lowInds + 1) - incidenceAngle(lowInds));
        acceptanceAngleDay = incidenceAngle(lowInds) + (0.5 - acceptancePercentageDay(lowInds))/slope;
end

acceptancePercentageColc = acceptancePercentageColc/acceptancePercentageColc(1);
lowInds = find(acceptancePercentageColc > 0.5);
if ~isempty(lowInds)
    % Take last
    lowInds = lowInds(end);
        slope = (acceptancePercentageColc(lowInds + 1) - acceptancePercentageColc(lowInds))/...
            (incidenceAngle(lowInds + 1) - incidenceAngle(lowInds));
        acceptanceAngleColc = incidenceAngle(lowInds) + (0.5 - acceptancePercentageColc(lowInds))/slope;
end

figure(spotF); set(gcf, 'position', [-1919 -149 1920 1104])
figure(rayFX); set(gcf, 'position', [-1919 -149 1920 1104])
figure(rayFY); set(gcf, 'position', [-1919 -149 1920 1104])
figure(tipF); set(gcf, 'position', [1 86 1680 869]) 

% Plot final results
sumF = figure; set(gcf, 'position', [1 86 1263 869]) 
subplot(3,3,1); hold on
plot(incidenceAngle, acceptancePercentageNight, 'r-', 'linewidth',2)
plot(acceptanceAngleNight, 0.5, 'ro','markersize',5, 'linewidth',2, 'MarkerFaceColor', 'r')
if plotColorsOnSummary
    for aAngle = 1:length(incidenceAngle)
        plot(incidenceAngle(aAngle), acceptancePercentageNight(aAngle), 'o', 'color', angleCols(aAngle,:))
    end
end

plot(incidenceAngle, acceptancePercentageDay, 'r:', 'linewidth',2)
plot(acceptanceAngleDay, 0.5, 'ro','markersize',5, 'linewidth',2, 'MarkerFaceColor', 'w')
if plotColorsOnSummary
    for aAngle = 1:length(incidenceAngle)
        plot(incidenceAngle(aAngle), acceptancePercentageDay(aAngle), 'o', 'color', angleCols(aAngle,:))
    end
end

plot(incidenceAngle, acceptancePercentageColc, '-', 'linewidth',2, 'color', [0.5 0 0.5])
plot(acceptanceAngleColc, 0.5, 'o','color', [0.5 0 0.5], 'markersize',5, 'linewidth',2, 'MarkerFaceColor', [0.5 0 0.5])
if plotColorsOnSummary
    for aAngle = 1:length(incidenceAngle)
        plot(incidenceAngle(aAngle), acceptancePercentageColc(aAngle), 'o', 'color', angleCols(aAngle,:))
    end
end
title('Acceptance function')
ylim([0 1]); 
ylabel('% entering receptor/passing cone tip'); 
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
if limitCOLCToExterior
    ylim([-25 50]); ylabel('Height (um)'); xlabel('Angle (deg)')
    set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', [-25 0 25 50],'XTick',[0 5 10 15 20]);
else
    ylim([-250 50]); ylabel('Height (um)'); xlabel('Angle (deg)')
    set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', [-200 -100 0 50],'XTick',[0 5 10 15 20]);
end

subplot(3,3,3); hold on
plot(incidenceAngle, colcRadius*1000, 'k-', 'linewidth',2)
if plotColorsOnSummary
    for aAngle = 1:length(incidenceAngle)
        plot(incidenceAngle(aAngle), colcRadius(aAngle)*1000, 'o', 'color', angleCols(aAngle,:))
    end
end
title('COLC radius')
ylim([0 100]); ylabel('Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', [0 25 50 75 100],'XTick',[0 5 10 15 20]);

subplot(3,3,5); hold on
plot(incidenceAngle, (focusXHeight - coneTipZ)*1000, 'k-', 'linewidth',2)
if plotColorsOnSummary
    for aAngle = 1:length(incidenceAngle)
        plot(incidenceAngle(aAngle), (focusXHeight(aAngle) - coneTipZ)*1000, 'o', 'color', angleCols(aAngle,:))
    end
end
title('X Focus Height from Tip')
line([0 20], [0 0], 'color', 'k', 'linewidth',2)
if limitCOLCToExterior
    ylim([-25 50]); ylabel('Height (um)'); xlabel('Angle (deg)')
    set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', [-25 0 25 50],'XTick',[0 5 10 15 20]);
else
    ylim([-250 50]); ylabel('Height (um)'); xlabel('Angle (deg)')
    set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', [-200 -100 0 50],'XTick',[0 5 10 15 20]);
end
subplot(3,3,6); hold on
plot(incidenceAngle, focusXRadius*1000, 'k-', 'linewidth',2)
if plotColorsOnSummary
    for aAngle = 1:length(incidenceAngle)
        plot(incidenceAngle(aAngle), focusXRadius(aAngle)*1000, 'o', 'color', angleCols(aAngle,:))
    end
end
title('X Focus Radius')
ylim([0 100]); ylabel('Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', [0 25 50 75 100],'XTick',[0 5 10 15 20]);

subplot(3,3,8); hold on
plot(incidenceAngle, (focusYHeight - coneTipZ)*1000, 'k-', 'linewidth',2)
if plotColorsOnSummary
    for aAngle = 1:length(incidenceAngle)
        plot(incidenceAngle(aAngle), (focusYHeight(aAngle)- coneTipZ)*1000, 'o', 'color', angleCols(aAngle,:))
    end
end
title('Y Focus Height from Tip')
line([0 20], [0 0], 'color', 'k', 'linewidth',2)
if limitCOLCToExterior
    ylim([-25 50]); ylabel('Height (um)'); xlabel('Angle (deg)')
    set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', [-25 0 25 50],'XTick',[0 5 10 15 20]);
else
    ylim([-250 50]); ylabel('Height (um)'); xlabel('Angle (deg)')
    set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', [-200 -100 0 50],'XTick',[0 5 10 15 20]);
end

subplot(3,3,9); hold on
plot(incidenceAngle, focusYRadius*1000, 'k-', 'linewidth',2)
if plotColorsOnSummary
    for aAngle = 1:length(incidenceAngle)
        plot(incidenceAngle(aAngle), focusYRadius(aAngle)*1000, 'o', 'color', angleCols(aAngle,:))
    end
end
title('Y Focus Radius')
ylim([0 100]); ylabel('Radius (um)'); xlabel('Angle (deg)')
set(gca,'TickDir','out', 'LineWidth', 1, 'FontSize', 20, 'YTick', [0 25 50 75 100],'XTick',[0 5 10 15 20]);

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