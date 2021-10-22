function [meanStretchedProfile, stdStretchedProfile, horizReference] = restretchProfile(profiles, numProfiles, xReference, ...
        originalRadius, originalStart, originalLength, referenceRadius, referenceLength, restretchLength, correctMissing, copyToEnd, fillToLimits)

    profileStretchArray = zeros(numProfiles, restretchLength);

    if copyToEnd | fillToLimits
            figure; subplot(1,3,1); hold on
    end

    centrePoint = zeros(numProfiles,1);

    for i = 1:numProfiles
       tempProfile = profiles(i,round(originalStart(i)):round(originalLength(i)));

       tempXReference = xReference(round(originalStart(i)):round(originalLength(i)));

       centrePoint(i) = mean(tempXReference(~isnan(tempProfile)));

       if copyToEnd | fillToLimits
            plot(tempProfile)
       end

       if correctMissing
            missingInds = find(isnan(tempProfile));
            goodInds = find(~isnan(tempProfile));

            tempProfile(missingInds) = interp1(tempXReference(goodInds), tempProfile(goodInds), ...
                tempXReference(missingInds), 'linear', NaN);
       end

       profileStretchArray(i,:) = interp1(tempXReference, tempProfile, ...
           tempXReference(1):((tempXReference(end)-tempXReference(1))/(restretchLength-1)):tempXReference(end), 'linear')/originalRadius(i);
   
    end

    % Scale back to radius
    profileStretchArray = profileStretchArray*referenceRadius;

    if copyToEnd | fillToLimits
        subplot(1,3,2)
        plot(profileStretchArray')
    end

    %%% Could include for CinC and Epicornea in case they are shifted away from end?
    if copyToEnd
        for i = 1:numProfiles
            profileStretchArray(i, copyToEnd:end) = nanmean(profileStretchArray(i, copyToEnd-1:copyToEnd+1)); 
        end
    end

    % Try centering by registering to most central
        % Using ICP registration with rotation and XY translation
        % Well aligned but would ideally limit to X translation
            %%% Check if that limit is possible

            %%% ICP is probably overkill, doing a PCA might be more sensible if primary axis ends up being on X axis
            % Would then need to normalize equal number of pts on curve.
    if fillToLimits
       [~, refProfile] = min(abs(centrePoint - mean(centrePoint))); 

       stretchedXReference = 1:restretchLength;

       tempInds = find(~isnan(profileStretchArray(refProfile,:)));

       tempRefCurve = [profileStretchArray(i,tempInds)', stretchedXReference(tempInds)', zeros(length(tempInds),1)];

       for i = 1:numProfiles
           tempInds = find(~isnan(profileStretchArray(i,:)));

           tempStretchedCurve = [profileStretchArray(i,tempInds)', stretchedXReference(tempInds)', zeros(length(tempInds),1)];

           [fittedCurve, M] = ICP_finite(tempRefCurve, tempStretchedCurve);

           %%% Plot with just X shift component - then replace that for use in averages
           M

           subplot(1,3,2); hold on
           plot(fittedCurve(:,2), fittedCurve(:,1), 'm')
       end
    end

    % Try recentering using by filling out to limits
    if fillToLimits & 0
       % Get limits
       minLimits = zeros(numProfiles,1);
       maxLimits = zeros(numProfiles,1);
       
       for i = 1:numProfiles
            tempInd = find(~isnan(profileStretchArray(i,:)));

            minLimits(i) = tempInd(1);
            maxLimits(i) = tempInd(end);
       end

       minInd = min(minLimits);
       maxInd = max(maxLimits);

       % Fill periperal value out to limits
       for i = 1:numProfiles
            profileStretchArray(i, minInd:minLimits(i)) = profileStretchArray(i,minLimits(i));
            profileStretchArray(i, maxLimits(i):maxInd) = profileStretchArray(i,maxLimits(i));
       end
    end

    horizReference = 0:(restretchLength-1)/(round(referenceLength)-1):(restretchLength-1);

    % Take average and then stretch to intended length
    meanStretchedProfile = nanmean(profileStretchArray);
    
    stdStretchedProfile = nanstd(profileStretchArray);

    % Test plot to check alignment
    if copyToEnd | fillToLimits
        subplot(1,3,2); hold on
        plot(meanStretchedProfile, 'r')

        subplot(1,3,3); hold on
        plot(profileStretchArray')
        plot(meanStretchedProfile, 'r')
    else
        figure; plot(profileStretchArray')
    end

    meanStretchedProfile = interp1((0:restretchLength-1), meanStretchedProfile, ...
        horizReference, 'linear');

    stdStretchedProfile = interp1((0:restretchLength-1), stdStretchedProfile, ...
        horizReference, 'linear');
end