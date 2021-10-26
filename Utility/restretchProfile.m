function [meanStretchedProfile, stdStretchedProfile, horizReference] = restretchProfile(profiles, numProfiles, xReference, ...
        originalRadius, originalStart, originalLength, referenceRadius, referenceLength, restretchLength, correctMissing, copyToBorder, alignProfiles)

    profileStretchArray = zeros(numProfiles, restretchLength);

    testPlot = 0;

    if alignProfiles & testPlot
         figure; subplot(1,3,1); hold on
    end

    centrePoint = zeros(numProfiles,1);

    for i = 1:numProfiles
       tempProfile = profiles(i,round(originalStart(i)):round(originalLength(i)));

       tempXReference = xReference(round(originalStart(i)):round(originalLength(i)));

       centrePoint(i) = mean(tempXReference(~isnan(tempProfile)));

       if alignProfiles & testPlot
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

    if alignProfiles & testPlot
        subplot(1,3,2)
        plot(profileStretchArray')
    end

    % Try centering by registering to most central
        % Using ICP registration for X translation
            %%% ICP is probably overkill, doing a PCA might be more sensible if primary axis ends up being on X axis
            % Would then need to normalize equal number of pts on curve.
    if alignProfiles
       [~, refProfile] = min(abs(centrePoint - mean(centrePoint))); 

       stretchedXReference = 1:restretchLength;

       tempInds = find(~isnan(profileStretchArray(refProfile,:)));

       tempRefCurve = [stretchedXReference(tempInds)', profileStretchArray(i,tempInds)', zeros(length(tempInds),1)];

       opts.Registration = 'X';
       opts.Verbose = 0;

       for i = 1:numProfiles
           tempInds = find(~isnan(profileStretchArray(i,:)));

           tempStretchedCurve = [stretchedXReference(tempInds)', profileStretchArray(i,tempInds)',  zeros(length(tempInds),1)];

        
           [fittedCurve, shiftX] = ICP_finite(tempRefCurve, tempStretchedCurve, opts);

           shiftX = -round(shiftX(1,4));

           if shiftX > 0
                profileStretchArray(i,1:end-shiftX) = profileStretchArray(i,shiftX+1:end);
                profileStretchArray(i,end-shiftX+1:end) = NaN;
           elseif shiftX < 0
                shiftX = abs(shiftX);
                profileStretchArray(i,shiftX+1:end) = profileStretchArray(i,1:end-shiftX);
                profileStretchArray(i,1:shiftX) = NaN;
           end

           subplot(1,3,2); hold on
           plot(fittedCurve(:,1), fittedCurve(:,2), 'c')
           plot(profileStretchArray(i,:), 'm')
       end
    end

    % Remove points without half number of labels - do for all
    for i = 1:restretchLength
        if sum(~isnan(profileStretchArray(:,i))) < numProfiles/2
            profileStretchArray(:,i) = NaN;
        end
    end

    horizReference = 0:(restretchLength-1)/(round(referenceLength)-1):(restretchLength-1);

    % Take average
    meanStretchedProfile = nanmean(profileStretchArray);
    
    stdStretchedProfile = nanstd(profileStretchArray);

    if alignProfiles & testPlot
        subplot(1,3,2); hold on
        plot(meanStretchedProfile, 'r')
    end

    % Correct borders if needed
    if copyToBorder == 1 % copy to end
        indsToFill = find(isnan(meanStretchedProfile));
        goodInds = find(~isnan(meanStretchedProfile));

        indsToFill(indsToFill < goodInds(end)) = [];

        meanStretchedProfile(indsToFill) = meanStretchedProfile(goodInds(end));
        stdStretchedProfile(indsToFill) = stdStretchedProfile(goodInds(end));
    elseif copyToBorder == -1 % copy to start
        indsToFill = find(isnan(meanStretchedProfile));
        goodInds = find(~isnan(meanStretchedProfile));

        indsToFill(indsToFill > goodInds(1)) = [];

        meanStretchedProfile(indsToFill) = meanStretchedProfile(goodInds(1));
        stdStretchedProfile(indsToFill) = stdStretchedProfile(goodInds(1));
    end

    % Test plot to check alignment
    if alignProfiles & testPlot

        subplot(1,3,3); hold on
        plot(profileStretchArray')
        plot(meanStretchedProfile, 'r')
    elseif testPlot
        figure; plot(profileStretchArray')
    end

    % Stretch to intended length
    meanStretchedProfile = interp1((0:restretchLength-1), meanStretchedProfile, ...
        horizReference, 'linear');
    stdStretchedProfile = interp1((0:restretchLength-1), stdStretchedProfile, ...
        horizReference, 'linear');
end