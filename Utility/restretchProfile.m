function [meanStretchedProfile, stdStretchedProfile, horizReference] = restretchProfile(profiles, numProfiles, xReference, ...
        originalRadius, originalStart, originalLength, referenceRadius, referenceLength, restretchLength, correctMissing)

    profileStretchArray = zeros(numProfiles, restretchLength);

    for i = 1:numProfiles
       tempProfile = profiles(i,round(originalStart(i)):round(originalLength(i)));

       tempXReference = xReference(round(originalStart(i)):round(originalLength(i)));

       if correctMissing
            missingInds = find(isnan(tempProfile));
            goodInds = find(~isnan(tempProfile));

            tempProfile(missingInds) = interp1(tempXReference(goodInds), tempProfile(goodInds), ...
                tempXReference(missingInds), 'linear', NaN);
       end

       profileStretchArray(i,:) = interp1(tempXReference, tempProfile, ...
           tempXReference(1):((tempXReference(end)-tempXReference(1))/(restretchLength-1)):tempXReference(end), 'linear')/originalRadius(i);
    end

    % Test plot to check alignment
%     figure; plot(profileStretchArray')

    horizReference = 0:(restretchLength-1)/(round(referenceLength)-1):(restretchLength-1);

    % Take average and then stretch to intended length
    meanStretchedProfile = nanmean(profileStretchArray)*referenceRadius;
    meanStretchedProfile = interp1((0:restretchLength-1), meanStretchedProfile, ...
        horizReference, 'linear');
    
    stdStretchedProfile = nanstd(profileStretchArray)*referenceRadius;
    stdStretchedProfile = interp1((0:restretchLength-1), stdStretchedProfile, ...
        horizReference, 'linear');

end