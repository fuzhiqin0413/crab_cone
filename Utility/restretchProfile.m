function [meanStretchedProfile, stdStretchedProfile] = restretchProfile(profiles, numProfiles, xReference, ...
        originalRadius, originalStart, originalLength, referenceRadius, referenceLength, restretchLength)

    profileStretchArray = zeros(numProfiles, restretchLength);

    for i = 1:numProfiles
       tempProfile = profiles(i,round(originalStart(i)):round(originalLength(i)));

       tempXReference = xReference(round(originalStart(i)):round(originalLength(i)));

       profileStretchArray(i,:) = interp1(tempXReference, tempProfile, ...
           tempXReference(1):(tempXReference(end)/(restretchLength-1)):tempXReference(end), 'linear')/originalRadius(i);
    end

    %%% Add in adjustment on radius...

    % Take average and then stretch to intended length
    meanStretchedProfile = mean(profileStretchArray)*referenceRadius;
    meanStretchedProfile = interp1((0:restretchLength-1), meanStretchedProfile, ...
        0:(restretchLength-1)/(round(referenceLength)-1):(restretchLength-1), 'linear');
    
    stdStretchedProfile = std(profileStretchArray)*referenceRadius;
    stdStretchedProfile = interp1((0:restretchLength-1), stdStretchedProfile, ...
        0:(restretchLength-1)/(round(referenceLength)-1):(restretchLength-1), 'linear');

end