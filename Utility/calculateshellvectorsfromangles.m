function [vecX, vecY, vecZ] = calculateshellvectorsfromangles(tempPhiPoint, tempThetaPoint)

if inclineVertical
    % Will not work if completely vertical..
    if tempThetaPoint == pi/2
        tempThetaPoint = 85/180*pi;
    elseif tempThetaPoint == -pi/2
        tempThetaPoint = -85/180*pi;
    end
end

vecX = cos(tempPhiPoint)*cos(tempThetaPoint);
vecY = sin(tempPhiPoint)*cos(tempThetaPoint);
vecZ = sin(tempThetaPoint);

% if pointing down in Z, flip
if vecZ < 0
    vecX = -vecX; vecY = -vecY; vecZ = -vecZ; 
end

% Check angle if vector is pointing towards center 
% (based on shorter distance at top of vector than base)

% Now use for flat as well
if 1; %vecZ ~= 0
    basedDist = sqrt((ringX(nearbyInds(pointToUse)) - volumeSize(1)/2)^2 + ...
        (ringY(nearbyInds(pointToUse)) - volumeSize(2)/2)^2);

    topDist = sqrt((ringX(nearbyInds(pointToUse)) + vecX - volumeSize(1)/2)^2 + ...
        (ringY(nearbyInds(pointToUse)) + vecY - volumeSize(2)/2)^2);

    if topDist > basedDist
       % Rotate by 180 degrees 

       vecX = -vecX;
       vecY = -vecY;
    end
end

end

