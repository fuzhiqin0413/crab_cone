function [average, standard, tempInds] = getProfileValue(heights, level, radius, angle, minAngleRange)
    tempInds = find(heights == level);
        
    average = NaN;
    standard = NaN;

    if ~isempty(tempInds)
        tempRadius = radius(tempInds);
        
        tempAngles = angle(tempInds);
        
        % Check angle range - aim for 3/4 of circle
        if max(tempAngles) - min(tempAngles) > minAngleRange

            average = mean(tempRadius);

            standard = std(tempRadius);
        else
            tempInds = [];
        end
    end

end

