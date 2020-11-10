function [lensRIVolume] = createluneberglens(radiusPixels, volumeSize)
    % As described in Nishidate 2011

    lensRIVolume = zeros(volumeSize, volumeSize, volumeSize)*NaN;
    
    % Iterate through x and y an z to assign RI based on pixels
    for i = 1:volumeSize
        
        for j = 1:volumeSize
            
            for k = 1:volumeSize
                x = (i-volumeSize/2);
                y = (j-volumeSize/2);
                z = (k-volumeSize/2);
                
                r = sqrt(x^2+y^2+z^2);

                if r < radius
                    lensRIVolume(i, j, k) = sqrt(2-(r/radiusPixels)^2);
                end
            end
        end
    end
    
end

