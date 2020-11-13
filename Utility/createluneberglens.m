function [lensRIVolume] = createluneberglens(radiusPixels, volumeSize)
    % As described in Nishidate 2011

    lensRIVolume = zeros(volumeSize)*NaN;
    
    % Iterate through x and y an z to assign RI based on pixels
    for i = 1:volumeSize(1)
        
        for j = 1:volumeSize(2)
            
            for k = 1:volumeSize(3)
                x = (i-volumeSize(1)/2);
                y = (j-volumeSize(2)/2);
                z = (k-volumeSize(3)/2);
                
                r = sqrt(x^2+y^2+z^2);

                if r < radiusPixels
                    lensRIVolume(i, j, k) = sqrt(2-(r/radiusPixels)^2);
                end
            end
        end
    end
end

