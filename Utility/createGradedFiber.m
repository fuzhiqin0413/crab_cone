function [lensRIVolume]  = createGradedFiber(radiusPixels, fiberLengthPixels, volumeSize)
    % As described in Nishidate 2011

    lensRIVolume = zeros(volumeSize)*NaN;
    
    if fiberLengthPixels > volumeSize(3)
       error('Fiber too long') 
    end
    
    % Centre fiber around middle of volume
    pixelsToFill = round(volumeSize(3)/2-fiberLengthPixels/2):round(volumeSize(3)/2+fiberLengthPixels/2);
     
    % Iterate through x and y to assign RI along z
    for i = 1:volumeSize(1)
        
        for j = 1:volumeSize(2)
            
            x = (i-volumeSize(1)/2);
            y = (j-volumeSize(2)/2);

            r = sqrt(x^2+y^2);

            if r < radiusPixels
                lensRIVolume(i,j,pixelsToFill) = sqrt(2.5-(r/radiusPixels)^2);
            end
        end
    end
end

