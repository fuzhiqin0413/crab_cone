function [lensRIVolume]  = creategradedfiber(radiusPixels, fiberLengthPixels, volumeSize, n0, alpha, voxelSize)
    % As described in Nishidate 2011 eqn. 30

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
                % Corresponds to eqn. 5.34 from Merchland 1978, where b = 1
                    % Note, 2.5 represents No^2, alpha is effectively 1/sqrt(2.5)
%                 lensRIVolume(i,j,pixelsToFill) = sqrt(2.5-(r/radiusPixels)^2);

                lensRIVolume(i,j,pixelsToFill) = sqrt(n0^2*(1-alpha^2*(r*voxelSize).^2));
            end
        end
    end
end

