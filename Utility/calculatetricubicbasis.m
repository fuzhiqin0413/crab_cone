function [p, dp] = calculatetricubicbasis(x, y, z, provide_dp)
    %tricubic interpolant basis vector 

    % As eqn 10 from nishdiate 2011 paper
    p = [1, x, y, z, x*y, y*z, z*x, x*y*z, ...
            x^2, y^2, z^2, x^2*y, x^2*z, y^2*x, y^2*z, z^2*x, z^2*y,...
            x^3, y^3, z^3]; % 1x20 - basis at point
    
    %most calls do not required dp, so just provide if required.    
    if provide_dp
        dp = zeros(3,20); %3x20 - derviatives of basis at point      

        %get derviatives
        dp(1,:) = [0, 1, 0, 0, y, 0, z, y*z,...
                    2*x, 0, 0, 2*x*y, 2*x*z, y^2, 0, z^2, 0,...
                    3*x^2, 0, 0]; % dx

        dp(2,:) = [0, 0, 1, 0, x, z, 0, x*z,...
                    0, 2*y, 0, x^2, 0, 2*y*x, 2*y*z, 0, z^2,...
                    0, 3*y^2, 0]; % dy

        dp(3,:) = [0, 0, 0, 1, 0, y, x, x*y,...
                    0, 0, 2*z, 0, x^2, 0, y^2, 2*z*x, 2*z*y,...
                    0, 0, 3*z^2]; %dz   
    else
       dp = NaN; 
    end
end

