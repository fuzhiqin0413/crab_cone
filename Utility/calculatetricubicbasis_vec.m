function [p] = calculatetricubicbasis_vec(x, y, z)
    %tricubic interpolant basis vector 
    %vectorized for all sampling coords at once.

    % As eqn 10 from nishdiate 2011 paper
    p = [ones(length(x),1), x, y, z, x.*y, y.*z, z.*x, x.*y.*z, ...
            x.^2, y.^2, z.^2, x.^2.*y, x.^2.*z, y.^2.*x, y.^2.*z, z.^2.*x, z.^2.*y,...
            x.^3, y.^3, z.^3]; % sx20 - basis at point
    
end

