function [p] = calculatetricubicbasis_vec(x, y, z)
    %tricubic interpolant basis vector 
    %vectorized for all sampling coords at once.

    % As eqn 10 from nishdiate 2011 paper
%     p = [ones(1,length(x)); x; y; z; x.*y; y.*z; z.*x; x.*y.*z; ...
%             x.^2; y.^2; z.^2; x.^2.*y; x.^2.*z; y.^2.*x; y.^2.*z; z.^2.*x; z.^2.*y;...
%             x.^3; y.^3; z.^3]; % 20xs - basis at point
       
    % Could be faster?
    p = ones(20,length(x));
    p(2,:) = x;
    p(3,:) = y; 
    p(4,:) = z; 
    p(5,:) = x.*y; 
    p(6,:) = y.*z; 
    p(7,:) = z.*x;
    p(8,:) = x.*y.*z;
    p(9,:) = x.^2; 
    p(10,:) = y.^2; 
    p(11,:) = z.^2; 
    p(12,:) = x.^2.*y; 
    p(13,:) = x.^2.*z; 
    p(14,:) = y.^2.*x;
    p(15,:) = y.^2.*z;
    p(16,:) = z.^2.*x;
    p(17,:) = z.^2.*y;
    p(18,:) = x.^3;
    p(19,:) = y.^3;
    p(20,:) = z.^3;
    
end

