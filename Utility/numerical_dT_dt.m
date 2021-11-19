function [dT_dt, n_of_X, dn_partial_dX] = numerical_dT_dt(X, vol_coords, vol_inds, lens_volume, fixedRI)
    % get change in ray vector based on change in refractive index and derivative as in nishidate 2011 paper
        
    % take sampling points first
    [sampling_inds, sampling_coords, point_dists] = getsamplingpoints(vol_coords, X);

    n_of_s = lens_volume(vol_inds(sampling_inds)); %nx1 vector of refractive indices at sampling points

    % As eqn 8
    W_of_X = zeros(size(sampling_coords,1), size(sampling_coords,1)); %nxn : diagonal matrix of weights per sampling point

    for i = 1:size(sampling_coords,1)
        %distance to sampling point
        d = point_dists(i);

        %As eqn 11, quartic spline weighting function
        W_of_X(i,i) = 1-6*d^2+8*d^3-3*d^4;
    end

    % As eqn 7   
    % A bit faster if done with single vectorized calc rather than per coord.
    P_of_X = calculatetricubicbasis_vec(sampling_coords(:,1), sampling_coords(:,2), sampling_coords(:,3)); %20xn : basis vector per sampling point

    % As eqn 5
    M_of_X = P_of_X*W_of_X*P_of_X'; %20x20

    % As eqn 6
    N_of_X = P_of_X*W_of_X; %20xn
        
    % Solve eqn 4
    % Supress warnings from backslash in calling function as they delay a lot
    a_of_X = M_of_X\(N_of_X*n_of_s); %20x1 : solves M*A=N*n for a, where N*n is 20x1  
    
    % Note that pinv is a lot slower than backslash
    % Also looks like error stops decreasing as function of delta around 10^-4
        % Same if choice occurs
%     a_of_X = pinv(M_of_X)*(N_of_X*n_of_s);

    [p_of_X, dp_partial_dX] = calculatetricubicbasis(X(1), X(2), X(3), 1); %20x1 and 20x3
    
    if isempty(fixedRI)
        % As eqn 3
        n_of_X = p_of_X*a_of_X; %1x1 : refractive index 
    else
       n_of_X = fixedRI; 
    end
        
    % As eqn 9
    dn_partial_dX = dp_partial_dX*a_of_X; %3x1 : spatial gradient of refractive index

    % As eqn 2
    dT_dt = n_of_X*dn_partial_dX; %3x1 vector
end