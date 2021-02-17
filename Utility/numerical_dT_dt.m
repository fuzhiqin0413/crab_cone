function [dT_dt, n_of_X, dn_partial_dX] = numerical_dT_dt(X, vol_coords, vol_inds, lens_volume)
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
    P_of_X = calculatetricubicbasis_vec(sampling_coords(:,1), sampling_coords(:,2), sampling_coords(:,3)); %nx20 : basis vector per sampling point

%     P_of_X = zeros(size(sampling_coords,1),20); %nx20 : basis vector per sampling point
% 
%     for i = 1:size(sampling_coords,1)
%         %place basis function for each sampling point
% 
%         P_of_X(i,:) = calculatetricubicbasis(sampling_coords(i,1), sampling_coords(i,2), sampling_coords(i,3), 0);
%     end

    % As eqn 5
    M_of_X = P_of_X'*W_of_X*P_of_X; %20x20

    % As eqn 6
    N_of_X = P_of_X'*W_of_X; %20xn

    % Solve eqn 4
    matA = M_of_X;
    
    matB = N_of_X*n_of_s;
    
    if rank(matB) >= 3
        % Change imported from anistropic
        a_of_X = matA\matB; %20x1 : solves M*A=N*n for a, where N*n is 20x1
    else
        % Backslash might end up degenerate, use alt method
        a_of_X = pinv(matA)*matB;
    end    
    
    [p_of_X, dp_partial_dX] = calculatetricubicbasis(X(1), X(2), X(3), 1); %20x1 and 20x3

    % As eqn 3
    n_of_X = p_of_X*a_of_X; %1x1 : refractive index 

    % As eqn 9
    dn_partial_dX = dp_partial_dX*a_of_X; %3x1 : spatial gradient of refractive index

    % As eqn 2
    dT_dt = n_of_X*dn_partial_dX; %3x1 vector
end