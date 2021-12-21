function [dt_ds] = numerical_dt_ds(X, T, vol_coords, vol_inds, lens_volume_seperated)
        % Get change in ray direction in discrete anistotropic medium based on christoffel symbols
            % merges method from Nishdate 2011 and 2013 papers
        
        %take sampling points first
        [sampling_inds, sampling_coords, point_dists] = getsamplingpoints(vol_coords, X);
        
        % Get epsilon values from cell array
        %epsilon_of_s = zeros(length(sampling_inds),3,3); %nx3x3 matrix of perimitivity at sampling points
        
        % Swaped 3x3 to come before point dimension for clarity 
        %%% Will use extra time on permutes, may be able to reduce later
        epsilon_of_s = zeros(3,3,length(sampling_inds)); %3x3xn matrix of permitivity at sampling points
        
        allEqualCells = zeros(3,3)*NaN;
        
        % Only unpack bottom triangal due to symmetry
        
%         for i = 1:3
%             for j = 1:i
%                 tempVolume = lens_volume_seperated{j,i};
%                 if numel(tempVolume) > 1
%                    epsilon_of_s(:,i,j) = tempVolume(vol_inds(sampling_inds)); 
%                 else
%                    % Just a single value, no need to unpack
%                    allEqualCells(i,j) = tempVolume;
%                 end
%             end
%         end

        for i = 1:3
            for j = 1:i
                tempVolume = lens_volume_seperated{i,j};
                if numel(tempVolume) > 1
                   epsilon_of_s(i,j,:) = tempVolume(vol_inds(sampling_inds)); 
                else
                   % Just a single value, no need to unpack
                   if ~isnan(tempVolume)
                       allEqualCells(i,j) = tempVolume;
                   end
                end
            end
        end

        %As eqn 8 (2011)
        W_of_X = zeros(size(sampling_coords,1), size(sampling_coords,1)); %nxn : diagonal matrix of weights per sampling point
        
        for i = 1:size(sampling_coords,1)
            %distance to sampling point
            d = point_dists(i);

            %As eqn 11 (2011), quartic spline weighting function
            W_of_X(i,i) = 1-6*d^2+8*d^3-3*d^4;
        end

        %As eqn 7 (2011)
        % A bit faster if done with single vectorized calc rather than per coord.
        P_of_X = calculatetricubicbasis_vec(sampling_coords(:,1), sampling_coords(:,2), sampling_coords(:,3)); %nx20 : basis vector per sampling point
        
        %As eqn 5 (2011)
        M_of_X = P_of_X'*W_of_X*P_of_X; %20x20
        
        % As eqn 6 (2011)
        N_of_X = P_of_X'*W_of_X; %20xn
        
        % Solve eqn 4 (2011)
        %a_of_X = zeros(20,3,3); %20x3x3
        a_of_X = zeros(3,3,20); %3x3x20
        
        % Only need to calculate the bottom triangle due to symmetry
%         for i = 1:3
%             for j = 1:i
%                 if isnan(allEqualCells(i,j))
%                     matA = M_of_X;
%                     
%                     matB = N_of_X*epsilon_of_s(:,i,j);
%                     
%                     % Just test matB, matA should be good...
%                     if rank(matB) >= 3 % && rank(matA) >= 3 
%                         a_of_X(:,i,j) = matA\matB; %20x1 : solves M*A=N*n for a, where N*n is 20x1
%                     else
%                         % if rank of either less than 3, use alt method
%                         a_of_X(:,i,j) = pinv(matA)*matB;
%                     end
%                 end
%             end
%         end

        for i = 1:3
            for j = 1:i
                if isnan(allEqualCells(i,j))
                    matA = M_of_X;
                    
                    matB = N_of_X*permute(epsilon_of_s(i,j,:), [3,1,2]);
                    
                    %%% Could just do if values go out of input range
                    warning('Maybe best not to use Pinv?')
                    % Just test matB, matA should be good...
                    if rank(matB) >= 3 % && rank(matA) >= 3 
                        a_of_X(i,j,:) = matA\matB; %20x1 : solves M*A=N*n for a, where N*n is 20x1
                    else
                        % if rank of either less than 3, use alt method
                        a_of_X(i,j,:) = pinv(matA)*matB;
                    end
                end
            end
        end

        % fill out top of triangle
        a_of_X(1,2,:) = a_of_X(2,1,:); a_of_X(1,3,:) = a_of_X(3,1,:);
        
        a_of_X(2,3,:) = a_of_X(3,2,:);
        
        [p_of_X, dp_partial_dX] = calculatetricubicbasis(X(1), X(2), X(3), 1); %20x1 and 20x3

%       % Eqv. eqn 3 (2011)
%       n_of_X = p_of_X'*a_of_X; %1x1 : refractive index 

        % Calculate approximated epsilon at point
        epsilon_of_X = zeros(3,3);
        
        % Just calculate bottom triangle
        for i = 1:3
            for j = 1:i
                if isnan(allEqualCells(i,j))
                    epsilon_of_X(i,j) = p_of_X*permute(a_of_X(i,j,:), [3,1,2]);
                else
                    epsilon_of_X(i,j) = allEqualCells(i,j);
                end
            end
        end
        
        % Fill top triangle before inversion
        epsilon_of_X(1,2) = epsilon_of_X(2,1); epsilon_of_X(1,3) = epsilon_of_X(3,1);
        
        epsilon_of_X(2,3) = epsilon_of_X(3,2);
        
        inv_epsilon_of_x = inv(epsilon_of_X);
        
%       % Eqv. eqn 9 (2011)
%       dn_partial_dX = dp_partial_dX'*a_of_X; %3x1 : spatial gradient of refractive index
        
        % Expand christoffel symbol value and calculate values
        gamma = zeros(3,3,3); % 3x3x3 tensor expansion
        
        % Only fill out bottom triangle for each i, as symbol jk = kj 
        % i placed on last dimension for clarity in symmetry
        
        for i = 1:3
            for j = 1:3
                for k = 1:j
                    
                    if inv_epsilon_of_x(i,1) ~= 0
                        term1= inv_epsilon_of_x(i,1)*(dp_partial_dX(k,:)*permute(a_of_X(1,j,:), [3,1,2]) + dp_partial_dX(j,:)*permute(a_of_X(1,k,:), [3,1,2]) - ...
                            dp_partial_dX(1,:)*permute(a_of_X(j,k,:), [3,1,2]));
                    else
                        term1 = 0;
                    end
                    
                    if inv_epsilon_of_x(i,2) ~= 0
                        term2 = inv_epsilon_of_x(i,2)*(dp_partial_dX(k,:)*permute(a_of_X(2,j,:), [3,1,2]) + dp_partial_dX(j,:)*permute(a_of_X(2,k,:), [3,1,2]) - ...
                            dp_partial_dX(2,:)*permute(a_of_X(j,k,:), [3,1,2]));
                    else
                        term2 = 0;
                    end
                    
                    if inv_epsilon_of_x(i,3) ~= 0
                        term3 = inv_epsilon_of_x(i,3)*(dp_partial_dX(k,:)*permute(a_of_X(3,j,:), [3,1,2]) + dp_partial_dX(j,:)*permute(a_of_X(3,k,:), [3,1,2]) - ...
                            dp_partial_dX(3,:)*permute(a_of_X(j,k,:), [3,1,2]));
                    else
                        term3 = 0;
                    end
                    
                    gamma(j,k,i) = 0.5*( term1 + term2 + term3 );
                end
            end
        end 
        
        % Fill out top triangle of symbols
        % Index of array swapped
        gamma(1,2,:) = gamma(2,1,:); gamma(1,3,:) = gamma(3,1,:);
        
        gamma(2,3,:) = gamma(3,2,:);
        
%     % Eqv. to eqn 2 (2011)
%       dT_dt = n_of_X*dn_partial_dX; %3x1 vector

        % Solve for change in dt
        dt_ds = zeros(3,1); %3x1 vector
        
        for i = 1:3
            % Expand eqn 32 (2013)
            dt_ds(i) = -(gamma(1,1,i)*T(1)^2 + gamma(1,2,i)*T(2)*T(1)+gamma(1,3,i)*T(3)*T(1)+...
                gamma(2,1,i)*T(1)*T(2) + gamma(2,2,i)*T(2)^2 + gamma(2,3,i)*T(3)*T(2)+...
                gamma(3,1,i)*T(1)*T(3) + gamma(3,2,i)*T(2)*T(3) + gamma(3,3,i)*T(3)^2);
            
        end  
end