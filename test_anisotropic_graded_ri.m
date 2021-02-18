% test algorithm from nishdate 2013 paper on axial gradient in slab

%%% Might be able to skip sort call if all are used... ?

%%% Can speed up 3x3 matrix calc due to symmetry

close all; 
clear;    
    
%create graded depth profile
depth = 10;
voxel_size = 0.0625;

vol_size = depth/voxel_size;

%%% Matrix is symetric, so can probably cut (at least) three terms
lens_volume = zeros(vol_size,vol_size,vol_size,3,3);

gradientType = 1;
% 0 is for eNegative, 1 is for axial
%%% Should pass to analytic function

if gradientType
    alpha = 0.5;
else
     a = 0.05;
     b = 0.2;
     c = 4; 
end

% Calculate refractive index
for i = 1:vol_size
    if gradientType
        lens_volume(:,:,i,3,3) = (1+alpha*sin(i*voxel_size))^2;
    else
        lens_volume(:,:,i,1,1) = a*((i*voxel_size)^2+c);
    
        lens_volume(:,:,i,2,2) = -lens_volume(:,:,i,1,1);

        lens_volume(:,:,i,3,3) = b^2*((i*voxel_size)^2+c);
    end
end

if gradientType
    lens_volume(:,:,:,1,1) = 1;

    lens_volume(:,:,:,2,2) = 1;
end

%%% Maybe be able to pack fewer given symmetry

%Seperate into cells for faster processing.
lens_volume_seperated = cell(3,3);

for i = 1:3
    for j = 1:3
        if all(lens_volume(:,:,:,i,j) == lens_volume(1,1,1,i,j))
            %If all are equal just store one value 
            lens_volume_seperated{i,j} = lens_volume(1,1,1,i,j);
        else
            lens_volume_seperated{i,j} = lens_volume(:,:,:,i,j);
        end
    end
end

[temp_x, temp_y, temp_z] = meshgrid(1:vol_size, 1:vol_size, 1:vol_size);
vol_inds = sub2ind([vol_size, vol_size, vol_size], temp_x(:), temp_y(:), temp_z(:));
vol_coords = ([temp_x(:), temp_y(:), temp_z(:)])*voxel_size;

z_steps = (1:vol_size)*voxel_size;

%set up initial rays.
% step_num = 20;
% temp_p = (0.5*depth/voxel_size+step_num):step_num:(2.5*depth/voxel_size-step_num);
% [temp_x, temp_y] = meshgrid(temp_p, temp_p);
ray_origins = [mean(vol_coords(:,1)), mean(vol_coords(:,2)), voxel_size]'/voxel_size; %3xn 

%trace rays by stepping across
figure; 
subplot(2,2,1);
imshow((permute(lens_volume(:,vol_size/2,:, 3, 3), [1 3 2])-1)/(max(lens_volume(:))-1));

%trace through
interp_type = '45';

for i = 1:size(ray_origins,2)
    %plot3(ray_origins(1,i), ray_origins(2,i), ray_origins(3,i),'.b');
    
    ray_t = [-0.3, 0.3, 1]'; %3x1 : ray will move from bottom to top
    
    %ray_t = [0.504, 0.3, 1]'; %3x1 : ray will move from bottom to top
    
    ray_x = ray_origins(:,i);
    
    xDirCos = ray_t(1);%/norm(ray_t);
    
    yDirCos = ray_t(2);%/norm(ray_t);
    
    % For slab, starts directly in gradient
    in_graded = 1;
    
    x = (ray_x)*voxel_size;
    t = ray_t;
    
    go = 1;

    last_nearby_x = x;

    % For unit matrix, radius 4 will take ~268 points
    nearby_inds = find(vol_coords(:,1) > x(1)-4*voxel_size & vol_coords(:,1) < x(1)+4*voxel_size & ...
        vol_coords(:,2) > x(2)-4*voxel_size & vol_coords(:,2) < x(2)+4*voxel_size & ...
        vol_coords(:,3) > x(3)-4*voxel_size & vol_coords(:,3) < x(3)+4*voxel_size);
    
    delta_s = 10^-4;
    
    ray_path = zeros(vol_size-1, 3)*NaN;
    
    current_step = 1;
    
    ray_path(current_step,:) = x;
    
    its = 1;
    
    while go
        if ~in_graded 
            ray_x = ray_x + initial_ray_t; 
        else

            x0 = x;
            t0 = t;

            [x, t, delta_s] = ray_interpolation(interp_type, x0, t0, delta_s, vol_coords(nearby_inds,:), ...
                vol_inds(nearby_inds), lens_volume_seperated, [1 1 1]*vol_size);
            
%             [x, t, delta_s] = ray_interpolation(interp_type, x0, t0, delta_s, [], ...
%                 [], [], [1 1 1]*vol_size);
            
            its = its + 1;
            
            %check if x has moved to include new voxels
            if norm(x-last_nearby_x) > voxel_size
                %if so, reselect nearby point
                
                last_nearby_x = x;
                
                nearby_inds = find(vol_coords(:,1) > x(1)-4*voxel_size & vol_coords(:,1) < x(1)+4*voxel_size & ...
                    vol_coords(:,2) > x(2)-4*voxel_size & vol_coords(:,2) < x(2)+4*voxel_size & ...
                    vol_coords(:,3) > x(3)-4*voxel_size & vol_coords(:,3) < x(3)+4*voxel_size);
                
                [x' t' delta_s*10^3]
                
                % record path for plotting later
                if x(3) >= z_steps(current_step+1)
                    current_step = current_step + 1;
                    
                    ray_path(current_step,:) = x;
                end
            end
            ray_x = x/voxel_size;
        end 
        if any(ray_x > vol_size) | any(ray_x < 1)
           go = 0; 
        end
    end
    
    ray_path(isnan(ray_path(:,1)),:) = [];
    
    error('Check these solutions, am a bit suspecious after GRIN paper')
    
    % Take analytic solution
    if gradientType
        % axial
        xPath = ray_path(1,1) + xDirCos./(1 + alpha*sin(z_steps(1))).* ...
            (z_steps - z_steps(1) - alpha*(cos(z_steps) - cos(z_steps(1))));

        yPath = ray_path(1,2) + yDirCos./(1 + alpha*sin(z_steps(1))).* ...
            (z_steps - z_steps(1) - alpha*(cos(z_steps) - cos(z_steps(1))));
    else
        % For e-negative
        xPath = ray_path(1,1) + b*xDirCos*(a*(z_steps(1)^2+c))/(a*sqrt(b^2*(z_steps(1)^2+c))).*...
            log((z_steps + sqrt(z_steps.^2+c))/(z_steps(1) + sqrt(z_steps(1).^2+c)));

        yPath = ray_path(1,2) + b*yDirCos*(a*(z_steps(1)^2+c))/(a*sqrt(b^2*(z_steps(1)^2+c))).*...
            log((z_steps + sqrt(z_steps.^2+c))/(z_steps(1) + sqrt(z_steps(1).^2+c)));
    end
     
    subplot(2,2,2);
    plot3(ray_path(:,1), ray_path(:,2), ray_path(:,3),'b'); hold on
    plot3(xPath, yPath, z_steps,'g');
    axis equal
    
    subplot(2,2,3); 
    plot(ray_path(:,3), ray_path(:,1), 'bx'); hold on
    plot(z_steps, xPath, 'g')
    axis equal
    
    subplot(2,2,4); 
    plot(ray_path(:,3), ray_path(:,2), 'bx'); hold on
    plot(z_steps, yPath, 'g')
    axis equal
    
    pause(0.01)
end

function [X, T, new_delta] = ray_interpolation(type, X0, T0, delta, vol_coords, vol_inds, lens_volume, vol_size)
    %%% Make this external function and pass function handle to evaluate.

    % RKF45 from Nishidate 2011 paper with variable step size 
    % using coeffs from Table 1

    error('Adapted treatment from *Upgrading ... for second order* paper, seems wrong')
    %%% Switch to a Nystrom method or system of first order equations with RKF4(5)
    % not Bettis 1973 - as abbreviated form of Nystom without T dependence
    % is assumed
    
    if 1
        k1 = numerical_dt_ds(X0, ...
            T0,...
            vol_coords, vol_inds, lens_volume, vol_size); 

        % New where both X and T are modified and delta shifted down

        if type ~= '1'
            k2 = numerical_dt_ds(X0 + delta*(T0/4), ...
                T0 +delta*(k1/4),...
                vol_coords, vol_inds, lens_volume, vol_size); 

            % Expanded eqns for X on terms below - all but 2 equal or less than 5e-16 off
            k3 = numerical_dt_ds(X0 + delta*T0*(3/32+9/32) + delta^2*(k1*(9/32*1/4)),  ...      
                T0 +delta*(k1*3/32 + k2*9/32),...
                vol_coords, vol_inds, lens_volume, vol_size); 

            k4 = numerical_dt_ds(X0+delta*T0*(1932/2197-7200/2197+7296/2197)+ ...
                delta^2*(k1*(-7200/2197*1/4 + 7296/2197*3/32) + k2*(7296/2197*9/32)),  ...
                T0 +delta*(k1*1932/2197 +k2*-7200/2197 +k3*7296/2197),...
                vol_coords, vol_inds, lens_volume, vol_size); 

            %%% Difference in K1 terms for X
            k5 = numerical_dt_ds(X0+delta*T0*(439/216 - 8 +3680/513 - 845/4104) + ...
                delta^2*(k1*(-8/4 + 3680/513*3/32 - 845/4104*1932/2197) +k2*(3680/513*9/32 -845/4104*-7200/2197) +k3*(-845/4104*7296/2197)), ...
                T0 +delta*(k1*439/216 +k2*-8 +k3*3680/513 +k4*-845/4104),...
                vol_coords, vol_inds, lens_volume, vol_size); 

            %%% Difference in K4 term for X
            k6 = numerical_dt_ds(X0+delta*T0*(-8/27 + 2 - 3544/2565 + 1859/4104 - 11/40)        + ...
                delta^2*(k1*(2*1/4 - 3544/2565*3/32 +1859/4104*1932/2197 -11/40*439/216) +k2*(-3544/2565*9/32 + 1859/4104*-7200/2197 - 11/40*-8) +k3*(1859/4104*7296/2197  - 11/40*3680/513) +k4*(-11/40*-845/4104)), ...
                T0 +delta*(k1*-8/27 +k2*2 +k3*-3544/2565 +k4*1859/4104 +k5*-11/40),...
                vol_coords, vol_inds, lens_volume, vol_size); 
        end
            
      %Old coeffs where only X is modified and delta kept up - kept to test
%     k2 = delta*numerical_dt_ds(X0+delta*T0/4        +delta*(k1/4), T0,...
%         vol_coords, vol_inds, lens_volume, vol_size); 
%     
%     k3 = delta*numerical_dt_ds(X0+delta*T0*3/8      +delta*(k1*3/32 +k2*9/32), T0,... 
%         vol_coords, vol_inds, lens_volume, vol_size); 
%     
%     k4 = delta*numerical_dt_ds(X0+delta*T0*12/13    +delta*(k1*1932/2197 +k2*-7200/2197 +k3*7296/2197), T0,...
%         vol_coords, vol_inds, lens_volume, vol_size); 
%     
%     k5 = delta*numerical_dt_ds(X0+delta*T0          +delta*(k1*439/216 +k2*-8 +k3*3680/513 +k4*-845/4104), T0,...
%         vol_coords, vol_inds, lens_volume, vol_size); 
%     
%     k6 = delta*numerical_dt_ds(X0+delta*T0/2        +delta*(k1*-8/27 +k2*2 +k3*-3544/2565 +k4*1859/4104 +k5*-11/40), T0,...
%         vol_coords, vol_inds, lens_volume, vol_size); 

    elseif 0
        k1 = analytic_dt_ds(X0, T0); 

        if type ~= '1'
            k2 = analytic_dt_ds(X0 + delta*(T0/4), T0 +delta*(k1/4)); 

            k3 = analytic_dt_ds(X0 + delta*T0*(3/32+9/32) + delta^2*(k1*(9/32*1/4)),  ...      
                T0 +delta*(k1*3/32 + k2*9/32)); 

            k4 = analytic_dt_ds(X0+delta*T0*(1932/2197-7200/2197+7296/2197)+ ...
                delta^2*(k1*(-7200/2197*1/4 + 7296/2197*3/32) + k2*(7296/2197*9/32)),  ...
                T0 +delta*(k1*1932/2197 +k2*-7200/2197 +k3*7296/2197)); 

            k5 = analytic_dt_ds(X0+delta*T0*(439/216 - 8 +3680/513 - 845/4104) + ...
                delta^2*(k1*(-8/4 + 3680/513*3/32 - 845/4104*1932/2197) +k2*(3680/513*9/32 -845/4104*-7200/2197) +k3*(-845/4104*7296/2197)), ...
                T0 +delta*(k1*439/216 +k2*-8 +k3*3680/513 +k4*-845/4104)); 

            k6 = analytic_dt_ds(X0+delta*T0*(-8/27 + 2 - 3544/2565 + 1859/4104 - 11/40)        + ...
                delta^2*(k1*(2*1/4 - 3544/2565*3/32 +1859/4104*1932/2197 -11/40*439/216) +k2*(-3544/2565*9/32 + 1859/4104*-7200/2197 - 11/40*-8) +k3*(1859/4104*7296/2197  - 11/40*3680/513) +k4*(-11/40*-845/4104)), ...
                T0 +delta*(k1*-8/27 +k2*2 +k3*-3544/2565 +k4*1859/4104 +k5*-11/40)); 
        end
        
        %Old coeffs where only X is modified and delta kept up - kept to test
%         k1 = delta*analytic_dt_ds(X0, T0); 
%         
%         k2 = delta*analytic_dt_ds(X0+delta*T0/4        +delta*(k1/4), T0); 
% 
%         k3 = delta*analytic_dt_ds(X0+delta*T0*3/8      +delta*(k1*3/32 +k2*9/32), T0); 
% 
%         k4 = delta*analytic_dt_ds(X0+delta*T0*12/13    +delta*(k1*1932/2197 +k2*-7200/2197 +k3*7296/2197), T0); 
% 
%         k5 = delta*analytic_dt_ds(X0+delta*T0          +delta*(k1*439/216 +k2*-8 +k3*3680/513 +k4*-845/4104), T0); 
% 
%         k6 = delta*analytic_dt_ds(X0+delta*T0/2        +delta*(k1*-8/27 +k2*2 +k3*-3544/2565 +k4*1859/4104 +k5*-11/40), T0); 
    end
    
    if type == '1'
        % Euler solver
        X = X0 + delta*T0;
        T = T0 + delta*k1;
        
        new_delta = delta;
        
    elseif type == '4'
        %%% Not tested (do against RK4 from RKF45), may need to extend to to modify T0
        
        % RK4 standard - works as Sharma
        % general description from Meggit and wiki is a confusing
        % explicit method on wiki and second-order eqn w/ Simpson's rule from Scaraborough - Numerical Mathmatical Analysis p 382 (ref- from Sharma) are helpful
            % can also derive w/ first order simultaneous eqns?
            % gives very similar results to first order eqns (0 1/2 1/2 1) instead of (0 1/8 1/8 1/2)
        k1 = analytic_dt_ds(X0, T0); 
        k2 = analytic_dt_ds(X0+delta*T0/2    +delta^2*k1/8, T0);
        k3 = analytic_dt_ds(X0+delta*T0/2    +delta^2*k2/8, T0);
        k4 = analytic_dt_ds(X0+delta*T0      +delta^2*k3/2, T0);
        X = X0 + delta*T0 + delta^2*(k1 + k2 + k3)/6; %b coeffs should add to 0.5
        T = T0 + delta*(k1+2*(k2+k3)+k4)/6; %b coeffs should add to 1
        
        new_delta = delta;
        
    elseif type == '45'      
        % Shifted delta down as in 5th order code
        T_4th = T0 + delta*(k1*25/216 +k3*1408/2565 +k4*2197/4104 + k5*-1/5);
        T_5th = T0 + delta*(k1*16/135 +k3*6656/12825 +k4*28561/56430 +k5*-9/50 +k6*2/55); %these b coeffs should sum to 1
        
        % from Upgrading Runge-Kutta-Fehlberg Method (RKFM) for Second Order Ordinary Differential Equations 
            % (corrected error in 5th order coeffs, k3c read 2856 instead of 28561)
            % in same manner as Scaraborough
        % 4th order coeffs
        k1C = 1408/2565*3/32 + 2197/4104*1932/2197 + -1/5*439/216;
        k2C = 1408/2565*9/32 + 2197/4104*-7200/2197 + -1/5*-8; %close to zero
        k3C = 2197/4104*7296/2197 + -1/5*3680/513;
        k4C = -1/5*-845/4104;
        X_4th = X0 + delta*T0 + delta^2*(k1*k1C + k2*k2C + k3*k3C + k4*k4C);
        
        %5th order coeffs
        k1C = 6656/12825*3/32 + 28561/56430*1932/2197 + -9/50*439/216 + 2/55*-8/27; 
        %k2C = 6656/12825*9/32 + 28561/56430*-7200/2197 + -9/50*-8 + 2/55*2; %close to zero
        k3C = 28561/56430*7296/2197 + -9/50*3680/513 + 2/55*-3544/2565;
        k4C = -9/50*-845/4104 + 2/55*1859/4104;
        k5C = 2/55*-11/40; 
        X_5th = X0 + delta*T0 + delta^2*(k1*k1C + k3*k3C + k4*k4C + k5*k5C); %these be coeffs should sum to 0.5
        
        % Get error - normalize position terms to meters
        % As eqn 16 (2011)
        er = max(abs([((X_5th')-(X_4th'))*10^3 T_5th'-T_4th']))/10^-12;
            
        %limit increase, if er goes to zero, will scale to inf.
         if er < 0.01
           er = 0.01; 
         end

        %lower bound (0.6 or 0.8 ???) should be scaled reltive to saftey factor (0.9)
            %so step does not decrease if it should increase
        if er < 0.8 | er > 1
           % As eqn 17 but divide by delta
           new_delta = 0.9*delta/er^0.2;
        else
           new_delta = delta; 
        end
        
        if isinf(new_delta) | isnan(er)
           error('Problem!') 
        end
        
        T = T_5th;
        X = X_5th;
        
    elseif type == '5'
        % Old version is when delta was included in k terms
        %T = T0 + k1*16/135 +k3*6656/12825 +k4*28561/56430 +k5*-9/50 +k6*2/55; %these b coeffs should sum to 1
        T = T0 + delta*(k1*16/135 +k3*6656/12825 +k4*28561/56430 +k5*-9/50 +k6*2/55); %these b coeffs should sum to 1
        
        %5th order coeffs - from paper above
        %TOC = (16/135 + 6656/12825 + 28561/56430 - 9/50 + 2/55); % Equals 1
        k1C = 6656/12825*3/32 + 28561/56430*1932/2197 + -9/50*439/216 + 2/55*-8/27;
        %k2C = 6656/12825*9/32 + 28561/56430*-7200/2197 + -9/50*-8 + 2/55*2; %close to zero
        k3C = 28561/56430*7296/2197 + -9/50*3680/513 + 2/55*-3544/2565;
        k4C = -9/50*-845/4104 + 2/55*1859/4104;
        k5C = 2/55*-11/40; 
        %X = X0 + delta*(T0+(k1*k1C + k3*k3C + k4*k4C + k5*k5C)); %these be coeffs should sum to 0.5
        X = X0 + delta*T0 + delta^2*(k1*k1C + k3*k3C + k4*k4C + k5*k5C);
        
        new_delta = delta;
    else
        error('Type not coded')
    end
end    

function [dt_ds] = analytic_dt_ds(X, T)

    % Should pass coefficients in for saftey...
    %%% Also pass in gradientType
    
    alpha = 0.5;

    a = 0.05;
    b = 0.2;
    c = 4;
    
    % Calculate epsilon and inverse
    epsilon_of_X = zeros(3,3);

%     epsilon_of_X(1,1) = 1;
%     
%     epsilon_of_X(2,2) = 1;
%     
%     epsilon_of_X(3,3) = (1+alpha*sin(X(3)))^2;
    
    % e-negative material
    epsilon_of_X(1,1) = a*(X(3)^2+c);
    
    epsilon_of_X(2,2) = -epsilon_of_X(1,1);
    
    epsilon_of_X(3,3) = b^2*(X(3)^2+c);
     
    inv_epsilon_of_x = inv(epsilon_of_X);
    
    % Eqv. eqn 9 (2011)
    
    % Expand christoffel symbol value and calculate values
    gamma = zeros(3,3,3); % 3x3x3 tensor expansion

    %%% Could fill out using symmetry (slightly different from epsilon)

    % Only has value when all are 3
    %gamma(3,3,3) = inv_epsilon_of_x(3,3)*(2*alpha*cos(X(3))*(alpha*sin(X(3))+1));
    
    % e-negative material
    gamma(1,1,3) = inv_epsilon_of_x(1,1)*2*a*X(3); gamma(1,3,1) = gamma(1,1,3);  
    
    gamma(2,2,3) = inv_epsilon_of_x(2,2)*-2*a*X(3); gamma(2,3,2) = gamma(2,2,3); 
    
    gamma(3,1,1) = inv_epsilon_of_x(3,3)*-2*a*X(3);
    
    gamma(3,2,2) = inv_epsilon_of_x(3,3)*2*a*X(3);
    
    gamma(3,3,3) = inv_epsilon_of_x(3,3)*2*b^2*X(3); 

    gamma = 0.5*gamma;
    
    % Eqv. to eqn 2 (2011)

    % Solve for change in dt
    dt_ds = zeros(3,1); %3x1 vector

    for i = 1:3
        % Expand eqn 32 (2013)
        dt_ds(i) = -(gamma(i,1,1)*T(1)^2 + gamma(i,2,1)*T(2)*T(1)+gamma(i,3,1)*T(3)*T(1)+...
            gamma(i,1,2)*T(1)*T(2) + gamma(i,2,2)*T(2)^2 + gamma(i,3,2)*T(3)*T(2)+...
            gamma(i,1,3)*T(1)*T(3) + gamma(i,2,3)*T(2)*T(3) + gamma(i,3,3)*T(3)^2);
        
        %dt_ds(i) = -(gamma(i,1,1)*T(1)^2 + gamma(i,2,2)*T(2)^2 + gamma(i,3,3)*T(3)^2);

    end  
end

function [dt_ds] = numerical_dt_ds(X, T, vol_coords, vol_inds, lens_volume_seperated, vol_size)
        %get change in refractive index based on christoffel symbols
            %merges method from Nishdate 2011 and 2013 papers
        
        %%% Indexes on off diagonals may be a mess here...
            
        %take sampling points first
        [sampling_inds, sampling_coords, point_dists] = getsamplingpoints(vol_coords, X);
        
        % Get epsilon values from cell array
        epsilon_of_s = zeros(length(sampling_inds),3,3);
        
        allEqualCells = zeros(3,3)*NaN;
        
        %%% Maybe able to unpack fewer given symmetry
        
        for i = 1:3
            for j =1:3
                tempVolume = lens_volume_seperated{i,j};
                if numel(tempVolume) > 1
                   epsilon_of_s(:,i,j) = tempVolume(vol_inds(sampling_inds));
                else
                   % Just a single value, no need to unpack
                   allEqualCells(i,j) = tempVolume;
                end
            end
        end

        % Hopefully faster to do using cell array 
%         [tempX, tempY, tempZ] = ind2sub(vol_size, vol_inds(sampling_inds));
%
%         for i = 1:length(sampling_inds)
%             epsilon_of_s(i,:,:) = lens_volume(tempX(i), tempY(i), tempZ(i), :, :); %nx3x3 matrix of refractive indices at sampling points
%         end
        
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

%         P_of_X = zeros(size(sampling_coords,1),20); %nx20 : basis vector per sampling point
%         
%         for i = 1:size(sampling_coords,1)
%             %place basis function for each sampling point
% 
%             P_of_X(i,:) = calculatetricubicbasis(sampling_coords(i,1), sampling_coords(i,2), sampling_coords(i,3), 0)';
%         end
        
        %As eqn 5 (2011)
        M_of_X = P_of_X'*W_of_X*P_of_X; %20x20
        
        % As eqn 6 (2011)
        N_of_X = P_of_X'*W_of_X; %20xn
        
        % Solve eqn 4 (2011) (extend to 3x3 matrix)
        a_of_X = zeros(20,3,3); %20x3x3
        
        %%% Only need to calcualte to the diagonal given symmetry
        
        for i = 1:3
            for j = 1:3
                if isnan(allEqualCells(i,j))
                    matA = M_of_X;
                    
                    matB = N_of_X*epsilon_of_s(:,i,j);
                    
                    warning('Maybe best not to use Pinv?')
                    % Just test matB, matA should be good...
                    if rank(matB) >= 3 % && rank(matA) >= 3 
                        a_of_X(:,i,j) = matA\matB; %20x1 : solves M*A=N*n for a, where N*n is 20x1
                    else
                        % if rank of either less than 3, use alt method
                        a_of_X(:,i,j) = pinv(matA)*matB;
                    end
                end
            end
        end
                
        [p_of_X, dp_partial_dX] = calculatetricubicbasis(X(1), X(2), X(3), 1); %20x1 and 20x3

%       % Eqv. eqn 3 (2011)
%       n_of_X = p_of_X'*a_of_X; %1x1 : refractive index 

        % Calculate approximated epsilon at point
        epsilon_of_X = zeros(3,3);
        
        %%% Can fill out using symmetry
        
        for i = 1:3
            for j = 1:3
                if isnan(allEqualCells(i,j))
                    epsilon_of_X(i,j) = p_of_X*a_of_X(:,i,j);
                else
                    epsilon_of_X(i,j) = allEqualCells(i,j);
                end
            end
        end
        
        inv_epsilon_of_x = inv(epsilon_of_X);
        
%       % Eqv. eqn 9 (2011)
%       dn_partial_dX = dp_partial_dX'*a_of_X; %3x1 : spatial gradient of refractive index
        
        % Expand christoffel symbol value and calculate values
        gamma = zeros(3,3,3); % 3x3x3 tensor expansion
        
        %%% Could fill out using symmetry (slightly different from epsilon )
        
        for i = 1:3
            for j = 1:3
                for k = 1:3
                    
                    if inv_epsilon_of_x(i,1) ~= 0
                        term1= inv_epsilon_of_x(i,1)*(dp_partial_dX(k,:)*a_of_X(:,1,j) + dp_partial_dX(j,:)*a_of_X(:,1,k) - dp_partial_dX(1,:)*a_of_X(:,j,k));
                    else
                        term1 = 0;
                    end
                    
                    if inv_epsilon_of_x(i,2) ~= 0
                        term2 = inv_epsilon_of_x(i,2)*(dp_partial_dX(k,:)*a_of_X(:,2,j) + dp_partial_dX(j,:)*a_of_X(:,2,k) - dp_partial_dX(2,:)*a_of_X(:,j,k));
                    else
                        term2 = 0;
                    end
                    
                    if inv_epsilon_of_x(i,3) ~= 0
                        term3 = inv_epsilon_of_x(i,3)*(dp_partial_dX(k,:)*a_of_X(:,3,j) + dp_partial_dX(j,:)*a_of_X(:,3,k) - dp_partial_dX(3,:)*a_of_X(:,j,k));
                    else
                        term3 = 0;
                    end
                    
                    gamma(i,j,k) = 0.5*( term1 + term2 + term3 );
                end
            end
        end 
        
%      % Eqv. to eqn 2 (2011)
%       dT_dt = n_of_X*dn_partial_dX; %3x1 vector

        % Solve for change in dt
        dt_ds = zeros(3,1); %3x1 vector
        
        for i = 1:3
            % Expand eqn 32 (2013)
            dt_ds(i) = -(gamma(i,1,1)*T(1)^2 + gamma(i,2,1)*T(2)*T(1)+gamma(i,3,1)*T(3)*T(1)+...
                gamma(i,1,2)*T(1)*T(2) + gamma(i,2,2)*T(2)^2 + gamma(i,3,2)*T(3)*T(2)+...
                gamma(i,1,3)*T(1)*T(3) + gamma(i,2,3)*T(2)*T(3) + gamma(i,3,3)*T(3)^2);
            
        end  
end