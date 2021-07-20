% Test algorithm from nishdate 2011 on luneberg lens

%%%to do:

%step size tends to decay rapidly...
    %still not exactly sure about step size  
    
%RK4 and RKF45 method seem to give similar accuracy    
    
%look for easy speed ups
    %when runing final RK5 method gradient does not need to be recomputed (?), huge save

close all; 
clear;    
    
radius = 1;
step_size = 0.0625; 

%create luneburg lens
vol_size = 3*radius/step_size;
lens_volume = ones(vol_size,vol_size,vol_size);

for i = 1:vol_size
    for j = 1:vol_size
        for k = 1:vol_size
            x = (i-vol_size/2)*step_size;
            y = (j-vol_size/2)*step_size;
            z = (k-vol_size/2)*step_size;
            r = sqrt(x^2+y^2+z^2);
            
            if r < radius
                lens_volume(i,j,k) = sqrt(2-(r/radius)^2);
            end
        end
    end
end

% create fiber
%%% note fully implement
% fiber_length = 10;
% 
% vol_size = 3*radius/step_size;
% lens_volume = ones(ceil([vol_size vol_size fiber_length*1.5/step_size]));
% 
% for i = 1:vol_size
%     for j = 1:vol_size
%        
%         x = (i-vol_size/2)*step_size;
%         y = (j-vol_size/2)*step_size;
%         r = sqrt(x^2+y^2);
% 
%         if r < radius
%             lens_volume(i,j,:) = sqrt(2.5-(r/radius)^2);
%         end
%     end
% end

[temp_x, temp_y, temp_z] = meshgrid(1:vol_size, 1:vol_size, 1:vol_size);
vol_inds = sub2ind([vol_size, vol_size, vol_size], temp_x(:), temp_y(:), temp_z(:));
vol_coords = ([temp_x(:), temp_y(:), temp_z(:)]-vol_size/2)*step_size;
%just keep points in lens
temp_inds = find(sqrt(vol_coords(:,1).^2 + vol_coords(:,2).^2 + vol_coords(:,3).^2) < radius);
vol_inds = vol_inds(temp_inds);
vol_coords = vol_coords(temp_inds,:);

%set up initial rays.
temp_p = (0.5*radius/step_size+5):5:(2.5*radius/step_size-5);
% [temp_y, temp_z] = meshgrid(temp_p, temp_p);

[temp_y, temp_z] = meshgrid(temp_p, mean(temp_p));

ray_origins = [ones(numel(temp_y),1), temp_y(:), temp_z(:)]'; %3xn 

%trace rays by stepping across
figure; 
subplot(1,3,1);
imshow((lens_volume(:,:,vol_size/2)-1)/(max(lens_volume(:))-1));

subplot(1,3,2);
hold on; axis equal

subplot(1,3,3);
hold on; axis equal

%trace through
interp_type = 1; 3;

for i = 1:size(ray_origins,2)
    %plot3(ray_origins(1,i), ray_origins(2,i), ray_origins(3,i),'.b');
    
    initial_ray_T = [1, 0, 0]'; %3x1 : ray will move from left to right
    ray_X = ray_origins(:,i);
    
    in_graded = 0;
    go = 1;
    delta_t = 0.001;
    
    while go
        if ~in_graded; ray_X = ray_X + initial_ray_T; end
        
        if sqrt((ray_X(1)-vol_size/2)^2+(ray_X(2)-vol_size/2)^2+(ray_X(3)-vol_size/2)^2)*step_size < radius
            
            if ~in_graded
                %%% Strictly speaking should test continuously as no interface to sphere
                
                %entry to graded
                in_graded = 1;
                ray_x_back = ray_X-initial_ray_T; %step back one
                
                %find intersection of sphere geometrically
                    %from https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection
                O = (ray_x_back - vol_size/2)*step_size;
                C = 0 ;
                L = C - O;
                t_ca = dot(L,initial_ray_T);
                d = sqrt(dot(L,L)-t_ca^2);
                t_hc = sqrt(radius^2-d^2);
                t_0 = t_ca - t_hc;
                
                X = O+initial_ray_T*t_0; %3x1 coordinates
                T = initial_ray_T; %3x1 vector
                
                % could calculate refraction as usual, but does not occur as grading starts from 1
                %surface_normal = ray_X - vol_size/2; surface_normal = surface_normal/norm(surface_normal);
            end

            X0 = X;
            T0 = T;
            
%             [X, T, delta_t] = ray_interpolation(interp_type, X0, T0, delta_t, vol_coords, vol_inds, lens_volume);

            [X, T, delta_t] = ray_interpolation('4S', 'iso', X0, T0, delta_t, vol_coords, vol_inds, lens_volume, 10^-12);
            
            if sqrt(X(1)^2+X(2)^2+X(3)^2) > radius
                %ray has exited lens on this step, need to iterate close to surface as in nishidate 2011 paper
                delta_t0 = delta_t;
                rad_current = sqrt(X(1)^2+X(2)^2+X(3)^2); %always larger outbound
                rad_back = sqrt(X0(1)^2+X0(2)^2+X0(3)^2);
                
                %%% Check in trace and test file
                error('Lambda should be rad - x0, not x1 - rad')
                
                lambda0 = (rad_current-radius)/(rad_current-rad_back);
                
                epsilon = sqrt(1.49e-8); %error term, suggest using sqrt of machine precision
                %for matlab sqrt(2.2204e-16), for paper sqrt(1.49e-8)
                
                its = 1;
                r_init = (1-sqrt(X(1)^2+X(2)^2+X(3)^2))*10^6;
                         
                %run loop until ray is barely moving forward
                % As algorithm 1 from paper
                while delta_t0*norm(T0) > epsilon && abs(lambda0-1) > 10^-3 % & lambda0 > 10^-3
                    SF = 1; %A from paper - saftey factor
                    
                    if lambda0 < 10^-3
                       T = T0; %%% included in loop on paper.
                       warning('Not sure on this');
                       break
                    end
                    
                    %this loop steps back from overshoot
                    while sqrt(X(1)^2+X(2)^2+X(3)^2) > radius
                        T = T0; %idea is that we keep using start vector, even though this isn't used...
                        
                        if (SF >= 0.1)
                           SF = SF - 0.1; 
                        else
                           SF = 0.1*SF;
                        end
                        
                        delta_t = SF*lambda0*delta_t0;
                        
                        %%% step ray forward with RK5 with fixed step
%                         [X, T] = ray_interpolation(4, X0, T0, delta_t, vol_coords, vol_inds, lens_volume);
                          [X, T] = ray_interpolation('5RKN', 'iso', X0, T0, delta_t, vol_coords, vol_inds, lens_volume, 10^-12);
                    end
                    
                    delta_t = (1-SF)*lambda0*delta_t0;
                    
                    %this loop steps forward to find new contact position from position just before overshoot.
                    while sqrt(X(1)^2+X(2)^2+X(3)^2) < radius
                        X0 = X; %think this is the case, we are stepping forward
                        T0 = T;
                        
%                         [X, T] = ray_interpolation(4, X0, T0, delta_t, vol_coords, vol_inds, lens_volume);
                        [X, T] = ray_interpolation('5RKN', 'iso', X0, T0, delta_t, vol_coords, vol_inds, lens_volume, 10^-12);
                    end
                    
                    %update parameters for next loop
                    rad_current = sqrt(X(1)^2+X(2)^2+X(3)^2); %always larger outbound
                    rad_back = sqrt(X0(1)^2+X0(2)^2+X0(3)^2);
                    lambda0 = (rad_current-radius)/(rad_current-rad_back);
                
                    delta_t0 = delta_t;
                    
                    its = its + 1;
                end
                
                in_graded = 0;
                initial_ray_T = T;
                
                subplot(1,3,3);
                plot(X(2)*10^6, X(3)*10^6,'xr');
                
                r_end = (1-sqrt(X(1)^2+X(2)^2+X(3)^2))*10^6;
                [its-1 r_init r_end]
            end
            
            
            ray_X = X/step_size + vol_size/2;
        end
        
        subplot(1,3,2);
        if in_graded
            plot3(ray_X(1), ray_X(2), ray_X(3),'.b');
        else
            plot3(ray_X(1), ray_X(2), ray_X(3),'.c');
        end
        
        if any(ray_X > vol_size)
           go = 0; 
        end
    end
    pause(0.01)
end
subplot(1,3,1);
title('Luneberg lens r = 1mm n = 1-1.41')

subplot(1,3,2);
% TR = SubdivideSphericalMesh(IcosahedronMesh,2);
% sphereCords = TR.X*radius/step_size+vol_size/2;
% plot3(sphereCords(:,1), sphereCords(:,2), sphereCords(:,3), '.r')
axis off

subplot(1,3,3);
ylim([-50 50]); xlim([-50 50]);
title('spot diagram (nm)')

% function [X, T, new_delta_t] = ray_interpolation(type, X0, T0, delta_t, vol_coords, vol_inds, lens_volume)
% 
%     if type == 1
%         % RK from Sharma - uses some abbreviations... is be forth order?
%         A = delta_t*numerical_dT_dt(X0, vol_coords, vol_inds, lens_volume); %3x1
%         B = delta_t*numerical_dT_dt(X0+delta_t/2*T0+delta_t/8*A, vol_coords, vol_inds, lens_volume); %3x1
%         C = delta_t*numerical_dT_dt(X0+delta_t*T0+delta_t/2*B, vol_coords, vol_inds, lens_volume); %3x1
% 
%         X = X0 + delta_t*(T0+(A+2*B)/6);
%         T = T0 + (A+4*B+C)/6;
%         
%         new_delta_t = delta_t;
%         
%     elseif type == 2
%         % RK4 standard - works as Sharma (but not exactly the same???)
%         % general description from Meggit and wiki is a confusing
%         % explicit method on wiki and second-order eqn w/ Simpson's rule from Scaraborough - Numerical Mathmatical Analysis p 382 (ref- from Sharma) are helpful
%             % can also derive w/ first order simultaneous eqns?
%             % gives very similar results to first order eqns (0 1/2 1/2 1) instead of (0 1/8 1/8 1/2)
%         k1 = delta_t*numerical_dT_dt(X0, vol_coords, vol_inds, lens_volume); 
%         k2 = delta_t*numerical_dT_dt(X0+delta_t*T0/2    +delta_t*k1/8, vol_coords, vol_inds, lens_volume);
%         k3 = delta_t*numerical_dT_dt(X0+delta_t*T0/2    +delta_t*k2/8, vol_coords, vol_inds, lens_volume);
%         k4 = delta_t*numerical_dT_dt(X0+delta_t*T0      +delta_t*k3/2, vol_coords, vol_inds, lens_volume);
%         X = X0 + delta_t*(T0+(k1 + k2 + k3)/6); %b coeffs should add to 0.5
%         T = T0 + (k1+2*(k2+k3)+k4)/6; %b coeffs should add to 1
%         
%         new_delta_t = delta_t;
%         
%     elseif type == 3
%         % RKF45 from Nishidate 2011 paper with variable step size 
%         % using coeffs from Table 1
%         
%         %%% Had note that I should change to anisotropic method, but I think that was actually wrong
%         
%         k1 = delta_t*numerical_dT_dt(X0, vol_coords, vol_inds, lens_volume); 
%         k2 = delta_t*numerical_dT_dt(X0+delta_t*T0/4        +delta_t*(k1/4), vol_coords, vol_inds, lens_volume); 
%         k3 = delta_t*numerical_dT_dt(X0+delta_t*T0*3/8      +delta_t*(k1*3/32 +k2*9/32), vol_coords, vol_inds, lens_volume); 
%         k4 = delta_t*numerical_dT_dt(X0+delta_t*T0*12/13    +delta_t*(k1*1932/2197 +k2*-7200/2197 +k3*7296/2197), vol_coords, vol_inds, lens_volume); 
%         k5 = delta_t*numerical_dT_dt(X0+delta_t*T0          +delta_t*(k1*439/216 +k2*-8 +k3*3680/513 +k4*-845/4104), vol_coords, vol_inds, lens_volume); 
%         k6 = delta_t*numerical_dT_dt(X0+delta_t*T0/2        +delta_t*(k1*-8/27 +k2*2 +k3*-3544/2565 +k4*1859/4104 +k5*-11/40), vol_coords, vol_inds, lens_volume); 
%         
%         T_4th = T0 +k1*25/216 +k3*1408/2565 +k4*2197/4104 + k5*-1/5;
%         T_5th = T0 +k1*16/135 +k3*6656/12825 +k4*28561/56430 +k5*-9/50 +k6*2/55; %these coeffs should sum to 1
%         
%         % from Upgrading Runge-Kutta-Fehlberg Method (RKFM) for Second Order Ordinary Differential Equations (note error in 5th order coeffs)
%                 % in same manner as Scaraborough
%         % 4th order coeffs
%         k1C = 1408/2565*3/32 + 2197/4104*1932/2197 + -1/5*439/216;
%         %k2C = 1408/2565*9/32 + 2197/4104*-7200/2197 + -1/5*-8 %close to zero
%         k3C = 2197/4104*7296/2197 + -1/5*3680/513;
%         k4C = -1/5*-845/4104;
%         X_4th = X0 + delta_t*(T0+(k1*k1C + k3*k3C + k4*k4C));
%         
%         %5th order coeffs
%         k1C = 6656/12825*3/32 + 28561/56430*1932/2197 + -9/50*439/216 + 2/55*-8/27;
%         %k2C = 6656/12825*9/32 + 28561/56430*-7200/2197 + -9/50*-8 + 2/55*2 %close to zero
%         k3C = 28561/56430*7296/2197 + -9/50*3680/513 + 2/55*-3544/2565;
%         k4C = -9/50*-845/4104 + 2/55*1859/4104;
%         k5C = 2/55*-11/40; 
%         X_5th = X0 + delta_t*(T0+(k1*k1C + k3*k3C + k4*k4C + k5*k5C)); %these coeffs should sum to 0.5
%         
%         % Get error - normalize position terms to meters
%         % As eqn 16
%         er = max(abs([((X_5th')-(X_4th'))*10^3 T_5th'-T_4th']))/10^-12;
%             
%         %limit increase, if er goes to zero, will scale to inf.
%          if er < 0.01
%            er = 0.01; 
%          end
% 
%         %lower bound (0.6 or 0.8 ???) should be scaled reltive to saftey factor (0.9)
%             %so step does not decrease if it should increase
%         if er < 0.8 | er > 1
%            % As eqn 17
%            new_delta_t = 0.9*delta_t*er^0.2;
%            % had er inverted - why?
%            %new_delta_t = 0.9*delta_t/er^0.2;
%         else
%            new_delta_t = delta_t; 
%         end
%         
%         if isinf(new_delta_t) | isnan(er)
%            error('Problem!') 
%         end
%         
%         T = T_5th;
%         X = X_5th;
%         
%     elseif type == 4
%         %RK5 from Nishidate 2011 paper with fixed stepsize
%         % using coeffs from table 5
%         k1 = delta_t*numerical_dT_dt(X0, vol_coords, vol_inds, lens_volume); 
%         k2 = delta_t*numerical_dT_dt(X0+delta_t*T0/4        +delta_t*(k1/4), vol_coords, vol_inds, lens_volume); 
%         k3 = delta_t*numerical_dT_dt(X0+delta_t*T0*3/8      +delta_t*(k1*3/32 +k2*9/32), vol_coords, vol_inds, lens_volume); 
%         k4 = delta_t*numerical_dT_dt(X0+delta_t*T0*12/13    +delta_t*(k1*1932/2197 +k2*-7200/2197 +k3*7296/2197), vol_coords, vol_inds, lens_volume); 
%         k5 = delta_t*numerical_dT_dt(X0+delta_t*T0          +delta_t*(k1*439/216 +k2*-8 +k3*3680/513 +k4*-845/4104), vol_coords, vol_inds, lens_volume); 
%         k6 = delta_t*numerical_dT_dt(X0+delta_t*T0/2        +delta_t*(k1*-8/27 +k2*2 +k3*-3544/2565 +k4*1859/4104 +k5*-11/40), vol_coords, vol_inds, lens_volume); 
%         
%         T = T0 +k1*16/135 +k3*6656/12825 +k4*28561/56430 +k5*-9/50 +k6*2/55; %these b coeffs should sum to 1
%         
%         %5th order coeffs
%         k1C = 6656/12825*3/32 + 28561/56430*1932/2197 + -9/50*439/216 + 2/55*-8/27;
%         %k2C = 6656/12825*9/32 + 28561/56430*-7200/2197 + -9/50*-8 + 2/55*2 %close to zero
%         k3C = 28561/56430*7296/2197 + -9/50*3680/513 + 2/55*-3544/2565;
%         k4C = -9/50*-845/4104 + 2/55*1859/4104;
%         k5C = 2/55*-11/40; 
%         X = X0 + delta_t*(T0+(k1*k1C + k3*k3C + k4*k4C + k5*k5C)); %these be coeffs should sum to 0.5
%         
%         new_delta_t = delta_t;
%     else
%         error('Type not coded')
%     end
% end    

% function [dT_dt, n_of_X, dn_partial_dX] = numerical_dT_dt(X, vol_coords, vol_inds, lens_volume)
%         % get change in ray vector based on change in refractive index and derivative as in nishidate 2011 paper
%         
%         % take sampling points first
%         [sampling_inds, sampling_coords, point_dists] = getsamplingpoints(vol_coords, X);
%         
%         n_of_s = lens_volume(vol_inds(sampling_inds)); %nx1 vector of refractive indices at sampling points
% 
%         % As eqn 8
%         W_of_X = zeros(size(sampling_coords,1), size(sampling_coords,1)); %nxn : diagonal matrix of weights per sampling point
%         
%         for i = 1:size(sampling_coords,1)
%             %distance to sampling point
%             d = point_dists(i);
% 
%             %As eqn 11, quartic spline weighting function
%             W_of_X(i,i) = 1-6*d^2+8*d^3-3*d^4;
%         end
% 
%         % As eqn 7
%         P_of_X = zeros(size(sampling_coords,1),20); %nx20 : basis vector per sampling point
%         
%         for i = 1:size(sampling_coords,1)
%             %place basis function for each sampling point
% 
%             P_of_X(i,:) = calculatetricubicbasis(sampling_coords(i,1), sampling_coords(i,2), sampling_coords(i,3), 0)';
%         end
% 
%         % As eqn 5
%         M_of_X = P_of_X'*W_of_X*P_of_X; %20x20
%         
%         % As eqn 6
%         N_of_X = P_of_X'*W_of_X; %20xn
%         
%         % Solve eqn 4
%         a_of_X = M_of_X\(N_of_X*n_of_s); %20x1 : solves M*A=N*n for a, where N*n is 20x1
% 
%         [p_of_X, dp_partial_dX] = calculatetricubicbasis(X(1), X(2), X(3), 1); %20x1 and 20x3
%         
%         % As eqn 3
%         n_of_X = p_of_X*a_of_X; %1x1 : refractive index 
%         
%         % As eqn 9
%         dn_partial_dX = dp_partial_dX*a_of_X; %3x1 : spatial gradient of refractive index
%         
%         % As eqn 2
%         dT_dt = n_of_X*dn_partial_dX; %3x1 vector
% end