function [X, T, new_delta] = ray_interpolation(type, RIType, X0, T0, delta, vol_coords, vol_inds, lens_volume)
    %%% Selectect interpolation and grin type with flags.

    % Copied from anistropic graded ri test

    % Mostly RKF45 from Nishidate 2011 paper with variable step size 
    % using coeffs from Table 1

    if strcmp(RIType, 'iso')
        k1 = delta*numerical_dT_dt(X0,...
            vol_coords, vol_inds, lens_volume); 

        if type ~= '1'
            k2 = delta*numerical_dT_dt(X0 + delta*(T0/4) + delta*(k1/4),...
                vol_coords, vol_inds, lens_volume); 

            k3 = delta*numerical_dT_dt(X0 + delta*T0*(3/8) + delta*(k1*3/32 +k2*9/32),...
                vol_coords, vol_inds, lens_volume); 

            k4 = delta*numerical_dT_dt(X0 + delta*T0*(12/13)+ ...
                delta*(k1*1932/2197 +k2*-7200/2197 +k3*7296/2197),...
                vol_coords, vol_inds, lens_volume); 

            k5 = delta*numerical_dT_dt(X0 + delta*T0 + ...
                delta*(k1*439/216 +k2*-8 +k3*3680/513 +k4*-845/4104),...
                vol_coords, vol_inds, lens_volume); 

            k6 = delta*numerical_dT_dt(X0 + delta*T0*(1/2) + ...
                delta*(k1*-8/27 +k2*2 +k3*-3544/2565 +k4*1859/4104 +k5*-11/40),...
                vol_coords, vol_inds, lens_volume); 

        end
    elseif strcmp(RIType, 'aniso')
        
        error('Adapted treatment from *Upgrading ... for second order* paper, seems wrong')
        % treat as anisotropic
        k1 = numerical_dt_ds(X0, ...
            T0,...
            vol_coords, vol_inds, lens_volume); 

        % New where both X and T are modified and delta shifted down
        if type ~= '1'
            negT = 1; % set to 0 for testing
            
            %%% Adjust to match iso
            k2 = numerical_dt_ds(X0 + delta*(T0/4), ...
                T0 +delta*(k1/4)*negT,...
                vol_coords, vol_inds, lens_volume); 

            % Expanded eqns for X on terms below - all but 2 equal or less than 5e-16 off
            k3 = numerical_dt_ds(X0 + delta*T0*(3/32+9/32) + delta^2*(k1*(9/32*1/4)),  ...      
                T0 +delta*(k1*3/32 + k2*9/32)*negT,...
                vol_coords, vol_inds, lens_volume); 

            k4 = numerical_dt_ds(X0 + delta*T0*(1932/2197-7200/2197+7296/2197) + ...
                delta^2*(k1*(-7200/2197*1/4 + 7296/2197*3/32) + k2*(7296/2197*9/32)),  ...
                T0 +delta*(k1*1932/2197 +k2*-7200/2197 +k3*7296/2197)*negT,...
                vol_coords, vol_inds, lens_volume); 

            %%% Difference in K1 terms for X
            k5 = numerical_dt_ds(X0 + delta*T0*(439/216 - 8 +3680/513 - 845/4104) + ...
                delta^2*(k1*(-8/4 + 3680/513*3/32 - 845/4104*1932/2197) +k2*(3680/513*9/32 -845/4104*-7200/2197) +k3*(-845/4104*7296/2197)), ...
                T0 +delta*(k1*439/216 +k2*-8 +k3*3680/513 +k4*-845/4104)*negT,...
                vol_coords, vol_inds, lens_volume); 

            %%% Difference in K4 term for X
            k6 = numerical_dt_ds(X0 + delta*T0*(-8/27 + 2 - 3544/2565 + 1859/4104 - 11/40) + ...
                delta^2*(k1*(2*1/4 - 3544/2565*3/32 +1859/4104*1932/2197 -11/40*439/216) +k2*(-3544/2565*9/32 + 1859/4104*-7200/2197 - 11/40*-8) +k3*(1859/4104*7296/2197  - 11/40*3680/513) +k4*(-11/40*-845/4104)), ...
                T0 +delta*(k1*-8/27 +k2*2 +k3*-3544/2565 +k4*1859/4104 +k5*-11/40)*negT,...
                vol_coords, vol_inds, lens_volume); 

        end
    else
       error('unknown RI type') 
    end
        
    
    if type == '1'
        % Euler solver
        X = X0 + delta*T0;
        T = T0 + k1;
        
        new_delta = delta;
        
    elseif type == '4'

        T = T0 + (k1*25/216 +k3*1408/2565 +k4*2197/4104 + k5*-1/5);

        % Not really sure these are correct
            %%% Would be good to get alternate reference to Upgrading paper
        k1C = 1408/2565*3/32 + 2197/4104*1932/2197 + -1/5*439/216;
        k2C = 1408/2565*9/32 + 2197/4104*-7200/2197 + -1/5*-8; %close to zero
        k3C = 2197/4104*7296/2197 + -1/5*3680/513;
        k4C = -1/5*-845/4104;
        X = X0 + delta*(T0 + k1*k1C + k2*k2C + k3*k3C + k4*k4C);
        
        % From Shwarma, generaly taken as reference.
%         A = delta*numerical_dT_dt(X0, vol_coords, vol_inds, lens_volume); 
%         B = delta*numerical_dT_dt(X0+delta*T0/2    +delta*A/8, vol_coords, vol_inds, lens_volume);
%         C = delta*numerical_dT_dt(X0+delta*T0    +delta*B/2, vol_coords, vol_inds, lens_volume);
%         
%         X2 = X0 + delta*(T0 + (A + 2*B)/6); %b coeffs should add to 0.5
%         T2 = T0 + (A + 4*B + C)/6;
%         
%         T - T2
%         
%         X - X2
        
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
        error('Interpolation type not coded')
    end
end    

