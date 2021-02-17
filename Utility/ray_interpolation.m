function [X, T, new_delta] = ray_interpolation(type, RIType, X0, T0, delta, vol_coords, vol_inds, lens_volume)
    %%% Selectect interpolation and grin type with flags.

    % Copied from anistropic graded ri test

    % Mostly RKF45 from Nishidate 2011 paper with variable step size 
    % using coeffs from Table 1

    new_delta = delta;
            
%     if strcmp(RIType, 'iso')
% 
%     elseif strcmp(RIType, 'aniso')
%         
%         error('Adapted treatment from *Upgrading ... for second order* paper, seems wrong')
%         % treat as anisotropic
%         k1 = numerical_dt_ds(X0, ...
%             T0,...
%             vol_coords, vol_inds, lens_volume); 
% 
%         % New where both X and T are modified and delta shifted down
%         if type ~= '1'
%             negT = 1; % set to 0 for testing
%             
%             %%% Adjust to match iso
%             k2 = numerical_dt_ds(X0 + delta*(T0/4), ...
%                 T0 +delta*(k1/4)*negT,...
%                 vol_coords, vol_inds, lens_volume); 
% 
%             % Expanded eqns for X on terms below - all but 2 equal or less than 5e-16 off
%             k3 = numerical_dt_ds(X0 + delta*T0*(3/32+9/32) + delta^2*(k1*(9/32*1/4)),  ...      
%                 T0 +delta*(k1*3/32 + k2*9/32)*negT,...
%                 vol_coords, vol_inds, lens_volume); 
% 
%             k4 = numerical_dt_ds(X0 + delta*T0*(1932/2197-7200/2197+7296/2197) + ...
%                 delta^2*(k1*(-7200/2197*1/4 + 7296/2197*3/32) + k2*(7296/2197*9/32)),  ...
%                 T0 +delta*(k1*1932/2197 +k2*-7200/2197 +k3*7296/2197)*negT,...
%                 vol_coords, vol_inds, lens_volume); 
% 
%             %%% Difference in K1 terms for X
%             k5 = numerical_dt_ds(X0 + delta*T0*(439/216 - 8 +3680/513 - 845/4104) + ...
%                 delta^2*(k1*(-8/4 + 3680/513*3/32 - 845/4104*1932/2197) +k2*(3680/513*9/32 -845/4104*-7200/2197) +k3*(-845/4104*7296/2197)), ...
%                 T0 +delta*(k1*439/216 +k2*-8 +k3*3680/513 +k4*-845/4104)*negT,...
%                 vol_coords, vol_inds, lens_volume); 
% 
%             %%% Difference in K4 term for X
%             k6 = numerical_dt_ds(X0 + delta*T0*(-8/27 + 2 - 3544/2565 + 1859/4104 - 11/40) + ...
%                 delta^2*(k1*(2*1/4 - 3544/2565*3/32 +1859/4104*1932/2197 -11/40*439/216) +k2*(-3544/2565*9/32 + 1859/4104*-7200/2197 - 11/40*-8) +k3*(1859/4104*7296/2197  - 11/40*3680/513) +k4*(-11/40*-845/4104)), ...
%                 T0 +delta*(k1*-8/27 +k2*2 +k3*-3544/2565 +k4*1859/4104 +k5*-11/40)*negT,...
%                 vol_coords, vol_inds, lens_volume); 
% 
%         end
%     else
%        error('unknown RI type') 
%     end
    
    if strcmp(RIType, 'aniso')
       error('Need to re formulate above correctly') 
    end
    
    switch type 
        case '1'
            % Euler solver
            X = X0 + delta*T0;
            T = T0 + delta*numerical_dT_dt(X0, vol_coords, vol_inds, lens_volume); 
        
        case '4S'
            %%% All three 4th order give similar results
            
            % From Shwarma, generaly taken as reference.
            %%% Weird note - slight numerical difference appears on T if terms not bracketed
            A = delta*numerical_dT_dt(X0, vol_coords, vol_inds, lens_volume); 
            B = delta*numerical_dT_dt(X0+delta*(T0/2+A/8), vol_coords, vol_inds, lens_volume);
            C = delta*numerical_dT_dt(X0+delta*(T0+B/2), vol_coords, vol_inds, lens_volume);

            X = X0 + delta*(T0 + (A + 2*B)/6); %b coeffs should add to 0.5
            T = T0 + (A + 4*B + C)/6;

            % RKN method from Kryzig - matches Shwarma
    %         k1N = 0.5*delta*numerical_dT_dt(X0, vol_coords, vol_inds, lens_volume); 
    %         k2N = 0.5*delta*numerical_dT_dt(X0+ 0.5*delta*(T0 + 0.5*k1N), vol_coords, vol_inds, lens_volume);
    %         k4N = 0.5*delta*numerical_dT_dt(X0+ delta*(T0 + B/2), vol_coords, vol_inds, lens_volume);
    %         
    %         XN = X0 + delta*(T0 + (k1N + 2*k2N)/3); %b coeffs should add to 0.5
    %         TN = T0 + (k1N + 4*k2N + k4N)/3;
        
        case '4F'
            % Original Felberg I was using, adapted for 2nd order

            k1 = delta*numerical_dT_dt(X0,...
                vol_coords, vol_inds, lens_volume); 

            k2 = delta*numerical_dT_dt(X0 + delta*(T0/4 + k1/4),...
                vol_coords, vol_inds, lens_volume); 

            k3 = delta*numerical_dT_dt(X0 + delta*(T0*3/8 + k1*3/32 +k2*9/32),...
                vol_coords, vol_inds, lens_volume); 

            k4 = delta*numerical_dT_dt(X0 + delta*(T0*12/13+ ...
                k1*1932/2197 +k2*-7200/2197 +k3*7296/2197),...
                vol_coords, vol_inds, lens_volume); 

            k5 = delta*numerical_dT_dt(X0 + delta*(T0 + ...
                k1*439/216 +k2*-8 +k3*3680/513 +k4*-845/4104),...
                vol_coords, vol_inds, lens_volume); 

%             k6 = delta*numerical_dT_dt(X0 + delta*(T0/2 + ...
%                 k1*-8/27 +k2*2 +k3*-3544/2565 +k4*1859/4104 +k5*-11/40),...
%                 vol_coords, vol_inds, lens_volume); 

            T = T0 + (k1*25/216 +k3*1408/2565 +k4*2197/4104 + k5*-1/5);

            % Not really sure these are correct
                %%% Would be good to get alternate reference to Upgrading paper
            k1C = 1408/2565*3/32 + 2197/4104*1932/2197 + -1/5*439/216;
            k2C = 1408/2565*9/32 + 2197/4104*-7200/2197 + -1/5*-8; %close to zero
            k3C = 2197/4104*7296/2197 + -1/5*3680/513;
            k4C = -1/5*-845/4104;
            X = X0 + delta*(T0 + k1*k1C + k2*k2C + k3*k3C + k4*k4C); 
        
        case '4B'
            % From Bettis 1973 FKN 4(5)
            k0B = numerical_dT_dt(X0, vol_coords, vol_inds, lens_volume); 
            k1B = numerical_dT_dt(X0 + delta*(T0/8 + delta*(k0B/128)), vol_coords, vol_inds, lens_volume); 
            k2B = numerical_dT_dt(X0 + delta*(T0/4 + delta*(k0B/96 + k1B/48)), vol_coords, vol_inds, lens_volume); 
            k3B = numerical_dT_dt(X0 + delta*(T0/2 + delta*(k0B/24 + k2B/12)), vol_coords, vol_inds, lens_volume); 
            k4B = numerical_dT_dt(X0 + delta*(T0*3/4 + delta*(k0B*9/128 + k2B*9/64 + k3B*9/128)), vol_coords, vol_inds, lens_volume); 
%             k5B = numerical_dT_dt(X0 + delta*(T0 + delta*(k0B*7/90 + k2B*4/15 + k3B/15 + k4B*4/45)), vol_coords, vol_inds, lens_volume); 

            X = X0 + delta*(T0 + delta*(k0B/6 + k3B/3)); 
            T = T0 + delta*(k2B*2/3 + k3B*-1/3 + k4B*2/3);
        
        case '45'    
            % Changed to using Bettis formula as I am more confident of the coefficients.
            k0B = numerical_dT_dt(X0, vol_coords, vol_inds, lens_volume); 
            k1B = numerical_dT_dt(X0 + delta*(T0/8 + delta*(k0B/128)), vol_coords, vol_inds, lens_volume); 
            k2B = numerical_dT_dt(X0 + delta*(T0/4 + delta*(k0B/96 + k1B/48)), vol_coords, vol_inds, lens_volume); 
            k3B = numerical_dT_dt(X0 + delta*(T0/2 + delta*(k0B/24 + k2B/12)), vol_coords, vol_inds, lens_volume); 
            k4B = numerical_dT_dt(X0 + delta*(T0*3/4 + delta*(k0B*9/128 + k2B*9/64 + k3B*9/128)), vol_coords, vol_inds, lens_volume); 
            % Can use final term as k0B in next step
            k5B = numerical_dT_dt(X0 + delta*(T0 + delta*(k0B*7/90 + k2B*4/15 + k3B/15 + k4B*4/45)), vol_coords, vol_inds, lens_volume); 
            
            X = X0 + delta*(T0 + delta*(k0B*7/90 + k2B*4/15 + k3B*1/15 + k4B*4/45)); 
            T = T0 + delta*(k0B*7/90 + k2B*16/45 + k3B*2/15 + k4B*16/45 + k5B*7/90);

            % These do not seem to work.
            %TEx = 4/15*(k0B*1/3 - k2B + k3B - k4B*1/3);
            %TEt = 7/15*(-k0B*1/6 + k2B*2/3 + k3B + k4B*2/3 - k5B*1/6);
            
            X4 = X0 + delta*(T0 + delta*(k0B/6 + k3B/3)); 
            T4 = T0 + delta*(k2B*2/3 + k3B*-1/3 + k4B*2/3);
            
            TEx = X4-X;
            TEt = T4-T;
            
            % Get error - normalize position terms to meters
            er = 10^-12/max(abs([TEx'*10^3 TEt']));
            
            if er < 0.8 | er > 1
               new_delta = 0.9*delta*er^0.2;
            else
               new_delta = delta; 
            end

            if isinf(new_delta) | isnan(er)
               error('Problem!') 
            end
        
        case '5B'    
            % Both of these end up very close the Shwarma
            k0B = numerical_dT_dt(X0, vol_coords, vol_inds, lens_volume); 
            k1B = numerical_dT_dt(X0 + delta*(T0/8 + delta*(k0B/128)), vol_coords, vol_inds, lens_volume); 
            k2B = numerical_dT_dt(X0 + delta*(T0/4 + delta*(k0B/96 + k1B/48)), vol_coords, vol_inds, lens_volume); 
            k3B = numerical_dT_dt(X0 + delta*(T0/2 + delta*(k0B/24 + k2B/12)), vol_coords, vol_inds, lens_volume); 
            k4B = numerical_dT_dt(X0 + delta*(T0*3/4 + delta*(k0B*9/128 + k2B*9/64 + k3B*9/128)), vol_coords, vol_inds, lens_volume); 
            k5B = numerical_dT_dt(X0 + delta*(T0 + delta*(k0B*7/90 + k2B*4/15 + k3B/15 + k4B*4/45)), vol_coords, vol_inds, lens_volume); 

            X = X0 + delta*(T0 + delta*(k0B*7/90 + k2B*4/15 + k3B*1/15 + k4B*4/45)); 
            T = T0 + delta*(k0B*7/90 + k2B*16/45 + k3B*2/15 + k4B*16/45 + k5B*7/90);
            
        case '5F'
           
            k1 = delta*numerical_dT_dt(X0,...
                vol_coords, vol_inds, lens_volume); 

            k2 = delta*numerical_dT_dt(X0 + delta*(T0/4 + k1/4),...
                vol_coords, vol_inds, lens_volume); 

            k3 = delta*numerical_dT_dt(X0 + delta*(T0*3/8 + k1*3/32 +k2*9/32),...
                vol_coords, vol_inds, lens_volume); 

            k4 = delta*numerical_dT_dt(X0 + delta*(T0*12/13+ ...
                k1*1932/2197 +k2*-7200/2197 +k3*7296/2197),...
                vol_coords, vol_inds, lens_volume); 

            k5 = delta*numerical_dT_dt(X0 + delta*(T0 + ...
                k1*439/216 +k2*-8 +k3*3680/513 +k4*-845/4104),...
                vol_coords, vol_inds, lens_volume); 

            k6 = delta*numerical_dT_dt(X0 + delta*(T0/2 + ...
                k1*-8/27 +k2*2 +k3*-3544/2565 +k4*1859/4104 +k5*-11/40),...
                vol_coords, vol_inds, lens_volume); 
            
            T = T0 + (k1*16/135 +k3*6656/12825 +k4*28561/56430 +k5*-9/50 +k6*2/55); %these b coeffs should sum to 1

            %5th order coeffs - from paper above
            k1C = 6656/12825*3/32 + 28561/56430*1932/2197 + -9/50*439/216 + 2/55*-8/27;
            k3C = 28561/56430*7296/2197 + -9/50*3680/513 + 2/55*-3544/2565;
            k4C = -9/50*-845/4104 + 2/55*1859/4104;
            k5C = 2/55*-11/40;
            
            X = X0 + delta*(T0+(k1*k1C + k3*k3C + k4*k4C + k5*k5C)); %these be coeffs should sum to 0.5

        otherwise 
            error('Interpolation type not coded')
    end
end    

