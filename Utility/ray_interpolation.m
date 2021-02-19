function [X, T, new_delta] = ray_interpolation(type, RIType, X0, T0, delta, vol_coords, vol_inds, lens_volume, tolerance)
    %%% Select interpolation and grin type with flags.

    new_delta = delta;
           
    % Calc function values for each interp type
    switch RIType
        case 'iso'
            f0 = numerical_dT_dt(X0, vol_coords, vol_inds, lens_volume);

            if strcmp(type, '4S')
                    % Sharma's 4th order - keep standard sytnax
                    %%% Weird note - slight numerical difference appears on T if terms not bracketed
                    A = delta*f0;
                    B = delta*numerical_dT_dt(X0+delta*(T0/2+A/8), vol_coords, vol_inds, lens_volume);
                    C = delta*numerical_dT_dt(X0+delta*(T0+B/2), vol_coords, vol_inds, lens_volume);
                    
            elseif strcmp(type, '4RKN') | strcmp(type, '5RKN') | strcmp(type, '45RKN')
                % Get 4th order fns - From Bettis 1973 RKN(4)5
                f1 = numerical_dT_dt(X0 + delta*(T0/8 + delta*(f0/128)), vol_coords, vol_inds, lens_volume); 
                f2 = numerical_dT_dt(X0 + delta*(T0/4 + delta*(f0/96 + f1/48)), vol_coords, vol_inds, lens_volume); 
                f3 = numerical_dT_dt(X0 + delta*(T0/2 + delta*(f0/24 + f2/12)), vol_coords, vol_inds, lens_volume); 
                f4 = numerical_dT_dt(X0 + delta*(T0*3/4 + delta*(f0*9/128 + f2*9/64 + f3*9/128)), vol_coords, vol_inds, lens_volume); 
            end
            
            if strcmp(type, '5RKN') | strcmp(type, '45RKN')
                %Also get 5th fns - Can apparently use this term as k0B in next step
                f5 = numerical_dT_dt(X0 + delta*(T0 + delta*(f0*7/90 + f2*4/15 + f3/15 + f4*4/45)), vol_coords, vol_inds, lens_volume); 

            end
        case 'aniso'
           error('Need to re formulate')    
        
        otherwise
            error('unknown RI type') 
    end
    
    % Then calculate values
    switch type 
        case '1'
            % Euler solver
            X = X0 + delta*T0;
            T = T0 + delta*f0; 
        
        case '4S'

            % From Sharma, generaly taken as reference.
            X = X0 + delta*(T0 + (A + 2*B)/6); %b coeffs should add to 0.5
            T = T0 + (A + 4*B + C)/6;
        
        case '4RKN'
            % From Bettis 1973 RKN(4)5
            X = X0 + delta*(T0 + delta*(f0/6 + f3/3)); 
            T = T0 + delta*(f2*2/3 + f3*-1/3 + f4*2/3);
        
        case '45RKN'    
            % Get 5th order values
            X = X0 + delta*(T0 + delta*(f0*7/90 + f2*4/15 + f3*1/15 + f4*4/45)); 
            T = T0 + delta*(f0*7/90 + f2*16/45 + f3*2/15 + f4*16/45 + f5*7/90);

            % These do not seem to work.
            %TEx = 4/15*(f0*1/3 - f2 + f3 - f4*1/3);
            %TEt = 7/15*(-f0*1/6 + f2*2/3 + f3 + f4*2/3 - f5*1/6);
            
            % So get 4th order values
            X4 = X0 + delta*(T0 + delta*(f0/6 + f3/3)); 
            T4 = T0 + delta*(f2*2/3 + f3*-1/3 + f4*2/3);
            
            % Then take difference
            TEx = X4-X;
            TEt = T4-T;
            
            % Kept calculation from Nishidate
            % Get error - normalize position terms to meters
            er = tolerance/max(abs([TEx'*10^3 TEt']));
            
            if er < 0.8 | er > 1
               new_delta = 0.9*delta*er^0.2;
            else
               new_delta = delta; 
            end

            if isinf(new_delta) | isnan(er) | new_delta == 0
               error('Problem!') 
            end
        
        case '5RKN'    

            X = X0 + delta*(T0 + delta*(f0*7/90 + f2*4/15 + f3*1/15 + f4*4/45)); 
            T = T0 + delta*(f0*7/90 + f2*16/45 + f3*2/15 + f4*16/45 + f5*7/90);

        otherwise 
            error('Interpolation type not coded')
    end
    
    % Seems like ray should be normalized, but adds lots of error.
%     T = T/norm(T);
end    

