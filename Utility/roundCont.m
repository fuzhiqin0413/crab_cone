function [roundedVector, roundDirection] = roundCont(vector, roundDirection)
    % Continous rounding for vectors
    
    if length(vector) ~= 3
       error('Vector length should be 3'); 
    end
    
    if isempty(roundDirection)
       % need to set rounding direction
       roundDirection = zeros(3,1);
       
       for i = 1:3
            test = round(vector(i)) - floor(vector(i));
            if test == 0
                % round down
                roundDirection(i) = -1;
                
            elseif test == 1
                % round up
                roundDirection(i) = 1;
            end
       end
    end
    
    roundedVector = zeros(3,1);
    
    for i = 1:3
        if roundDirection(i) == 1
            roundedVector(i) = ceil(vector(i));
        elseif roundDirection(i) == -1
            roundedVector(i) = floor(vector(i));
        end
    end
    
end

