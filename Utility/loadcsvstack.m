function [ arrayOut ] = loadcsvstack( address, loadSmall)
% Copy tiff stack saved by amira into volume. 
    %colour image will be converted to greyscale, but unique values for
    %colours are not guaranteed
% address - is folder to load
% loadSmall - flag, if set will arrayOut will be in uint8 format, 
              %otherwise, will be double
                   
    currentDirectory = pwd;          
              
    % Find which files to load.
    cd(address); directoryContents = dir; includedFiles = [];
    
    % natural langueg sort then remove from cells
    directoryNames = natsortfiles({directoryContents.name});
    
    %directoryNames = {directoryContents.name};

    for i = 1:length(directoryContents)
        if strfind(directoryNames{i}, 'csv')   
            
            includedFiles = [ includedFiles i ];
        end
    end
    
    % Load each file into an array.
    counter = 1; nImages = length(includedFiles);
    
    for  i = 1:nImages
        
        % Initalize array given dimensions of first image.
        if i == 1
            
            [tempImage] = readmatrix(directoryNames{includedFiles(i)}, 'TreatAsMissing', '-');
            
            % Make square
            size2Use = max(size(tempImage));
            
            %%% Note cropped dimension is offset up from bottom
            
            if loadSmall
                
                arrayOut = zeros(size2Use, size2Use, nImages, 'uint8')*NaN;
                
                arrayOut(size2Use-size(tempImage,1)+1:end,size2Use-size(tempImage,2)+1:end,counter) = uint8(tempImage);
                
            else
                
                arrayOut = zeros(size2Use, size2Use, nImages)*NaN;
                
                arrayOut(size2Use-size(tempImage,1)+1:end,size2Use-size(tempImage,2)+1:end,counter) = tempImage;
                
            end
            
        else
            
            tempImage = readmatrix(directoryNames{includedFiles(i)}, 'TreatAsMissing', '-');
            
            if loadSmall
                
                arrayOut(size2Use-size(tempImage,1)+1:end,size2Use-size(tempImage,2)+1:end,counter) = uint8(tempImage);
                
            else
                
                arrayOut(size2Use-size(tempImage,1)+1:end,size2Use-size(tempImage,2)+1:end,counter) = tempImage;
                
            end
            
        end
        
        counter = counter + 1;
    
        pause(0.01)
    end
    
    cd(currentDirectory)
    
end

