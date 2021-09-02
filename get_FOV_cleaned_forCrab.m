%main file for doing raytracing
%split into three main sections
%1 loading and preprocessing
%2 initilizing rays origins and targets
%3 actual ray tracing

close all

%load files and preprocess volume.
if 1
    clear;
    
    outPutFile = 'fineTracerOutputMed.mat';
    
    preLoadCoarseRays = 0; %can do a test with few points on sphere and then refine (looks in folderToSave)
    coarseRayFile = '';
    
    saveOutPut = 0;
    folderToSave = '';

    folderToGet = '/Users/gavintaylor/Documents/Company/Client Projects/Cones MPI/Data/Sync CT/Amira project/Matlab/Inital_1_cone'; 
    sphericalAxis = [];
    
    pixelScale = 2*0.16;
    dsP = 30; %down sampling for plotting surfaces
    dsFact = 1; %resize input volumes if desired
    
    cd(folderToGet)
    dirContents = dir; 
    includedFiles = [];
    for i = 1:length(dirContents)
        if strfind(dirContents(i).name, 'tif')     
            includedFiles = [ includedFiles i ];
        end
    end

    numberOfImages = length(includedFiles);

    counter = 1;
    for  i = 1:floor(length(includedFiles))
        if strfind(dirContents(includedFiles(i)).name, 'tif')
            if i == 1
                warning('can this be updated to _w_by_h version? check following usage...');
                [temp, map] = imread(dirContents(includedFiles(i)).name, 'TIFF');
                imageStack = zeros(size(temp,1), size(temp,2), numberOfImages);
                imageStack(:,:,counter) = temp;
                counter = counter + 1;
            else
                imageStack(:,:,counter) = imread(dirContents(includedFiles(i)).name);
                counter = counter + 1;
            end
        end
    end

    %get volume indicies - these are surfaces, no need to find later
    ToAirInds = find(imageStack == 8);
    ToInterInds = find(imageStack == 9);
    ToCorneaInds = find(imageStack == 10);
    ConeInds = find(imageStack == 7);
    
    %small close to remove mesh problems
    borderVolume = zeros(size(imageStack));
    borderVolume(ToCorneaInds) = 1; borderVolume(ToInterInds) = 2; borderVolume(ToAirInds) = 3;
    coneVolume = zeros(size(imageStack));
    coneVolume(ConeInds) = 1;
    
    %resize if dsFact set
    newSize = round(size(imageStack)/dsFact);
%    generalVolume = round(resize(generalVolume, newSize));

    %ToAirInds = find(generalVolume == 3);
    %ToInterInds = find(generalVolume == 2);
    %ToCorneaInds = find(generalVolume == 1);

    %get subscripts of surfaces
    [corneaSurfX, corneaSurfY, corneaSurfZ] = ind2sub(newSize, ToCorneaInds);
    [interSurfX, interSurfY, interSurfZ] = ind2sub(newSize, ToInterInds);
    [airSurfX, airSurfY, airSurfZ] = ind2sub(newSize, ToAirInds);
    
    corneaSurf = [corneaSurfX, corneaSurfY, corneaSurfZ];
    interconeSurf = [interSurfX, interSurfY, interSurfZ];
    internalSurf = [airSurfX, airSurfY, airSurfZ];
    fullConeSurf = [[airSurfX; interSurfX], [airSurfY; interSurfY] , [airSurfZ; interSurfZ]];
    
    %find max range of labels
    xFullRange = [min([corneaSurfX', interSurfX', airSurfX'])-1 max([corneaSurfX', interSurfX', airSurfX'])+1];
    yFullRange = [min([corneaSurfY', interSurfY', airSurfY'])-1 max([corneaSurfY', interSurfY', airSurfY'])+1];
    zFullRange = [min([corneaSurfZ', interSurfZ', airSurfZ'])-1 max([corneaSurfZ', interSurfZ', airSurfZ'])+1];

    maxD = max([xFullRange(2) yFullRange(2) zFullRange(2)]);

    corneaMean = mean(corneaSurf);
    coneMean = mean(fullConeSurf);
    
    %get axes
    corneaAxis = pca([corneaSurfX-corneaMean(1), corneaSurfY-corneaMean(2), corneaSurfZ-corneaMean(3)]);
    corneaAxis = corneaAxis(:,3);
    coneAxes = pca([[airSurfX; interSurfX]-coneMean(1), [airSurfY; interSurfY]-coneMean(2), [airSurfZ; interSurfZ]-coneMean(3)]);
    coneAxis = coneAxes(:,1);
    
    %get limits along major axis
    rerToAxis = matrix2rotatevectors([0 0 1]',coneAxis'); %  
    rerToCone = matrix2rotatevectors(coneAxis',[0 0 1]');
    tempRot = [[airSurfX; interSurfX]-coneMean(1), [airSurfY; interSurfY]-coneMean(2), [airSurfZ; interSurfZ]-coneMean(3)]*rerToAxis;
    %figure; plot3(tempRot(:,1), tempRot(:,2), tempRot(:,3), '.'); axis equal
    minZ = min(floor(tempRot(:,3)));
    maxZ = max(ceil(tempRot(:,3)));
    
    %% make sliding plane for cone
    planeRange = 200; %%%should define range dynamically
    planeSubs = zeros((planeRange*2+1)^2,3);
    planeImageRef = zeros((planeRange*2+1)^2,2);
    lineSubs = [(-planeRange:planeRange)'*coneAxes(1,2), (-planeRange:planeRange)'*coneAxes(2,2), (-planeRange:planeRange)'*coneAxes(3,2)];
    ind = 1;
    for i = -planeRange:planeRange
        planeSubs(ind:ind+planeRange*2,:) = round([coneAxes(1,3)*i+lineSubs(:,1), coneAxes(2,3)*i+lineSubs(:,2), coneAxes(3,3)*i+lineSubs(:,3) ]);
        planeImageRef(ind:ind+planeRange*2,:) = [(i+planeRange+1)*ones(planeRange*2+1,1), (-planeRange:planeRange)'+planeRange+1];
        ind = ind+planeRange*2+1;
    end
    [planeSubs, tempInds] = unique((planeSubs),'rows');
    planeImageRef = planeImageRef(tempInds,:);
    %% make gradient map 
    gradMap = zeros(newSize);
    numDiscSteps = 20;
    lastZ = minZ;
    firstZ = 0;
    bordersInVolume = borderVolume*0;
    
    %method of using distance from center line
%     for i = minZ:maxZ
%        testPlane = round([coneMean(1)+planeSubs(:,1)+coneAxis(1)*i, coneMean(2)+planeSubs(:,2)+coneAxis(2)*i, coneMean(3)+planeSubs(:,3)+coneAxis(3)*i]);
%        planeCenter = [coneMean(1)+coneAxis(1)*i, coneMean(2)+coneAxis(2)*i, coneMean(3)+coneAxis(3)*i];
%        
%        toUseInds = find(testPlane(:,1) >= 1 & testPlane(:,1) <= newSize(1) & testPlane(:,2) >= 1 & testPlane(:,2) <= newSize(2) ...
%            & testPlane(:,3) >= 1 & testPlane(:,3) <= newSize(3));
%        testInds = sub2ind(newSize, testPlane(toUseInds,1), testPlane(toUseInds,2), testPlane(toUseInds,3));
%        %don't assign values if we are at the level of cornea
%        if sum(borderVolume(testInds) == 1) == 0 
%%% add cone size check if used later
%            indsInCone = find(coneVolume(testInds));
%            coneInPlane = testPlane(toUseInds(indsInCone),:);
%            distToCenter = sqrt((coneInPlane(:,1)-planeCenter(1)).^2+(coneInPlane(:,2)-planeCenter(2)).^2+(coneInPlane(:,3)-planeCenter(3)).^2);
%            distToCenter = distToCenter-min(distToCenter);
%            distToCenter = ceil(distToCenter/max(distToCenter)*numDiscSteps);
%            gradMap(testInds(indsInCone)) = distToCenter;
%            lastZ = i;
%        else
%           break; 
%        end
%     end

    %method of using dist from border, more useful later on
    for i = minZ:maxZ
       testPlane = round([coneMean(1)+planeSubs(:,1)+coneAxis(1)*i, coneMean(2)+planeSubs(:,2)+coneAxis(2)*i, coneMean(3)+planeSubs(:,3)+coneAxis(3)*i]);

       toUseInds = find(testPlane(:,1) >= 1 & testPlane(:,1) <= newSize(1) & testPlane(:,2) >= 1 & testPlane(:,2) <= newSize(2) ...
           & testPlane(:,3) >= 1 & testPlane(:,3) <= newSize(3));
       testInds = sub2ind(newSize, testPlane(toUseInds,1), testPlane(toUseInds,2), testPlane(toUseInds,3));
       %don't assign values if we are at the level of cornea and need to a few inds on border
       if sum(borderVolume(testInds) == 1) == 0 
           %check we are in a decet size portion of cone, otherwise it can
           %be weird
           if sum(coneVolume(testInds)) > 10000
               %[sum(borderVolume(testInds) > 1)  sum(coneVolume(testInds))]

               if firstZ == 0
                  firstZ = i; 
               end
               
               %create image with border cones included
               imForDist = zeros(planeRange*2+1,planeRange*2+1);

               %get cone border inds in this slice
               indsOnBorder = find(borderVolume(testInds) > 1);
               indsInCone = find(coneVolume(testInds));
               
               %label inds on borer for later interp
               bordersInVolume(testInds(indsOnBorder)) = 1;
               
               %place into image - need to be careful of image flipping...
               imInds = sub2ind([planeRange*2+1,planeRange*2+1], planeImageRef(toUseInds(indsOnBorder),1), planeImageRef(toUseInds(indsOnBorder),2));
               imForDist(imInds) = 1;
               
%                figure; imshow(imForDist);
%                figure; hold on
%                plot(planeImageRef(toUseInds(indsOnBorder),1), planeImageRef(toUseInds(indsOnBorder),2),'.b')
%                plot(planeImageRef(toUseInds(indsInCone),1), planeImageRef(toUseInds(indsInCone),2),'.r')
               
               %get dist map, and then dists associated with cone
               imForDist = bwdist(imForDist);
               imIndsForCones = sub2ind([planeRange*2+1,planeRange*2+1], planeImageRef(toUseInds(indsInCone),1), planeImageRef(toUseInds(indsInCone),2));
               distToBorders = -imForDist(imIndsForCones);

%                coneInPlane = testPlane(toUseInds(indsInCone),:);
%                distToCenter = sqrt((coneInPlane(:,1)-planeCenter(1)).^2+(coneInPlane(:,2)-planeCenter(2)).^2+(coneInPlane(:,3)-planeCenter(3)).^2);
               distToBorders = distToBorders-min(distToBorders);
               distToBorders = ceil(distToBorders/max(distToBorders)*numDiscSteps);
               gradMap(testInds(indsInCone)) = distToBorders;
               lastZ = i;
           end
       else
          break; 
       end
    end
    
    %close is better than scattered interp
    %all values can't be put into interp because of memory limit, and
    %putting in limited number of points picks weird neighbours.
    gradMap = imclose(gradMap, strel('sphere',1));
    
    %get values from first and last Z planes to interpolate
    %last
    testPlaneLast = round([coneMean(1)+planeSubs(:,1)+coneAxis(1)*lastZ, coneMean(2)+planeSubs(:,2)+coneAxis(2)*lastZ, coneMean(3)+planeSubs(:,3)+coneAxis(3)*lastZ]);
    toUseIndsLast = find(testPlaneLast(:,1) >= 1 & testPlaneLast(:,1) <= newSize(1) & testPlaneLast(:,2) >= 1 & testPlaneLast(:,2) <= newSize(2) ...
           & testPlaneLast(:,3) >= 1 & testPlaneLast(:,3) <= newSize(3));
    
    testIndsLast = sub2ind(newSize, testPlaneLast(toUseIndsLast,1), testPlaneLast(toUseIndsLast,2), testPlaneLast(toUseIndsLast,3));
    indsInConeLast = find(coneVolume(testIndsLast));
    
    %first
    testPlaneFirst = round([coneMean(1)+planeSubs(:,1)+coneAxis(1)*firstZ, coneMean(2)+planeSubs(:,2)+coneAxis(2)*firstZ, coneMean(3)+planeSubs(:,3)+coneAxis(3)*firstZ]);
    toUseIndsFirst = find(testPlaneFirst(:,1) >= 1 & testPlaneFirst(:,1) <= newSize(1) & testPlaneFirst(:,2) >= 1 & testPlaneFirst(:,2) <= newSize(2) ...
           & testPlaneFirst(:,3) >= 1 & testPlaneFirst(:,3) <= newSize(3));
    
    testIndsFirst = sub2ind(newSize, testPlaneFirst(toUseIndsFirst,1), testPlaneFirst(toUseIndsFirst,2), testPlaneFirst(toUseIndsFirst,3));
    indsInConeFirst = find(coneVolume(testIndsFirst));
    
    toInterpInds = find(gradMap == 0 & coneVolume == 1);
    [interpX, interpY, interpZ] = ind2sub(newSize, toInterpInds);

    %border inds used previously, don't use all or interp is too large
    %good to include these in interp as some border points missed in close
    bordersUsedInds = find(bordersInVolume == 1);
    [tempX, tempY, tempZ] = ind2sub(newSize, bordersUsedInds(1:100:length(bordersUsedInds)));
    
    %use nearest neighbour interp to carry values before cone in cone forward
    toInterpSubs = [testPlaneLast(toUseIndsLast(indsInConeLast),:)', testPlaneFirst(toUseIndsFirst(indsInConeFirst),:)', [tempX, tempY, tempZ]']';
    toInterpVals = [gradMap(testIndsLast(indsInConeLast))' gradMap(testIndsFirst(indsInConeFirst))' ones(1,length(tempX))*numDiscSteps]';
    
    %figure; plot3(toInterpSubs(:,1), toInterpSubs(:,2), toInterpSubs(:,3),'.'); hold on
    %plot3(interpX, interpY, interpZ,'r.');
    
    interpFn = scatteredInterpolant(toInterpSubs(:,1), toInterpSubs(:,2), toInterpSubs(:,3), toInterpVals, 'nearest', 'nearest');
    interpVals = interpFn(interpX, interpY, interpZ);
    gradMap(toInterpInds) = interpVals;
  
    %filter a bit to clean up borders
    gradMap = medfilt3(gradMap, [5,5,5]).*coneVolume;
    figure; imshow(reshape(gradMap(:,round(coneMean(3)),:),[newSize([1, 3])])/max(gradMap(:))); hold on



 %toInterpInds = find(gradMap(:,:,round(coneMean(3))) == 0 & coneVolume(:,:,round(coneMean(3))) == 1);
 %[interpX, interpY] = ind2sub(newSize(1:2), toInterpInds);
 %plot(interpY,interpX,'rx')
unique(gradMap(:))'
    %%
    %now get interfaces - places subscripts in cells
    surfaceOfGradient = cell(numDiscSteps-1,1);
    %place interfaces in border volume

    borderVolume = zeros(size(imageStack));

    cols = lines(numDiscSteps);
    %go until second highest value, highest is on outside...
    for i =  1:numDiscSteps - 1
        tempMap = gradMap*0;
        %take outer interface
        inds2Use = find(gradMap <= i & gradMap > 0);
        tempMap(inds2Use) = 1;
        tempMap = imdilate(tempMap, strel('sphere',1))-tempMap; %should just leave border
        inds2Use = find(tempMap);
        
        %remove cornea inds
        inds2Use = setdiff(inds2Use,ToCorneaInds);
        
        %code internal borders with negative values
        borderVolume(inds2Use) = -i;
        [tempX, tempY, tempZ] = ind2sub(newSize, inds2Use);
        surfaceOfGradient{i} = [tempX, tempY, tempZ];
        
        %[tempX, tempY, tempZ] = ind2sub(newSize, inds2Use);
        %plot3(tempX(1:dsP:length(tempX)), tempY(1:dsP:length(tempX)), tempZ(1:dsP:length(tempX)), '.', 'color', cols(i,:))
        
        %figure; imshow(tempMap(:,:,round(coneMean(3))));
    end
    borderVolume(ToCorneaInds) = 1; borderVolume(ToInterInds) = 2; borderVolume(ToAirInds) = 3;
    figure; imshow((reshape(borderVolume(:,round(coneMean(3)),:),[newSize([1, 3])])+5)/8);
    
    clear tempMap
%%  place refractive indexs in grad index
    %RIProfile = 1.5; %homoy
    RIProfile = fliplr(1.45:0.1/(numDiscSteps-1):1.55); %graded, last is furtherst from center
    %RIProfile = fliplr([1.8 1.6 ones(1,numDiscSteps-2)*1.5]); %homg with border . internal borders present but should be ignored
    %RIProfile = fliplr([1.8 1.6 1.45:0.1/(numDiscSteps-3):1.55]); %graded with border
    
    RIMap = gradMap;
    
    
    if length(RIProfile) ~= numDiscSteps
        error(' did not update RI profile length')  
    end
    
    for i = 1:numDiscSteps
       RIMap(gradMap == i) = RIProfile(i); 
    end
    RIMap(imageStack == 1) = 1.35*(1.4/1.28); %set cornea (from slice w/o water so scaled up such that intercone matches)
    RIMap(imageStack == 3) = 1.4; %set intercone
    RIMap(RIMap == 0) = 1.34; %set remainder to viterous body...
    
    %set interfaces to zero
    RIMap(borderVolume ~= 0) = 0;
    
    rnger = max(newSize)*0.2;
    RIMapTemp = imclose(RIMap, strel('sphere',2));    
    figure; imshow((reshape(RIMapTemp(:,round(coneMean(3)),:),[newSize([1, 3])])-1.3)/0.3); hold on
        line( [-5 3]*corneaAxis(3)*rnger+corneaMean(3), [-5 3]*corneaAxis(1)*rnger+corneaMean(1),'color', 'm','linewidth',1);
    line( [-3 3.5]*coneAxis(3)*rnger+coneMean(3),[-3 3.5]*coneAxis(1)*rnger+coneMean(1) , 'color', 'g','linewidth',1);
    cH = colorbar;
    set(cH, 'ticks', 0:1/3:1, 'TickLabels', {'1.3', '1.4', '1.5', '1.6'})
    
    figure; imshow(((permute(RIMapTemp(round(coneMean(1)),:,:), [2 3 1]))-1.3)/0.3); hold on
      %% plot to check it worked
    figure;
    plot3(corneaSurfX(1:dsP:length(corneaSurfX)), corneaSurfY(1:dsP:length(corneaSurfX)), corneaSurfZ(1:dsP:length(corneaSurfX)), 'm.'); hold on; 
    plot3(interSurfX(1:dsP:length(interSurfX)), interSurfY(1:dsP:length(interSurfX)), interSurfZ(1:dsP:length(interSurfX)), 'y.'); 
    plot3(airSurfX(1:dsP/2:length(airSurfX)), airSurfY(1:dsP/2:length(airSurfX)), airSurfZ(1:dsP/2:length(airSurfX)), 'g.');  title('Surfaces found and lens axis');

    line([-3.6 1]*corneaAxis(1)*rnger+corneaMean(1),  [-3.6 1]*corneaAxis(2)*rnger+corneaMean(2), [-3.6 1]*corneaAxis(3)*rnger+corneaMean(3), 'color', 'm','linewidth',2);
    line([-3 2.8]*coneAxis(1)*rnger+coneMean(1),  [-3 2.8]*coneAxis(2)*rnger+coneMean(2), [-3 2.8]*coneAxis(3)*rnger+coneMean(3), 'color', 'g','linewidth',2);
    axis equal; 
    %axis off; set(gcf, 'Color','k', 'InvertHardcopy', 'off');
    %[caz,cel] = view(gca);
    %view(-0.85,-90)
    
    angle = atan2(norm(cross(corneaAxis,coneAxis)), dot(corneaAxis,coneAxis))/pi*180
    
%     tempInds = find(gradMap == 1);
%     [tempX, tempY, tempZ] = ind2sub(newSize, tempInds);
%     plot3(tempX(1:dsP:length(tempX)), tempY(1:dsP:length(tempX)), tempZ(1:dsP:length(tempX)), 'c.')
%     tempInds = find(gradMap == 6);Ga
%     [tempX, tempY, tempZ] = ind2sub(newSize, tempInds);
%     plot3(tempX(1:dsP:length(tempX)), tempY(1:dsP:length(tempX)), tempZ(1:dsP:length(tempX)), 'w.')
    
    %plot planes
    tempPlane = [coneMean(1)+planeSubs(:,1), coneMean(2)+planeSubs(:,2), coneMean(3)+planeSubs(:,3)];
    %plot3(tempPlane(:,1), tempPlane(:,2), tempPlane(:,3), 'k.'); hold on; 
    %plot3(tempPlane(:,1)+coneAxis(1)*maxZ, tempPlane(:,2)+coneAxis(2)*maxZ, tempPlane(:,3)+coneAxis(3)*maxZ, 'r.'); hold on; 
    %plot3(tempPlane(:,1)+coneAxis(1)*minZ, tempPlane(:,2)+coneAxis(2)*minZ, tempPlane(:,3)+coneAxis(3)*minZ, 'r.'); hold on; 
end

%preallocate a grid to project rays onto for each point source
if 1
    %parameters
    rDensity = 3;   %subdivision coeff influences points on sphere: 1 - 42, 2 - 162, 3 - 642, 4 - 2562, 5 - 10242
    raySpacing = 5; %spacing for ray intercepts
    pointDistance = 10^6;
    createSphere = 0;
    
    %create sphere with equally distributed points, set to zero to project from a single point along lens axis
    if createSphere
        TR = SubdivideSphericalMesh(IcosahedronMesh,rDensity);
        sphereCords = TR.X*pointDistance;
        sphereCords = [sphereCords(:,1)+lensCenter(1), sphereCords(:,2)+lensCenter(2), sphereCords(:,3)+lensCenter(3)];

        %if coarse solution already tested figure out bounds of field of view
            %note: this run should then have a higher rDensity value than the coarse run
        if preLoadCoarseRays
            cd(folderToSave);
            coarseData = load(coarseRayFile);

            %find the mean minimum distance between the coarse points
            minDists = zeros(size(coarseData.sphereCords,1),1);
            for i = 1:size(coarseData.sphereCords,1)
                rayOrig = coarseData.rayOriginParams{i,2};
                minDs = sqrt((coarseData.sphereCords(:,1)-rayOrig(1)).^2 + (coarseData.sphereCords(:,2)-rayOrig(2)).^2 + ...
                    (coarseData.sphereCords(:,3)-rayOrig(3)).^2);
                minDs(i) = [];
                minDists(i) = min(minDs);
            end
            dRnge = mean(minDists)*1.5;

            %find good ray origins
            toKeep = [];
            for i = 1:size(coarseData.rayOriginParams,1)
                if coarseData.rayOriginParams{i,3}
                   toKeep = [toKeep i]; 
                end
            end
            coarseData.sphereCords = coarseData.sphereCords(toKeep,:);

            %now find points on our new sphere within dRnge of the good points
            %on the coarse sphere
            toKeep = [];
            for i = 1:size(sphereCords,1);
                minDs = sqrt((coarseData.sphereCords(:,1)-sphereCords(i,1)).^2 + (coarseData.sphereCords(:,2)-sphereCords(i,2)).^2 + ...
                    (coarseData.sphereCords(:,3)-sphereCords(i,3)).^2);
                if min(minDs) < dRnge
                    toKeep = [toKeep i];
                end
            end

            %plot to test
            figure; hold on
            plot3(sphereCords(:,1), sphereCords(:,2), sphereCords(:,3), '.');
            sphereCords = sphereCords(toKeep,:); 
            plot3(coarseData.sphereCords(:,1), coarseData.sphereCords(:,2), coarseData.sphereCords(:,3), 'rx');
            plot3(sphereCords(:,1), sphereCords(:,2), sphereCords(:,3), 'og'); title(sprintf('%i new rays to sample', size(sphereCords,1)));
        end
        
        %set to zero to make grid across lens from each point source
        %set to one to just test a line of points across lens
        radialArrangment = 0; %probably leave in 0 for this mode
        
    else 
        %for single origin - mostly used for debug
        radialArrangment = 1; %may want to set to 1 to test grid in this mode
        
        rayRadiuses = (-1:0.15:1)*200; %defines radius at which points are put
        %be careful with input radiuses as there is no check they are on surface
        
        oneLineSet = 0; %draws points along a single line
        quadrantSet = 0; %draws line of points given angular spacing across just one quarant
        %if neither set will draw lines of points at angular spacing across entire surface
        angularSpacing = 30/180*pi; 
        
        sphereCords = coneAxis'*pointDistance; %set a single point at infinity on lens axis
        %sphereCords = corneaAxis'*pointDistance;
        
        sphereCords = [sphereCords(:,1)+coneMean(1), sphereCords(:,2)+coneMean(2), sphereCords(:,3)+coneMean(3)];
    end
 
    rayOriginParams = cell(size(sphereCords,1),3);
    %holds almost all of the data from the ray tracer, x index is ray origin
    %y 1 is an array of rayObject structures
    %y 2 is the coords of the point
    %y 3 is the number of rays from this point entering the retina
    
    for i = 1:size(sphereCords,1)
        rayOriginParams{i,2} = sphereCords(i,:); 
    
        ptsGridX = []; ptsGridY = [];

        %choose how to set points depending on flags
        if oneLineSet
            angles = 90/180*pi;
        else if quadrantSet
                angles = 0:angularSpacing:2*pi/4;
            else if ~oneLineSet & ~quadrantSet
                angles = 0:angularSpacing:2*pi; 
        end;end;end

        %draw points on grid
        xPts = sin(angles);
        yPts = cos(angles);
        for j = 1:length(rayRadiuses)
            ptsGridX = [ptsGridX xPts*rayRadiuses(j)];
            ptsGridY = [ptsGridY yPts*rayRadiuses(j)];
        end
        ptsGridX = (ptsGridX)'; ptsGridY = (ptsGridY)';
        ptsGridZ = zeros(length(ptsGridY),1);

        %calculate starting parameters of rays and save into rayObject structures
        ptsGridNew = [ptsGridX ptsGridY ptsGridZ]*rerToCone;
        ptsGridNew = [ptsGridNew(:,1)+coneMean(1) ptsGridNew(:,2)+coneMean(2) ptsGridNew(:,3)+coneMean(3)];
        clear raysFromPoint;
        raysFromPoint(size(ptsGridNew,1),:) = rayObject;

        for j = 1:size(ptsGridNew,1)
            raysFromPoint(j).segmentStartPoints(1,:) = sphereCords(i,:);
            vec = [ptsGridNew(j,1)-sphereCords(i,1) ptsGridNew(j,2)-sphereCords(i,2) ptsGridNew(j,3)-sphereCords(i,3)];
            raysFromPoint(j).segmentVectors(1,:) = vec./sqrt(vec(1)^2+vec(2)^2+vec(3)^2);
        end
        rayOriginParams{i,1} = raysFromPoint;

        rayOriginParams{i,3} = size(ptsGridNew,1);
        
        if rem(i,10) == 0
            display(sprintf('rays computed for %i of %i ray origins', i, size(sphereCords,1)));
        end
    end
    rayOriginParamsOld = rayOriginParams;
else
    %if not computed this time, use the previous ones
    rayOriginParams = rayOriginParamsOld;
end

%ray tracing section (finally)
if 1
    jumpForward = 3; %really needed to avoid getting stuck...
    
    surfaceFitSize = 30; %put somewhere and modify relevent file...
    
    %test plots
    debugPlot = 0;

    %set up parameters for line drawing in volume
    grid3D.nx = xFullRange(2);
    grid3D.ny = yFullRange(2);
    grid3D.nz = zFullRange(2);
    grid3D.minBound = [1 1 1]';
    grid3D.maxBound = [xFullRange(2), yFullRange(2), zFullRange(2)]';

    %process for each ray origin
    for i = 1:size(sphereCords,1); 
        raysFromPoint = rayOriginParams{i,1};

        if debugPlot
            figure; hold on
            plot3(corneaSurfX(1:dsP:length(corneaSurfX)), corneaSurfY(1:dsP:length(corneaSurfX)), corneaSurfZ(1:dsP:length(corneaSurfX)), 'b.'); hold on; 
            plot3(interSurfX(1:dsP:length(interSurfX)), interSurfY(1:dsP:length(interSurfX)), interSurfZ(1:dsP:length(interSurfX)), 'c.'); 
            plot3(airSurfX(1:dsP/2:length(airSurfX)), airSurfY(1:dsP/2:length(airSurfX)), airSurfZ(1:dsP/2:length(airSurfX)), 'g.'); 
            %scatter3(raysFromPoint(1).segmentStartPoints(1,1), raysFromPoint(1).segmentStartPoints(1,2), raysFromPoint(1).segmentStartPoints(1,3), 50, 'k');
            axis equal
        end

        toKeep = [];
        
        %now process for each ray from specific origin

        for j = 1:rayOriginParams{i,3} 
            goodRay = 1; % assume it is good until proven bad

            %test making one border vol per ray origin - dodgy hack
            tempBorderVol = borderVolume;
            
%             line([raysFromPoint(j).segmentVectors(1,1)*1000 0]+raysFromPoint(j).segmentStartPoints(1,1),...
%                     [raysFromPoint(j).segmentVectors(1,2)*1000 0]+raysFromPoint(j).segmentStartPoints(1,2), ...
%                     [raysFromPoint(j).segmentVectors(1,3)*1000 0]+raysFromPoint(j).segmentStartPoints(1,3))

            %draw a line given previously calculated vector from rayorigin to intersect on projected surface of the lens
            [xCords, yCords, zCords] = amanatidesWooAlgorithm_preAlloc(raysFromPoint(j).segmentStartPoints(1,:), raysFromPoint(j).segmentVectors(1,:), grid3D, 0); %takes this in volume coordinates
            
            %check ray is in volume (not guaranteed) 
            if ~isnan(xCords)
                %get indices where ray intersects lens
                rayInds = sub2ind(newSize, xCords, yCords, zCords); 
                corneaIntersects = find(tempBorderVol(rayInds) == 1);
                coneIntersect = find(tempBorderVol(rayInds) > 1);
                
                %check cornea is intersected
                if isempty(corneaIntersects)
                    goodRay = 0;
                else
                   %if cone is intersected, check cornea is first
                   if ~isempty(coneIntersect)
                      if coneIntersect(1) < corneaIntersects(1)
                         goodRay = 0; 
                      end
                   end
                end

                if goodRay
                    %the ray has passed the cornea, its good to rock and roll.
                    if debugPlot
                        plot3(xCords(1:corneaIntersects(1)), yCords(1:corneaIntersects(1)), zCords(1:corneaIntersects(1)), 'k');  
                        plot3(xCords(1), yCords(1), zCords(1), 'ro');  
                    end
                    
                    %refract first time

                    %take last ray ind above zero before cornea inersect to get outer RI
                    tempInd = find(RIMap(rayInds(1:corneaIntersects(1))) > 0); %ok, its always cornea...
                    R1 = RIMap(rayInds(tempInd(end)));
                    %and 1st ind above zero after cornea intersect for inner RI
                    tempInd = find(RIMap(rayInds(corneaIntersects(1):end)) > 0);
                    R2 = RIMap(rayInds(tempInd(1)+corneaIntersects(1)-1));
                 
                    rayInt = [xCords(corneaIntersects(1)), yCords(corneaIntersects(1)), zCords(corneaIntersects(1))];
                    
                    %remove other cornea intersects on the initial ray from this line... 
                    %%% still needed, but should be obselet based on correct refrecative index position
                    if length(corneaIntersects) > 1
                        tempDiff = diff(corneaIntersects);
                        tempStep = find(tempDiff > 1);
                        if ~isempty(tempStep)
                            %step, wipe up until first val greater than 1 (offset 1 given diff)
                            tempBorderVol(rayInds(1:corneaIntersects(tempStep(1)))) = 0;
                        else
                            %no steps blank all
                            tempBorderVol(rayInds(corneaIntersects)) = 0;
                        end
                    else
                        tempBorderVol(rayInds(corneaIntersects)) = 0;
                    end
                    
                    %calculate angle of refraction at interface
                    [interceptPoint, outRay, transmission] = AngleOfRefractedRay(rayInt, raysFromPoint(j).segmentStartPoints(1,:), raysFromPoint(j).segmentVectors(1,:),...
                        corneaSurf, R1, R2, NaN, -raysFromPoint(j).segmentVectors(1,:)', surfaceFitSize);
                    
                    %save values for next segment of the ray
                    raysFromPoint(j).rayStrength = raysFromPoint(j).rayStrength*transmission;
                    raysFromPoint(j).segmentStartPoints(2,:) = interceptPoint;
                    raysFromPoint(j).segmentVectors(2,:) = outRay;
                    raysFromPoint(j).currentSegment = 2;
                    raysFromPoint(j).RIs(:,2) = [R1, R2]; 
                end
            else
               goodRay = 0; 
            end
            
            %continue in loop until ray leaves cones
            while goodRay
                segNum = raysFromPoint(j).currentSegment;    
                %next ray segment, draw a new line
                %need to stepforward a bit so last point not interersected.
                [xCords, yCords, zCords] = amanatidesWooAlgorithm_preAlloc(raysFromPoint(j).segmentStartPoints(segNum,:)+raysFromPoint(j).segmentVectors(segNum,:)*0, ... %jumpForward
                        raysFromPoint(j).segmentVectors(segNum,:), grid3D, 0); %takes this in volume coordinates
   
                rayInds = sub2ind(newSize, xCords, yCords, zCords); 
                intersects = find(tempBorderVol(rayInds) ~= 0);
                
                %take r1 from previous steps R2 to continue refraction
                R1 = raysFromPoint(j).RIs(2,segNum);
                
                %test that intersect is followed by different RI value to R1
                while ~isempty(intersects)
                    %get first RI afer value
                    tempInd = find(RIMap(rayInds(intersects(1):end)) > 0);
                    R2 = RIMap(rayInds(tempInd(1)+intersects(1)-1));
                    if R1 == R2
                       %remove first set of continous intersects
                       tempDiff = diff(intersects);
                       tempStep = find(tempDiff > 1);
                       if ~isempty(tempStep)
                           %remove up to first step in diff
                           intersects(1:tempStep(1)) = [];
                       else
                           intersects = [];
                       end
                    else
                       %we are ok
                       break;
                   end
                end
                
                if ~isempty(intersects)
                    if debugPlot
                        plot3(xCords(1:intersects(1)), yCords(1:intersects(1)), zCords(1:intersects(1)), 'm');  
                        plot3(xCords(1), yCords(1), zCords(1), 'rx');  
                    end
                    
                    %and 1st ind above zero after intersect for inner RI
                    tempInd = find(RIMap(rayInds(intersects(1):end)) > 0);
                    R2 = RIMap(rayInds(tempInd(1)+intersects(1)-1));
                    
                    %figure; imshow((borderVolume(:,:,round(zCords(intersects(1))))+5)/8);
                    %hold on; plot(yCords(intersects(1)),xCords(intersects(1)),'rx');
                    %figure; plot(RIMap(rayInds(intersects(1):end))')
                    
                    if R1 == R2
                       error('imaginary fucking interface') 
                    end
                    
                    rayInt = [xCords(intersects(1)), yCords(intersects(1)), zCords(intersects(1))];
                    
                    %now need to determine which surface to use
                    if tempBorderVol(rayInds(intersects(1))) < 0
                        %on inner gradient, get out of cell array
                        surface2Use = surfaceOfGradient{-tempBorderVol(rayInds(intersects(1)))};
                    elseif tempBorderVol(rayInds(intersects(1))) > 1
                        %refraction on cone surface
                        surface2Use = fullConeSurf;
                    elseif tempBorderVol(rayInds(intersects(1))) == 1
                        %probably shouldn't get ack to cornea
                        surface2Use = corneaSurf;
                    else
                       error('no surf selected'); 
                    end
                    
                     %remove extra voxels in this intersect
                    if length(intersects) > 1
                        tempDiff = diff(intersects);
                        tempStep = find(tempDiff > 1);
                        if ~isempty(tempStep)
                            %step, wipe up until first val greater than 1 (offset 1 given diff)
                            tempBorderVol(rayInds(1:intersects(tempStep(1)))) = 0;
                        else
                            %no steps blank all
                            tempBorderVol(rayInds(intersects)) = 0;
                        end
                    else
                        tempBorderVol(rayInds(intersects)) = 0;
                    end
                    
                    %calculate angle of refraction at interface
                    [interceptPoint, outRay, transmission] = AngleOfRefractedRay(rayInt, raysFromPoint(j).segmentStartPoints(segNum,:), raysFromPoint(j).segmentVectors(segNum,:),...
                        surface2Use, R1, R2, NaN, -raysFromPoint(j).segmentVectors(segNum,:)', surfaceFitSize);

                    %save values for next segment of the ray
                    raysFromPoint(j).rayStrength = raysFromPoint(j).rayStrength*transmission;
                    raysFromPoint(j).segmentStartPoints(segNum+1,:) = interceptPoint;
                    raysFromPoint(j).segmentVectors(segNum+1,:) = outRay;
                    raysFromPoint(j).currentSegment = segNum+1;
                    raysFromPoint(j).RIs(:,segNum+1) = [R1, R2]; 
                    
                else
                   %no intersects, the ray has left the cone
                   goodRay = 0;
                   
                   if debugPlot & ~isempty(xCords)
                        plot3(xCords, yCords, zCords, 'k');  
                        plot3(xCords(1), yCords(1), zCords(1), 'ro');  
                   end
                end
                    
            end
            
            %list the good rays
            if raysFromPoint(j).currentSegment > 1
                toKeep = [toKeep j];  
            end
        end

        %only keep good rays
        if ~isempty(toKeep) 
            rayOriginParams{i,1} = raysFromPoint(toKeep) ;
            rayOriginParams{i,3} = length(toKeep);
            if debugPlot; title('kept!'); end
        else
            rayOriginParams{i,1} = [];
            rayOriginParams{i,3} = 0;
        end

        if debugPlot; close all; end %but a breakpoint here to view each debug plot
        if rem(i,10) == 0
            display(sprintf('%i of %i ray origins tested', i, size(sphereCords,1)));
        end
    end
end
%%
%plot ray outputs, only really useful if single point used
if ~createSphere
    rayStepSize = 0.1;
    display1stCrossing = 0; %tries to calculate where rays cross the optical axis

    ier = size(sphereCords,1);

    figure; hold on; axis equal
    rer = matrix2rotatevectors([0 -1 0]',[0 -1 0]'); %coneAxis');
    tempCornea = [corneaSurfX-corneaMean(1), corneaSurfY-corneaMean(2), corneaSurfZ-corneaMean(3)]*rer;

    yOffset = 0; min(tempCornea(:,2));
    plot3(tempCornea(1:dsP:length(corneaSurfX),1), tempCornea(1:dsP:length(corneaSurfX),2)-yOffset, tempCornea(1:dsP:length(corneaSurfX),3), 'm.');

    tempInter = [interSurfX-corneaMean(1), interSurfY-corneaMean(2), interSurfZ-corneaMean(3)]*rer;
    plot3(tempInter(1:dsP:length(interSurfX),1), tempInter(1:dsP:length(interSurfX),2)-yOffset, tempInter(1:dsP:length(interSurfX),3), 'y.');

    tempInner = [airSurfX-corneaMean(1), airSurfY-corneaMean(2), airSurfZ-corneaMean(3)]*rer;
    plot3(tempInner(1:dsP:length(airSurfX),1), tempInner(1:dsP:length(airSurfX),2)-yOffset, tempInner(1:dsP:length(airSurfX),3), 'g.');

    dBehind = 1000;

    if rayOriginParams{ier,3} > 0
        raysFromPoint = rayOriginParams{ier,1};
        cols = jet(rayOriginParams{ier,3});

        xCrossing = zeros(rayOriginParams{ier,3},1);

        for i = 1:rayOriginParams{ier,3}
            
            %1st segment
            distToLens = sqrt((raysFromPoint(i).segmentStartPoints(1,1)-corneaMean(1)).^2+(raysFromPoint(i).segmentStartPoints(1,2)-corneaMean(2)).^2+(raysFromPoint(i).segmentStartPoints(1,3)-corneaMean(3)).^2);

            raypath = ((distToLens-maxD):rayStepSize:distToLens)'*raysFromPoint(i).segmentVectors(1,:);
            raypath = ([raypath(:,1)+raysFromPoint(i).segmentStartPoints(1,1),  raypath(:,2)+raysFromPoint(i).segmentStartPoints(1,2),...
                raypath(:,3)+raysFromPoint(i).segmentStartPoints(1,3)]);

            ptDist = sqrt((raysFromPoint(i).segmentStartPoints(2,1)-raypath(:,1)).^2+(raysFromPoint(i).segmentStartPoints(2,2)-raypath(:,2)).^2+(raysFromPoint(i).segmentStartPoints(2,3)-raypath(:,3)).^2);

            [blah, ind] = min(ptDist);
            raypath = raypath(ind-1000:ind,:);
            tempRay = [raypath(:,1)-corneaMean(1), raypath(:,2)-corneaMean(2), raypath(:,3)-corneaMean(3)]*rer;
            plot3(tempRay(:,1), tempRay(:,2)-yOffset, tempRay(:,3), 'color', 'w');

            for j = 2:raysFromPoint(i).currentSegment-1
                raypath = (1:rayStepSize:maxD)'*raysFromPoint(i).segmentVectors(j,:);
                
                raypath = ([raypath(:,1)+raysFromPoint(i).segmentStartPoints(j,1),  raypath(:,2)+raysFromPoint(i).segmentStartPoints(j,2),...
                    raypath(:,3)+raysFromPoint(i).segmentStartPoints(j,3)]);

                actualDist = sqrt((raysFromPoint(i).segmentStartPoints(j,1)-raysFromPoint(i).segmentStartPoints(j+1,1)).^2+(raysFromPoint(i).segmentStartPoints(j,2)-raysFromPoint(i).segmentStartPoints(j+1,2)).^2+(raysFromPoint(i).segmentStartPoints(j,3)-raysFromPoint(i).segmentStartPoints(j+1,3)).^2);
                ptDist = sqrt((raysFromPoint(i).segmentStartPoints(j,1)-raypath(:,1)).^2+(raysFromPoint(i).segmentStartPoints(j,2)-raypath(:,2)).^2+(raysFromPoint(i).segmentStartPoints(j,3)-raypath(:,3)).^2);

                inds = find(ptDist < actualDist);
                raypath = raypath(inds,:);
                
                tempRay = [raypath(:,1)-corneaMean(1), raypath(:,2)-corneaMean(2), raypath(:,3)-corneaMean(3)]*rer;
                plot3(tempRay(:,1), tempRay(:,2)-yOffset, tempRay(:,3), 'color', 'r');
            end
            
            %last segment
            raypath = (1:rayStepSize:maxD)'*raysFromPoint(i).segmentVectors(end,:);
            raypath = ([raypath(:,1)+raysFromPoint(i).segmentStartPoints(end,1),  raypath(:,2)+raysFromPoint(i).segmentStartPoints(end,2),...
                raypath(:,3)+raysFromPoint(i).segmentStartPoints(end,3)]);

            inds = find(raypath(:,1) >= raysFromPoint(i).segmentStartPoints(end,1));
            raypath = raypath(1:dBehind,:);
            tempRay = [raypath(:,1)-corneaMean(1), raypath(:,2)-corneaMean(2), raypath(:,3)-corneaMean(3)]*rer;
            
            %test last outer RI to determine exit location
            if raysFromPoint(i).RIs(2,end) == 1.4
                plot3(tempRay(:,1), tempRay(:,2)-yOffset, tempRay(:,3), 'color', 'c');
            else
                plot3(tempRay(:,1), tempRay(:,2)-yOffset, tempRay(:,3), 'color', 'w');
            end
        end
    end
    axis equal; axis off; set(gcf, 'Color','k', 'InvertHardcopy', 'off');
    %[caz,cel] = view(gca);
    line([-3.6 1.3]*corneaAxis(1)*rnger,  [-3.6 1.3]*corneaAxis(2)*rnger, [-3.6 1.3]*corneaAxis(3)*rnger, 'color', 'm','linewidth',2);
    line([-4 3.1]*coneAxis(1)*rnger+coneMean(1)-corneaMean(1),  [-4 3.1]*coneAxis(2)*rnger+coneMean(2)-corneaMean(2), [-4 3.1]*coneAxis(3)*rnger+coneMean(3)-corneaMean(3), 'color', 'g','linewidth',2);
 
    %%% remove all plot offsets!    
    view(-0.85,-90)
end
%% save  - ideally all paremters related to the calculation...
if saveOutPut
    cd(folderToSave);
    save(outPutFile, 'rayOriginParams', 'folderToGet', 'airIndex', 'lensIndex', 'vbIndex', 'surfaceFitSize', 'sphereCords');
end
