function [interceptPoint, outRay, transmission ] = AngleOfRefractedRay(coarseIntPt, originPt, originRay, surface, rIn, rOut, polAngle, approximateNormal, surfaceFitSize)
    %originPt and originRay should be in x,y,z coordinate vector
    %coarseIntPt is the surface voxel that is intersected.
    %polAngle is not used currently.
    %approximateNormal is used to specify the direction of ray propogation as normal direction may be flipped
    %surfaceFitSize describes how many points are incorperated into surface. 
        %This needs tuning. Large value may smooth too much, small value may be too rough.
    
    %%%% !!! need to update as in compound eye analysis, border and constraints
    
    debugPlots = 0; scler = 25;
    
    %get all points within some range
    surfInds = find(sqrt((surface(:,1)-coarseIntPt(1)).^2 + (surface(:,2)-coarseIntPt(2)).^2 + (surface(:,3)-coarseIntPt(3)).^2) < surfaceFitSize);
    %[~, orignalPointInd] = min(sqrt((surface(surfInds,1)-coarseIntPt(1)).^2 + (surface(surfInds,2)-coarseIntPt(2)).^2 + (surface(surfInds,3)-coarseIntPt(3)).^2));
    
    %use pca to align 3rd principle axis with z axis - should make it flat (enough)
    %surfPts = [surface(surfInds,1)-coarseIntPt(1),surface(surfInds,2)-coarseIntPt(2),surface(surfInds,3)-coarseIntPt(3)];
    centM = mean(surface(surfInds,:));
    surfPts = surface(surfInds,:)-centM;
    pcaVecs = pca(surfPts);

    rer = vectorRotationFromAtoB([0 0 1]',pcaVecs(:,3)');
    rer2 = vectorRotationFromAtoB(pcaVecs(:,3)',[0 0 1]'); 
    approximateNormal = approximateNormal'*rer;
    
    %test if approx normal points up, if so, rotate down
    %note here approx normal is copy of ray, so indicates ray must come from above...
%     approximateNormalOrig = approximateNormal;
%     approximateNormal = approximateNormalOrig'*rer;
%     if approximateNormal(3) < 0
%         rer = vectorRotationFromAtoB([0 0 -1]',pcaVecs(:,3)');
%         rer2 = vectorRotationFromAtoB(pcaVecs(:,3)',[0 0 -1]');
%         approximateNormal = approximateNormalOrig'*rer;
%     else
%         rer2 = vectorRotationFromAtoB(pcaVecs(:,3)',[0 0 1]');  
%     end
    
    newPts = surfPts*rer;
    
    %shift base to zero as cf tool can't translate horizontally
%     [~, minimaInd] = min(sqrt(newPts(:,1)^2+newPts(:,2)^2));
%    originOffset = newPts(minimaInd,:);
%    newPts = newPts - originOffset; 

    opts = fitoptions( 'Method', 'LinearLeastSquares' );
    
%   can't force direction because direction is unsure.
%   often produces opposite signs on sqaures without other constraint - not ideal.
%   also x/y choice is often weird
    %opts.Lower = [-Inf -Inf -Inf 0 0 0]; %
    %opts.Upper = [Inf Inf Inf Inf Inf Inf];
    
    %fit surface to rotated points - use cftool to tests different fns (poly22 usually works well)
    fittedPoly = fit(newPts(:,1:2), newPts(:,3),'poly22', opts);
    
    %test plots
    if debugPlots == 1;
        close all; figure;
        subplot(1,3,1);
        plot3(surface(surfInds,1),surface(surfInds,2),surface(surfInds,3), '.'); hold on
        scatter3(coarseIntPt(:,1), coarseIntPt(:,2), coarseIntPt(:,3),50,'g'); 
        line([coarseIntPt(1) coarseIntPt(1)+pcaVecs(1,1)], [coarseIntPt(2) coarseIntPt(2)+pcaVecs(2,1)], [coarseIntPt(3) coarseIntPt(3)+pcaVecs(3,1)], 'color', 'b')
        line([coarseIntPt(1) coarseIntPt(1)+pcaVecs(1,2)], [coarseIntPt(2) coarseIntPt(2)+pcaVecs(2,2)], [coarseIntPt(3) coarseIntPt(3)+pcaVecs(3,2)], 'color', 'g')
        line([coarseIntPt(1) coarseIntPt(1)+pcaVecs(1,3)], [coarseIntPt(2) coarseIntPt(2)+pcaVecs(2,3)], [coarseIntPt(3) coarseIntPt(3)+pcaVecs(3,3)], 'color', 'r')
        subplot(1,3,2);
        plot3(newPts(:,1), newPts(:,2), newPts(:,3),'.'); hold on
        plot(fittedPoly);
    end
    
    rotatedOriginPt = (originPt-centM)*rer;
    rotOriginRay = (originRay)*rer;
    
    %test plots
    if debugPlots == 1;
        subplot(1,3,1); 
        line([-originRay(1) 0]*scler+coarseIntPt(1), [-originRay(2) 0]*scler+coarseIntPt(2), [-originRay(3) 0]*scler+coarseIntPt(3),'color','g');
        subplot(1,3,2);
        line([-rotOriginRay(1) 0]*scler, [-rotOriginRay(2) 0]*scler, [-rotOriginRay(3) 0]*scler,'color','g');
        line([rotOriginRay(1) 0]*scler, [rotOriginRay(2) 0]*scler, [rotOriginRay(3) 0]*scler,'color','g');
    end
    
    %function to solve to calculate fine intersection of ray with surface
    function F = solveTheFn(in,poly,pt,ang)
        F = [in(3) - poly(in(1),in(2)) ;
            in(1) - pt(1) - in(4)*ang(1) ;
            in(2) - pt(2) - in(4)*ang(2) ;
            in(3) - pt(3) - in(4)*ang(3)];
    end
    options = optimoptions('fsolve','Display','none','TolFun',1e-6);
    f = @(x)solveTheFn(x,fittedPoly,rotatedOriginPt,rotOriginRay);
    
    aproxD = sqrt((originPt(1)-coarseIntPt(1))^2+(originPt(2)-coarseIntPt(2))^2+(originPt(3)-coarseIntPt(3))^2);
    
    %use optimizer to find fine intersection
    [fineInt] = fsolve(f,[0 0 0 aproxD],options);
    
%     %rotate intercept point back
     interceptPoint = (fineInt(1:3))*rer2+centM;
     %warn if big seperation
     if sqrt((interceptPoint(1)-coarseIntPt(1))^2+(interceptPoint(2)-coarseIntPt(2))^2+(interceptPoint(3)-coarseIntPt(3))^2) > 5
        %sqrt((interceptPoint(1)-coarseIntPt(1))^2+(interceptPoint(2)-coarseIntPt(2))^2+(interceptPoint(3)-coarseIntPt(3))^2)
        %warning('check point seperation');
     end
    %test plots
    if debugPlots == 1;
        subplot(1,3,1); axis equal
        scatter3(interceptPoint(1),interceptPoint(2),interceptPoint(3),5,'r');
        subplot(1,3,2); axis equal
        scatter3(fineInt(1),fineInt(2),fineInt(3),5,'r');
        line([-rotOriginRay(1)*scler 0]+fineInt(1), [-rotOriginRay(2)*scler 0]+fineInt(2), [-rotOriginRay(3)*scler 0]+fineInt(3),'color','m');
        scatter3(0,0,0,50,'g');
    end
    
    %get normal of surface
    [dX, dY] = differentiate(fittedPoly, fineInt(1), fineInt(2));
    if isnan(dX)
        dX = 0;
    end
    if isnan(dY)
        dY = 0;
    end
    
    N = [dX, dY, -1]; 
    N = N/sqrt(N(1)^2+N(2)^2+N(3)^2);
    
    %update this from compound eye
    
    %flip normal to same side as ray propogation if they are more than 90deg apart
    if sqrt((approximateNormal(1)-N(1))^2+(approximateNormal(2)-N(2))^2+(approximateNormal(3)-N(3))^2) > sqrt(2)
        N = N*-1;
    end
    
    %test plots
    if debugPlots == 1;
        subplot(1,3,3);
        plot(fittedPoly); hold on; axis equal
        scatter3(fineInt(1),fineInt(2),fineInt(3),5,'r'); xlabel('x'); ylabel('y'); zlabel('z');
        line([0 N(1)]*scler+fineInt(1), [0 N(2)]*scler+fineInt(2), [0 N(3)]*scler+fineInt(3),'color','b');
        line([-rotOriginRay(1) 0]*scler+fineInt(1), [-rotOriginRay(2) 0]*scler+fineInt(2), [-rotOriginRay(3) 0]*scler+fineInt(3),'color','g');
        %line([-N(1) 0]*scler+fineInt(1), [-N(2) 0]*scler+fineInt(2), [-N(3) 0]*scler+fineInt(3),'color','k');
        line([approximateNormal(1) 0]*scler+fineInt(1), [approximateNormal(2) 0]*scler+fineInt(2), [approximateNormal(3) 0]*scler+fineInt(3),'color','c');
    end
    
    %check if angle within 90 deg of normal (should be solved in previous step...)
    if sqrt((-rotOriginRay(1)-N(1))^2+(-rotOriginRay(2)-N(2))^2+(-rotOriginRay(3)-N(3))^2) < sqrt(2)
        %check if TIR occurs
        nRatio = rIn/rOut;
        cosI = -dot(N, rotOriginRay);
        sinT2 = nRatio^2*(1-cosI^2);
        cosT = sqrt(1-sinT2);
        
        %from Bram de Greve doc, same results
        if sinT2 < 1
           
            
        %if 1-((rIn/rOut)^2)*dot(cross(N,rotOriginRay),cross(N,rotOriginRay)) > 0
            %not sure where this is from???
            %outRay = rIn/rOut*cross(N,cross(-N,rotOriginRay))-N*sqrt(1-((rIn/rOut)^2)*dot(cross(N,rotOriginRay),cross(N,rotOriginRay)));
            if isnan(polAngle)
                    r0rth = (rIn*cosI-rOut*cosT)/(rIn*cosI+rOut*cosT);
                    rPar = (rOut*cosI-rIn*cosT)/(rOut*cosI+rIn*cosT);
            else
                error('not coded') 
            end
            
            % Squaring should really be done above...
            %check if transmission dominates
            if (r0rth*r0rth+rPar*rPar)/2 < 0.5
                % If not, do refraction
                outRay =  nRatio*rotOriginRay + (nRatio*cosI-cosT)*N;
                %determine transmission for unpolarized light
                transmission = 1-(r0rth*r0rth+rPar*rPar)/2;    
            else
                % If not, do reflection
                outRay = rotOriginRay-2*dot(N,rotOriginRay)*N;
                transmission = (r0rth*r0rth+rPar*rPar)/2;    
            end
        else
            %total internal reflection has occured
            outRay = rotOriginRay-2*dot(N,rotOriginRay)*N;
            %all light is reflected, so simple
            transmission = 1;
        end
    else
        outRay = [NaN NaN NaN];
        transmission = 0;
        warning('Problem')
    end
    
    %test plots
    if debugPlots == 1;
        subplot(1,3,3); 
        line([0 outRay(1)]*scler+fineInt(1), [0 outRay(2)]*scler+fineInt(2), [0 outRay(3)]*scler+fineInt(3),'color','r');
        %line([0 -outRay(1)]*scler+fineInt(1), [0 -outRay(2)]*scler+fineInt(2), [0 -outRay(3)]*scler+fineInt(3),'color','r');
    end   
    
    %rotate the ray back
    outRay = outRay*rer2;
    
    %test plots
    if debugPlots == 1;
        subplot(1,3,1); 
        line([0 outRay(1)]*scler+interceptPoint(1), [0 outRay(2)]*scler+interceptPoint(2), [0 outRay(3)]*scler+interceptPoint(3),'color','r');
        close all;
    end
end